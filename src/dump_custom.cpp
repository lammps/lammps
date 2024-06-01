// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "dump_custom.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_store_atom.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

// customize by adding keyword
// also customize compute_property_atom.cpp

enum{ID,MOL,PROC,PROCP1,TYPE,ELEMENT,MASS,
     X,Y,Z,XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q,MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,
     COMPUTE,FIX,VARIABLE,IVEC,DVEC,IARRAY,DARRAY};
enum{LT,LE,GT,GE,EQ,NEQ,XOR};

static constexpr int ONEFIELD = 32;
static constexpr int DELTA = 1048576;

/* ---------------------------------------------------------------------- */

DumpCustom::DumpCustom(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg), idregion(nullptr), thresh_array(nullptr), thresh_op(nullptr),
    thresh_value(nullptr), thresh_last(nullptr), thresh_fix(nullptr), thresh_fixID(nullptr),
    thresh_first(nullptr), earg(nullptr), vtype(nullptr), vformat(nullptr), columns(nullptr),
    columns_default(nullptr), choose(nullptr), dchoose(nullptr), clist(nullptr),
    field2index(nullptr), argindex(nullptr), id_compute(nullptr), compute(nullptr), id_fix(nullptr),
    fix(nullptr), id_variable(nullptr), variable(nullptr), vbuf(nullptr), id_custom(nullptr),
    custom(nullptr), custom_flag(nullptr), typenames(nullptr), header_choice(nullptr),
    pack_choice(nullptr)
{
  if (narg == 5) error->all(FLERR,"No dump {} arguments specified", style);

  clearstep = 1;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal dump {} command: output frequency must be > 0", style);

  // expand args if any have wildcard character "*"
  // ok to include trailing optional args,
  //   so long as they do not have "*" between square brackets
  // nfield may be shrunk below if extra optional args exist

  expand = 0;
  nfield = nargnew = utils::expand_args(FLERR,narg-5,&arg[5],1,earg,lmp);
  if (earg != &arg[5]) expand = 1;

  // allocate field vectors

  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];
  memory->create(field2index,nfield,"dump:field2index");
  memory->create(argindex,nfield,"dump:argindex");

  buffer_allow = 1;
  buffer_flag = 1;

  triclinic_general = 0;
  nthresh = 0;
  nthreshlast = 0;

  // computes, fixes, variables which the dump accesses

  ncompute = 0;
  nfix = 0;
  nvariable = 0;
  ncustom = 0;

  // process attributes
  // ioptional = start of additional optional args in expanded args

  ioptional = parse_fields(nfield,earg);

  if (ioptional < nfield &&
      strcmp(style,"image") != 0 && strcmp(style,"movie") != 0)
    error->all(FLERR,"Invalid attribute {} in dump {} command",earg[ioptional],style);

  // noptional = # of optional args
  // reset nfield to subtract off optional args
  // reset ioptional to what it would be in original arg list
  // only dump image and dump movie styles process optional args,
  //   they do not use expanded earg list

  int noptional = nfield - ioptional;
  nfield -= noptional;
  size_one = nfield;
  ioptional = narg - noptional;

  // atom selection arrays

  maxlocal = 0;

  // default element name for all types = C

  ntypes = atom->ntypes;
  typenames = new char*[ntypes+1];
  for (int itype = 1; itype <= ntypes; itype++)
    typenames[itype] = utils::strdup("C");

  // setup format strings

  vformat = new char*[nfield];
  std::string cols;

  cols.clear();
  for (int i = 0; i < nfield; i++) {
    if (vtype[i] == Dump::INT) cols += "%d ";
    else if (vtype[i] == Dump::DOUBLE) cols += "%g ";
    else if (vtype[i] == Dump::STRING) cols += "%s ";
    else if (vtype[i] == Dump::BIGINT) cols += BIGINT_FORMAT " ";
    vformat[i] = nullptr;
  }
  if (nfield > 0) cols.resize(cols.size()-1);
  format_default = utils::strdup(cols);

  format_column_user = new char*[nfield];
  for (int i = 0; i < nfield; i++) format_column_user[i] = nullptr;

  // setup column string

  cols.clear();
  keyword_user.resize(nfield);
  for (int iarg = 0; iarg < nfield; iarg++) {
    key2col[earg[iarg]] = iarg;
    keyword_user[iarg].clear();
    if (cols.size()) cols += " ";
    cols += earg[iarg];
  }
  columns_default = utils::strdup(cols);
}

/* ---------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
  // if wildcard expansion occurred, free earg memory from expand_args()
  // could not do in constructor, b/c some derived classes process earg

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  delete[] pack_choice;
  delete[] vtype;
  memory->destroy(field2index);
  memory->destroy(argindex);

  delete[] idregion;
  memory->destroy(thresh_array);
  memory->destroy(thresh_op);
  memory->destroy(thresh_value);
  memory->destroy(thresh_last);

  // check nfix in case all fixes have already been deleted

  for (int i = 0; i < nthreshlast; i++) {
    if (modify->nfix) modify->delete_fix(thresh_fixID[i]);
    delete[] thresh_fixID[i];
  }
  memory->sfree(thresh_fix);
  memory->sfree(thresh_fixID);
  memory->destroy(thresh_first);

  for (int i = 0; i < ncompute; i++) delete[] id_compute[i];
  memory->sfree(id_compute);
  delete[] compute;

  for (int i = 0; i < nfix; i++) delete[] id_fix[i];
  memory->sfree(id_fix);
  delete[] fix;

  for (int i = 0; i < nvariable; i++) delete[] id_variable[i];
  memory->sfree(id_variable);
  delete[] variable;
  for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
  delete[] vbuf;

  for (int i = 0; i < ncustom; i++) delete[] id_custom[i];
  memory->sfree(id_custom);
  memory->sfree(custom);
  memory->sfree(custom_flag);
  memory->destroy(choose);
  memory->destroy(dchoose);
  memory->destroy(clist);

  for (int i = 1; i <= ntypes; i++) delete[] typenames[i];
  delete[] typenames;

  if (vformat) {
    for (int i = 0; i < nfield; i++) delete[] vformat[i];
    delete[] vformat;
  }

  if (format_column_user) {
    for (int i = 0; i < nfield; i++) delete[] format_column_user[i];
    delete[] format_column_user;
  }

  delete[] columns_default;
  delete[] columns;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::init_style()
{
  // assemble ITEMS: column string from defaults and user values

  delete[] columns;
  std::string combined;
  int icol = 0;
  for (const auto &item : utils::split_words(columns_default)) {
    if (combined.size()) combined += " ";
    if (keyword_user[icol].size()) combined += keyword_user[icol];
    else combined += item;
    ++icol;
  }
  columns = utils::strdup(combined);

  // format = copy of default or user-specified line format

  delete[] format;
  if (format_line_user) format = utils::strdup(format_line_user);
  else format = utils::strdup(format_default);

  // tokenize the format string and add space at end of each format element
  // if user-specified int/float format exists, use it instead
  // if user-specified column format exists, use it instead
  // lo priority = line, medium priority = int/float, hi priority = column

  auto words = utils::split_words(format);
  if ((int) words.size() < nfield)
    error->all(FLERR,"Dump_modify format line is too short: {}", format);

  int i=0;
  for (const auto &word : words) {
    if (i >= nfield) break;
    delete[] vformat[i];

    if (format_column_user[i])
      vformat[i] = utils::strdup(std::string(format_column_user[i]) + " ");
    else if (vtype[i] == Dump::INT && format_int_user)
      vformat[i] = utils::strdup(std::string(format_int_user) + " ");
    else if (vtype[i] == Dump::DOUBLE && format_float_user)
      vformat[i] = utils::strdup(std::string(format_float_user) + " ");
    else if (vtype[i] == Dump::BIGINT && format_bigint_user)
      vformat[i] = utils::strdup(std::string(format_bigint_user) + " ");
    else vformat[i] = utils::strdup(word + " ");

    // remove trailing blank on last column's format
    if (i == nfield-1) vformat[i][strlen(vformat[i])-1] = '\0';

    ++i;
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  if (binary && domain->triclinic == 0)
    header_choice = &DumpCustom::header_binary;
  else if (binary && triclinic_general == 1)
    header_choice = &DumpCustom::header_binary_triclinic_general;
  else if (binary && domain->triclinic == 1)
    header_choice = &DumpCustom::header_binary_triclinic;
  else if (!binary && domain->triclinic == 0)
    header_choice = &DumpCustom::header_item;
  else if (!binary && triclinic_general == 1)
    header_choice = &DumpCustom::header_item_triclinic_general;
  else if (!binary && domain->triclinic == 1)
    header_choice = &DumpCustom::header_item_triclinic;

  if (binary) write_choice = &DumpCustom::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpCustom::write_string;
  else write_choice = &DumpCustom::write_lines;

  // find current ptr for each compute,fix,variable and custom atom property
  // check that fix frequency is acceptable

  for (i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i]) error->all(FLERR,"Could not find dump {} compute ID {}",style,id_compute[i]);
  }

  for (i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i]) error->all(FLERR,"Could not find dump {} fix ID {}", style, id_fix[i]);
    if (nevery % fix[i]->peratom_freq)
      error->all(FLERR,"Dump {} and fix not computed at compatible times", style);
  }

  for (i = 0; i < nvariable; i++) {
    int ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump {} variable name {}", style, id_variable[i]);
    variable[i] = ivariable;
  }

  int icustom,flag,cols;
  for (int i = 0; i < ncustom; i++) {
    icustom = atom->find_custom(id_custom[i],flag,cols);
    if (icustom < 0)
      error->all(FLERR,"Could not find dump {} atom property name", style);
    custom[i] = icustom;
    if (!flag && !cols) custom_flag[i] = IVEC;
    else if (flag && !cols) custom_flag[i] = DVEC;
    else if (!flag && cols) custom_flag[i] = IARRAY;
    else if (flag && cols) custom_flag[i] = DARRAY;
  }

  // check validity of region

  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR,"Region {} for dump {} does not exist", idregion, style);

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_header(bigint ndump)
{
  if (!header_choice) error->all(FLERR, "Must not use 'run pre no' after creating a new dump");

  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::format_magic_string_binary()
{
  // use negative ntimestep as marker for new format
  bigint fmtlen = strlen(MAGIC_STRING);
  bigint marker = -fmtlen;
  fwrite(&marker, sizeof(bigint), 1, fp);
  fwrite(MAGIC_STRING, sizeof(char), fmtlen, fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::format_endian_binary()
{
  int endian = ENDIAN;
  fwrite(&endian, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::format_revision_binary()
{
  int revision = FORMAT_REVISION;
  fwrite(&revision, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_unit_style_binary()
{
  int len = 0;
  if (unit_flag && !unit_count) {
    ++unit_count;
    len = strlen(update->unit_style);
    fwrite(&len, sizeof(int), 1, fp);
    fwrite(update->unit_style, sizeof(char), len, fp);
  } else {
    fwrite(&len, sizeof(int), 1, fp);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_columns_binary()
{
  int len = strlen(columns);
  fwrite(&len, sizeof(int), 1, fp);
  fwrite(columns, sizeof(char), len, fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_time_binary()
{
  char flag = time_flag ? 1 : 0;
  fwrite(&flag, sizeof(char), 1, fp);

  if (time_flag) {
    double t = compute_time();
    fwrite(&t, sizeof(double), 1, fp);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_format_binary()
{
  format_magic_string_binary();
  format_endian_binary();
  format_revision_binary();
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary(bigint ndump)
{
  header_format_binary();

  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);

  header_unit_style_binary();
  header_time_binary();
  header_columns_binary();

  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary_triclinic(bigint ndump)
{
  header_format_binary();

  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&boxxy,sizeof(double),1,fp);
  fwrite(&boxxz,sizeof(double),1,fp);
  fwrite(&boxyz,sizeof(double),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);

  header_unit_style_binary();
  header_time_binary();
  header_columns_binary();

  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary_triclinic_general(bigint ndump)
{
  header_format_binary();

  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  int triclinic_general_flag = 2;
  fwrite(&triclinic_general_flag,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
  fwrite(domain->avec,3*sizeof(double),1,fp);
  fwrite(domain->bvec,3*sizeof(double),1,fp);
  fwrite(domain->cvec,3*sizeof(double),1,fp);
  fwrite(domain->boxlo,3*sizeof(double),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);

  header_unit_style_binary();
  header_time_binary();
  header_columns_binary();

  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item(bigint ndump)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fmt::print(fp,"ITEM: UNITS\n{}\n",update->unit_style);
  }
  if (time_flag) fmt::print(fp,"ITEM: TIME\n{:.16}\n",compute_time());

  fmt::print(fp,"ITEM: TIMESTEP\n{}\n"
             "ITEM: NUMBER OF ATOMS\n{}\n",
             update->ntimestep, ndump);

  fmt::print(fp,"ITEM: BOX BOUNDS {}\n"
             "{:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e}\n",
             boundstr,boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi);

  fmt::print(fp,"ITEM: ATOMS {}\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item_triclinic(bigint ndump)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fmt::print(fp,"ITEM: UNITS\n{}\n",update->unit_style);
  }
  if (time_flag) fmt::print(fp,"ITEM: TIME\n{:.16}\n",compute_time());

  fmt::print(fp,"ITEM: TIMESTEP\n{}\n"
             "ITEM: NUMBER OF ATOMS\n{}\n",
             update->ntimestep, ndump);

  fmt::print(fp,"ITEM: BOX BOUNDS xy xz yz {}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e}\n",
             boundstr,boxxlo,boxxhi,boxxy,boxylo,boxyhi,boxxz,boxzlo,boxzhi,boxyz);

  fmt::print(fp,"ITEM: ATOMS {}\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item_triclinic_general(bigint ndump)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fmt::print(fp,"ITEM: UNITS\n{}\n",update->unit_style);
  }
  if (time_flag) fmt::print(fp,"ITEM: TIME\n{:.16}\n",compute_time());

  fmt::print(fp,"ITEM: TIMESTEP\n{}\nITEM: NUMBER OF ATOMS\n{}\n", update->ntimestep, ndump);

  fmt::print(fp,"ITEM: BOX BOUNDS abc origin {}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e} {:>1.16e}\n",
             boundstr,
             domain->avec[0],domain->avec[1],domain->avec[2],domain->boxlo[0],
             domain->bvec[0],domain->bvec[1],domain->bvec[2],domain->boxlo[1],
             domain->cvec[0],domain->cvec[1],domain->cvec[2],domain->boxlo[2]);

  fmt::print(fp,"ITEM: ATOMS {}\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count()
{
  int i;

  // grow choose and variable vbuf arrays if needed

  const int nlocal = atom->nlocal;
  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;

    memory->destroy(choose);
    memory->destroy(dchoose);
    memory->destroy(clist);
    memory->create(choose,maxlocal,"dump:choose");
    memory->create(dchoose,maxlocal,"dump:dchoose");
    memory->create(clist,maxlocal,"dump:clist");

    for (i = 0; i < nvariable; i++) {
      memory->destroy(vbuf[i]);
      memory->create(vbuf[i],maxlocal,"dump:vbuf");
    }
  }

  // invoke Computes for per-atom quantities
  // cannot invoke before first run, otherwise invoke if necessary

  if (ncompute) {
    for (i = 0; i < ncompute; i++) {
      if (!compute[i]->is_initialized())
        error->all(FLERR,"Dump compute ID {} cannot be invoked before initialization by a run",
          compute[i]->id);
      if (!(compute[i]->invoked_flag & Compute::INVOKED_PERATOM)) {
        compute[i]->compute_peratom();
        compute[i]->invoked_flag |= Compute::INVOKED_PERATOM;
      }
    }
  }

  // evaluate atom-style Variables for per-atom quantities

  if (nvariable)
    for (i = 0; i < nvariable; i++)
      input->variable->compute_atom(variable[i],igroup,vbuf[i],1,0);

  // choose all local atoms for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if not in group

  if (igroup) {
    int *mask = atom->mask;
    for (i = 0; i < nlocal; i++)
      if (!(mask[i] & groupbit))
        choose[i] = 0;
  }

  // un-choose if not in region

  if (idregion) {
    auto region = domain->get_region_by_id(idregion);
    region->prematch();
    double **x = atom->x;
    for (i = 0; i < nlocal; i++)
      if (choose[i] && region->match(x[i][0],x[i][1],x[i][2]) == 0)
        choose[i] = 0;
  }

  // un-choose if any threshold criterion isn't met

  if (nthresh) {
    double *ptr,*ptrhold;
    double *values;
    double value;
    int nstride,lastflag;

    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      // customize by adding to if statement

      if (thresh_array[ithresh] == ID) {
        tagint *tag = atom->tag;
        for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MOL) {
        if (!atom->molecule_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        tagint *molecule = atom->molecule;
        for (i = 0; i < nlocal; i++) dchoose[i] = molecule[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == PROC) {
        for (i = 0; i < nlocal; i++) dchoose[i] = me;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == PROCP1) {
        for (i = 0; i < nlocal; i++) dchoose[i] = me;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == TYPE) {
        int *type = atom->type;
        for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ELEMENT) {
        int *type = atom->type;
        for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MASS) {
        if (atom->rmass) {
          ptr = atom->rmass;
          nstride = 1;
        } else {
          double *mass = atom->mass;
          int *type = atom->type;
          for (i = 0; i < nlocal; i++) dchoose[i] = mass[type[i]];
          ptr = dchoose;
          nstride = 1;
        }

      } else if (thresh_array[ithresh] == X) {
        ptr = &atom->x[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == Y) {
        ptr = &atom->x[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == Z) {
        ptr = &atom->x[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == XS) {
        double **x = atom->x;
        double boxxlo = domain->boxlo[0];
        double invxprd = 1.0/domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][0] - boxxlo) * invxprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YS) {
        double **x = atom->x;
        double boxylo = domain->boxlo[1];
        double invyprd = 1.0/domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][1] - boxylo) * invyprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZS) {
        double **x = atom->x;
        double boxzlo = domain->boxlo[2];
        double invzprd = 1.0/domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][2] - boxzlo) * invzprd;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
            h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
            h_inv[3]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double xprd = domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][0] + ((image[i] & IMGMASK) - IMGMAX) * xprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double yprd = domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][1] +
            ((image[i] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double zprd = domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][2] + ((image[i] >> IMG2BITS) - IMGMAX) * zprd;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *h = domain->h;
        int xbox,ybox,zbox;
        for (i = 0; i < nlocal; i++) {
          xbox = (image[i] & IMGMASK) - IMGMAX;
          ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        }
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *h = domain->h;
        int ybox,zbox;
        for (i = 0; i < nlocal; i++) {
          ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][1] + h[1]*ybox + h[3]*zbox;
        }
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *h = domain->h;
        int zbox;
        for (i = 0; i < nlocal; i++) {
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][2] + h[2]*zbox;
        }
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double boxxlo = domain->boxlo[0];
        double invxprd = 1.0/domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][0] - boxxlo) * invxprd +
            (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == YSU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double boxylo = domain->boxlo[1];
        double invyprd = 1.0/domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] =
            (x[i][1] - boxylo) * invyprd +
            (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == ZSU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double boxzlo = domain->boxlo[2];
        double invzprd = 1.0/domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][2] - boxzlo) * invzprd +
            (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
            h_inv[5]*(x[i][1]-boxlo[1]) +
            h_inv[4]*(x[i][2]-boxlo[2]) +
            (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YSUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
            h_inv[3]*(x[i][2]-boxlo[2]) +
            (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZSUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]) +
            (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == IX) {
        imageint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == IY) {
        imageint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == IZ) {
        imageint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == VX) {
        ptr = &atom->v[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == VY) {
        ptr = &atom->v[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == VZ) {
        ptr = &atom->v[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == FX) {
        ptr = &atom->f[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == FY) {
        ptr = &atom->f[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == FZ) {
        ptr = &atom->f[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == Q) {
        if (!atom->q_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = atom->q;
        nstride = 1;
      } else if (thresh_array[ithresh] == MUX) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][0];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUY) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][1];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUZ) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][2];
        nstride = 4;
      } else if (thresh_array[ithresh] == MU) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][3];
        nstride = 4;

      } else if (thresh_array[ithresh] == RADIUS) {
        if (!atom->radius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->radius;
        nstride = 1;
      } else if (thresh_array[ithresh] == DIAMETER) {
        if (!atom->radius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        double *radius = atom->radius;
        for (i = 0; i < nlocal; i++) dchoose[i] = 2.0*radius[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == OMEGAX) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAY) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAZ) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMX) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMY) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMZ) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQX) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQZ) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == COMPUTE) {
        i = nfield + ithresh;
        if (argindex[i] == 0) {
          ptr = compute[field2index[i]]->vector_atom;
          nstride = 1;
        } else {
          ptr = &compute[field2index[i]]->array_atom[0][argindex[i]-1];
          nstride = compute[field2index[i]]->size_peratom_cols;
        }

      } else if (thresh_array[ithresh] == FIX) {
        i = nfield + ithresh;
        if (argindex[i] == 0) {
          ptr = fix[field2index[i]]->vector_atom;
          nstride = 1;
        } else {
          ptr = &fix[field2index[i]]->array_atom[0][argindex[i]-1];
          nstride = fix[field2index[i]]->size_peratom_cols;
        }

      } else if (thresh_array[ithresh] == VARIABLE) {
        i = nfield + ithresh;
        ptr = vbuf[field2index[i]];
        nstride = 1;

      } else if (thresh_array[ithresh] == IVEC) {
        i = nfield + ithresh;
        int iwhich = custom[field2index[i]];
        int *ivector = atom->ivector[iwhich];
        for (i = 0; i < nlocal; i++)
          dchoose[i] = ivector[i];
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == DVEC) {
        i = nfield + ithresh;
        int iwhich = custom[field2index[i]];
        ptr = atom->dvector[iwhich];
        nstride = 1;

      } else if (thresh_array[ithresh] == IARRAY) {
        i = nfield + ithresh;
        int iwhich = custom[field2index[i]];
        int **iarray = atom->iarray[iwhich];
        int icol = argindex[i] - 1;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = iarray[i][icol];
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == DARRAY) {
        i = nfield + ithresh;
        int iwhich = custom[field2index[i]];
        double **darray = atom->darray[iwhich];
        ptr = &darray[0][argindex[i]-1];
        nstride = atom->dcols[iwhich];
      }

      // unselect atoms that don't meet threshold criterion
      // compare to single value or values stored in threshfix
      // copy ptr attribute into thresh_fix if this is first comparison

      if (thresh_last[ithresh] < 0) {
        lastflag = 0;
        value = thresh_value[ithresh];
      } else {
        lastflag = 1;
        int ilast = thresh_last[ithresh];
        values = thresh_fix[ilast]->vstore;
        ptrhold = ptr;
        if (thresh_first[ilast]) {
          thresh_first[ilast] = 0;
          for (i = 0; i < nlocal; i++, ptr += nstride) values[i] = *ptr;
          ptr = ptrhold;
        }
      }

      if (thresh_op[ithresh] == LT) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr >= values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr >= value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == LE) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr > values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr > value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == GT) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr <= values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr <= value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == GE) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr < values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr < value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == EQ) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr != values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr != value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == NEQ) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr == values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr == value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == XOR) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if ((choose[i] && *ptr == 0.0 && values[i] == 0.0) ||
                (*ptr != 0.0 && values[i] != 0.0))
              choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if ((choose[i] && *ptr == 0.0 && value == 0.0) ||
                (*ptr != 0.0 && value != 0.0))
              choose[i] = 0;
        }
      }

      // update values stored in threshfix

      if (lastflag) {
        ptr = ptrhold;
        for (i = 0; i < nlocal; i++, ptr += nstride) values[i] = *ptr;
      }
    }
  }

  // compress choose flags into clist
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack(tagint *ids)
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);

  if (ids) {
    tagint *tag = atom->tag;
    for (int i = 0; i < nchoose; i++)
      ids[i] = tag[clist[i]];
  }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCustom::convert_string(int n, double *mybuf)
{
  int i,j;

  int offset = 0;
  int m = 0;
  for (i = 0; i < n; i++) {
    if (offset + nfield*ONEFIELD > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    for (j = 0; j < nfield; j++) {
      if (vtype[j] == Dump::INT)
        offset += sprintf(&sbuf[offset],vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE)
        offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
      else if (vtype[j] == Dump::STRING)
        offset += sprintf(&sbuf[offset],vformat[j],typenames[(int) mybuf[m]]);
      else if (vtype[j] == Dump::BIGINT)
        offset += sprintf(&sbuf[offset],vformat[j],
                          static_cast<bigint> (mybuf[m]));
      m++;
    }
    offset += sprintf(&sbuf[offset],"\n");
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_string(int n, double *mybuf)
{
  if (mybuf)
    fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_lines(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < nfield; j++) {
      if (vtype[j] == Dump::INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
      else if (vtype[j] == Dump::STRING)
        fprintf(fp,vformat[j],typenames[(int) mybuf[m]]);
      else if (vtype[j] == Dump::BIGINT)
        fprintf(fp,vformat[j],static_cast<bigint> (mybuf[m]));
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpCustom::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  has_id = 0;

  for (int iarg = 0; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_id;
      if (sizeof(tagint) == sizeof(smallint)) vtype[iarg] = Dump::INT;
      else vtype[iarg] = Dump::BIGINT;
      has_id = 1;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[iarg] = &DumpCustom::pack_molecule;
      if (sizeof(tagint) == sizeof(smallint)) vtype[iarg] = Dump::INT;
      else vtype[iarg] = Dump::BIGINT;
    } else if (strcmp(arg[iarg],"proc") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_proc;
      vtype[iarg] = Dump::INT;
    } else if (strcmp(arg[iarg],"procp1") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_procp1;
      vtype[iarg] = Dump::INT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_type;
      vtype[iarg] = Dump::INT;
    } else if (strcmp(arg[iarg],"element") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_type;
      vtype[iarg] = Dump::STRING;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_mass;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"x") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_x_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_x;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_y_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_y;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_z_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_z;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_xs_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_xs;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_ys_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_ys;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_zs_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_zs;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_xu_triclinic_general;
      else if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_xu_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_xu;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_yu_triclinic_general;
      else if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_yu_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_yu;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_zu_triclinic_general;
      else if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_zu_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_zu;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"xsu") == 0) {
      if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_xsu_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_xsu;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"ysu") == 0) {
      if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_ysu_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_ysu;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"zsu") == 0) {
      if (domain->triclinic) pack_choice[iarg] = &DumpCustom::pack_zsu_triclinic;
      else pack_choice[iarg] = &DumpCustom::pack_zsu;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"ix") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_ix;
      vtype[iarg] = Dump::INT;
    } else if (strcmp(arg[iarg],"iy") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_iy;
      vtype[iarg] = Dump::INT;
    } else if (strcmp(arg[iarg],"iz") == 0) {
      pack_choice[iarg] = &DumpCustom::pack_iz;
      vtype[iarg] = Dump::INT;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_vx_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_vx;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_vy_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_vy;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_vz_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_vz;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"fx") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_fx_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_fx;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_fy_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_fy;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_fz_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_fz;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[iarg] = &DumpCustom::pack_q;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_mux_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_mux;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_muy_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_muy;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_muz_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_muz;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[iarg] = &DumpCustom::pack_mu;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[iarg] = &DumpCustom::pack_radius;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[iarg] = &DumpCustom::pack_diameter;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_omegax_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_omegax;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_omegay_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_omegay;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_omegaz_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_omegaz;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_angmomx_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_angmomx;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_angmomy_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_angmomy;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_angmomz_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_angmomz;
      vtype[iarg] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_tqx_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_tqx;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_tqy_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_tqy;
      vtype[iarg] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      if (triclinic_general) pack_choice[iarg] = &DumpCustom::pack_tqz_triclinic_general;
      else pack_choice[iarg] = &DumpCustom::pack_tqz;
      vtype[iarg] = Dump::DOUBLE;

    // compute or fix or variable or custom vector/array

    } else {
      int n,flag,cols;
      ArgInfo argi(arg[iarg], ArgInfo::COMPUTE | ArgInfo::FIX | ArgInfo::VARIABLE |
                   ArgInfo::DNAME | ArgInfo::INAME);
      argindex[iarg] = argi.get_index1();
      auto name = argi.get_name();
      Compute *icompute = nullptr;
      Fix *ifix = nullptr;

      switch (argi.get_type()) {

      case ArgInfo::UNKNOWN:
        error->all(FLERR,"Invalid attribute {} in dump {} command",arg[iarg],style);
        break;

      case ArgInfo::NONE:
        // ignore because this may be a valid argument for a derived dump style class
        return iarg;
        break;

      // compute value = c_ID
      // if no trailing [], then arg is set to 0, else arg is int between []

      case ArgInfo::COMPUTE:
        pack_choice[iarg] = &DumpCustom::pack_compute;
        vtype[iarg] = Dump::DOUBLE;

        icompute = modify->get_compute_by_id(name);
        if (!icompute) error->all(FLERR,"Could not find dump {} compute ID: {}", style, name);
        if (icompute->peratom_flag == 0)
          error->all(FLERR,"Dump {} compute {} does not compute per-atom info", style, name);
        if (argi.get_dim() == 0 && icompute->size_peratom_cols > 0)
          error->all(FLERR,"Dump {} compute {} does not calculate per-atom vector", style, name);
        if (argi.get_dim() > 0 && icompute->size_peratom_cols == 0)
          error->all(FLERR,"Dump {} compute {} does not calculate per-atom array", style, name);
        if (argi.get_dim() > 0 && argi.get_index1() > icompute->size_peratom_cols)
          error->all(FLERR,"Dump {} compute {} vector is accessed out-of-range", style, name);

        field2index[iarg] = add_compute(name);
        break;

      // fix value = f_ID
      // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::FIX:
        pack_choice[iarg] = &DumpCustom::pack_fix;
        vtype[iarg] = Dump::DOUBLE;

        ifix = modify->get_fix_by_id(name);
        if (!ifix) error->all(FLERR,"Could not find dump {} fix ID: {}", style, name);
        if (ifix->peratom_flag == 0)
          error->all(FLERR,"Dump {} fix {} does not compute per-atom info", style, name);
        if (argi.get_dim() == 0 && ifix->size_peratom_cols > 0)
          error->all(FLERR,"Dump {} fix {} does not compute per-atom vector", style, name);
        if (argi.get_dim() > 0 && ifix->size_peratom_cols == 0)
          error->all(FLERR,"Dump {} fix {} does not compute per-atom array", style, name);
        if (argi.get_dim() > 0 && argi.get_index1() > ifix->size_peratom_cols)
          error->all(FLERR,"Dump {} fix {} vector is accessed out-of-range", style, name);

        field2index[iarg] = add_fix(name);
        break;

      // variable value = v_name

      case ArgInfo::VARIABLE:
        pack_choice[iarg] = &DumpCustom::pack_variable;
        vtype[iarg] = Dump::DOUBLE;

        n = input->variable->find(name);
        if (n < 0) error->all(FLERR,"Could not find dump {} variable name {}", style, name);
        if (input->variable->atomstyle(n) == 0)
          error->all(FLERR,"Dump {} variable {} is not atom-style variable", style, name);

        field2index[iarg] = add_variable(name);
        break;

      // custom per-atom floating point vector or array = d_ID d2_ID

      case ArgInfo::DNAME:
        pack_choice[iarg] = &DumpCustom::pack_custom;
        vtype[iarg] = Dump::DOUBLE;

        n = atom->find_custom(name,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", name);
        if (argindex[iarg] == 0) {
          if (!flag || cols)
            error->all(FLERR,"Property double vector {} for dump {} does not exist", name, style);
        } else {
          if (!flag || !cols)
            error->all(FLERR,"Property double array {} for dump {} does not exist", name, style);
          if (argindex[iarg] > atom->dcols[n])
            error->all(FLERR,"Dump {} property array {} is accessed out-of-range", style, name);
        }

        field2index[iarg] = add_custom(name,1);
        break;

      // custom per-atom integer vector or array = i_ID or i2_ID

      case ArgInfo::INAME:
        pack_choice[iarg] = &DumpCustom::pack_custom;
        vtype[iarg] = Dump::INT;

        n = atom->find_custom(name,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", name);
        if (argindex[iarg] == 0) {
          if (flag || cols)
            error->all(FLERR,"Property integer vector {} for dump {} does not exist", name, style);
        } else {
          if (flag || !cols)
            error->all(FLERR,"Property integer array {} for dump {} does not exist", name, style);
          if (argindex[iarg] > atom->icols[n])
            error->all(FLERR,"Dump {} property array {} is accessed out-of-range", style, name);
        }

        field2index[iarg] = add_custom(name,0);
        break;

      // no match

      default:
        return iarg;
        break;
      }
    }
  }

  return narg;
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_compute(const char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete[] compute;
  compute = new Compute*[ncompute+1];

  id_compute[ncompute] = utils::strdup(id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_fix(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete[] fix;
  fix = new Fix*[nfix+1];

  id_fix[nfix] = utils::strdup(id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_variable(const char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;

  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                     "dump:id_variable");
  delete[] variable;
  variable = new int[nvariable+1];
  delete[] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = nullptr;

  id_variable[nvariable] = utils::strdup(id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   add custom atom property to list used by dump
   return index of where this property is in Atom class custom lists
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_custom(const char *id, int flag)
{
  int icustom;
  for (icustom = 0; icustom < ncustom; icustom++)
    if (strcmp(id,id_custom[icustom]) == 0) break;
  if (icustom < ncustom) return icustom;

  id_custom = (char **) memory->srealloc(id_custom,(ncustom+1)*sizeof(char *),"dump:id_custom");
  custom = (int *) memory->srealloc(custom,(ncustom+1)*sizeof(int),"dump:custom");
  custom_flag = (int *) memory->srealloc(custom_flag,(ncustom+1)*sizeof(int),"dump:custom_flag");

  id_custom[ncustom] = utils::strdup(id);
  custom_flag[ncustom] = flag;
  ncustom++;

  return ncustom-1;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      delete[] idregion;
      idregion = nullptr;
    } else {
      delete[] idregion;
      if (!domain->get_region_by_id(arg[1]))
        error->all(FLERR,"Dump_modify region {} does not exist", arg[1]);
      idregion = utils::strdup(arg[1]);
    }
    return 2;
  }

  if (strcmp(arg[0],"triclinic/general") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    triclinic_general = utils::logical(FLERR,arg[1],false,lmp);
    if (triclinic_general && !domain->triclinic_general)
      error->all(FLERR,"Dump_modify triclinic/general cannot be used "
                 "if simulation box is not general triclinic");
    return 2;
  }

  if (strcmp(arg[0],"format") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify format", error);

    if (strcmp(arg[1],"none") == 0) {
      // just clear format_column_user allocated by this dump child class
      for (int i = 0; i < nfield; i++) {
        delete[] format_column_user[i];
        format_column_user[i] = nullptr;
      }
      return 2;
    }

    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify format", error);

    if (strcmp(arg[1],"int") == 0) {
      delete[] format_int_user;
      format_int_user = utils::strdup(arg[2]);
      delete[] format_bigint_user;
      int n = strlen(format_int_user) + 8;
      format_bigint_user = new char[n];
      // replace "d" in format_int_user with bigint format specifier
      // use of &str[1] removes leading '%' from BIGINT_FORMAT string
      char *ptr = strchr(format_int_user,'d');
      if (ptr == nullptr)
        error->all(FLERR,"Dump_modify int format does not contain d character");
      char str[8];
      sprintf(str,"%s",BIGINT_FORMAT);
      *ptr = '\0';
      sprintf(format_bigint_user,"%s%s%s",format_int_user,&str[1],ptr+1);
      *ptr = 'd';

    } else if (strcmp(arg[1],"float") == 0) {
      delete[] format_float_user;
      format_float_user = utils::strdup(arg[2]);

    } else {
      int i = utils::inumeric(FLERR,arg[1],false,lmp) - 1;
      if (i < 0 || i >= nfield)
        error->all(FLERR,"Unknown dump_modify format ID keyword: {}", arg[1]);
      delete[] format_column_user[i];
      format_column_user[i] = utils::strdup(arg[2]);
    }
    return 3;
  }

  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR,"Number of dump_modify element names does not match number of atom types");

    for (int i = 1; i <= ntypes; i++) delete[] typenames[i];
    delete[] typenames;
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = utils::strdup(arg[itype]);
    }
    return ntypes+1;
  }

  if (strcmp(arg[0],"refresh") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify refresh", error);
    ArgInfo argi(arg[1],ArgInfo::COMPUTE);
    if ((argi.get_type() != ArgInfo::COMPUTE) || (argi.get_dim() != 0))
      error->all(FLERR,"Illegal dump_modify command");
    if (refreshflag) error->all(FLERR,"Dump_modify can only have one refresh");

    refreshflag = 1;
    idrefresh = argi.copy_name();
    return 2;
  }

  if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify thresh", error);
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
        memory->destroy(thresh_array);
        memory->destroy(thresh_op);
        memory->destroy(thresh_value);
        thresh_array = nullptr;
        thresh_op = nullptr;
        thresh_value = nullptr;
        thresh_last = nullptr;
        for (int i = 0; i < nthreshlast; i++) {
          modify->delete_fix(thresh_fixID[i]);
          delete[] thresh_fixID[i];
        }
        thresh_fix = nullptr;
        thresh_fixID = nullptr;
        thresh_first = nullptr;
      }
      nthresh = nthreshlast = 0;
      return 2;
    }

    if (narg < 4) utils::missing_cmd_args(FLERR, "dump_modify thresh", error);

    // grow threshold arrays

    memory->grow(thresh_array,nthresh+1,"dump:thresh_array");
    memory->grow(thresh_op,(nthresh+1),"dump:thresh_op");
    memory->grow(thresh_value,(nthresh+1),"dump:thresh_value");
    memory->grow(thresh_last,(nthresh+1),"dump:thresh_last");

    // set attribute type of threshold
    // customize by adding to if statement

    if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
    else if (strcmp(arg[1],"mol") == 0) thresh_array[nthresh] = MOL;
    else if (strcmp(arg[1],"proc") == 0) thresh_array[nthresh] = PROC;
    else if (strcmp(arg[1],"procp1") == 0) thresh_array[nthresh] = PROCP1;
    else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;
    else if (strcmp(arg[1],"mass") == 0) thresh_array[nthresh] = MASS;

    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XS;
    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSTRI;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YS;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSTRI;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZS;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSTRI;

    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XU;
    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XUTRI;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YU;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YUTRI;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZU;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZUTRI;

    else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XSU;
    else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSUTRI;
    else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YSU;
    else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSUTRI;
    else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZSU;
    else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSUTRI;

    else if (strcmp(arg[1],"ix") == 0) thresh_array[nthresh] = IX;
    else if (strcmp(arg[1],"iy") == 0) thresh_array[nthresh] = IY;
    else if (strcmp(arg[1],"iz") == 0) thresh_array[nthresh] = IZ;
    else if (strcmp(arg[1],"vx") == 0) thresh_array[nthresh] = VX;
    else if (strcmp(arg[1],"vy") == 0) thresh_array[nthresh] = VY;
    else if (strcmp(arg[1],"vz") == 0) thresh_array[nthresh] = VZ;
    else if (strcmp(arg[1],"fx") == 0) thresh_array[nthresh] = FX;
    else if (strcmp(arg[1],"fy") == 0) thresh_array[nthresh] = FY;
    else if (strcmp(arg[1],"fz") == 0) thresh_array[nthresh] = FZ;

    else if (strcmp(arg[1],"q") == 0) thresh_array[nthresh] = Q;
    else if (strcmp(arg[1],"mux") == 0) thresh_array[nthresh] = MUX;
    else if (strcmp(arg[1],"muy") == 0) thresh_array[nthresh] = MUY;
    else if (strcmp(arg[1],"muz") == 0) thresh_array[nthresh] = MUZ;
    else if (strcmp(arg[1],"mu") == 0) thresh_array[nthresh] = MU;

    else if (strcmp(arg[1],"radius") == 0) thresh_array[nthresh] = RADIUS;
    else if (strcmp(arg[1],"diameter") == 0) thresh_array[nthresh] = DIAMETER;
    else if (strcmp(arg[1],"omegax") == 0) thresh_array[nthresh] = OMEGAX;
    else if (strcmp(arg[1],"omegay") == 0) thresh_array[nthresh] = OMEGAY;
    else if (strcmp(arg[1],"omegaz") == 0) thresh_array[nthresh] = OMEGAZ;
    else if (strcmp(arg[1],"angmomx") == 0) thresh_array[nthresh] = ANGMOMX;
    else if (strcmp(arg[1],"angmomy") == 0) thresh_array[nthresh] = ANGMOMY;
    else if (strcmp(arg[1],"angmomz") == 0) thresh_array[nthresh] = ANGMOMZ;
    else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
    else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
    else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;

    // compute or fix or variable or custom vector/array
    // must grow field2index and argindex arrays, since access is beyond nfield

    else {
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");

      int n,flag,cols;
      ArgInfo argi(arg[1], ArgInfo::COMPUTE | ArgInfo::FIX | ArgInfo::VARIABLE |
                   ArgInfo::DNAME | ArgInfo::INAME);
      argindex[nfield+nthresh] = argi.get_index1();
      auto name = argi.get_name();
      Compute *icompute = nullptr;
      Fix *ifix = nullptr;

      switch (argi.get_type()) {

      case ArgInfo::UNKNOWN:
        error->all(FLERR,"Invalid attribute in dump modify command");
        break;

      // compute value = c_ID
      // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::COMPUTE:
        thresh_array[nthresh] = COMPUTE;

        icompute = modify->get_compute_by_id(name);
        if (!icompute) error->all(FLERR,"Could not find dump modify compute ID {}",name);
        if (icompute->peratom_flag == 0)
          error->all(FLERR,"Dump modify compute ID {} does not compute per-atom info",name);
        if (argi.get_dim() == 0 && icompute->size_peratom_cols > 0)
          error->all(FLERR,"Dump modify compute ID {} does not compute per-atom vector",name);
        if (argi.get_index1() > 0 && icompute->size_peratom_cols == 0)
          error->all(FLERR,"Dump modify compute ID {} does not compute per-atom array",name);
        if (argi.get_index1() > 0 && argi.get_index1() > icompute->size_peratom_cols)
          error->all(FLERR,"Dump modify compute ID {} vector is not large enough",name);

        field2index[nfield+nthresh] = add_compute(name);
        break;

      // fix value = f_ID
      // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::FIX:
        thresh_array[nthresh] = FIX;

        ifix = modify->get_fix_by_id(name);
        if (!ifix) error->all(FLERR,"Could not find dump modify fix ID: {}",name);

        if (ifix->peratom_flag == 0)
          error->all(FLERR,"Dump modify fix ID {} does not compute per-atom info",name);
        if (argi.get_dim() == 0 && ifix->size_peratom_cols > 0)
          error->all(FLERR,"Dump modify fix ID {} does not compute per-atom vector",name);
        if (argi.get_index1() > 0 && ifix->size_peratom_cols == 0)
          error->all(FLERR,"Dump modify fix ID {} does not compute per-atom array",name);
        if (argi.get_index1() > 0 && argi.get_index1() > ifix->size_peratom_cols)
          error->all(FLERR,"Dump modify fix ID {} vector is not large enough",name);

        field2index[nfield+nthresh] = add_fix(name);
        break;

      // variable value = v_ID

      case ArgInfo::VARIABLE:
        thresh_array[nthresh] = VARIABLE;
        n = input->variable->find(name);
        if (n < 0) error->all(FLERR,"Could not find dump modify variable name: {}", name);
        if (input->variable->atomstyle(n) == 0)
          error->all(FLERR,"Dump modify variable {} is not atom-style variable", name);

        field2index[nfield+nthresh] = add_variable(name);
        break;

      // custom per atom floating point vector or array

      case ArgInfo::DNAME:
        n = atom->find_custom(name,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", name);
        if (argindex[nfield+nthresh] == 0) {
          if (!flag || cols)
            error->all(FLERR,"Property double vector {} for dump {} does not exist", name, style);
          thresh_array[nthresh] = DVEC;
        } else {
          if (!flag || !cols)
            error->all(FLERR,"Property double array {} for dump {} does not exist", name, style);
          if (argindex[nfield+nthresh] > atom->dcols[n])
            error->all(FLERR,"Dump {} property array {} is accessed out-of-range", style, name);
          thresh_array[nthresh] = DARRAY;
        }

        field2index[nfield+nthresh] = add_custom(name,thresh_array[nthresh]);
        break;

      // custom per atom integer vector or array

      case ArgInfo::INAME:
        n = atom->find_custom(name,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", name);
        if (argindex[nfield+nthresh] == 0) {
          if (flag || cols)
            error->all(FLERR,"Property integer vector {} for dump {} does not exist", name, style);
          thresh_array[nthresh] = IVEC;
        } else {
          if (flag || !cols)
            error->all(FLERR,"Property integer array {} for dump {} does not exist", name, style);
          if (argindex[nfield+nthresh] > atom->icols[n])
            error->all(FLERR,"Dump {} property array {} is accessed out-of-range", style, name);
          thresh_array[nthresh] = IARRAY;
        }

        field2index[nfield+nthresh] = add_custom(name,thresh_array[nthresh]);
        break;

      // no match

      default:
        error->all(FLERR,"Invalid dump_modify thresh attribute: {}", name);
        break;
      }
    }

    // set operation type of threshold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else if (strcmp(arg[2],"|^") == 0) thresh_op[nthresh] = XOR;
    else error->all(FLERR,"Invalid dump_modify thresh operator");

    // set threshold value as number or special LAST keyword
    // create FixStore to hold LAST values, should work with restart
    // id = dump-ID + nthreshlast + DUMP_STORE, fix group = dump group

    if (strcmp(arg[3],"LAST") != 0) {
      thresh_value[nthresh] = utils::numeric(FLERR,arg[3],false,lmp);
      thresh_last[nthresh] = -1;
    } else {
      thresh_fix = (FixStoreAtom **)
        memory->srealloc(thresh_fix,(nthreshlast+1)*sizeof(FixStoreAtom *),"dump:thresh_fix");
      thresh_fixID = (char **)
        memory->srealloc(thresh_fixID,(nthreshlast+1)*sizeof(char *),"dump:thresh_fixID");
      memory->grow(thresh_first,(nthreshlast+1),"dump:thresh_first");

      std::string threshid = fmt::format("{}{}_DUMP_STORE",id,nthreshlast);
      thresh_fixID[nthreshlast] = utils::strdup(threshid);
      threshid += fmt::format(" {} STORE/ATOM 1 0 0 1", group->names[igroup]);
      thresh_fix[nthreshlast] = dynamic_cast<FixStoreAtom *>(modify->add_fix(threshid));

      thresh_last[nthreshlast] = nthreshlast;
      thresh_first[nthreshlast] = 1;
      nthreshlast++;
    }

    nthresh++;
    return 4;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

double DumpCustom::memory_usage()
{
  double bytes = Dump::memory_usage();
  bytes += memory->usage(choose,maxlocal);
  bytes += memory->usage(dchoose,maxlocal);
  bytes += memory->usage(clist,maxlocal);
  bytes += memory->usage(vbuf,nvariable,maxlocal);
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpCustom::pack_compute(int n)
{
  double *vector = compute[field2index[n]]->vector_atom;
  double **array = compute[field2index[n]]->array_atom;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_atom;
  double **array = fix[field2index[n]]->array_atom;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = vector[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_custom(int n)
{
  int flag = custom_flag[field2index[n]];
  int iwhich = custom[field2index[n]];
  int index = argindex[n];

  if (flag == IVEC) {
    int *ivector = atom->ivector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = ivector[clist[i]];
      n += size_one;
    }
  } else if (flag == DVEC) {
    double *dvector = atom->dvector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = dvector[clist[i]];
      n += size_one;
    }
  } else if (flag == IARRAY) {
    index--;
    int **iarray = atom->iarray[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = iarray[clist[i]][index];
      n += size_one;
    }
  } else if (flag == DARRAY) {
    index--;
    double **darray = atom->darray[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = darray[clist[i]][index];
      n += size_one;
    }
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump custom can output
   the atom property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpCustom::pack_id(int n)
{
  tagint *tag = atom->tag;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = tag[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_molecule(int n)
{
  tagint *molecule = atom->molecule;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = molecule[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_proc(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = me;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_procp1(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = me+1;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_type(int n)
{
  int *type = atom->type;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = type[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_mass(int n)
{
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  if (rmass) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = rmass[clist[i]];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = mass[type[clist[i]]];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_x(int n)
{
  double **x = atom->x;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = x[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_y(int n)
{
  double **x = atom->x;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = x[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_z(int n)
{
  double **x = atom->x;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = x[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_x_triclinic_general(int n)
{
  double **x = atom->x;
  double xtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_coords(x[clist[i]],xtri);
    buf[n] = xtri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_y_triclinic_general(int n)
{
  double **x = atom->x;
  double xtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_coords(x[clist[i]],xtri);
    buf[n] = xtri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_z_triclinic_general(int n)
{
  double **x = atom->x;
  double xtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_coords(x[clist[i]],xtri);
    buf[n] = xtri[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xs(int n)
{
  double **x = atom->x;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (x[clist[i]][0] - boxxlo) * invxprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ys(int n)
{
  double **x = atom->x;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (x[clist[i]][1] - boxylo) * invyprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zs(int n)
{
  double **x = atom->x;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (x[clist[i]][2] - boxzlo) * invzprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xs_triclinic(int n)
{
  int j;
  double **x = atom->x;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[0]*(x[j][0]-boxlo[0]) + h_inv[5]*(x[j][1]-boxlo[1]) +
      h_inv[4]*(x[j][2]-boxlo[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ys_triclinic(int n)
{
  int j;
  double **x = atom->x;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[1]*(x[j][1]-boxlo[1]) + h_inv[3]*(x[j][2]-boxlo[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zs_triclinic(int n)
{
  double **x = atom->x;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = h_inv[2]*(x[clist[i]][2]-boxlo[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xu(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double xprd = domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = x[j][0] + ((image[j] & IMGMASK) - IMGMAX) * xprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_yu(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double yprd = domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = x[j][1] + ((image[j] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zu(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double zprd = domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = x[j][2] + ((image[j] >> IMG2BITS) - IMGMAX) * zprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *h = domain->h;
  int xbox,ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    xbox = (image[j] & IMGMASK) - IMGMAX;
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    buf[n] = x[j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_yu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *h = domain->h;
  int ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    buf[n] = x[j][1] + h[1]*ybox + h[3]*zbox;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *h = domain->h;
  int zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    buf[n] = x[j][2] + h[2]*zbox;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xu_triclinic_general(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *h = domain->h;
  double xu[3];
  int xbox,ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    xbox = (image[j] & IMGMASK) - IMGMAX;
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    xu[0] = x[j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    xu[1] = x[j][1] + h[1]*ybox + h[3]*zbox;
    xu[2] = x[j][2] + h[2]*zbox;
    domain->restricted_to_general_coords(xu);
    buf[n] = xu[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_yu_triclinic_general(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *h = domain->h;
  double xu[3];
  int xbox,ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    xbox = (image[j] & IMGMASK) - IMGMAX;
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    xu[0] = x[j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    xu[1] = x[j][1] + h[1]*ybox + h[3]*zbox;
    xu[2] = x[j][2] + h[2]*zbox;
    domain->restricted_to_general_coords(xu);
    buf[n] = xu[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zu_triclinic_general(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *h = domain->h;
  double xu[3];
  int xbox,ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    xbox = (image[j] & IMGMASK) - IMGMAX;
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    xu[0] = x[j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    xu[1] = x[j][1] + h[1]*ybox + h[3]*zbox;
    xu[2] = x[j][2] + h[2]*zbox;
    domain->restricted_to_general_coords(xu);
    buf[n] = xu[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xsu(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = (x[j][0] - boxxlo) * invxprd + (image[j] & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ysu(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = (x[j][1] - boxylo) * invyprd + (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zsu(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = (x[j][2] - boxzlo) * invzprd + (image[j] >> IMG2BITS) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xsu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[0]*(x[j][0]-boxlo[0]) + h_inv[5]*(x[j][1]-boxlo[1]) +
      h_inv[4]*(x[j][2]-boxlo[2]) + (image[j] & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ysu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[1]*(x[j][1]-boxlo[1]) + h_inv[3]*(x[j][2]-boxlo[2]) +
      (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zsu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  imageint *image = atom->image;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[2]*(x[j][2]-boxlo[2]) + (image[j] >> IMG2BITS) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ix(int n)
{
  imageint *image = atom->image;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (image[clist[i]] & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_iy(int n)
{
  imageint *image = atom->image;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (image[clist[i]] >> IMGBITS & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_iz(int n)
{
  imageint *image = atom->image;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (image[clist[i]] >> IMG2BITS) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vx(int n)
{
  double **v = atom->v;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = v[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vy(int n)
{
  double **v = atom->v;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = v[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vz(int n)
{
  double **v = atom->v;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = v[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vx_triclinic_general(int n)
{
  double **v = atom->v;
  double vtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(v[clist[i]],vtri);
    buf[n] = vtri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vy_triclinic_general(int n)
{
  double **v = atom->v;
  double vtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(v[clist[i]],vtri);
    buf[n] = vtri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vz_triclinic_general(int n)
{
  double **v = atom->v;
  double vtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(v[clist[i]],vtri);
    buf[n] = vtri[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fx(int n)
{
  double **f = atom->f;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = f[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fy(int n)
{
  double **f = atom->f;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = f[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fz(int n)
{
  double **f = atom->f;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = f[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fx_triclinic_general(int n)
{
  double **f = atom->f;
  double ftri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(f[clist[i]],ftri);
    buf[n] = ftri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fy_triclinic_general(int n)
{
  double **f = atom->f;
  double ftri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(f[clist[i]],ftri);
    buf[n] = ftri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fz_triclinic_general(int n)
{
  double **f = atom->f;
  double ftri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(f[clist[i]],ftri);
    buf[n] = ftri[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_q(int n)
{
  double *q = atom->q;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = q[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_mux(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_muy(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_muz(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_mu(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][3];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_mux_triclinic_general(int n)
{
  double **mu = atom->mu;
  double mutri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(mu[clist[i]],mutri);
    buf[n] = mutri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_muy_triclinic_general(int n)
{
  double **mu = atom->mu;
  double mutri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(mu[clist[i]],mutri);
    buf[n] = mutri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_muz_triclinic_general(int n)
{
  double **mu = atom->mu;
  double mutri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(mu[clist[i]],mutri);
    buf[n] = mutri[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_radius(int n)
{
  double *radius = atom->radius;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = radius[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_diameter(int n)
{
  double *radius = atom->radius;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = 2.0*radius[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegax(int n)
{
  double **omega = atom->omega;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = omega[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegay(int n)
{
  double **omega = atom->omega;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = omega[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegaz(int n)
{
  double **omega = atom->omega;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = omega[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegax_triclinic_general(int n)
{
  double **omega = atom->omega;
  double omegatri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(omega[clist[i]],omegatri);
    buf[n] = omegatri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegay_triclinic_general(int n)
{
  double **omega = atom->omega;
  double omegatri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(omega[clist[i]],omegatri);
    buf[n] = omegatri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegaz_triclinic_general(int n)
{
  double **omega = atom->omega;
  double omegatri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(omega[clist[i]],omegatri);
    buf[n] = omegatri[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomx(int n)
{
  double **angmom = atom->angmom;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = angmom[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomy(int n)
{
  double **angmom = atom->angmom;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = angmom[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomz(int n)
{
  double **angmom = atom->angmom;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = angmom[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomx_triclinic_general(int n)
{
  double **angmom = atom->angmom;
  double angmomtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(angmom[clist[i]],angmomtri);
    buf[n] = angmomtri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomy_triclinic_general(int n)
{
  double **angmom = atom->angmom;
  double angmomtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(angmom[clist[i]],angmomtri);
    buf[n] = angmomtri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomz_triclinic_general(int n)
{
  double **angmom = atom->angmom;
  double angmomtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(angmom[clist[i]],angmomtri);
    buf[n] = angmomtri[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqx(int n)
{
  double **torque = atom->torque;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = torque[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqy(int n)
{
  double **torque = atom->torque;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = torque[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqz(int n)
{
  double **torque = atom->torque;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = torque[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqx_triclinic_general(int n)
{
  double **torque = atom->torque;
  double tqtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(torque[clist[i]],tqtri);
    buf[n] = tqtri[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqy_triclinic_general(int n)
{
  double **torque = atom->torque;
  double tqtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(torque[clist[i]],tqtri);
    buf[n] = tqtri[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqz_triclinic_general(int n)
{
  double **torque = atom->torque;
  double tqtri[3];

  for (int i = 0; i < nchoose; i++) {
    domain->restricted_to_general_vector(torque[clist[i]],tqtri);
    buf[n] = tqtri[2];
    n += size_one;
  }
}
