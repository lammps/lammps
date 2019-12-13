/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "dump_custom.h"
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_store.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "variable.h"
#include "utils.h"

using namespace LAMMPS_NS;

// customize by adding keyword
// also customize compute_atom_property.cpp

enum{ID,MOL,PROC,PROCP1,TYPE,ELEMENT,MASS,
     X,Y,Z,XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q,MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,
     COMPUTE,FIX,VARIABLE,INAME,DNAME};
enum{LT,LE,GT,GE,EQ,NEQ,XOR};

#define INVOKED_PERATOM 8
#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCustom::DumpCustom(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg),
  idregion(NULL), thresh_array(NULL), thresh_op(NULL), thresh_value(NULL),
  thresh_last(NULL), thresh_fix(NULL), thresh_fixID(NULL), thresh_first(NULL),
  earg(NULL), vtype(NULL), vformat(NULL), columns(NULL), choose(NULL),
  dchoose(NULL), clist(NULL), field2index(NULL), argindex(NULL), id_compute(NULL),
  compute(NULL), id_fix(NULL), fix(NULL), id_variable(NULL), variable(NULL),
  vbuf(NULL), id_custom(NULL), flag_custom(NULL), typenames(NULL),
  pack_choice(NULL)
{
  if (narg == 5) error->all(FLERR,"No dump custom arguments specified");

  clearstep = 1;

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal dump custom command");

  // expand args if any have wildcard character "*"
  // ok to include trailing optional args,
  //   so long as they do not have "*" between square brackets
  // nfield may be shrunk below if extra optional args exist

  expand = 0;
  nfield = nargnew = input->expand_args(narg-5,&arg[5],1,earg);
  if (earg != &arg[5]) expand = 1;

  // allocate field vectors

  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];
  memory->create(field2index,nfield,"dump:field2index");
  memory->create(argindex,nfield,"dump:argindex");

  buffer_allow = 1;
  buffer_flag = 1;
  iregion = -1;
  idregion = NULL;

  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;
  thresh_last = NULL;

  nthreshlast = 0;
  thresh_fix = NULL;
  thresh_fixID = NULL;
  thresh_first = NULL;

  // computes, fixes, variables which the dump accesses

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;

  nvariable = 0;
  id_variable = NULL;
  variable = NULL;
  vbuf = NULL;

  ncustom = 0;
  id_custom = NULL;
  flag_custom = NULL;

  // process attributes
  // ioptional = start of additional optional args in expanded args

  ioptional = parse_fields(nfield,earg);

  if (ioptional < nfield &&
      strcmp(style,"image") != 0 && strcmp(style,"movie") != 0)
    error->all(FLERR,"Invalid attribute in dump custom command");

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
  choose = NULL;
  dchoose = NULL;
  clist = NULL;

  // default element name for all types = C

  ntypes = atom->ntypes;
  typenames = new char*[ntypes+1];
  for (int itype = 1; itype <= ntypes; itype++) {
    typenames[itype] = new char[2];
    strcpy(typenames[itype],"C");
  }

  // setup format strings

  vformat = new char*[size_one];

  format_default = new char[4*size_one+1];
  format_default[0] = '\0';

  for (int i = 0; i < size_one; i++) {
    if (vtype[i] == Dump::INT) strcat(format_default,"%d ");
    else if (vtype[i] == Dump::DOUBLE) strcat(format_default,"%g ");
    else if (vtype[i] == Dump::STRING) strcat(format_default,"%s ");
    else if (vtype[i] == Dump::BIGINT) strcat(format_default,BIGINT_FORMAT " ");
    vformat[i] = NULL;
  }

  format_column_user = new char*[size_one];
  for (int i = 0; i < size_one; i++) format_column_user[i] = NULL;

  // setup column string

  int n = 0;
  for (int iarg = 0; iarg < nfield; iarg++) n += strlen(earg[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 0; iarg < nfield; iarg++) {
    strcat(columns,earg[iarg]);
    strcat(columns," ");
  }
}

/* ---------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
  // if wildcard expansion occurred, free earg memory from expand_args()
  // could not do in constructor, b/c some derived classes process earg

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  delete [] pack_choice;
  delete [] vtype;
  memory->destroy(field2index);
  memory->destroy(argindex);

  delete [] idregion;
  memory->destroy(thresh_array);
  memory->destroy(thresh_op);
  memory->destroy(thresh_value);
  memory->destroy(thresh_last);

  // check nfix in case all fixes have already been deleted

  for (int i = 0; i < nthreshlast; i++) {
    if (modify->nfix) modify->delete_fix(thresh_fixID[i]);
    delete [] thresh_fixID[i];
  }
  memory->sfree(thresh_fix);
  memory->sfree(thresh_fixID);
  memory->destroy(thresh_first);

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  memory->sfree(id_fix);
  delete [] fix;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  memory->sfree(id_variable);
  delete [] variable;
  for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
  delete [] vbuf;

  for (int i = 0; i < ncustom; i++) delete [] id_custom[i];
  memory->sfree(id_custom);
  delete [] flag_custom;

  memory->destroy(choose);
  memory->destroy(dchoose);
  memory->destroy(clist);

  for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
  delete [] typenames;

  if (vformat) {
    for (int i = 0; i < size_one; i++) delete [] vformat[i];
    delete [] vformat;
  }

  if (format_column_user) {
    for (int i = 0; i < size_one; i++) delete [] format_column_user[i];
    delete [] format_column_user;
  }

  delete [] columns;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::init_style()
{
  // format = copy of default or user-specified line format

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;
  char *r_token;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element
  // if user-specified int/float format exists, use it instead
  // if user-specified column format exists, use it instead
  // lo priority = line, medium priority = int/float, hi priority = column

  char *ptr;
  for (int i = 0; i < size_one; i++) {
    r_token = format;
    if (i == 0) ptr = utils::strtok_r(r_token," \0",&r_token);
    else ptr = utils::strtok_r(NULL," \0",&r_token);
    if (ptr == NULL) error->all(FLERR,"Dump_modify format line is too short");
    delete [] vformat[i];

    if (format_column_user[i]) {
      vformat[i] = new char[strlen(format_column_user[i]) + 2];
      strcpy(vformat[i],format_column_user[i]);
    } else if (vtype[i] == Dump::INT && format_int_user) {
      vformat[i] = new char[strlen(format_int_user) + 2];
      strcpy(vformat[i],format_int_user);
    } else if (vtype[i] == Dump::DOUBLE && format_float_user) {
      vformat[i] = new char[strlen(format_float_user) + 2];
      strcpy(vformat[i],format_float_user);
    } else if (vtype[i] == Dump::BIGINT && format_bigint_user) {
      vformat[i] = new char[strlen(format_bigint_user) + 2];
      strcpy(vformat[i],format_bigint_user);
    } else {
      vformat[i] = new char[strlen(ptr) + 2];
      strcpy(vformat[i],ptr);
    }

    vformat[i] = strcat(vformat[i]," ");
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  if (binary && domain->triclinic == 0)
    header_choice = &DumpCustom::header_binary;
  else if (binary && domain->triclinic == 1)
    header_choice = &DumpCustom::header_binary_triclinic;
  else if (!binary && domain->triclinic == 0)
    header_choice = &DumpCustom::header_item;
  else if (!binary && domain->triclinic == 1)
    header_choice = &DumpCustom::header_item_triclinic;

  if (binary) write_choice = &DumpCustom::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpCustom::write_string;
  else write_choice = &DumpCustom::write_lines;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find dump custom compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump custom fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->peratom_freq)
      error->all(FLERR,"Dump custom and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump custom variable name");
    variable[i] = ivariable;
  }

  int icustom;
  for (int i = 0; i < ncustom; i++) {
    icustom = atom->find_custom(id_custom[i],flag_custom[i]);
    if (icustom < 0)
      error->all(FLERR,"Could not find custom per-atom property ID");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for dump custom does not exist");
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary(bigint ndump)
{
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
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary_triclinic(bigint ndump)
{
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
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item(bigint ndump)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fprintf(fp,"ITEM: UNITS\n%s\n",update->unit_style);
  }
  if (time_flag) fprintf(fp,"ITEM: TIME\n%.16g\n",compute_time());

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%-1.16e %-1.16e\n",boxxlo,boxxhi);
  fprintf(fp,"%-1.16e %-1.16e\n",boxylo,boxyhi);
  fprintf(fp,"%-1.16e %-1.16e\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item_triclinic(bigint ndump)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fprintf(fp,"ITEM: UNITS\n%s\n",update->unit_style);
  }
  if (time_flag) fprintf(fp,"ITEM: TIME\n%.16g\n",compute_time());

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
  fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxxlo,boxxhi,boxxy);
  fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxylo,boxyhi,boxxz);
  fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxzlo,boxzhi,boxyz);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count()
{
  int i;

  // grow choose and variable vbuf arrays if needed

  int nlocal = atom->nlocal;
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
  // only if within a run or minimize
  // else require that computes are current
  // this prevents a compute from being invoked by the WriteDump class

  if (ncompute) {
    if (update->whichflag == 0) {
      for (i = 0; i < ncompute; i++)
        if (compute[i]->invoked_peratom != update->ntimestep)
          error->all(FLERR,"Compute used in dump between runs is not current");
    } else {
      for (i = 0; i < ncompute; i++) {
        if (!(compute[i]->invoked_flag & INVOKED_PERATOM)) {
          compute[i]->compute_peratom();
          compute[i]->invoked_flag |= INVOKED_PERATOM;
        }
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

  if (iregion >= 0) {
    Region *region = domain->regions[iregion];
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
    int nlocal = atom->nlocal;

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
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->q;
        nstride = 1;
      } else if (thresh_array[ithresh] == MUX) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][0];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUY) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][1];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUZ) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][2];
        nstride = 4;
      } else if (thresh_array[ithresh] == MU) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
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

      } else if (thresh_array[ithresh] == DNAME) {
        int iwhich,tmp;
        i = nfield + ithresh;
        iwhich = atom->find_custom(id_custom[field2index[i]],tmp);
        ptr = atom->dvector[iwhich];
        nstride = 1;

      } else if (thresh_array[ithresh] == INAME) {
        int iwhich,tmp;
        i = nfield + ithresh;
        iwhich = atom->find_custom(id_custom[field2index[i]],tmp);

        int *ivector = atom->ivector[iwhich];
        for (i = 0; i < nlocal; i++)
          dchoose[i] = ivector[i];
        ptr = dchoose;
        nstride = 1;
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
    if (offset + size_one*ONEFIELD > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    for (j = 0; j < size_one; j++) {
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
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_lines(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
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

  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpCustom::pack_id;
      if (sizeof(tagint) == sizeof(smallint)) vtype[i] = Dump::INT;
      else vtype[i] = Dump::BIGINT;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_molecule;
      if (sizeof(tagint) == sizeof(smallint)) vtype[i] = Dump::INT;
      else vtype[i] = Dump::BIGINT;
    } else if (strcmp(arg[iarg],"proc") == 0) {
      pack_choice[i] = &DumpCustom::pack_proc;
      vtype[i] = Dump::INT;
    } else if (strcmp(arg[iarg],"procp1") == 0) {
      pack_choice[i] = &DumpCustom::pack_procp1;
      vtype[i] = Dump::INT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &DumpCustom::pack_type;
      vtype[i] = Dump::INT;
    } else if (strcmp(arg[iarg],"element") == 0) {
      pack_choice[i] = &DumpCustom::pack_type;
      vtype[i] = Dump::STRING;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[i] = &DumpCustom::pack_mass;
      vtype[i] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &DumpCustom::pack_x;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &DumpCustom::pack_y;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &DumpCustom::pack_z;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_xs_triclinic;
      else pack_choice[i] = &DumpCustom::pack_xs;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_ys_triclinic;
      else pack_choice[i] = &DumpCustom::pack_ys;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_zs_triclinic;
      else pack_choice[i] = &DumpCustom::pack_zs;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_xu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_xu;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_yu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_yu;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_zu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_zu;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"xsu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_xsu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_xsu;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"ysu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_ysu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_ysu;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"zsu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_zsu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_zsu;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"ix") == 0) {
      pack_choice[i] = &DumpCustom::pack_ix;
      vtype[i] = Dump::INT;
    } else if (strcmp(arg[iarg],"iy") == 0) {
      pack_choice[i] = &DumpCustom::pack_iy;
      vtype[i] = Dump::INT;
    } else if (strcmp(arg[iarg],"iz") == 0) {
      pack_choice[i] = &DumpCustom::pack_iz;
      vtype[i] = Dump::INT;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[i] = &DumpCustom::pack_vx;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[i] = &DumpCustom::pack_vy;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[i] = &DumpCustom::pack_vz;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"fx") == 0) {
      pack_choice[i] = &DumpCustom::pack_fx;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      pack_choice[i] = &DumpCustom::pack_fy;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      pack_choice[i] = &DumpCustom::pack_fz;
      vtype[i] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_q;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_mux;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_muy;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_muz;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_mu;
      vtype[i] = Dump::DOUBLE;

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_radius;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_diameter;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_omegax;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_omegay;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_omegaz;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_angmomx;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_angmomy;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_angmomz;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqx;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqy;
      vtype[i] = Dump::DOUBLE;
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqz;
      vtype[i] = Dump::DOUBLE;

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_compute;
      vtype[i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump custom command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump custom compute ID");
      if (modify->compute[n]->peratom_flag == 0)
        error->all(FLERR,"Dump custom compute does not compute per-atom info");
      if (argindex[i] == 0 && modify->compute[n]->size_peratom_cols > 0)
        error->all(FLERR,
                   "Dump custom compute does not calculate per-atom vector");
      if (argindex[i] > 0 && modify->compute[n]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Dump custom compute does not calculate per-atom array");
      if (argindex[i] > 0 &&
          argindex[i] > modify->compute[n]->size_peratom_cols)
        error->all(FLERR,"Dump custom compute vector is accessed out-of-range");

      field2index[i] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_fix;
      vtype[i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump custom command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump custom fix ID");
      if (modify->fix[n]->peratom_flag == 0)
        error->all(FLERR,"Dump custom fix does not compute per-atom info");
      if (argindex[i] == 0 && modify->fix[n]->size_peratom_cols > 0)
        error->all(FLERR,"Dump custom fix does not compute per-atom vector");
      if (argindex[i] > 0 && modify->fix[n]->size_peratom_cols == 0)
        error->all(FLERR,"Dump custom fix does not compute per-atom array");
      if (argindex[i] > 0 &&
          argindex[i] > modify->fix[n]->size_peratom_cols)
        error->all(FLERR,"Dump custom fix vector is accessed out-of-range");

      field2index[i] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_variable;
      vtype[i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump custom variable name");
      if (input->variable->atomstyle(n) == 0)
        error->all(FLERR,"Dump custom variable is not atom-style variable");

      field2index[i] = add_variable(suffix);
      delete [] suffix;

    // custom per-atom floating point value = d_ID

    } else if (strncmp(arg[iarg],"d_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_custom;
      vtype[i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      argindex[i] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if (n < 0)
        error->all(FLERR,"Could not find custom per-atom property ID");

      if (tmp != 1)
        error->all(FLERR,"Custom per-atom property ID is not floating point");

      field2index[i] = add_custom(suffix,1);
      delete [] suffix;

    // custom per-atom integer value = i_ID

    } else if (strncmp(arg[iarg],"i_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_custom;
      vtype[i] = Dump::INT;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      argindex[i] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if (n < 0)
        error->all(FLERR,"Could not find custom per-atom property ID");

      if (tmp != 0)
        error->all(FLERR,"Custom per-atom property ID is not integer");

      field2index[i] = add_custom(suffix,0);
      delete [] suffix;

    } else return iarg;
  }

  return narg;
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_variable(char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;

  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                     "dump:id_variable");
  delete [] variable;
  variable = new int[nvariable+1];
  delete [] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   add custom atom property to list used by dump
   return index of where this property is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_custom(char *id, int flag)
{
  int icustom;
  for (icustom = 0; icustom < ncustom; icustom++)
    if ((strcmp(id,id_custom[icustom]) == 0)
        && (flag == flag_custom[icustom])) break;
  if (icustom < ncustom) return icustom;

  id_custom = (char **)
    memory->srealloc(id_custom,(ncustom+1)*sizeof(char *),"dump:id_custom");
  flag_custom = (int *)
    memory->srealloc(flag_custom,(ncustom+1)*sizeof(int),"dump:flag_custom");

  int n = strlen(id) + 1;
  id_custom[ncustom] = new char[n];
  strcpy(id_custom[ncustom],id);
  flag_custom[ncustom] = flag;

  ncustom++;
  return ncustom-1;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = domain->find_region(arg[1]);
      if (iregion == -1)
        error->all(FLERR,"Dump_modify region ID does not exist");
      delete [] idregion;
      int n = strlen(arg[1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[1]);
    }
    return 2;
  }

  if (strcmp(arg[0],"format") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");

    if (strcmp(arg[1],"none") == 0) {
      // just clear format_column_user allocated by this dump child class
      for (int i = 0; i < size_one; i++) {
        delete [] format_column_user[i];
        format_column_user[i] = NULL;
      }
      return 2;
    }

    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");

    if (strcmp(arg[1],"int") == 0) {
      delete [] format_int_user;
      int n = strlen(arg[2]) + 1;
      format_int_user = new char[n];
      strcpy(format_int_user,arg[2]);
      delete [] format_bigint_user;
      n = strlen(format_int_user) + 8;
      format_bigint_user = new char[n];
      // replace "d" in format_int_user with bigint format specifier
      // use of &str[1] removes leading '%' from BIGINT_FORMAT string
      char *ptr = strchr(format_int_user,'d');
      if (ptr == NULL)
        error->all(FLERR,
                   "Dump_modify int format does not contain d character");
      char str[8];
      sprintf(str,"%s",BIGINT_FORMAT);
      *ptr = '\0';
      sprintf(format_bigint_user,"%s%s%s",format_int_user,&str[1],ptr+1);
      *ptr = 'd';

    } else if (strcmp(arg[1],"float") == 0) {
      delete [] format_float_user;
      int n = strlen(arg[2]) + 1;
      format_float_user = new char[n];
      strcpy(format_float_user,arg[2]);

    } else {
      int i = force->inumeric(FLERR,arg[1]) - 1;
      if (i < 0 || i >= size_one)
        error->all(FLERR,"Illegal dump_modify command");
      if (format_column_user[i]) delete [] format_column_user[i];
      int n = strlen(arg[2]) + 1;
      format_column_user[i] = new char[n];
      strcpy(format_column_user[i],arg[2]);
    }
    return 3;
  }

  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR,"Dump_modify element names do not match atom types");

    for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
    delete [] typenames;
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      int n = strlen(arg[itype]) + 1;
      typenames[itype] = new char[n];
      strcpy(typenames[itype],arg[itype]);
    }
    return ntypes+1;
  }

  if (strcmp(arg[0],"refresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strncmp(arg[1],"c_",2) != 0)
      error->all(FLERR,"Illegal dump_modify command");
    if (refreshflag) error->all(FLERR,"Dump modify can only have one refresh");

    refreshflag = 1;
    int n = strlen(arg[1]);
    refresh = new char[n];
    strcpy(refresh,&arg[1][2]);
    return 2;
  }

  if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
        memory->destroy(thresh_array);
        memory->destroy(thresh_op);
        memory->destroy(thresh_value);
        thresh_array = NULL;
        thresh_op = NULL;
        thresh_value = NULL;
        thresh_last = NULL;
        for (int i = 0; i < nthreshlast; i++) {
          modify->delete_fix(thresh_fixID[i]);
          delete [] thresh_fixID[i];
        }
        thresh_fix = NULL;
        thresh_fixID = NULL;
        thresh_first = NULL;
      }
      nthresh = nthreshlast = 0;
      return 2;
    }

    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");

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

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2index and argindex arrays, since access is beyond nfield

    else if (strncmp(arg[1],"c_",2) == 0) {
      thresh_array[nthresh] = COMPUTE;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump modify command");
        argindex[nfield+nthresh] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nfield+nthresh] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify compute ID");

      if (modify->compute[n]->peratom_flag == 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom info");
      if (argindex[nfield+nthresh] == 0 &&
          modify->compute[n]->size_peratom_cols > 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom vector");
      if (argindex[nfield+nthresh] > 0 &&
          modify->compute[n]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom array");
      if (argindex[nfield+nthresh] > 0 &&
          argindex[nfield+nthresh] > modify->compute[n]->size_peratom_cols)
        error->all(FLERR,"Dump modify compute ID vector is not large enough");

      field2index[nfield+nthresh] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"f_",2) == 0) {
      thresh_array[nthresh] = FIX;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump modify command");
        argindex[nfield+nthresh] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nfield+nthresh] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify fix ID");

      if (modify->fix[n]->peratom_flag == 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom info");
      if (argindex[nfield+nthresh] == 0 &&
          modify->fix[n]->size_peratom_cols > 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom vector");
      if (argindex[nfield+nthresh] > 0 &&
          modify->fix[n]->size_peratom_cols == 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom array");
      if (argindex[nfield+nthresh] > 0 &&
          argindex[nfield+nthresh] > modify->fix[n]->size_peratom_cols)
        error->all(FLERR,"Dump modify fix ID vector is not large enough");

      field2index[nfield+nthresh] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_ID
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"v_",2) == 0) {
      thresh_array[nthresh] = VARIABLE;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      argindex[nfield+nthresh] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify variable name");
      if (input->variable->atomstyle(n) == 0)
        error->all(FLERR,"Dump modify variable is not atom-style variable");

      field2index[nfield+nthresh] = add_variable(suffix);
      delete [] suffix;

    // custom per atom floating point value = d_ID
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"d_",2) == 0) {
      thresh_array[nthresh] = DNAME;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
      argindex[nfield+nthresh] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if ((n < 0) || (tmp != 1))
        error->all(FLERR,"Could not find dump modify "
                   "custom atom floating point property ID");

      field2index[nfield+nthresh] = add_custom(suffix,1);
      delete [] suffix;

    // custom per atom integer value = i_ID
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"i_",2) == 0) {
      thresh_array[nthresh] = INAME;
      memory->grow(field2index,nfield+nthresh+1,"dump:field2index");
      memory->grow(argindex,nfield+nthresh+1,"dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
      argindex[nfield+nthresh] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if ((n < 0) || (tmp != 0))
        error->all(FLERR,"Could not find dump modify "
                   "custom atom integer property ID");

      field2index[nfield+nthresh] = add_custom(suffix,0);
      delete [] suffix;

    } else error->all(FLERR,"Invalid dump_modify thresh attribute");

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
      thresh_value[nthresh] = force->numeric(FLERR,arg[3]);
      thresh_last[nthresh] = -1;
    } else {
      thresh_fix = (FixStore **)
        memory->srealloc(thresh_fix,(nthreshlast+1)*sizeof(FixStore *),
                         "dump:thresh_fix");
      thresh_fixID = (char **)
        memory->srealloc(thresh_fixID,(nthreshlast+1)*sizeof(char *),
                         "dump:thresh_fixID");
      memory->grow(thresh_first,(nthreshlast+1),"dump:thresh_first");

      int n = strlen(id) + strlen("_DUMP_STORE") + 8;
      thresh_fixID[nthreshlast] = new char[n];
      strcpy(thresh_fixID[nthreshlast],id);
      sprintf(&thresh_fixID[nthreshlast][strlen(id)],"%d",nthreshlast);
      strcat(thresh_fixID[nthreshlast],"_DUMP_STORE");

      char **newarg = new char*[6];
      newarg[0] = thresh_fixID[nthreshlast];
      newarg[1] = group->names[igroup];
      newarg[2] = (char *) "STORE";
      newarg[3] = (char *) "peratom";
      newarg[4] = (char *) "1";
      newarg[5] = (char *) "1";
      modify->add_fix(6,newarg);
      thresh_fix[nthreshlast] = (FixStore *) modify->fix[modify->nfix-1];
      delete [] newarg;

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

bigint DumpCustom::memory_usage()
{
  bigint bytes = Dump::memory_usage();
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

  int index = field2index[n];

  if (flag_custom[index] == 0) { // integer
    int iwhich,tmp;
    iwhich = atom->find_custom(id_custom[index],tmp);

    int *ivector = atom->ivector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = ivector[clist[i]];
      n += size_one;
    }
  } else if (flag_custom[index] == 1) { // double
    int iwhich,tmp;
    iwhich = atom->find_custom(id_custom[index],tmp);

    double *dvector = atom->dvector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = dvector[clist[i]];
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
