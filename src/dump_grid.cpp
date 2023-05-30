// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "dump_grid.h"

#include "arg_info.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "grid2d.h"
#include "grid3d.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

// customize by adding keyword
// also customize compute_atom_property.cpp

enum {COMPUTE,FIX};

#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpGrid::DumpGrid(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg), idregion(nullptr), earg(nullptr), vtype(nullptr),
  vformat(nullptr), columns(nullptr), columns_default(nullptr),
  field2index(nullptr), field2grid(nullptr), field2data(nullptr),
  argindex(nullptr), id_compute(nullptr), id_fix(nullptr), pack_choice(nullptr)
{
  if (narg == 5) error->all(FLERR,"No dump grid arguments specified");

  clearstep = 1;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal dump grid nevery value: {}", nevery);

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
  field2source = new int[nfield];
  field2index = new int[nfield];
  field2grid = new int[nfield];
  field2data = new int[nfield];
  argindex = new int[nfield];

  buffer_allow = 1;
  buffer_flag = 1;

  dimension = domain->dimension;

  // for 2d, set nzgrid = 1 for dump grid and grid/vtk files

  if (dimension == 2) nzgrid = 1;

  // computes and fixes which the dump accesses

  ncompute = 0;
  nfix = 0;

  // process attributes
  // ioptional = start of additional optional args in expanded args

  ioptional = parse_fields(nfield,earg);

  if (ioptional < nfield &&
      strcmp(style,"image") != 0 && strcmp(style,"movie") != 0)
    error->all(FLERR,"Invalid attribute {} in dump {} command", earg[ioptional],style);

  // noptional = # of optional args
  // reset nfield to subtract off optional args
  // reset ioptional to what it would be in original arg list
  // only dump image and dump movie styles process optional args,
  //   they do not use expanded earg list

  int noptional = nfield - ioptional;
  nfield -= noptional;
  size_one = nfield;
  ioptional = narg - noptional;

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
  cols.resize(cols.size()-1);
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

DumpGrid::~DumpGrid()
{
  // if wildcard expansion occurred, free earg memory from expand_args()
  // could not do in constructor, b/c some derived classes process earg

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  delete[] pack_choice;
  delete[] vtype;
  delete[] field2source;
  delete[] field2index;
  delete[] field2grid;
  delete[] field2data;
  delete[] argindex;

  delete[] idregion;

  for (int i = 0; i < ncompute; i++) delete[] id_compute[i];
  memory->sfree(id_compute);

  for (int i = 0; i < nfix; i++) delete[] id_fix[i];
  memory->sfree(id_fix);

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

void DumpGrid::init_style()
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
    error->all(FLERR,"Dump_modify format line is too short");

  int i=0;
  for (const auto &word : words) {
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
    header_choice = &DumpGrid::header_binary;
  else if (binary && domain->triclinic == 1)
    header_choice = &DumpGrid::header_binary_triclinic;
  else if (!binary && domain->triclinic == 0)
    header_choice = &DumpGrid::header_item;
  else if (!binary && domain->triclinic == 1)
    header_choice = &DumpGrid::header_item_triclinic;

  if (binary) write_choice = &DumpGrid::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpGrid::write_string;
  else write_choice = &DumpGrid::write_lines;

  // find current ptr for each compute and fix

  for (i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i])
      error->all(FLERR,"Could not find dump grid compute ID {}",id_compute[i]);
  }

  for (i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i])
      error->all(FLERR,"Could not find dump grid fix ID {}", id_fix[i]);
  }

  // check that grid sizes for all fields are the same

  Grid2d *grid2d = nullptr;
  Grid3d *grid3d = nullptr;
  int nxtmp,nytmp,nztmp;
  for (int i = 0; i < nfield; i++) {
     if (dimension == 2) {
       if (field2source[i] == COMPUTE)
         grid2d = (Grid2d *) compute[field2index[i]]->get_grid_by_index(field2grid[i]);
       else
        grid2d = (Grid2d *) fix[field2index[i]]->get_grid_by_index(field2grid[i]);
      if (i == 0) grid2d->get_size(nxgrid,nygrid);
      else {
        grid2d->get_size(nxtmp,nytmp);
        if ((nxtmp != nxgrid) || (nytmp != nygrid))
          error->all(FLERR,"Dump grid field grid sizes do not match");
      }

    } else {
      if (field2source[i] == COMPUTE)
        grid3d = (Grid3d *) compute[field2index[i]]->get_grid_by_index(field2grid[i]);
      else
        grid3d = (Grid3d *) fix[field2index[i]]->get_grid_by_index(field2grid[i]);
      if (i == 0) grid3d->get_size(nxgrid,nygrid,nzgrid);
      else {
        grid3d->get_size(nxtmp,nytmp,nztmp);
        if ((nxtmp != nxgrid) || (nytmp != nygrid) || (nztmp != nzgrid))
          error->all(FLERR,"Dump grid field grid sizes do not match");
      }
    }
  }

  // check validity of region

  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR,"Region {} for dump grid does not exist", idregion);

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_unit_style_binary()
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

void DumpGrid::header_columns_binary()
{
  int len = strlen(columns);
  fwrite(&len, sizeof(int), 1, fp);
  fwrite(columns, sizeof(char), len, fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_time_binary()
{
  char flag = time_flag ? 1 : 0;
  fwrite(&flag, sizeof(char), 1, fp);

  if (time_flag) {
    double t = compute_time();
    fwrite(&t, sizeof(double), 1, fp);
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_format_binary()
{
  format_magic_string_binary();
  format_endian_binary();
  format_revision_binary();
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_binary(bigint ndump)
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
  fwrite(&domain->dimension,sizeof(int),1,fp);
  fwrite(&nxgrid,sizeof(int),1,fp);
  fwrite(&nygrid,sizeof(int),1,fp);
  fwrite(&nzgrid,sizeof(int),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);

  header_unit_style_binary();
  header_time_binary();
  header_columns_binary();

  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_binary_triclinic(bigint ndump)
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
  fwrite(&domain->dimension,sizeof(int),1,fp);
  fwrite(&nxgrid,sizeof(int),1,fp);
  fwrite(&nygrid,sizeof(int),1,fp);
  fwrite(&nzgrid,sizeof(int),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);

  header_unit_style_binary();
  header_time_binary();
  header_columns_binary();

  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_item(bigint /*ndump*/)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fmt::print(fp,"ITEM: UNITS\n{}\n",update->unit_style);
  }
  if (time_flag) fmt::print(fp,"ITEM: TIME\n{:.16}\n",compute_time());

  fmt::print(fp,"ITEM: TIMESTEP\n{}\n",update->ntimestep);
  fmt::print(fp,"ITEM: BOX BOUNDS {}\n"
             "{:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e}\n",
             boundstr,boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi);
  fmt::print(fp,"ITEM: DIMENSION\n{}\n",domain->dimension);
  fmt::print(fp,"ITEM: GRID SIZE nx ny nz\n{} {} {}\n",nxgrid,nygrid,nzgrid);
  fmt::print(fp,"ITEM: GRID CELLS {}\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::header_item_triclinic(bigint /*ndump*/)
{
  if (unit_flag && !unit_count) {
    ++unit_count;
    fmt::print(fp,"ITEM: UNITS\n{}\n",update->unit_style);
  }
  if (time_flag) fmt::print(fp,"ITEM: TIME\n{:.16}\n",compute_time());

  fmt::print(fp,"ITEM: TIMESTEP\n{}\n",update->ntimestep);
  fmt::print(fp,"ITEM: BOX BOUNDS xy xz yz {}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e}\n"
             "{:>1.16e} {:>1.16e} {:>1.16e}\n",
             boundstr,boxxlo,boxxhi,boxxy,boxylo,boxyhi,boxxz,boxzlo,boxzhi,boxyz);
  fmt::print(fp,"ITEM: DIMENSION\n{}\n",domain->dimension);
  fmt::print(fp,"ITEM: GRID SIZE nx ny nz\n{} {} {}\n",nxgrid,nygrid,nzgrid);
  fmt::print(fp,"ITEM: GRID CELLS {}\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::format_magic_string_binary()
{
  // use negative ntimestep as marker for new format
  bigint fmtlen = strlen(MAGIC_STRING);
  bigint marker = -fmtlen;
  fwrite(&marker, sizeof(bigint), 1, fp);
  fwrite(MAGIC_STRING, sizeof(char), fmtlen, fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::format_endian_binary()
{
  int endian = ENDIAN;
  fwrite(&endian, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::format_revision_binary()
{
  int revision = FORMAT_REVISION;
  fwrite(&revision, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

int DumpGrid::count()
{
  int i;

  // set current size for portion of grid on each proc
  // may change between dump snapshots due to load balancing

  Grid2d *grid2d = nullptr;
  Grid3d *grid3d = nullptr;

  if (dimension == 2) {
    if (field2source[0] == COMPUTE)
      grid2d = (Grid2d *) compute[field2index[0]]->get_grid_by_index(field2grid[0]);
    else if (field2source[0] == FIX)
      grid2d = (Grid2d *) fix[field2index[0]]->get_grid_by_index(field2grid[0]);
    else error->all(FLERR, "Unsupported grid data source type {}", field2source[0]);
    grid2d->get_bounds_owned(nxlo_in,nxhi_in,nylo_in,nyhi_in);
  } else {
    if (field2source[0] == COMPUTE)
      grid3d = (Grid3d *) compute[field2index[0]]->get_grid_by_index(field2grid[0]);
    else if (field2source[0] == FIX)
      grid3d = (Grid3d *) fix[field2index[0]]->get_grid_by_index(field2grid[0]);
    else error->all(FLERR, "Unsupported grid data source type {}", field2source[0]);
    grid3d->get_bounds_owned(nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in);
  }

  // invoke Computes for per-grid quantities
  // cannot invoke before first run, otherwise invoke if necessary

  if (ncompute) {
    if (update->first_update == 0)
      error->all(FLERR,"Dump compute cannot be invoked before first run");
    for (i = 0; i < ncompute; i++) {
      if (!(compute[i]->invoked_flag & Compute::INVOKED_PERGRID)) {
        compute[i]->compute_pergrid();
        compute[i]->invoked_flag |= Compute::INVOKED_PERGRID;
      }
    }
  }

  // return count of grid points I own

  int ngrid;
  if (dimension == 2)
    ngrid = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1);
  else
    ngrid = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) * (nzhi_in-nzlo_in+1);

  return ngrid;
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack(tagint *ids)
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);

  // ids = list of my grid IDs

  if (ids) {
    int m = 0;
    if (dimension == 2) {
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            ids[m++] = iy*nxgrid + ix + 1;
    } else if (dimension == 3) {
        for (int iz = nzlo_in; iz <= nzhi_in; iz++)
          for (int iy = nylo_in; iy <= nyhi_in; iy++)
            for (int ix = nxlo_in; ix <= nxhi_in; ix++)
              ids[m++] = iz*nygrid*nxgrid + iy*nxgrid + ix + 1;
    }
  }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpGrid::convert_string(int n, double *mybuf)
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
      else if (vtype[j] == Dump::BIGINT)
        offset += sprintf(&sbuf[offset],vformat[j], static_cast<bigint> (mybuf[m]));
      m++;
    }
    offset += sprintf(&sbuf[offset],"\n");
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_string(int n, double *mybuf)
{
  if (mybuf)
    fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpGrid::write_lines(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < nfield; j++) {
      if (vtype[j] == Dump::INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
      else if (vtype[j] == Dump::BIGINT) fprintf(fp,vformat[j],static_cast<bigint> (mybuf[m]));
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpGrid::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  for (int iarg = 0; iarg < narg; iarg++) {

    char *id;
    int igrid,idata,index;
    int iflag =
      utils::check_grid_reference((char *) "Dump grid",
                                  arg[iarg],nevery,id,igrid,idata,index,lmp);

    // arg is not a valid Grid reference
    // assume it's an additional dump grid option and return

    if (iflag < 0) return iarg;

    // grid reference is to a compute or fix

    if (iflag == ArgInfo::COMPUTE) {
      auto icompute = lmp->modify->get_compute_by_id(id);
      field2index[iarg] = add_compute(id,icompute);
      field2source[iarg] = COMPUTE;
    } else if (iflag == ArgInfo::FIX) {
      auto ifix = modify->get_fix_by_id(id);
      field2index[iarg] = add_fix(id,ifix);
      field2source[iarg] = FIX;
    }

    delete [] id;
    argindex[iarg] = index;
    vtype[iarg] = Dump::DOUBLE;
    field2grid[iarg] = igrid;
    field2data[iarg] = idata;

    if (dimension == 2) pack_choice[iarg] = &DumpGrid::pack_grid2d;
    else pack_choice[iarg] = &DumpGrid::pack_grid3d;
  }

  return narg;
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpGrid::add_compute(const std::string &id, Compute *cptr)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (id == id_compute[icompute]) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **) memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  id_compute[ncompute] = utils::strdup(id);
  compute.push_back(cptr);

  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpGrid::add_fix(const std::string &id, Fix *fptr)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (id == id_fix[ifix]) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **) memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  id_fix[nfix] = utils::strdup(id);
  fix.push_back(fptr);

  nfix++;
  return nfix-1;
}

/* ---------------------------------------------------------------------- */

int DumpGrid::modify_param(int narg, char **arg)
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

  if (strcmp(arg[0],"format") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");

    if (strcmp(arg[1],"none") == 0) {
      // just clear format_column_user allocated by this dump child class
      for (int i = 0; i < nfield; i++) {
        delete[] format_column_user[i];
        format_column_user[i] = nullptr;
      }
      return 2;
    }

    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");

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
        error->all(FLERR,"Illegal dump_modify command");
      delete[] format_column_user[i];
      format_column_user[i] = utils::strdup(arg[2]);
    }
    return 3;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

double DumpGrid::memory_usage()
{
  double bytes = Dump::memory_usage();
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of 2d and 3d grid data
   from either compute or fix
------------------------------------------------------------------------- */

void DumpGrid::pack_grid2d(int n)
{
  int index = argindex[n];

  if (index == 0) {
    double **vec2d;
    if (field2source[n] == COMPUTE)
      vec2d = (double **) compute[field2index[n]]->get_griddata_by_index(field2data[n]);
    else if (field2source[n] == FIX)
      vec2d = (double **) fix[field2index[n]]->get_griddata_by_index(field2data[n]);
    else error->all(FLERR, "Unsupported grid data source type {}", field2source[n]);
    for (int iy = nylo_in; iy <= nyhi_in; iy++)
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf[n] = vec2d[iy][ix];
        n += size_one;
      }
  } else {
    double ***array2d;
    if (field2source[n] == COMPUTE)
      array2d = (double ***) compute[field2index[n]]->get_griddata_by_index(field2data[n]);
    else if (field2source[n] == FIX)
      array2d = (double ***) fix[field2index[n]]->get_griddata_by_index(field2data[n]);
    else error->all(FLERR, "Unsupported grid data source type {}", field2source[n]);
    index--;
    for (int iy = nylo_in; iy <= nyhi_in; iy++)
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf[n] = array2d[iy][ix][index];
        n += size_one;
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpGrid::pack_grid3d(int n)
{
  int index = argindex[n];

  if (index == 0) {
    double ***vec3d;
    if (field2source[n] == COMPUTE)
      vec3d = (double ***) compute[field2index[n]]->get_griddata_by_index(field2data[n]);
    else if (field2source[n] == FIX)
      vec3d = (double ***) fix[field2index[n]]->get_griddata_by_index(field2data[n]);
    else error->all(FLERR, "Unsupported grid data source type {}", field2source[n]);
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
          buf[n] = vec3d[iz][iy][ix];
          n += size_one;
        }
  } else {
    double ****array3d;
    if (field2source[n] == COMPUTE)
      array3d = (double ****) compute[field2index[n]]->get_griddata_by_index(field2data[n]);
    else if (field2source[n] == FIX)
      array3d = (double ****) fix[field2index[n]]->get_griddata_by_index(field2data[n]);
    else error->all(FLERR, "Unsupported grid data source type {}", field2source[n]);
    index--;
    for (int iz = nzlo_in; iz <= nzhi_in; iz++)
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
          buf[n] = array3d[iz][iy][ix][index];
          n += size_one;
        }
  }
}
