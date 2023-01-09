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

#include "dump_xyz.h"

#include "atom.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

#define ONELINE 128
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpXYZ::DumpXYZ(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  typenames(nullptr)
{
  if (narg != 5) error->all(FLERR,"Illegal dump xyz command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump xyz filename");

  size_one = 5;

  buffer_allow = 1;
  buffer_flag = 1;
  sort_flag = 1;
  sortcol = 0;

  delete[] format_default;

  format_default = utils::strdup("%s %g %g %g");

  ntypes = atom->ntypes;
  typenames = nullptr;
}

/* ---------------------------------------------------------------------- */

DumpXYZ::~DumpXYZ()
{
  delete[] format_default;
  format_default = nullptr;

  if (typenames) {
    for (int i = 1; i <= ntypes; i++)
      delete [] typenames[i];
    delete [] typenames;
    typenames = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::init_style()
{
  // format = copy of default or user-specified line format

  delete [] format;

  if (format_line_user)
    format = utils::strdup(fmt::format("{}\n", format_line_user));
  else
    format = utils::strdup(fmt::format("{}\n", format_default));

  // initialize typenames array to be backward compatible by default
  // a 32-bit int can be maximally 10 digits plus sign

  if (typenames == nullptr) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[12];
      sprintf(typenames[itype],"%d",itype);
    }
  }

  // setup function ptr

  if (buffer_flag == 1) write_choice = &DumpXYZ::write_string;
  else write_choice = &DumpXYZ::write_lines;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpXYZ::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR, "Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++)
        delete [] typenames[i];

      delete [] typenames;
      typenames = nullptr;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = utils::strdup(arg[itype]);
    }

    return ntypes+1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_header(bigint n)
{
  if (me == 0) {
    auto header = fmt::format("{}\n Atoms. Timestep: {}", n, update->ntimestep);
    if (time_flag) header += fmt::format(" Time: {:.6f}", compute_time());
    header += "\n";
    fmt::print(fp, header);
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (ids) ids[n++] = tag[i];
    }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpXYZ::convert_string(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
                      typenames[static_cast<int> (mybuf[m+1])],
                      mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_string(int n, double *mybuf)
{
  if (mybuf)
    fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            typenames[static_cast<int> (mybuf[m+1])],
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }
}
