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

#include "string.h"
#include "dump_xyz.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpXYZ::DumpXYZ(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump xyz command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump xyz filename");

  size_one = 5;
  sort_flag = 1;
  sortcol = 0;

  char *str = (char *) "%d %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_header(bigint n)
{
  if (me == 0) {
    fprintf(fp,BIGINT_FORMAT "\n",n);
    fprintf(fp,"Atoms\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpXYZ::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::pack(int *ids)
{
  int m,n;

  int *tag = atom->tag;
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

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_data(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
	    static_cast<int> (mybuf[m+1]),mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }
}
