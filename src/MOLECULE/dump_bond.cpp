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
#include "dump_bond.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpBond::DumpBond(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal dump bond command");
  if (atom->molecular == 0)
    error->all("Cannot use dump bond with non-molecular system");

  size_one = 3;

  char *str = (char *) "%d %d %d %d";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);
}

/* ---------------------------------------------------------------------- */

void DumpBond::init()
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

void DumpBond::write_header(int ndump)
{
  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",update->ntimestep);
    fprintf(fp,"ITEM: NUMBER OF BONDS\n");
    fprintf(fp,"%d\n",ndump);
    fprintf(fp,"ITEM: BONDS\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpBond::count()
{
  index = 0;

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int i,j,k;

  int m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (j = 0; j < num_bond[i]; j++) {
      k = atom->map(bond_atom[i][j]);
      if (k >= 0 && !(mask[k] & groupbit)) continue;
      if (bond_type[i][j] == 0) continue;
      m++;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int DumpBond::pack()
{
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int i,j,k,type,iatom;

  int m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (j = 0; j < num_bond[i]; j++) {
      iatom = bond_atom[i][j];
      k = atom->map(iatom);
      if (k >= 0 && !(mask[k] & groupbit)) continue;
      type = bond_type[i][j];
      if (type == 0) continue;
      buf[m++] = type;
      buf[m++] = tag[i];
      buf[m++] = iatom;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void DumpBond::write_data(int n, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    index++;
    fprintf(fp,format,
	    index, static_cast<int> (buf[m]),
	    static_cast<int> (buf[m+1]), static_cast<int> (buf[m+2]));
    m += size_one;
  }
}
