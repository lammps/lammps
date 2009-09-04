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
  if (narg != 5) error->all("Illegal dump xyz command");
  if (binary || multiproc) error->all("Invalid dump xyz filename");

  size_one = 5;

  char *str = (char *) "%d %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  // allocate global array for atom coords if group is all

  if (igroup == 0) {
    natoms = static_cast<int> (atom->natoms);
    if (natoms <= 0) error->all("Invalid natoms for dump xyz");
    if (atom->tag_consecutive() == 0)
      error->all("Atom IDs must be consecutive for dump xyz");
    types = (int *) memory->smalloc(natoms*sizeof(int),"dump:types");
    coords = (float *) memory->smalloc(3*natoms*sizeof(float),"dump:coords");
  }

  ntotal = 0;
}

/* ---------------------------------------------------------------------- */

DumpXYZ::~DumpXYZ()
{
  if (igroup == 0) {
    memory->sfree(types);
    memory->sfree(coords);
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::init()
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

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf and global coords array
------------------------------------------------------------------------- */

double DumpXYZ::memory_usage()
{
  double bytes = maxbuf * sizeof(double);
  if (igroup == 0) {
    bytes += natoms * sizeof(int);
    bytes += 3*natoms * sizeof(float);
  }
  return bytes;
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_header(int n)
{
  // all procs realloc types & coords if necessary

  if (igroup == 0 && n != natoms) {
    memory->sfree(types);
    memory->sfree(coords);
    if (atom->tag_consecutive() == 0)
      error->all("Atom IDs must be consecutive for dump xyz");
    natoms = n;
    types = (int *) memory->smalloc(natoms*sizeof(int),"dump:types");
    coords = (float *) memory->smalloc(3*natoms*sizeof(float),"dump:coords");
  }

  // only proc 0 writes header

  if (me == 0) {
    fprintf(fp,"%d\n",n);
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

int DumpXYZ::pack()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
    }

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_data(int n, double *buf)
{
  // for group = all, spread buf atom coords into global arrays
  // if last chunk of atoms in this snapshot, write global arrays to file

  if (igroup == 0) {
    int j,tag;

    int m = 0;
    for (int i = 0; i < n; i++) {
      tag = static_cast<int> (buf[m]) - 1;
      types[tag] = static_cast<int> (buf[m+1]);
      j = 3*tag;
      coords[j++] = buf[m+2];
      coords[j++] = buf[m+3];
      coords[j] = buf[m+4];
      m += size_one;
    }

    ntotal += n;
    if (ntotal == natoms) {
      write_frame();
      ntotal = 0;
    }

  // for group != all, write unsorted type,x,y,z to file

  } else {
    int m = 0;
    for (int i = 0; i < n; i++) {
      fprintf(fp,format,
	      static_cast<int> (buf[m+1]),buf[m+2],buf[m+3],buf[m+4]);
      m += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZ::write_frame()
{
  int m = 0;
  for (int i = 0; i < natoms; i++) {
    fprintf(fp,format,types[i],coords[m],coords[m+1],coords[m+2]);
    m += 3;
  }
}
