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

#include "stdlib.h"
#include "string.h"
#include "atom_vec_hybrid.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecHybrid::AtomVecHybrid(LAMMPS *lmp, int narg, char **arg) :
  AtomVec(lmp, narg, arg)
{
  int i,k;

  if (narg < 1) error->all("Illegal atom_style command");

  // create sub-styles

  nstyles = narg;
  styles = new AtomVec*[nstyles];
  keywords = new char*[nstyles];

  for (i = 0; i < narg; i++) {
    for (k = 0; k < i; k++)
      if (strcmp(arg[i],keywords[k]) == 0) 
	error->all("Atom style hybrid cannot use same atom style twice");
    if (strcmp(arg[i],"hybrid") == 0) 
      error->all("Atom style hybrid cannot have hybrid as an argument");
    styles[i] = atom->new_avec(arg[i],0,NULL);
    keywords[i] = new char[strlen(arg[i])+1];
    strcpy(keywords[i],arg[i]);
  }

  // hybrid settings are MAX or MIN of sub-style settings
  // hybrid sizes are minimial values plus extra values for each sub-style

  size_comm = 3;
  size_reverse = 3;
  size_border = 6;
  size_data_atom = 5;
  size_data_vel = 4;
  xcol_data = 3;

  for (k = 0; k < nstyles; k++) {
    molecular = MAX(molecular,styles[k]->molecular);
    bonds_allow = MAX(bonds_allow,styles[k]->bonds_allow);
    angles_allow = MAX(angles_allow,styles[k]->angles_allow);
    dihedrals_allow = MAX(dihedrals_allow,styles[k]->dihedrals_allow);
    impropers_allow = MAX(impropers_allow,styles[k]->impropers_allow);
    mass_type = MAX(mass_type,styles[k]->mass_type);
    shape_type = MAX(shape_type,styles[k]->shape_type);
    dipole_type = MAX(dipole_type,styles[k]->dipole_type);
    comm_x_only = MIN(comm_x_only,styles[k]->comm_x_only);
    comm_f_only = MIN(comm_f_only,styles[k]->comm_f_only);
    ghost_velocity = MAX(ghost_velocity,styles[k]->ghost_velocity);

    size_comm += styles[k]->size_comm - 3;
    size_reverse += styles[k]->size_reverse - 3;
    size_border += styles[k]->size_border - 6;
    size_data_atom += styles[k]->size_data_atom - 5;
    size_data_vel += styles[k]->size_data_vel - 4;
  }
}

/* ---------------------------------------------------------------------- */

AtomVecHybrid::~AtomVecHybrid()
{
  for (int k = 0; k < nstyles; k++) delete styles[k];
  delete [] styles;
  for (int k = 0; k < nstyles; k++) delete [] keywords[k];
  delete [] keywords;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void AtomVecHybrid::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;

  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;
  for (int k = 0; k < nstyles; k++) styles[k]->grow(nmax);
  atom->nextra_grow = tmp;

  // pointers for arrays used directly by hybrid style

  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::reset_special()
{
  for (int k = 0; k < nstyles; k++) styles[k]->reset_special();
}

/* ----------------------------------------------------------------------
   copy array values for all sub-styles
------------------------------------------------------------------------- */

void AtomVecHybrid::copy(int i, int j)
{
  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;
  for (int k = 0; k < nstyles; k++) styles[k]->copy(i,j);
  atom->nextra_grow = tmp;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j);
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_comm(int n, int *list, double *buf,
			     int pbc_flag, int *pbc)
{
  int i,j,k,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      for (k = 0; k < nstyles; k++)
	m += styles[k]->pack_comm_one(j,&buf[m]);
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      for (k = 0; k < nstyles; k++)
	m += styles[k]->pack_comm_one(j,&buf[m]);
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::unpack_comm(int n, int first, double *buf)
{
  int i,k,last;

  int m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    for (k = 0; k < nstyles; k++)
      m += styles[k]->unpack_comm_one(i,&buf[m]);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_reverse(int n, int first, double *buf)
{
  int i,k,last;

  int m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    for (k = 0; k < nstyles; k++)
      m += styles[k]->pack_reverse_one(i,&buf[m]);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    for (k = 0; k < nstyles; k++)
      m += styles[k]->unpack_reverse_one(j,&buf[m]);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_border(int n, int *list, double *buf,
			       int pbc_flag, int *pbc)
{
  int i,j,k,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      for (k = 0; k < nstyles; k++)
	m += styles[k]->pack_border_one(j,&buf[m]);
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      for (k = 0; k < nstyles; k++)
	m += styles[k]->pack_border_one(j,&buf[m]);
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::unpack_border(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    for (k = 0; k < nstyles; k++)
      m += styles[k]->unpack_border_one(i,&buf[m]);
  }
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   pack each sub-style one after the other
------------------------------------------------------------------------- */

int AtomVecHybrid::pack_exchange(int i, double *buf)
{
  int k;

  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;

  int m = 0;
  for (k = 0; k < nstyles; k++) 
    m += styles[k]->pack_exchange(i,&buf[m]);

  atom->nextra_grow = tmp;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for single atom received from another proc
   unpack each sub-style one after the other
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

int AtomVecHybrid::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;

  int m = 0;
  for (int k = 0; k < nstyles; k++) {
    m += styles[k]->unpack_exchange(&buf[m]);
    atom->nlocal--;
  }

  atom->nextra_grow = tmp;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      m += modify->fix[atom->extra_grow[iextra]]->
	unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecHybrid::size_restart()
{
  int tmp = atom->nextra_restart;
  atom->nextra_restart = 0;

  int n = 0;
  for (int k = 0; k < nstyles; k++)
    n += styles[k]->size_restart();

  atom->nextra_restart = tmp;

  int nlocal = atom->nlocal;
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++) 
      for (int i = 0; i < nlocal; i++)
	n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   pack each sub-style one after the other
------------------------------------------------------------------------- */

int AtomVecHybrid::pack_restart(int i, double *buf)
{
  int tmp = atom->nextra_restart;
  atom->nextra_restart = 0;

  int m = 0;
  for (int k = 0; k < nstyles; k++)
    m += styles[k]->pack_restart(i,&buf[m]);

  atom->nextra_restart = tmp;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++) 
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
   unpack each sub-style one after the other
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

int AtomVecHybrid::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      atom->extra = memory->grow_2d_double_array(atom->extra,nmax,
						 atom->nextra_store,
						 "atom:extra");
  }

  int tmp = atom->nextra_store;
  atom->nextra_store = 0;

  int m = 0;
  for (int k = 0; k < nstyles; k++) {
    m += styles[k]->unpack_restart(&buf[m]);
    atom->nlocal--;
  }
  atom->nextra_store = tmp;

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   create each sub-style one after the other
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

void AtomVecHybrid::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  for (int k = 0; k < nstyles; k++) {
    styles[k]->create_atom(itype,coord);
    atom->nlocal--;
  }
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

void AtomVecHybrid::data_atom(double *coord, int imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one("Invalid atom ID in Atoms section of data file");

  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one("Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;
  mask[nlocal] = 1;

  // each sub-style parses sub-style specific values

  int m = 5;
  for (int k = 0; k < nstyles; k++) 
    m += styles[k]->data_atom_hybrid(nlocal,&values[m]);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecHybrid::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);

  // each sub-style parses sub-style specific values

  int n = 3;
  for (int k = 0; k < nstyles; k++) 
    n += styles[k]->data_vel_hybrid(m,&values[n]);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double AtomVecHybrid::memory_usage()
{
  double bytes = 0.0;
  for (int k = 0; k < nstyles; k++) bytes += styles[k]->memory_usage();
  return bytes;
}
