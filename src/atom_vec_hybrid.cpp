
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
#include "atom_vec_hybrid.h"
#include "atom.h"
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
  int i,m;

  if (narg < 1) error->all("Illegal atom_style command");

  // create sub-styles

  nstyles = narg;
  styles = new AtomVec*[nstyles];
  keywords = new char*[nstyles];

  for (i = 0; i < narg; i++) {
    for (m = 0; m < i; m++)
      if (strcmp(arg[i],keywords[m]) == 0) 
	error->all("Atom style hybrid cannot use same atom style twice");
    if (strcmp(arg[i],"hybrid") == 0) 
      error->all("Atom style hybrid cannot have hybrid as an argument");
    styles[i] = atom->new_avec(arg[i],0,NULL);
    keywords[i] = new char[strlen(arg[i])+1];
    strcpy(keywords[i],arg[i]);
  }

  // hybrid settings are MAX or MIN of sub-style settings
  // size_border has +1 for hybrid[] value that is also communicated

  for (m = 0; m < nstyles; m++) {
    molecular = MAX(molecular,styles[m]->molecular);
    bonds_allow = MAX(bonds_allow,styles[m]->bonds_allow);
    angles_allow = MAX(angles_allow,styles[m]->angles_allow);
    dihedrals_allow = MAX(dihedrals_allow,styles[m]->dihedrals_allow);
    impropers_allow = MAX(impropers_allow,styles[m]->impropers_allow);
    mass_type = MAX(mass_type,styles[m]->mass_type);
    dipole_type = MAX(dipole_type,styles[m]->dipole_type);
    comm_x_only = MIN(comm_x_only,styles[m]->comm_x_only);
    comm_f_only = MIN(comm_f_only,styles[m]->comm_f_only);
    size_comm = MAX(size_comm,styles[m]->size_comm);
    size_reverse = MAX(size_reverse,styles[m]->size_reverse);
    size_border = MAX(size_border,styles[m]->size_border) + 1;
  }
}

/* ---------------------------------------------------------------------- */

AtomVecHybrid::~AtomVecHybrid()
{
  for (int m = 0; m < nstyles; m++) delete styles[m];
  delete [] styles;
  for (int m = 0; m < nstyles; m++) delete [] keywords[m];
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
  for (int m = 0; m < nstyles; m++) {
    atom->nextra_grow = 0;
    styles[m]->grow(nmax);
    atom->nextra_grow = tmp;
  }

  // pointers for arrays used directly by hybrid style

  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;

  hybrid = atom->hybrid = (int *) 
    memory->srealloc(atom->hybrid,nmax*sizeof(int),"atom:hybrid");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::reset_ptrs()
{
  for (int m = 0; m < nstyles; m++) styles[m]->reset_ptrs();

  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;

  hybrid = atom->hybrid;
}

/* ----------------------------------------------------------------------
   copy array values based on hybrid style of atom i
   zero auxiliary arrays for all other styles before copy
------------------------------------------------------------------------- */

void AtomVecHybrid::copy(int i, int j)
{
  int ihybrid = hybrid[i];
  for (int m = 0; m < nstyles; m++)
    if (m != ihybrid) styles[m]->zero_owned(j);
  styles[ihybrid]->copy(i,j);
  hybrid[j] = ihybrid;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_comm(int n, int *list, double *buf,
			     int pbc_flag, double *pbc_dist)
{
  int i,j,m;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      m += styles[hybrid[j]]->pack_comm_one(j,&buf[m]);
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + pbc_dist[0];
      buf[m++] = x[j][1] + pbc_dist[1];
      buf[m++] = x[j][2] + pbc_dist[2];
      m += styles[hybrid[j]]->pack_comm_one(j,&buf[m]);
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    m += styles[hybrid[i]]->unpack_comm_one(i,&buf[m]);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    m += styles[hybrid[i]]->pack_reverse_one(i,&buf[m]);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    m += styles[hybrid[j]]->unpack_reverse_one(j,&buf[m]);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecHybrid::pack_border(int n, int *list, double *buf,
			       int pbc_flag, double *pbc_dist)
{
  int i,j,m;

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
      buf[m++] = hybrid[j];
      m += styles[hybrid[j]]->pack_border_one(j,&buf[m]);
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + pbc_dist[0];
      buf[m++] = x[j][1] + pbc_dist[1];
      buf[m++] = x[j][2] + pbc_dist[2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = hybrid[j];
      m += styles[hybrid[j]]->pack_border_one(j,&buf[m]);
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack border ghost atom data
   zero auxiliary ghost arrays for all styles before unpack
   grow() is called in zero_ghost() and here (in case zero_ghost is no-op)
------------------------------------------------------------------------- */

void AtomVecHybrid::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  for (m = 0; m < nstyles; m++) styles[m]->zero_ghost(n,first);

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
    hybrid[i] = static_cast<int> (buf[m++]);
    m += styles[hybrid[i]]->unpack_border_one(i,&buf[m]);
  }
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   sub-style does packing
   append hybrid[i] and increment count stored in buf[0]
------------------------------------------------------------------------- */

int AtomVecHybrid::pack_exchange(int i, double *buf)
{
  int m = styles[hybrid[i]]->pack_exchange(i,buf);
  buf[m++] = hybrid[i];
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for single atom received from another proc
   sub-style does unpacking
   grow() must occur here so arrays for all sub-styles are grown
   extract hybrid[nlocal] from end of buf
------------------------------------------------------------------------- */

int AtomVecHybrid::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = static_cast<int> (buf[0]);
  hybrid[nlocal] = static_cast<int> (buf[m-1]);
  int tmp = styles[hybrid[nlocal]]->unpack_exchange(buf);
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecHybrid::size_restart()
{
  int i;
  int nlocal = atom->nlocal;

  int n = 0;
  for (i = 0; i < nlocal; i++)
    n += styles[hybrid[i]]->size_restart_one(i);

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++) 
      for (i = 0; i < nlocal; i++)
	n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   sub-style does packing
   append hybrid[i] and increment count stored in buf[0]
------------------------------------------------------------------------- */

int AtomVecHybrid::pack_restart(int i, double *buf)
{
  int m = styles[hybrid[i]]->pack_restart(i,buf);
  buf[m++] = hybrid[i];
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
   sub-style does unpacking
   grow() must occur here so arrays for all sub-styles are grown
   zero auxiliary arrays for all other styles before unpack
   extract hybrid[nlocal] from end of buf
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

  int m = static_cast<int> (buf[0]);
  int ihybrid = static_cast<int> (buf[m-1]);
  for (int m = 0; m < nstyles; m++)
    if (m != ihybrid) styles[m]->zero_owned(nlocal);
  hybrid[nlocal] = ihybrid;

  // size of extra unpack in sub-style includes end-of-buf entry of hybrid
  // avoid this by resetting buf[0] to one less

  buf[0] = m-1;
  int tmp = styles[ihybrid]->unpack_restart(buf);
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord for ihybrid style
   sub-style does create
   grow() must occur here so arrays for all sub-styles are grown
   zero auxiliary arrays for all other styles before create
------------------------------------------------------------------------- */

void AtomVecHybrid::create_atom(int itype, double *coord, int ihybrid)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  for (int m = 0; m < nstyles; m++)
    if (m != ihybrid) styles[m]->zero_owned(nlocal);
  hybrid[nlocal] = ihybrid;
  styles[ihybrid]->create_atom(itype,coord,0);
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   zero auxiliary arrays for all other styles before unpack
   sub-style will increment nlocal
------------------------------------------------------------------------- */

void AtomVecHybrid::data_atom(double *coord, int imagetmp, char **values,
			      int ihybrid)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  for (int m = 0; m < nstyles; m++)
    if (m != ihybrid) styles[m]->zero_owned(nlocal);
  hybrid[nlocal] = ihybrid;
  styles[ihybrid]->data_atom(coord,imagetmp,values,0);
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecHybrid::data_vel(int m, char *line, int ihybrid)
{
  styles[ihybrid]->data_vel(m,line,0);
}

/* ----------------------------------------------------------------------
   set data file parameters for ihybrid sub-style
------------------------------------------------------------------------- */

void AtomVecHybrid::data_params(int ihybrid)
{
  size_data_atom = styles[ihybrid]->size_data_atom;
  size_data_vel = styles[ihybrid]->size_data_vel;
  xcol_data = styles[ihybrid]->xcol_data;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

int AtomVecHybrid::memory_usage()
{
  int bytes = 0;
  for (int m = 0; m < nstyles; m++)
    bytes += styles[m]->memory_usage();
  if (atom->memcheck("hybrid")) bytes += nmax * sizeof(int);
  return bytes;
}
