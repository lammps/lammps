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

#include "fix_minimize.h"
#include "atom.h"
#include "domain.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMinimize::FixMinimize(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvector(0), peratom(NULL), vectors(NULL)
{
  // register callback to this fix from Atom class
  // don't perform initial allocation here, must wait until add_vector()

  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixMinimize::~FixMinimize()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored data

  memory->destroy(peratom);
  if (vectors) {
    for (int m = 0; m < nvector; m++) memory->destroy(vectors[m]);
    memory->sfree(vectors);
  }
}

/* ---------------------------------------------------------------------- */

int FixMinimize::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with N elements per atom
------------------------------------------------------------------------- */

void FixMinimize::add_vector(int n)
{
  memory->grow(peratom,nvector+1,"minimize:peratom");
  peratom[nvector] = n;

  vectors = (double **)
    memory->srealloc(vectors,(nvector+1)*sizeof(double *),"minimize:vectors");
  memory->create(vectors[nvector],atom->nmax*n,"minimize:vector");

  int ntotal = n*atom->nlocal;
  for (int i = 0; i < ntotal; i++) vectors[nvector][i] = 0.0;
  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

double *FixMinimize::request_vector(int m)
{
  return vectors[m];
}

/* ----------------------------------------------------------------------
   store box size at beginning of line search
------------------------------------------------------------------------- */

void FixMinimize::store_box()
{
  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];
  boxhi[0] = domain->boxhi[0];
  boxhi[1] = domain->boxhi[1];
  boxhi[2] = domain->boxhi[2];
}

/* ----------------------------------------------------------------------
   reset x0 for atoms that moved across PBC via reneighboring in line search
   x0 = 1st vector
   must do minimum_image using original box stored at beginning of line search
   swap & set_global_box() change to original box, then restore current box
------------------------------------------------------------------------- */

void FixMinimize::reset_coords()
{
  box_swap();
  domain->set_global_box();

  double **x = atom->x;
  double *x0 = vectors[0];
  int nlocal = atom->nlocal;
  double dx,dy,dz,dx0,dy0,dz0;

  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    dx = dx0 = x[i][0] - x0[n];
    dy = dy0 = x[i][1] - x0[n+1];
    dz = dz0 = x[i][2] - x0[n+2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dx0) x0[n] = x[i][0] - dx;
    if (dy != dy0) x0[n+1] = x[i][1] - dy;
    if (dz != dz0) x0[n+2] = x[i][2] - dz;
    n += 3;
  }

  box_swap();
  domain->set_global_box();
}

/* ----------------------------------------------------------------------
   swap current box size with stored box size
------------------------------------------------------------------------- */

void FixMinimize::box_swap()
{
  double tmp;

  tmp = boxlo[0];
  boxlo[0] = domain->boxlo[0];
  domain->boxlo[0] = tmp;
  tmp = boxlo[1];
  boxlo[1] = domain->boxlo[1];
  domain->boxlo[1] = tmp;
  tmp = boxlo[2];
  boxlo[2] = domain->boxlo[2];
  domain->boxlo[2] = tmp;

  tmp = boxhi[0];
  boxhi[0] = domain->boxhi[0];
  domain->boxhi[0] = tmp;
  tmp = boxhi[1];
  boxhi[1] = domain->boxhi[1];
  domain->boxhi[1] = tmp;
  tmp = boxhi[2];
  boxhi[2] = domain->boxhi[2];
  domain->boxhi[2] = tmp;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixMinimize::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvector; m++)
    bytes += atom->nmax*peratom[m]*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixMinimize::grow_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
    memory->grow(vectors[m],peratom[m]*nmax,"minimize:vector");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixMinimize::copy_arrays(int i, int j, int /*delflag*/)
{
  int m,iper,nper,ni,nj;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    nj = nper*j;
    for (iper = 0; iper < nper; iper++) vectors[m][nj++] = vectors[m][ni++];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixMinimize::pack_exchange(int i, double *buf)
{
  int m,iper,nper,ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    for (iper = 0; iper < nper; iper++) buf[n++] = vectors[m][ni++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixMinimize::unpack_exchange(int nlocal, double *buf)
{
  int m,iper,nper,ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*nlocal;
    for (iper = 0; iper < nper; iper++) vectors[m][ni++] = buf[n++];
  }
  return n;
}
