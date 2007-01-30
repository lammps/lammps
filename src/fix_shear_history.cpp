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
#include "stdio.h"
#include "fix_shear_history.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXTOUCH 15

/* ---------------------------------------------------------------------- */

FixShearHistory::FixShearHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  restart_peratom = 1;

  // perform initial allocation of atom-based arrays
  // register with atom class

  npartner = NULL;
  partner = NULL;
  shearpartner = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;
}

/* ---------------------------------------------------------------------- */

FixShearHistory::~FixShearHistory()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->sfree(npartner);
  memory->destroy_2d_int_array(partner);
  memory->destroy_3d_double_array(shearpartner);
}

/* ---------------------------------------------------------------------- */

int FixShearHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::init()
{
  if (atom->tag_enable == 0) 
    error->all("Pair style granular with history requires atoms have IDs");
}

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   so can be exchanged with atoms
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange()
{
  int i,j,k,m;

  // zero npartners for all current atoms

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) npartner[i] = 0;

  // copy shear info from neighbor list atoms to atom arrays
  // nlocal = nlocal_neighbor = nlocal when neighbor list last built,
  //   which might be pre-insert on this step

  int numneigh;
  int *neighs,*touch;
  double *firstshear,*shear;
  int *tag = atom->tag;
  nlocal = neighbor->nlocal_neighbor;

  for (i = 0; i < nlocal; i++) {
    neighs = neighbor->firstneigh[i];
    touch = neighbor->firsttouch[i];
    firstshear = neighbor->firstshear[i];
    numneigh = neighbor->numneigh[i];
    for (k = 0; k < numneigh; k++) {
      if (touch[k]) {
	shear = &firstshear[3*k];
	j = neighs[k];
	if (npartner[i] < MAXTOUCH) {
	  m = npartner[i];
	  partner[i][m] = tag[j];
	  shearpartner[i][m][0] = shear[0];
	  shearpartner[i][m][1] = shear[1];
	  shearpartner[i][m][2] = shear[2];
	}
	npartner[i]++;
	if (j < nlocal) {
	  if (npartner[j] < MAXTOUCH) {
	    m = npartner[j];
	    partner[j][m] = tag[i];
	    shearpartner[j][m][0] = -shear[0];
	    shearpartner[j][m][1] = -shear[1];
	    shearpartner[j][m][2] = -shear[2];
	  }
	  npartner[j]++;
	}
      }
    }
  }

  // test for too many touching neighbors

  int flag = 0;
  for (i = 0; i < nlocal; i++)
    if (npartner[i] >= MAXTOUCH) flag = 1;
  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all("Too many touching neighbors - boost MAXTOUCH");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

int FixShearHistory::memory_usage()
{
  int nmax = atom->nmax;
  int bytes = nmax * sizeof(int);
  bytes += nmax*MAXTOUCH * sizeof(int);
  bytes += nmax*MAXTOUCH*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::grow_arrays(int nmax)
{
  npartner = (int *) memory->srealloc(npartner,nmax*sizeof(int),
				      "shear_history:npartner");
  partner = memory->grow_2d_int_array(partner,nmax,MAXTOUCH,
				      "shear_history:partner");
  shearpartner = 
    memory->grow_3d_double_array(shearpartner,nmax,MAXTOUCH,3,
				 "shear_history:shearpartner");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::copy_arrays(int i, int j)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    shearpartner[j][m][0] = shearpartner[i][m][0];
    shearpartner[j][m][1] = shearpartner[i][m][1];
    shearpartner[j][m][2] = shearpartner[i][m][2];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    buf[m++] = shearpartner[i][n][0];
    buf[m++] = shearpartner[i][n][1];
    buf[m++] = shearpartner[i][n][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (buf[m++]);
    shearpartner[nlocal][n][0] = buf[m++];
    shearpartner[nlocal][n][1] = buf[m++];
    shearpartner[nlocal][n][2] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixShearHistory::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 4*npartner[i] + 2;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    buf[m++] = shearpartner[i][n][0];
    buf[m++] = shearpartner[i][n][1];
    buf[m++] = shearpartner[i][n][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixShearHistory::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    shearpartner[nlocal][n][0] = extra[nlocal][m++];
    shearpartner[nlocal][n][1] = extra[nlocal][m++];
    shearpartner[nlocal][n][2] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixShearHistory::maxsize_restart()
{
  return 4*MAXTOUCH + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixShearHistory::size_restart(int nlocal)
{
  return 4*npartner[nlocal] + 2;
}
