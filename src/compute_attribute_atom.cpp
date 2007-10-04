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
#include "compute_attribute_atom.h"
#include "atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{X,Y,Z,XU,YU,ZU,VX,VY,VZ,FX,FY,FZ,XYZ,V,F};

/* --------------------------------------------------------------------- */

ComputeAttributeAtom::ComputeAttributeAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute ke/atom command");

  peratom_flag = 1;
  size_peratom = 0;
  allocate = 1;

  if (strcmp(arg[3],"x") == 0) which = X;
  else if (strcmp(arg[3],"y") == 0) which = Y;
  else if (strcmp(arg[3],"z") == 0) which = Z;
  else if (strcmp(arg[3],"xu") == 0) which = XU;
  else if (strcmp(arg[3],"yu") == 0) which = YU;
  else if (strcmp(arg[3],"zu") == 0) which = ZU;
  else if (strcmp(arg[3],"vx") == 0) which = VX;
  else if (strcmp(arg[3],"vy") == 0) which = VY;
  else if (strcmp(arg[3],"vz") == 0) which = VZ;
  else if (strcmp(arg[3],"fx") == 0) which = FX;
  else if (strcmp(arg[3],"fy") == 0) which = FY;
  else if (strcmp(arg[3],"fz") == 0) which = FZ;

  else if (strcmp(arg[3],"xyz") == 0) {
    which = XYZ;
    size_peratom = 3;
    allocate = 0;
  } else if (strcmp(arg[3],"v") == 0) {
    which = V;
    size_peratom = 3;
    allocate = 0;
  } else if (strcmp(arg[3],"f") == 0) {
    which = F;
    size_peratom = 3;
    allocate = 0;
  } else error->all("Illegal compute attribute/atom command");

  nmax = 0;
  s_attribute = NULL;
  v_attribute = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeAttributeAtom::~ComputeAttributeAtom()
{
  if (allocate) {
    memory->sfree(s_attribute);
    memory->destroy_2d_double_array(v_attribute);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeAttributeAtom::compute_peratom()
{
  // grow attribute array if necessary

  if (allocate && atom->nlocal > nmax) {
    if (size_peratom == 0) {
      memory->sfree(s_attribute);
      nmax = atom->nmax;
      s_attribute = (double *) 
	memory->smalloc(nmax*sizeof(double),
			"compute/attribute/atom:s_attribute");
    } else {
      memory->destroy_2d_double_array(v_attribute);
      nmax = atom->nmax;
      v_attribute =  
	memory->create_2d_double_array(nmax,size_peratom,
				       "compute/attribute/atom:v_attribute");
    }
  }

  // fill attribute vector with appropriate atom value
  // or simply set pointer to exisitng atom vector

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  if (which == X) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = x[i][0];
      else s_attribute[i] = 0.0;
  } else if (which == Y) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = x[i][1];
      else s_attribute[i] = 0.0;
  } else if (which == Z) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = x[i][2];
      else s_attribute[i] = 0.0;

  } else if (which == XU) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	s_attribute[i] = x[i][0] + ((image[i] & 1023) - 512) * xprd;
      else s_attribute[i] = 0.0;
  } else if (which == YU) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	s_attribute[i] = x[i][1] + ((image[i] >> 10 & 1023) - 512) * yprd;
      else s_attribute[i] = 0.0;
  } else if (which == ZU) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	s_attribute[i] = x[i][2] + ((image[i] >> 20) - 512) * zprd;
      else s_attribute[i] = 0.0;

  } else if (which == VX) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = v[i][0];
      else s_attribute[i] = 0.0;
  } else if (which == VY) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = v[i][1];
      else s_attribute[i] = 0.0;
  } else if (which == VZ) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = v[i][2];
      else s_attribute[i] = 0.0;

  } else if (which == FX) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = f[i][0];
      else s_attribute[i] = 0.0;
  } else if (which == FY) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = f[i][1];
      else s_attribute[i] = 0.0;
  } else if (which == FZ) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) s_attribute[i] = f[i][2];
      else s_attribute[i] = 0.0;

  } else if (which == XYZ) v_attribute = x;
  else if (which == V) v_attribute = v;
  else if (which == F) v_attribute = f;

  // set appropriate compute ptr to local array

  if (size_peratom == 0) scalar_atom = s_attribute;
  else vector_atom = v_attribute;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeAttributeAtom::memory_usage()
{
  double bytes = 0.0;
  if (allocate) {
    if (size_peratom == 0) bytes = nmax * sizeof(double);
    else bytes = size_peratom * nmax * sizeof(double);
  }
  return bytes;
}
