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

/* ----------------------------------------------------------------------
   Contributing authors:
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
----------------------------------------------------------------------- */

#include "compute_rheo_property_atom.h"

#include "atom.h"
#include "atom_vec.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_surface.h"
#include "compute_rheo_vshift.h"
#include "error.h"
#include "fix_rheo.h"
#include "fix_rheo_thermal.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace RHEO_NS;

/* ---------------------------------------------------------------------- */

ComputeRHEOPropertyAtom::ComputeRHEOPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), fix_rheo(nullptr), fix_thermal(nullptr),  compute_interface(nullptr),
  compute_kernel(nullptr), compute_surface(nullptr), compute_vshift(nullptr),
  index(nullptr), pack_choice(nullptr)
{
  if (narg < 4)  utils::missing_cmd_args(FLERR, "compute property/atom", error);

  peratom_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  thermal_flag = interface_flag = surface_flag = shift_flag = 0;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];
  index = new int[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"phase") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_phase;
    } else if (strcmp(arg[iarg],"chi") == 0) {
      interface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_chi;
    } else if (strcmp(arg[iarg],"surface") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface;
    } else if (strcmp(arg[iarg],"surface/r") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_r;
    } else if (strcmp(arg[iarg],"surface/divr") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_divr;
    } else if (strcmp(arg[iarg],"surface/nx") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_nx;
    } else if (strcmp(arg[iarg],"surface/ny") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_ny;
    } else if (strcmp(arg[iarg],"surface/nz") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_nz;
    } else if (strcmp(arg[iarg],"coordination") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_coordination;
    } else if (strcmp(arg[iarg],"cv") == 0) {
      thermal_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_cv;
    } else if (strcmp(arg[iarg],"shift/vx") == 0) {
      shift_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_shift_vx;
    } else if (strcmp(arg[iarg],"shift/vy") == 0) {
      shift_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_shift_vy;
    } else if (strcmp(arg[iarg],"shift/vz") == 0) {
      shift_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_shift_vz;
    } else {
      index[i] = atom->avec->property_atom(arg[iarg]);
      if (index[i] < 0)
        error->all(FLERR,
                   "Invalid keyword {} for atom style {} in compute rheo/property/atom command ",
                   atom->get_style(), arg[iarg]);
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_atom_style;

      if (strcmp(arg[iarg],"temperature") == 0) thermal_flag = 1;
      if (strcmp(arg[iarg],"heatflow") == 0) thermal_flag = 1;
      if (strcmp(arg[iarg],"conductivity") == 0) thermal_flag = 1;
    }
  }

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOPropertyAtom::~ComputeRHEOPropertyAtom()
{
  delete[] pack_choice;
  delete[] index;
  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/viscosity");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  if (interface_flag && !(fix_rheo->interface_flag))
    error->all(FLERR, "Cannot request interfacial property without corresponding option in fix rheo");
  if (surface_flag && !(fix_rheo->surface_flag))
    error->all(FLERR, "Cannot request surface property without corresponding option in fix rheo");
  if (shift_flag && !(fix_rheo->shift_flag))
    error->all(FLERR, "Cannot request velocity shifting property without corresponding option in fix rheo");
  if (thermal_flag && !(fix_rheo->thermal_flag))
    error->all(FLERR, "Cannot request thermal property without fix rheo/thermal");

  compute_interface = fix_rheo->compute_interface;
  compute_kernel = fix_rheo->compute_kernel;
  compute_surface = fix_rheo->compute_surface;
  compute_vshift = fix_rheo->compute_vshift;

  if (thermal_flag) {
    fixes = modify->get_fix_by_style("rheo/thermal");
    fix_thermal = dynamic_cast<FixRHEOThermal *>(fixes[0]);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow vector or array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    if (nvalues == 1) {
      memory->destroy(vector_atom);
      memory->create(vector_atom,nmax,"rheo/property/atom:vector");
    } else {
      memory->destroy(array_atom);
      memory->create(array_atom,nmax,nvalues,"rheo/property/atom:array");
    }
  }

  // fill vector or array with per-atom values

  if (nvalues == 1) {
    buf = vector_atom;
    (this->*pack_choice[0])(0);
  } else {
    if (nmax) buf = &array_atom[0][0];
    else buf = nullptr;
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOPropertyAtom::memory_usage()
{
  double bytes = (double)nmax * nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute rheo/property/atom can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_phase(int n)
{
  int *status = atom->status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = (status[i] & PHASECHECK);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_chi(int n)
{
  double *chi = compute_interface->chi;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = chi[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface(int n)
{
  int *status = atom->status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = (status[i] & SURFACECHECK);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface_r(int n)
{
  double *rsurface = compute_surface->rsurface;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = rsurface[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface_divr(int n)
{
  double *divr = compute_surface->divr;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = divr[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface_nx(int n)
{
  double **nsurface = compute_surface->nsurface;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = nsurface[i][0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface_ny(int n)
{
  double **nsurface = compute_surface->nsurface;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = nsurface[i][1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}


/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface_nz(int n)
{
  double **nsurface = compute_surface->nsurface;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = nsurface[i][2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_coordination(int n)
{
  int *coordination = compute_kernel->coordination;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = coordination[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_cv(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = fix_thermal->calc_cv(i);
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_shift_vx(int n)
{
  double **vshift = compute_vshift->vshift;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = vshift[i][0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_shift_vy(int n)
{
  double **vshift = compute_vshift->vshift;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = vshift[i][1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_shift_vz(int n)
{
  double **vshift = compute_vshift->vshift;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = vshift[i][2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_atom_style(int n)
{
  atom->avec->pack_property_atom(index[n],&buf[n],nvalues,groupbit);
}
