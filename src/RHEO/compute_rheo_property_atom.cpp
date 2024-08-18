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
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "compute_rheo_property_atom.h"

#include "atom.h"
#include "atom_vec.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_surface.h"
#include "compute_rheo_vshift.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "fix_rheo_oxidation.h"
#include "fix_rheo_pressure.h"
#include "fix_rheo_thermal.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace RHEO_NS;

/* ---------------------------------------------------------------------- */

ComputeRHEOPropertyAtom::ComputeRHEOPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), avec_index(nullptr), col_index(nullptr), col_t_index(nullptr),
    buf(nullptr), pack_choice(nullptr), fix_rheo(nullptr), fix_pressure(nullptr),
    fix_thermal(nullptr), compute_interface(nullptr), compute_kernel(nullptr),
    compute_surface(nullptr), compute_vshift(nullptr), compute_grad(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "compute property/atom", error);

  peratom_flag = 1;
  int dim = domain->dimension;

  // Determine number of values
  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "surface/n/*") == 0) {
      nvalues += dim;
    } else if (strcmp(arg[iarg], "shift/v/*") == 0) {
      nvalues += dim;
    } else if (strcmp(arg[iarg], "grad/v/*") == 0) {
      nvalues += dim * dim;
    } else if (strcmp(arg[iarg], "stress/v/*") == 0) {
      nvalues += dim * dim;
    } else if (strcmp(arg[iarg], "stress/t/*") == 0) {
      nvalues += dim * dim;
    } else {
      nvalues += 1;
    }
  }

  if (nvalues == 1)
    size_peratom_cols = 0;
  else
    size_peratom_cols = nvalues;

  pressure_flag = thermal_flag = interface_flag = 0;
  surface_flag = shift_flag = shell_flag = 0;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];
  avec_index = new int[nvalues];
  col_index = new int[nvalues];
  col_t_index = new int[nvalues];

  int i = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "phase") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_phase;
    } else if (strcmp(arg[iarg], "status") == 0) {
      // Short hand for "rheo_status"
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_status;
    } else if (strcmp(arg[iarg], "chi") == 0) {
      interface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_chi;
    } else if (strcmp(arg[iarg], "surface") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface;
    } else if (strcmp(arg[iarg], "surface/r") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_r;
    } else if (strcmp(arg[iarg], "surface/divr") == 0) {
      surface_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_surface_divr;
    } else if (strcmp(arg[iarg], "coordination") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_coordination;
    } else if (strcmp(arg[iarg], "pressure") == 0) {
      pressure_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_pressure;
    } else if (strcmp(arg[iarg], "cv") == 0) {
      thermal_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_cv;
    } else if (utils::strmatch(arg[iarg], "^surface/n/")) {
      surface_flag = 1;
      i += add_vector_component(arg[iarg], i, &ComputeRHEOPropertyAtom::pack_surface_n) - 1;
    } else if (utils::strmatch(arg[iarg], "^shift/v/")) {
      shift_flag = 1;
      i += add_vector_component(arg[iarg], i, &ComputeRHEOPropertyAtom::pack_shift_v) - 1;
    } else if (utils::strmatch(arg[iarg], "^grad/v/")) {
      i += add_tensor_component(arg[iarg], i, &ComputeRHEOPropertyAtom::pack_gradv) - 1;
    } else if (utils::strmatch(arg[iarg], "^stress/v/")) {
      i += add_tensor_component(arg[iarg], i, &ComputeRHEOPropertyAtom::pack_viscous_stress) - 1;
    } else if (utils::strmatch(arg[iarg], "^stress/t/")) {
      i += add_tensor_component(arg[iarg], i, &ComputeRHEOPropertyAtom::pack_total_stress) - 1;
    } else if (strcmp(arg[iarg], "energy") == 0) {
      avec_index[i] = atom->avec->property_atom("esph");
      if (avec_index[i] < 0)
        error->all(FLERR,
                   "Invalid keyword {} for atom style {} in compute rheo/property/atom command ",
                   arg[iarg], atom->get_style());
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_atom_style;
      thermal_flag = 1;
    } else if (strcmp(arg[iarg], "nbond/shell") == 0) {
      shell_flag = 1;
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_nbond_shell;
    } else {
      avec_index[i] = atom->avec->property_atom(arg[iarg]);
      if (avec_index[i] < 0)
        error->all(FLERR,
                   "Invalid keyword {} for atom style {} in compute rheo/property/atom command ",
                   arg[iarg], atom->get_style());
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_atom_style;

      if (strcmp(arg[iarg], "temperature") == 0) thermal_flag = 1;
      if (strcmp(arg[iarg], "heatflow") == 0) thermal_flag = 1;
      if (strcmp(arg[iarg], "conductivity") == 0) thermal_flag = 1;
    }
    i++;
  }

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOPropertyAtom::~ComputeRHEOPropertyAtom()
{
  delete[] pack_choice;
  delete[] avec_index;
  delete[] col_index;
  delete[] col_t_index;
  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0)
    error->all(FLERR, "Need to define fix rheo to use compute rheo/property/atom");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  if (interface_flag && !(fix_rheo->interface_flag))
    error->all(FLERR,
               "Cannot request interfacial property without corresponding option in fix rheo");
  if (surface_flag && !(fix_rheo->surface_flag))
    error->all(FLERR, "Cannot request surface property without corresponding option in fix rheo");
  if (shift_flag && !(fix_rheo->shift_flag))
    error->all(
        FLERR,
        "Cannot request velocity shifting property without corresponding option in fix rheo");
  if (thermal_flag && !(fix_rheo->thermal_flag))
    error->all(FLERR, "Cannot request thermal property without fix rheo/thermal");

  compute_interface = fix_rheo->compute_interface;
  compute_kernel = fix_rheo->compute_kernel;
  compute_surface = fix_rheo->compute_surface;
  compute_vshift = fix_rheo->compute_vshift;
  compute_grad = fix_rheo->compute_grad;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::setup()
{
  if (thermal_flag) {
    auto fixes = modify->get_fix_by_style("rheo/thermal");
    fix_thermal = dynamic_cast<FixRHEOThermal *>(fixes[0]);
  }

  if (pressure_flag) {
    auto fixes = modify->get_fix_by_style("rheo/pressure");
    fix_pressure = dynamic_cast<FixRHEOPressure *>(fixes[0]);
  }

  if (shell_flag) {
    auto fixes = modify->get_fix_by_style("rheo/oxidation");
    fix_oxidation = dynamic_cast<FixRHEOOxidation *>(fixes[0]);
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
      memory->create(vector_atom, nmax, "rheo/property/atom:vector");
    } else {
      memory->destroy(array_atom);
      memory->create(array_atom, nmax, nvalues, "rheo/property/atom:array");
    }
  }

  // fill vector or array with per-atom values

  if (nvalues == 1) {
    buf = vector_atom;
    (this->*pack_choice[0])(0);
  } else {
    if (nmax)
      buf = &array_atom[0][0];
    else
      buf = nullptr;
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOPropertyAtom::memory_usage()
{
  double bytes = (double) nmax * nvalues * sizeof(double);
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
  int *status = atom->rheo_status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = (status[i] & PHASECHECK);
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_status(int n)
{
  int *status = atom->rheo_status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = status[i];
    else
      buf[n] = 0.0;
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
    if (mask[i] & groupbit)
      buf[n] = chi[i];
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface(int n)
{
  int *status = atom->rheo_status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double label;
  for (int i = 0; i < nlocal; i++) {
    label = 0;
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_SURFACE) label = 1.0;
      if (status[i] & STATUS_SPLASH) label = 2.0;
    }
    buf[n] = label;
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
    if (mask[i] & groupbit)
      buf[n] = rsurface[i];
    else
      buf[n] = 0.0;
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
    if (mask[i] & groupbit)
      buf[n] = divr[i];
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_surface_n(int n)
{
  double **nsurface = compute_surface->nsurface;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int index = col_index[n];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = nsurface[i][index];
    else
      buf[n] = 0.0;
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
    if (mask[i] & groupbit)
      buf[n] = coordination[i];
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_cv(int n)
{
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = fix_thermal->calc_cv(type[i]);
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_pressure(int n)
{
  int *type = atom->type;
  int *mask = atom->mask;
  double *rho = atom->rho;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = fix_pressure->calc_pressure(rho[i], type[i]);
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_viscous_stress(int n)
{
  double **gradv = compute_grad->gradv;
  double *viscosity = atom->viscosity;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int index = col_index[n];
  int dim = domain->dimension;
  int a = index / dim;
  int b = index % dim;
  int index_transpose = b * dim + a;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = viscosity[i] * (gradv[i][index] + gradv[i][index_transpose]);
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_total_stress(int n)
{
  double **gradv = compute_grad->gradv;
  double *viscosity = atom->viscosity;
  double *rho = atom->rho;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int index = col_index[n];
  int dim = domain->dimension;
  int a = index / dim;
  int b = index % dim;
  int index_transpose = b * dim + a;
  double p;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (index == index_transpose)
        p = fix_pressure->calc_pressure(rho[i], type[i]);
      else
        p = 0.0;
      buf[n] = viscosity[i] * (gradv[i][index] + gradv[i][index_transpose]) + p;
    } else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_nbond_shell(int n)
{
  int *nbond = fix_oxidation->nbond;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = nbond[i];
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_shift_v(int n)
{
  double **vshift = compute_vshift->vshift;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int index = col_index[n];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = vshift[i][index];
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_gradv(int n)
{
  double **gradv = compute_grad->gradv;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int index = col_index[n];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      buf[n] = gradv[i][index];
    else
      buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_atom_style(int n)
{
  atom->avec->pack_property_atom(avec_index[n], &buf[n], nvalues, groupbit);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOPropertyAtom::add_tensor_component(char *option, int i, FnPtrPack pack_function)
{
  int shift;
  int dim = domain->dimension;
  if (((std::string) option).back() == '*') {
    for (int a = 0; a < dim; a++) {
      for (int b = 0; b < dim; b++) {
        pack_choice[i + a * dim + b] = pack_function;
        col_index[i + a * dim + b] = a * dim + b;
      }
    }
    shift = dim * dim;
  } else {
    int index = -1;
    int dim_error = 0;

    if (utils::strmatch(option, "xx$")) {
      index = 0;
    } else if (utils::strmatch(option, "xy$")) {
      index = 1;
    } else if (utils::strmatch(option, "xz$")) {
      index = 2;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "yx$")) {
      if (dim == 2)
        index = 2;
      else
        index = 3;
    } else if (utils::strmatch(option, "yy$")) {
      if (dim == 2)
        index = 3;
      else
        index = 4;
    } else if (utils::strmatch(option, "yz$")) {
      index = 5;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zx$")) {
      index = 6;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zy$")) {
      index = 7;
      if (dim == 2) dim_error = 1;
    } else if (utils::strmatch(option, "zz$")) {
      index = 8;
      if (dim == 2) dim_error = 1;
    } else {
      error->all(FLERR, "Invalid compute rheo/property/atom property {}", option);
    }

    if (dim_error)
      error->all(FLERR, "Invalid compute rheo/property/atom property {} in 2D", option);

    pack_choice[i] = pack_function;
    col_index[i] = index;
    shift = 1;
  }

  return shift;
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOPropertyAtom::add_vector_component(char *option, int i, FnPtrPack pack_function)
{
  int shift;
  int dim = domain->dimension;
  if (((std::string) option).back() == '*') {
    for (int a = 0; a < dim; a++) {
      pack_choice[i + a] = pack_function;
      col_index[i + a] = a;
    }
    shift = dim;
  } else {
    int index = -1;
    if (utils::strmatch(option, "x$")) {
      index = 0;
    } else if (utils::strmatch(option, "y$")) {
      index = 1;
    } else if (utils::strmatch(option, "z$")) {
      if (dim == 2)
        error->all(FLERR, "Invalid compute rheo/property/atom property {} in 2D", option);
      index = 2;
    } else {
      error->all(FLERR, "Invalid compute rheo/property/atom property {}", option);
    }

    pack_choice[i] = pack_function;
    col_index[i] = index;
    shift = 1;
  }

  return shift;
}
