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
   Contributing author: Shifeng Ke (Zhejiang University)
------------------------------------------------------------------------- */

#include "compute_fep_ta.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum { X, Y, Z };

/* ---------------------------------------------------------------------- */

ComputeFEPTA::ComputeFEPTA(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR, "Illegal number of arguments in compute fep/ta");

  scalar_flag = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  vector = new double[size_vector];

  fepinitflag = 0;    // avoid init to run entirely when called by write_data

  temp_fep = utils::numeric(FLERR, arg[3], false, lmp);

  if (strcmp(arg[4], "xy") == 0) {
    tan_axis1 = X;
    tan_axis2 = Y;
    norm_axis = Z;
  } else if (strcmp(arg[4], "xz") == 0) {
    tan_axis1 = X;
    tan_axis2 = Z;
    norm_axis = Y;
  } else if (strcmp(arg[4], "yz") == 0) {
    tan_axis1 = Y;
    tan_axis2 = Z;
    norm_axis = X;
  } else
    error->all(FLERR, "Illegal arguments in compute fep/ta");

  scale_factor = utils::numeric(FLERR, arg[5], false, lmp);

  // optional keywords

  tailflag = 0;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "tail") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal optional keyword in compute fep/ta");
      tailflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal optional keyword in compute fep/ta");
  }

  // allocate space for position, force, energy, virial arrays

  x_orig = nullptr;
  f_orig = nullptr;
  peatom_orig = keatom_orig = nullptr;
  pvatom_orig = kvatom_orig = nullptr;

  allocate_storage();

  fixgpu = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeFEPTA::~ComputeFEPTA()
{
  delete[] vector;

  deallocate_storage();
}

/* ---------------------------------------------------------------------- */

void ComputeFEPTA::init()
{
  if (!fepinitflag)    // avoid init to run entirely when called by write_data
    fepinitflag = 1;
  else
    return;

  // setup and error checks

  if (domain->dimension == 2) { error->all(FLERR, "Cannot compute fep/ta in 2d simulation"); }

  if (tailflag) {
    if (force->pair->tail_flag == 0)
      error->all(FLERR,
                 "Compute fep/ta tail when pair style does not "
                 "compute tail corrections");
  }

  // detect if package gpu is present

  fixgpu = modify->get_fix_by_id("package_gpu");

  if (comm->me == 0)
    utils::logmesg(lmp,
                   "FEP/TA settings ...\n"
                   "  temperature = {:f}\n"
                   "  scale factor = {:f}\n"
                   "  tail {}\n",
                   temp_fep, scale_factor, tailflag ? "yes" : "no");
}

/* ---------------------------------------------------------------------- */

void ComputeFEPTA::compute_vector()
{
  double pe0, pe1;

  eflag = 1;
  vflag = 0;

  invoked_vector = update->ntimestep;

  if (atom->nmax > nmax) {    // reallocate working arrays if necessary
    deallocate_storage();
    allocate_storage();
  }

  backup_xfev();    // backup position, force, energy, virial array values
  backup_box();     // backup box size

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(eflag, vflag);
    timer->stamp(Timer::PAIR);
  }

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
    timer->stamp(Timer::BOND);
  }

  if (force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(eflag, vflag);
    timer->stamp(Timer::KSPACE);
  }

  // accumulate force/energy/virial from /gpu pair styles
  // this is required as to empty the answer queue,
  // otherwise the force compute on the GPU in the next step would be incorrect
  if (fixgpu) fixgpu->post_force(vflag);

  pe0 = compute_pe();

  change_box();

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(eflag, vflag);
    timer->stamp(Timer::PAIR);
  }

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
    timer->stamp(Timer::BOND);
  }

  if (force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(eflag, vflag);
    timer->stamp(Timer::KSPACE);
  }

  // accumulate force/energy/virial from /gpu pair styles
  // this is required as to empty the answer queue,
  // otherwise the force compute on the GPU in the next step would be incorrect
  if (fixgpu) fixgpu->post_force(vflag);

  pe1 = compute_pe();

  restore_xfev();    // restore position, force, energy, virial array values
  restore_box();     // restore box size

  vector[0] = pe1 - pe0;
  vector[1] = exp(-(pe1 - pe0) / (force->boltz * temp_fep));
  vector[2] = area_orig * (scale_factor - 1.0);
}

/* ----------------------------------------------------------------------
   obtain potential energy from lammps accumulators
------------------------------------------------------------------------- */

double ComputeFEPTA::compute_pe()
{
  double eng, eng_potential;

  eng = 0.0;
  if (force->pair) eng = force->pair->eng_vdwl + force->pair->eng_coul;

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) eng += force->bond->energy;
    if (force->angle) eng += force->angle->energy;
    if (force->dihedral) eng += force->dihedral->energy;
    if (force->improper) eng += force->improper->energy;
  }

  MPI_Allreduce(&eng, &eng_potential, 1, MPI_DOUBLE, MPI_SUM, world);

  if (tailflag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    eng_potential += force->pair->etail / volume;
  }

  if (force->kspace) eng_potential += force->kspace->energy;

  return eng_potential;
}

/* ----------------------------------------------------------------------
   apply changes to box
------------------------------------------------------------------------- */

void ComputeFEPTA::change_box()
{
  int i;
  double **x = atom->x;
  int natom = atom->nlocal + atom->nghost;

  for (i = 0; i < natom; i++) domain->x2lamda(x[i], x[i]);

  domain->boxhi[tan_axis1] *= sqrt(scale_factor);
  domain->boxlo[tan_axis1] *= sqrt(scale_factor);
  domain->boxhi[tan_axis2] *= sqrt(scale_factor);
  domain->boxlo[tan_axis2] *= sqrt(scale_factor);
  domain->boxhi[norm_axis] /= scale_factor;
  domain->boxlo[norm_axis] /= scale_factor;

  domain->set_global_box();
  domain->set_local_box();

  // remap atom position
  for (i = 0; i < natom; i++) domain->lamda2x(x[i], x[i]);

  if (force->kspace) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   backup box size
------------------------------------------------------------------------- */

void ComputeFEPTA::backup_box()
{
  for (int i = 0; i < domain->dimension; i++) {
    boxhi_orig[i] = domain->boxhi[i];
    boxlo_orig[i] = domain->boxlo[i];
  }

  area_orig = domain->prd[tan_axis1] * domain->prd[tan_axis2];
}

/* ----------------------------------------------------------------------
   restore box size to original values
------------------------------------------------------------------------- */

void ComputeFEPTA::restore_box()
{
  for (int i = 0; i < domain->dimension; i++) {
    domain->boxhi[i] = boxhi_orig[i];
    domain->boxlo[i] = boxlo_orig[i];
  }

  domain->set_global_box();
  domain->set_local_box();

  if (force->kspace) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   manage storage for position, force, energy, virial arrays
------------------------------------------------------------------------- */

void ComputeFEPTA::allocate_storage()
{
  nmax = atom->nmax;
  memory->create(x_orig, nmax, 3, "fep:x_orig");
  memory->create(f_orig, nmax, 3, "fep:f_orig");
  memory->create(peatom_orig, nmax, "fep:peatom_orig");
  memory->create(pvatom_orig, nmax, 6, "fep:pvatom_orig");
  if (force->kspace) {
    memory->create(keatom_orig, nmax, "fep:keatom_orig");
    memory->create(kvatom_orig, nmax, 6, "fep:kvatom_orig");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFEPTA::deallocate_storage()
{
  memory->destroy(x_orig);
  memory->destroy(f_orig);
  memory->destroy(peatom_orig);
  memory->destroy(pvatom_orig);
  memory->destroy(keatom_orig);
  memory->destroy(kvatom_orig);

  x_orig = nullptr;
  f_orig = nullptr;
  peatom_orig = keatom_orig = nullptr;
  pvatom_orig = kvatom_orig = nullptr;
}

/* ----------------------------------------------------------------------
   backup and restore arrays with position, force, energy, virial
------------------------------------------------------------------------- */

void ComputeFEPTA::backup_xfev()
{
  int i;

  int natom = atom->nlocal + atom->nghost;

  double **x = atom->x;
  for (i = 0; i < natom; i++) {
    x_orig[i][0] = x[i][0];
    x_orig[i][1] = x[i][1];
    x_orig[i][2] = x[i][2];
  }

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f_orig[i][0] = f[i][0];
    f_orig[i][1] = f[i][1];
    f_orig[i][2] = f[i][2];
  }

  eng_vdwl_orig = force->pair->eng_vdwl;
  eng_coul_orig = force->pair->eng_coul;

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) eng_bond_orig = force->bond->energy;
    if (force->angle) eng_angle_orig = force->angle->energy;
    if (force->dihedral) eng_dihedral_orig = force->dihedral->energy;
    if (force->improper) eng_improper_orig = force->improper->energy;
  }

  pvirial_orig[0] = force->pair->virial[0];
  pvirial_orig[1] = force->pair->virial[1];
  pvirial_orig[2] = force->pair->virial[2];
  pvirial_orig[3] = force->pair->virial[3];
  pvirial_orig[4] = force->pair->virial[4];
  pvirial_orig[5] = force->pair->virial[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++) peatom_orig[i] = peatom[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom_orig[i][0] = pvatom[i][0];
      pvatom_orig[i][1] = pvatom[i][1];
      pvatom_orig[i][2] = pvatom[i][2];
      pvatom_orig[i][3] = pvatom[i][3];
      pvatom_orig[i][4] = pvatom[i][4];
      pvatom_orig[i][5] = pvatom[i][5];
    }
  }

  if (force->kspace) {
    energy_orig = force->kspace->energy;
    kvirial_orig[0] = force->kspace->virial[0];
    kvirial_orig[1] = force->kspace->virial[1];
    kvirial_orig[2] = force->kspace->virial[2];
    kvirial_orig[3] = force->kspace->virial[3];
    kvirial_orig[4] = force->kspace->virial[4];
    kvirial_orig[5] = force->kspace->virial[5];

    if (update->eflag_atom) {
      double *keatom = force->kspace->eatom;
      for (i = 0; i < natom; i++) keatom_orig[i] = keatom[i];
    }
    if (update->vflag_atom) {
      double **kvatom = force->kspace->vatom;
      for (i = 0; i < natom; i++) {
        kvatom_orig[i][0] = kvatom[i][0];
        kvatom_orig[i][1] = kvatom[i][1];
        kvatom_orig[i][2] = kvatom[i][2];
        kvatom_orig[i][3] = kvatom[i][3];
        kvatom_orig[i][4] = kvatom[i][4];
        kvatom_orig[i][5] = kvatom[i][5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFEPTA::restore_xfev()
{
  int i;

  int natom = atom->nlocal + atom->nghost;

  double **x = atom->x;
  for (i = 0; i < natom; i++) {
    x[i][0] = x_orig[i][0];
    x[i][1] = x_orig[i][1];
    x[i][2] = x_orig[i][2];
  }

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f[i][0] = f_orig[i][0];
    f[i][1] = f_orig[i][1];
    f[i][2] = f_orig[i][2];
  }

  force->pair->eng_vdwl = eng_vdwl_orig;
  force->pair->eng_coul = eng_coul_orig;

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->energy = eng_bond_orig;
    if (force->angle) force->angle->energy = eng_angle_orig;
    if (force->dihedral) force->dihedral->energy = eng_dihedral_orig;
    if (force->improper) force->improper->energy = eng_improper_orig;
  }

  force->pair->virial[0] = pvirial_orig[0];
  force->pair->virial[1] = pvirial_orig[1];
  force->pair->virial[2] = pvirial_orig[2];
  force->pair->virial[3] = pvirial_orig[3];
  force->pair->virial[4] = pvirial_orig[4];
  force->pair->virial[5] = pvirial_orig[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++) peatom[i] = peatom_orig[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom[i][0] = pvatom_orig[i][0];
      pvatom[i][1] = pvatom_orig[i][1];
      pvatom[i][2] = pvatom_orig[i][2];
      pvatom[i][3] = pvatom_orig[i][3];
      pvatom[i][4] = pvatom_orig[i][4];
      pvatom[i][5] = pvatom_orig[i][5];
    }
  }

  if (force->kspace) {
    force->kspace->energy = energy_orig;
    force->kspace->virial[0] = kvirial_orig[0];
    force->kspace->virial[1] = kvirial_orig[1];
    force->kspace->virial[2] = kvirial_orig[2];
    force->kspace->virial[3] = kvirial_orig[3];
    force->kspace->virial[4] = kvirial_orig[4];
    force->kspace->virial[5] = kvirial_orig[5];

    if (update->eflag_atom) {
      double *keatom = force->kspace->eatom;
      for (i = 0; i < natom; i++) keatom[i] = keatom_orig[i];
    }
    if (update->vflag_atom) {
      double **kvatom = force->kspace->vatom;
      for (i = 0; i < natom; i++) {
        kvatom[i][0] = kvatom_orig[i][0];
        kvatom[i][1] = kvatom_orig[i][1];
        kvatom[i][2] = kvatom_orig[i][2];
        kvatom[i][3] = kvatom_orig[i][3];
        kvatom[i][4] = kvatom_orig[i][4];
        kvatom[i][5] = kvatom_orig[i][5];
      }
    }
  }
}
