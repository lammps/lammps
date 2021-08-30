/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Rene Halver (JSC)
------------------------------------------------------------------------- */

#include "scafacos.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "pair_hybrid.h"

#include <cstdlib>
#include <cstring>

// ScaFaCoS library

#include "fcs.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Scafacos::Scafacos(LAMMPS *lmp) : KSpace(lmp)
{
  me = comm->me;
  initialized = 0;

  maxatom = 0;
  xpbc = nullptr;
  epot = nullptr;
  efield = nullptr;
  fcs = nullptr;
}

/* ---------------------------------------------------------------------- */

void Scafacos::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal scafacos command");

  method = utils::strdup(arg[0]);
  tolerance = utils::numeric(FLERR, arg[1], false, lmp);

  // optional ScaFaCoS library setting defaults
  // choose the correct default tolerance type for chosen method
  // throw an error if a not yet supported solver is chosen
  if (strcmp(method, "fmm") == 0) {
    tolerance_type = FCS_TOLERANCE_TYPE_ENERGY;
    fmm_tuning_flag = 0;
  } else if (strcmp(method, "p3m") == 0 || strcmp(method, "p2nfft") == 0 ||
             strcmp(method, "ewald") == 0) {
    tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
  } else if (strcmp(method, "direct") == 0) {
    ;    // direct summation has no tolerance type
  } else {
    error->all(FLERR, "Unsupported ScaFaCoS method");
  }
}

/* ---------------------------------------------------------------------- */

Scafacos::~Scafacos()
{
  delete[] method;

  memory->destroy(xpbc);
  memory->destroy(epot);
  memory->destroy(efield);

  // clean up of the ScaFaCoS handle and internal arrays
  fcs_destroy((FCS) fcs);
}

/* ---------------------------------------------------------------------- */

void Scafacos::init()
{
  // error checks
  if (me == 0) {
    utils::logmesg(lmp, "Setting up ScaFaCoS with solver {} ...\n", method);

    if (strcmp(method, "p3m") == 0)
      error->warning(FLERR, "Virial computation for P3M not available");

    if (strcmp(method, "ewald") == 0)
      error->warning(FLERR, "Virial computation for Ewald not available");
  }

  if (!atom->q_flag) error->all(FLERR, "Kspace style requires atom attribute q");

  if (domain->dimension == 2) error->all(FLERR, "Cannot use ScaFaCoS with 2d simulation");

  if (domain->triclinic) error->all(FLERR, "Cannot use ScaFaCoS with triclinic domain yet");

  if (atom->natoms > INT_MAX && sizeof(int) != 8)
    error->all(FLERR, "ScaFaCoS atom count exceeds 2B");

  if ((atom->molecular != Atom::ATOMIC) && (atom->nbonds + atom->nangles + atom->ndihedrals) > 0) {
    int flag = 0;
    if ((force->special_coul[1] == 1.0) && (force->special_coul[2] == 1.0) &&
        (force->special_coul[3] == 1.0))
      ++flag;

    PairHybrid *ph = reinterpret_cast<PairHybrid *>(force->pair_match("^hybrid", 0));
    if (ph) {
      for (int isub = 0; isub < ph->nstyles; ++isub)
        if (force->pair_match("coul/exclude", 0, isub)) ++flag;
    } else {
      if (force->pair_match("coul/exclude", 0)) ++flag;
    }
    if (!flag)
      error->all(FLERR,
                 "Must use pair style coul/exclude or 'special_bonds coul 1.0 1.0 1.0'"
                 "for molecular charged systems with ScaFaCos");
  }

  FCSResult result;

  // one-time initialization of ScaFaCoS

  qqrd2e = force->qqrd2e;

  if (!initialized) {
    result = fcs_init((FCS *) &fcs, method, world);
    check_result((void *) &result);

    setup_handle();

    // using other methods lead to termination of the program,
    // since they have no tolerance tuning methods
    if (strcmp(method, "fmm") == 0 || strcmp(method, "p3m") == 0 || strcmp(method, "p2nfft") == 0 ||
        strcmp(method, "ewald") == 0) {
      result = fcs_set_tolerance((FCS) fcs, tolerance_type, tolerance);
      check_result((void *) &result);
    }

    double **x = atom->x;
    double *q = atom->q;
    int nlocal = atom->nlocal;

    if (strcmp(method, "fmm") == 0) {
      if (fmm_tuning_flag == 1)
        fcs_fmm_set_internal_tuning((FCS) fcs, FCS_FMM_INHOMOGENOUS_SYSTEM);
      else
        fcs_fmm_set_internal_tuning((FCS) fcs, FCS_FMM_HOMOGENOUS_SYSTEM);
    }

    // for the FMM at least one particle is required per process
    if (strcmp(method, "fmm") == 0) {
      int empty = (nlocal == 0) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_SUM, world);
      if (empty > 0)
        fcs_set_redistribute((FCS) fcs, 1);
      else
        fcs_set_redistribute((FCS) fcs, 0);
    }

    result = fcs_tune((FCS) fcs, nlocal, &x[0][0], q);
    check_result((void *) &result);
    // more useful here, since the parameters should be tuned now
    if (me == 0) fcs_print_parameters((FCS) fcs);
  }

  initialized = 1;
}

/* ---------------------------------------------------------------------- */

void Scafacos::compute(int eflag, int vflag)
{
  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  const double qscale = qqrd2e;
  FCSResult result;

  // for the FMM at least one particle is required per process
  if (strcmp(method, "fmm")) {
    int empty = (nlocal == 0) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_SUM, world);
    if (empty > 0)
      fcs_set_redistribute((FCS) fcs, 1);
    else
      fcs_set_redistribute((FCS) fcs, 0);
  }

  ev_init(eflag, vflag);

  // grow xpbc, epot, efield if necessary

  if (nlocal > maxatom || maxatom == 0) {
    memory->destroy(xpbc);
    memory->destroy(epot);
    memory->destroy(efield);
    maxatom = atom->nmax;
    memory->create(xpbc, 3 * maxatom, "scafacos:xpbc");
    memory->create(epot, maxatom, "scafacos:epot");
    memory->create(efield, maxatom, 3, "scafacos:efield");
  }

  if (vflag_global) {

    // for P3M or Ewald we cannot compute the virial. skip it.

    if ((strcmp(method, "p3m") != 0) && (strcmp(method, "ewald") != 0)) {
      result = fcs_set_compute_virial((FCS) fcs, 1);
      check_result((void *) &result);
    }
  }

  // pack coords into xpbc and apply PBC
  memcpy(xpbc, &x[0][0], 3 * nlocal * sizeof(double));

  if (domain->xperiodic || domain->yperiodic || domain->zperiodic) {
    int j = 0;
    for (int i = 0; i < nlocal; i++) {
      domain->remap(&xpbc[j]);
      j += 3;
    }
  }

  // if simulation box has changed, call fcs_tune()

  if (box_has_changed()) {
    setup_handle();
    result = fcs_tune((FCS) fcs, nlocal, xpbc, q);
    check_result((void *) &result);
  }

  // invoke ScaFaCoS solver

  result = fcs_run((FCS) fcs, nlocal, xpbc, q, &efield[0][0], epot);
  check_result((void *) &result);

  // extract virial

  if (vflag_global) {

    // for P3M or Ewald we cannot compute the virial. skip it.

    if ((strcmp(method, "p3m") != 0) && (strcmp(method, "ewald") != 0)) {
      result = fcs_get_virial((FCS) fcs, virial_int);
      check_result((void *) &result);

      virial[0] = virial_int[0];
      virial[1] = virial_int[1];
      virial[2] = virial_int[2];
      virial[3] = virial_int[4];
      virial[4] = virial_int[5];
      virial[5] = virial_int[8];
    }
  }

  // apply Efield to each particle
  // accumulate total energy

  double **f = atom->f;

  double qone;
  double myeng = 0.0;

  for (int i = 0; i < nlocal; i++) {
    qone = q[i] * qscale;
    f[i][0] += qone * efield[i][0];
    f[i][1] += qone * efield[i][1];
    f[i][2] += qone * efield[i][2];
    myeng += 0.5 * qone * epot[i];
  }

  if (eflag_atom) {
    for (int i = 0; i < nlocal; i++) eatom[i] = 0.5 * qscale * q[i] * epot[i];
  }

  MPI_Allreduce(&myeng, &energy, 1, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

int Scafacos::modify_param(int narg, char **arg)
{
  // add any Scafacos options here you want to expose to LAMMPS
  // syntax: kspace_modify scafacos keyword value1 value2 ...
  //   keyword = tolerance
  //     value1 = energy, energy_rel, etc
  // everyone of these should have a default, so user doesn't need to set

  if (strcmp(arg[0], "scafacos") != 0) return 0;

  if (strcmp(arg[1], "tolerance") == 0) {
    if (narg < 3) error->all(FLERR, "Illegal kspace_modify command (tolerance)");
    if (strcmp(arg[2], "energy") == 0)
      tolerance_type = FCS_TOLERANCE_TYPE_ENERGY;
    else if (strcmp(arg[2], "energy_rel") == 0)
      tolerance_type = FCS_TOLERANCE_TYPE_ENERGY_REL;
    else if (strcmp(arg[2], "field") == 0)
      tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
    else if (strcmp(arg[2], "field_rel") == 0)
      tolerance_type = FCS_TOLERANCE_TYPE_FIELD_REL;
    else if (strcmp(arg[2], "potential") == 0)
      tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL;
    else if (strcmp(arg[2], "potential_rel") == 0)
      tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL_REL;
    else
      error->all(FLERR, "Illegal kspace_modify command (tolerance argument)");
    // check if method is compatatible to chosen tolerance type
    if ((strcmp(method, "fmm") == 0 &&
         (tolerance_type != FCS_TOLERANCE_TYPE_ENERGY &&
          tolerance_type != FCS_TOLERANCE_TYPE_ENERGY_REL)) ||
        (strcmp(method, "p2nfft") == 0 &&
         (tolerance_type != FCS_TOLERANCE_TYPE_FIELD &&
          tolerance_type != FCS_TOLERANCE_TYPE_POTENTIAL)) ||
        (strcmp(method, "p3m") == 0 && (tolerance_type != FCS_TOLERANCE_TYPE_FIELD)) ||
        (strcmp(method, "ewald") == 0 && (tolerance_type != FCS_TOLERANCE_TYPE_FIELD)))
      error->all(FLERR, "Illegal kspace_modify command \
                          (invalid tolerance / method combination)");
    return 3;
  }

  //   keyword = fmm_inhomogen_tuning
  //     value1 = 0, 1
  //       0 -> homogenous system (default)
  //       1 -> inhomogenous system (more internal tuning is provided (sequential!))
  if (strcmp(arg[1], "fmm_tuning") == 0) {
    if (me == 0) utils::logmesg(lmp, "ScaFaCoS setting fmm inhomogen tuning ...\n");

    if (narg < 3) error->all(FLERR, "Illegal kspace_modify command (fmm_tuning)");
    fmm_tuning_flag = utils::inumeric(FLERR, arg[2], false, lmp);
    return 3;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double Scafacos::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) maxatom * sizeof(double);
  bytes += (double) 3 * maxatom * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
    setup of ScaFaCoS handle with common parameters
------------------------------------------------------------------------- */

void Scafacos::setup_handle()
{
  FCSResult result;

  // store simulation box params

  // setup periodicity
  old_periodicity[0] = domain->xperiodic;
  old_periodicity[1] = domain->yperiodic;
  old_periodicity[2] = domain->zperiodic;

  // setup box origin (lower left front corner of the system)
  old_origin[0] = domain->boxlo[0];
  old_origin[1] = domain->boxlo[1];
  old_origin[2] = domain->boxlo[2];

  // setup box vectors (base vectors of the system box)
  old_box_x[0] = domain->prd[0];
  old_box_x[1] = old_box_x[2] = 0.0;
  old_box_y[1] = domain->prd[1];
  old_box_y[0] = old_box_y[2] = 0.0;
  old_box_z[2] = domain->prd[2];
  old_box_z[1] = old_box_z[0] = 0.0;

  // setup number of atoms in the system
  old_natoms = atom->natoms;

  // store parameters to ScaFaCoS handle
  result = fcs_set_box_a((FCS) fcs, old_box_x);
  check_result((void *) &result);

  result = fcs_set_box_b((FCS) fcs, old_box_y);
  check_result((void *) &result);

  result = fcs_set_box_c((FCS) fcs, old_box_z);
  check_result((void *) &result);

  result = fcs_set_box_origin((FCS) fcs, old_origin);
  check_result((void *) &result);

  result = fcs_set_periodicity((FCS) fcs, old_periodicity);
  check_result((void *) &result);

  result = fcs_set_total_particles((FCS) fcs, old_natoms);
  check_result((void *) &result);

  // allow ScaFaCoS to calculate the near field computations for now
  // TODO: allow the delegation of the near field computations
  //       within LAMMPS
  //       (near_field_flag = 1 -> enables the internal near field calcs
  //                          0 -> disables the internal near field calcs
  int near_field_flag = 1;
  result = fcs_set_near_field_flag((FCS) fcs, near_field_flag);
  check_result((void *) &result);
}

/* ----------------------------------------------------------------------
    check if box parameters changed, requiring a new call to fcs_tune
------------------------------------------------------------------------- */

bool Scafacos::box_has_changed()
{
  int *periodicity = domain->periodicity;
  double *prd = domain->prd;

  bool changed = (periodicity[0] != old_periodicity[0]) || (periodicity[1] != old_periodicity[1]) ||
      (periodicity[2] != old_periodicity[2]) || (domain->boundary[0][0] != old_origin[0]) ||
      (domain->boundary[1][0] != old_origin[1]) || (domain->boundary[2][0] != old_origin[2]) ||
      (prd[0] != old_box_x[0]) || (prd[1] != old_box_y[1]) || (prd[2] != old_box_z[2]) ||
      (atom->natoms != old_natoms);

  return changed;
}

/* ----------------------------------------------------------------------
   check ScaFaCoS result for error condition
------------------------------------------------------------------------- */

void Scafacos::check_result(void *result_p)
{
  FCSResult result = *(FCSResult *) result_p;
  if (!result) return;

  error->one(FLERR, "ScaFaCoS: {}\n{}\n", fcs_result_get_function(result),
             fcs_result_get_message(result));
}
