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
   Contributing author: Agilio Padua (ENS de Lyon & CNRS)
------------------------------------------------------------------------- */

#include "compute_mbar.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "timer.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { PAIR, ATOM };
enum { CHARGE };

/* ----------------------------------------------------------------------
   number of args consumed by a grid spec: 1 for a vector-style variable
   (v_name), 3 for the bare numeric form "lo hi n"
------------------------------------------------------------------------- */

static int grid_nargs(const char *arg)
{
  return utils::strmatch(arg, "^v_") ? 1 : 3;
}

/* ---------------------------------------------------------------------- */

ComputeMBAR::ComputeMBAR(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR, "Illegal number of arguments in compute mbar");

  scalar_flag = 0;
  vector_flag = 1;
  extvector = 0;

  const int ntypes = atom->ntypes;
  mbarinitflag = 0;    // avoid init to run entirely when called by write_data

  temp_mbar = utils::numeric(FLERR, arg[3], false, lmp);

  // each perturbation supplies its own grid (a vector-style variable) holding
  // the absolute values that the perturbed parameter takes at each state. The
  // number of states is the common length of these grids.

  int iarg = 4;
  const int pertstart = iarg;

  npert = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "pair") == 0) {
      if (iarg + 6 > narg) error->all(FLERR, "Illegal pair attribute in compute mbar");
      int g = grid_nargs(arg[iarg + 5]);
      if (iarg + 5 + g > narg) error->all(FLERR, "Illegal pair attribute in compute mbar");
      npert++;
      iarg += 5 + g;
    } else if (strcmp(arg[iarg], "atom") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal atom attribute in compute mbar");
      int g = grid_nargs(arg[iarg + 3]);
      if (iarg + 3 + g > narg) error->all(FLERR, "Illegal atom attribute in compute mbar");
      npert++;
      iarg += 3 + g;
    } else
      break;
  }

  if (npert == 0) error->all(FLERR, "Illegal syntax in compute mbar");
  perturb = new Perturb[npert];
  for (int m = 0; m < npert; m++) {
    perturb[m].gridname = nullptr;
    perturb[m].grid = nullptr;
  }

  // parse perturbation keywords and resolve each grid variable

  npert = 0;
  nlambda = 0;
  chgflag = 0;

  iarg = pertstart;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "pair") == 0) {
      perturb[npert].which = PAIR;
      perturb[npert].pstyle = utils::strdup(arg[iarg + 1]);
      perturb[npert].pparam = utils::strdup(arg[iarg + 2]);
      utils::bounds(FLERR, arg[iarg + 3], 1, ntypes, perturb[npert].ilo, perturb[npert].ihi, error);
      utils::bounds(FLERR, arg[iarg + 4], 1, ntypes, perturb[npert].jlo, perturb[npert].jhi, error);
      {
        int g = grid_nargs(arg[iarg + 5]);
        set_grid(&arg[iarg + 5], g, npert);
        iarg += 5 + g;
      }
      npert++;
    } else if (strcmp(arg[iarg], "atom") == 0) {
      perturb[npert].which = ATOM;
      if (strcmp(arg[iarg + 1], "charge") == 0) {
        perturb[npert].aparam = CHARGE;
        chgflag = 1;
      } else
        error->all(FLERR, "Illegal atom argument in compute mbar");
      utils::bounds(FLERR, arg[iarg + 2], 1, ntypes, perturb[npert].ilo, perturb[npert].ihi, error);
      {
        int g = grid_nargs(arg[iarg + 3]);
        set_grid(&arg[iarg + 3], g, npert);
        iarg += 3 + g;
      }
      npert++;
    } else
      break;
  }

  size_vector = nlambda;
  vector = new double[size_vector];

  // optional keywords

  tailflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "tail") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal optional keyword in compute mbar");
      tailflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal optional keyword in compute mbar");
  }

  // allocate pair style arrays

  for (int m = 0; m < npert; m++) {
    if (perturb[m].which == PAIR)
      memory->create(perturb[m].array_orig, ntypes + 1, ntypes + 1, "mbar:array_orig");
  }

  // charge, force, energy, virial per-atom arrays are allocated in init()

  f_orig = nullptr;
  q_orig = nullptr;
  peatom_orig = keatom_orig = nullptr;
  pvatom_orig = kvatom_orig = nullptr;

  fixgpu = nullptr;
}

/* ----------------------------------------------------------------------
   resolve the grid spec for perturbation m into perturb[m].grid, the
   absolute values that the perturbed parameter takes at each state. The grid
   is given either as a single vector-style variable (nargs == 1, arg[0] is
   v_name) or as the bare numeric form "lo hi n" (nargs == 3) producing n
   equally spaced values from lo to hi. All grids must share the same length,
   which is the number of states to sample, nlambda.
------------------------------------------------------------------------- */

void ComputeMBAR::set_grid(char **arg, int nargs, int m)
{
  double *grid;
  double *linear = nullptr;
  int n;

  if (nargs == 1) {    // single vector-style variable v_name

    if (!utils::strmatch(arg[0], "^v_"))
      error->all(FLERR, "Grid for compute mbar perturbation must be a vector-style variable");

    perturb[m].gridname = utils::strdup(arg[0] + 2);
    int gridvar = input->variable->find(perturb[m].gridname);
    if (gridvar < 0)
      error->all(FLERR, "Variable name {} for compute mbar does not exist", perturb[m].gridname);
    if (!input->variable->vectorstyle(gridvar))
      error->all(FLERR, "Variable {} for compute mbar must be vector style", perturb[m].gridname);

    n = input->variable->compute_vector(gridvar, &grid);
    if (n == 0) error->all(FLERR, "No grid values in compute mbar");

  } else {    // bare numeric form: lo hi n

    double lo = utils::numeric(FLERR, arg[0], false, lmp);
    double hi = utils::numeric(FLERR, arg[1], false, lmp);
    n = utils::inumeric(FLERR, arg[2], false, lmp);
    if (n < 2) error->all(FLERR, "Number of states in compute mbar grid must be >= 2");

    linear = new double[n];
    for (int k = 0; k < n; k++) linear[k] = lo + (hi - lo) * k / (n - 1);
    grid = linear;
  }

  if (nlambda == 0)
    nlambda = n;
  else if (n != nlambda)
    error->all(FLERR, "All perturbation grids in compute mbar must have the same length");

  perturb[m].grid = new double[n];
  for (int k = 0; k < n; k++) perturb[m].grid[k] = grid[k];

  delete[] linear;
}

/* ---------------------------------------------------------------------- */

ComputeMBAR::~ComputeMBAR()
{
  delete[] vector;

  for (int m = 0; m < npert; m++) {
    if (perturb[m].which == PAIR) {
      delete[] perturb[m].pstyle;
      delete[] perturb[m].pparam;
      memory->destroy(perturb[m].array_orig);
    }
    delete[] perturb[m].gridname;
    delete[] perturb[m].grid;
  }
  delete[] perturb;

  deallocate_storage();
}

/* ---------------------------------------------------------------------- */

void ComputeMBAR::init()
{
  int i, j;

  if (!mbarinitflag)    // avoid init to run entirely when called by write_data
    mbarinitflag = 1;
  else
    return;

  // setup and error checks

  pairflag = 0;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];

    if (force->pair == nullptr) error->all(FLERR, "compute mbar pair requires pair interactions");

    if (pert->which == PAIR) {
      pairflag = 1;
      Pair *pair = nullptr;
      if (lmp->suffix_enable) {
        if (lmp->suffix)
          pair = force->pair_match(fmt::format("{}/{}", pert->pstyle, lmp->suffix), 1);
        if ((pair == nullptr) && lmp->suffix2)
          pair = force->pair_match(fmt::format("{}/{}", pert->pstyle, lmp->suffix2), 1);
      }

      if (pair == nullptr) pair = force->pair_match(pert->pstyle, 1);
      if (pair == nullptr) error->all(FLERR, "Compute mbar pair style {} not found", pert->pstyle);

      void *ptr = pair->extract(pert->pparam, pert->pdim);
      if (ptr == nullptr)
        error->all(FLERR, "Compute mbar pair style {} param {} not supported", pert->pstyle,
                   pert->pparam);

      pert->array = (double **) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if ((strcmp(force->pair_style, "hybrid") == 0 ||
           strcmp(force->pair_style, "hybrid/overlay") == 0)) {
        auto *pair = dynamic_cast<PairHybrid *>(force->pair);
        for (i = pert->ilo; i <= pert->ihi; i++)
          for (j = MAX(pert->jlo, i); j <= pert->jhi; j++)
            if (!pair->check_ijtype(i, j, pert->pstyle))
              error->all(FLERR,
                         "compute mbar type pair range is not valid for "
                         "pair hybrid sub-style");
      }

    } else if (pert->which == ATOM) {
      if (pert->aparam == CHARGE) {
        if (!atom->q_flag) error->all(FLERR, "compute mbar requires atom attribute charge");
      }
    }
  }

  if (tailflag) {
    if (force->pair->tail_flag == 0)
      error->all(FLERR,
                 "Compute mbar tail when pair style does not "
                 "compute tail corrections");
  }

  // allocate per-atom storage now that force->kspace is non-null

  allocate_storage();

  // detect if package gpu is present

  fixgpu = modify->get_fix_by_id("package_gpu");

  if (comm->me == 0) {
    auto mesg = fmt::format("MBAR settings ...\n  temperature = {:f}\n", temp_mbar);
    mesg += fmt::format("  tail {}\n", (tailflag ? "yes" : "no"));
    mesg += fmt::format("  states = {}\n", nlambda);
    for (int m = 0; m < npert; m++) {
      Perturb *pert = &perturb[m];
      if (pert->which == PAIR)
        mesg += fmt::format("  pair {} {} {}-{} {}-{} grid:", pert->pstyle, pert->pparam, pert->ilo,
                            pert->ihi, pert->jlo, pert->jhi);
      else if (pert->which == ATOM)
        mesg += fmt::format("  atom charge {}-{} grid:", pert->ilo, pert->ihi);
      for (int k = 0; k < nlambda; k++) mesg += fmt::format(" {:g}", pert->grid[k]);
      mesg += "\n";
    }
    utils::logmesg(lmp, mesg);
  }
}

/* ----------------------------------------------------------------------
   for the current configuration, evaluate the reduced potential energy
   U(lambda_k)/kT at each lambda state, producing one row of the u_kn
   matrix to be used by an external MBAR solver such as pymbar
------------------------------------------------------------------------- */

void ComputeMBAR::compute_vector()
{
  // flag that we only need to compute the global energy
  int eflag = ENERGY_GLOBAL | ENERGY_ONLY;
  int vflag = VIRIAL_NONE;

  invoked_vector = update->ntimestep;

  if (atom->nmax > nmax) {    // reallocate working arrays if necessary
    deallocate_storage();
    allocate_storage();
  }

  backup_qfev();      // backup charge, force, energy, virial array values
  backup_params();    // backup pair parameters

  const double kT = force->boltz * temp_mbar;

  for (int k = 0; k < nlambda; k++) {

    // compute with perturbation parameters at each lambda state k

    perturb_params(k);

    timer->stamp();
    if (force->pair && force->pair->compute_flag) {
      force->pair->compute(eflag, vflag);
      timer->stamp(Timer::PAIR);
    }
    if (chgflag && force->kspace && force->kspace->compute_flag) {
      force->kspace->compute(eflag, vflag);
      timer->stamp(Timer::KSPACE);
    }

    // accumulate force/energy/virial from /gpu pair styles
    // this is required as to empty the answer queue,
    // otherwise the force compute on the GPU in the next step would be incorrect
    if (fixgpu) fixgpu->post_force(vflag);

    vector[k] = compute_epair() / kT;

    restore_qfev();      // restore charge, force, energy, virial array values
    restore_params();    // restore pair parameters
  }
}

/* ----------------------------------------------------------------------
   obtain pair energy from lammps accumulators
------------------------------------------------------------------------- */

double ComputeMBAR::compute_epair()
{
  double eng, eng_pair;

  eng = 0.0;
  if (force->pair) eng = force->pair->eng_vdwl + force->pair->eng_coul;
  MPI_Allreduce(&eng, &eng_pair, 1, MPI_DOUBLE, MPI_SUM, world);

  if (tailflag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    eng_pair += force->pair->etail / volume;
  }

  if (chgflag && force->kspace) eng_pair += force->kspace->energy;

  return eng_pair;
}

/* ----------------------------------------------------------------------
   set pair and atom parameters to their absolute values for state k,
   taken from each perturbation's own grid
------------------------------------------------------------------------- */

void ComputeMBAR::perturb_params(int k)
{
  int i, j;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];
    const double value = pert->grid[k];

    if (pert->which == PAIR) {    // set the perturbed parameter to its state value
      // Each grid holds the absolute value the parameter takes at every state,
      // independent of the configuration's sampling state. The parameter may be
      // a soft-core activation (lambda) parameter, or any other pair parameter
      // such as epsilon or sigma; the user supplies the appropriate grid.
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo, i); j <= pert->jhi; j++)
          pert->array[i][j] = value;

    } else if (pert->which == ATOM) {

      if (pert->aparam == CHARGE) {    // set charges to their state value
        int *atype = atom->type;
        double *q = atom->q;
        int *mask = atom->mask;
        int natom = atom->nlocal + atom->nghost;

        for (i = 0; i < natom; i++)
          if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
            if (mask[i] & groupbit) q[i] = value;
      }
    }
  }

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  if (pairflag) force->pair->reinit();

  // reset KSpace charges if charges have changed

  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}

/* ----------------------------------------------------------------------
   backup pair parameters
------------------------------------------------------------------------- */

void ComputeMBAR::backup_params()
{
  int i, j;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];
    if (pert->which == PAIR) {
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo, i); j <= pert->jhi; j++) pert->array_orig[i][j] = pert->array[i][j];
    }
  }
}

/* ----------------------------------------------------------------------
   restore pair parameters to original values
------------------------------------------------------------------------- */

void ComputeMBAR::restore_params()
{
  int i, j;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];
    if (pert->which == PAIR) {
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo, i); j <= pert->jhi; j++) pert->array[i][j] = pert->array_orig[i][j];
    }
  }

  if (pairflag) force->pair->reinit();

  // reset KSpace charges if charges have changed

  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}

/* ----------------------------------------------------------------------
   manage storage for charge, force, energy, virial arrays
------------------------------------------------------------------------- */

void ComputeMBAR::allocate_storage()
{
  nmax = atom->nmax;
  memory->create(f_orig, nmax, 3, "mbar:f_orig");
  memory->create(peatom_orig, nmax, "mbar:peatom_orig");
  memory->create(pvatom_orig, nmax, 6, "mbar:pvatom_orig");
  if (chgflag) {
    memory->create(q_orig, nmax, "mbar:q_orig");
    if (force->kspace) {
      memory->create(keatom_orig, nmax, "mbar:keatom_orig");
      memory->create(kvatom_orig, nmax, 6, "mbar:kvatom_orig");
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeMBAR::deallocate_storage()
{
  memory->destroy(f_orig);
  memory->destroy(peatom_orig);
  memory->destroy(pvatom_orig);
  memory->destroy(q_orig);
  memory->destroy(keatom_orig);
  memory->destroy(kvatom_orig);

  f_orig = nullptr;
  q_orig = nullptr;
  peatom_orig = keatom_orig = nullptr;
  pvatom_orig = kvatom_orig = nullptr;
}

/* ----------------------------------------------------------------------
   backup and restore arrays with charge, force, energy, virial
------------------------------------------------------------------------- */

void ComputeMBAR::backup_qfev()
{
  int i;

  int nall = atom->nlocal + atom->nghost;
  int natom = atom->nlocal;
  if (force->newton || (force->kspace && force->kspace->tip4pflag)) natom += atom->nghost;

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f_orig[i][0] = f[i][0];
    f_orig[i][1] = f[i][1];
    f_orig[i][2] = f[i][2];
  }

  eng_vdwl_orig = force->pair->eng_vdwl;
  eng_coul_orig = force->pair->eng_coul;

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

  if (chgflag) {
    double *q = atom->q;
    for (i = 0; i < nall; i++) q_orig[i] = q[i];

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
}

/* ---------------------------------------------------------------------- */

void ComputeMBAR::restore_qfev()
{
  int i;

  int nall = atom->nlocal + atom->nghost;
  int natom = atom->nlocal;
  if (force->newton || (force->kspace && force->kspace->tip4pflag)) natom += atom->nghost;

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f[i][0] = f_orig[i][0];
    f[i][1] = f_orig[i][1];
    f[i][2] = f_orig[i][2];
  }

  force->pair->eng_vdwl = eng_vdwl_orig;
  force->pair->eng_coul = eng_coul_orig;

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

  if (chgflag) {
    double *q = atom->q;
    for (i = 0; i < nall; i++) q[i] = q_orig[i];

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
}
