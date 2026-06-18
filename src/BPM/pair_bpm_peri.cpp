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
   Contributing author: Claude Opus 4.8 (Anthropic), under the direction of
   Joel Clemmer (SNL) and Axel Kohlmeyer (Temple U).

   Short-range contact pair style for the BPM peridynamics model. It supplies
   the repulsive contact force between non-bonded near pairs (bonded pairs are
   excluded via the 1-2 special weight and handled by bond_style bpm/peri).
   Companion to bond_style bpm/peri. Derived from the contact term in the PERI
   package pair_peri_pmb.cpp (Mike Parks, SNL).
------------------------------------------------------------------------- */

#include "pair_bpm_peri.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "lattice.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBPMPeri::PairBPMPeri(LAMMPS *_lmp) : Pair(_lmp), kspring(nullptr), cut(nullptr), index_vfrac(-1)
{
  writedata = 1;
  single_enable = 1;
}

/* ---------------------------------------------------------------------- */

PairBPMPeri::~PairBPMPeri()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(kspring);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairBPMPeri::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double fxtmp, fytmp, fztmp;
  double r, rsq, dr, rk, vfrac_eff;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double *special_lj = force->special_lj;
  double *vfrac = atom->dvector[index_vfrac];

  const double lc = domain->lattice->xlattice;
  // contact onset distance. Legacy uses min(0.9*r0_ref, 1.35*lc); for non-bonded
  // pairs (reference separation exceeds the horizon) this is always 1.35*lc, so
  // the dropped reference positions (x0) are not needed here.
  const double d_ij = 1.35 * lc;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    // accumulate the force on atom i in registers across the inner loop and
    // write it back once, avoiding repeated indexed loads/stores of f[i][*]
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // skip bonded pairs (excluded via the 1-2 special weight) and 1-3/1-4
      if (special_lj[sbmask(j)] == 0.0) continue;
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;
      r = sqrt(rsq);
      if (r >= d_ij) continue;

      // short-range repulsion (the 15 constant is from the EMU Theory Manual).
      // The mean nodal volume keeps the contact Newton-third-law balanced; for
      // uniform nodal volume this matches legacy's vfrac[j] form.
      dr = r - d_ij;
      vfrac_eff = 0.5 * (vfrac[i] + vfrac[j]);
      rk = 15.0 * kspring[itype][jtype] * vfrac_eff * (dr / cut[itype][jtype]);
      fpair = (r > 0.0) ? -(rk / r) : 0.0;

      fxtmp += delx * fpair;
      fytmp += dely * fpair;
      fztmp += delz * fpair;
      if (newton_pair || (j < nlocal)) {
        f[j][0] -= delx * fpair;
        f[j][1] -= dely * fpair;
        f[j][2] -= delz * fpair;
      }

      if (eflag) evdwl = 0.5 * rk * dr;
      // the peridynamic stress integrates over both nodal volumes, so the virial
      // carries one more vfrac factor than the (one-vfrac) contact force/energy
      // (matches legacy PERI's fpair*vfrac[i]); consistent with bond_style bpm/peri
      if (evflag)
        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair * vfrac_eff, delx, dely, delz);
    }

    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBPMPeri::allocate()
{
  allocated = 1;
  const int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(kspring, np1, np1, "pair:kspring");
  memory->create(cut, np1, np1, "pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBPMPeri::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style bpm/peri command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   pair_coeff <i> <j> <contact stiffness> <horizon>
------------------------------------------------------------------------- */

void PairBPMPeri::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  double cut_one = utils::numeric(FLERR, arg[3], false, lmp);
  if (cut_one <= 0.0) error->all(FLERR, 3, "Invalid horizon value {} for pair style bpm/peri", cut_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      kspring[i][j] = k_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBPMPeri::init_style()
{
  // the contact onset distance needs a uniform cubic lattice spacing
  if (!domain->lattice)
    error->all(FLERR, Error::NOLASTLINE, "Pair style bpm/peri requires a lattice be defined");
  if ((domain->lattice->xlattice != domain->lattice->ylattice) ||
      (domain->lattice->xlattice != domain->lattice->zlattice))
    error->all(FLERR, Error::NOLASTLINE,
               "Pair style bpm/peri requires equal lattice spacing in x, y, and z");

  // per-atom nodal volume (vfrac), supplied by the user via fix property/atom
  int flag, cols;
  index_vfrac = atom->find_custom("vfrac", flag, cols);
  if ((index_vfrac < 0) || (flag != 1) || (cols != 0))
    error->all(FLERR, Error::NOLASTLINE,
               "Pair style bpm/peri requires a per-atom vfrac property; add "
               "'fix <ID> all property/atom d_vfrac ghost yes'");

  // handshake: reject a non-bpm/peri bond style. The check tolerates force->bond
  // being unset because this runs during create_bonds (before bond_style in the
  // standard BPM deck order); the requirement is enforced once a bond style is
  // attached.
  if (force->bond && !utils::strmatch(force->bond_style, "^bpm/peri"))
    error->all(FLERR, Error::NOLASTLINE, "Pair style bpm/peri requires bond style bpm/peri");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBPMPeri::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair_coeff for pair style bpm/peri must be set explicitly");

  kspring[j][i] = kspring[i][j];
  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBPMPeri::write_restart(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&kspring[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBPMPeri::read_restart(FILE *fp)
{
  allocate();

  int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &kspring[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&kspring[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ---------------------------------------------------------------------- */

double PairBPMPeri::single(int i, int j, int itype, int jtype, double rsq,
                           double /*factor_coul*/, double /*factor_lj*/, double &fforce)
{
  double *vfrac = atom->dvector[index_vfrac];
  const double d_ij = 1.35 * domain->lattice->xlattice;

  double r = sqrt(rsq);
  if (r >= d_ij) {
    fforce = 0.0;
    return 0.0;
  }

  double dr = r - d_ij;
  double vfrac_eff = 0.5 * (vfrac[i] + vfrac[j]);
  double rk = 15.0 * kspring[itype][jtype] * vfrac_eff * (dr / cut[itype][jtype]);
  fforce = (r > 0.0) ? -(rk / r) : 0.0;

  return 0.5 * rk * dr;
}
