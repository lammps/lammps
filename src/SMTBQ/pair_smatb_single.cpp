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
   Contributing author: Daniele Rapetti (iximiel@gmail.com)
------------------------------------------------------------------------- */

#include "pair_smatb_single.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSMATBSingle::PairSMATBSingle(LAMMPS *_lmp) :
    Pair(_lmp), nmax(0), on_eb(nullptr), r0(0), p(0), A(0), q(0), QSI(0), cutOffStart(0),
    cutOffEnd(0), cutOffEnd2(0), a3(0), a4(0), a5(0), x3(0), x4(0), x5(0)
{
  single_enable = 0;              // 1 if single() routine exists
  restartinfo = 1;                // 1 if pair style writes restart info
  respa_enable = 0;               // 1 if inner/middle/outer rRESPA routines
  one_coeff = 0;                  // 1 if allows only one coeff * * call
  manybody_flag = 1;              // 1 if a manybody potential
  no_virial_fdotr_compute = 0;    // 1 if does not invoke virial_fdotr_compute()
  writedata = 1;                  // 1 if writes coeffs to data file
  ghostneigh = 0;                 // 1 if pair style needs neighbors of ghosts

  // set comm size needed by this Pair
  comm_forward = 1;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

PairSMATBSingle::~PairSMATBSingle()
{
  if (copymode) { return; }
  memory->destroy(on_eb);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::compute(int eflag, int vflag)
{
  int i, j, ii, jj, jnum;
  double xtmp, ytmp, ztmp, del[3], fpair;
  double dijsq, dij;
  double espo, aexpp, qsiexpq, eb_i, Fb, Fr;
  double polyval, polyval2, polyval3, polyval4, polyval5;

  if (eflag || vflag) {
    ev_setup(eflag, vflag);
    eng_vdwl = 0;
  } else {
    evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
  }

  // grow on_eb array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(on_eb, nmax, "pair_smatb:on_eb");
  }

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int newton_pair = force->newton_pair;

  // zero out on_eb
  memset(on_eb, 0, nall * sizeof(double));

  int inum = list->inum;
  int *ilist = list->ilist;
  int *jlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // FIRST LOOP: CALCULATES the squared bonding energy and accumulate it in on_eb for each atom
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;
      del[0] = xtmp - x[j][0];
      del[1] = ytmp - x[j][1];
      del[2] = ztmp - x[j][2];
      dijsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];

      if (dijsq < cutOffEnd2) {
        dij = sqrt(dijsq);
        if (dij < cutOffStart) {
          qsiexpq = (QSI * QSI) * exp(2.0 * q * (1.0 - dij / r0));
        } else {
          polyval = dij - cutOffEnd;
          polyval3 = polyval * polyval * polyval;
          polyval4 = polyval3 * polyval;
          polyval5 = polyval4 * polyval;
          qsiexpq = x5 * polyval5 + x4 * polyval4 + x3 * polyval3;
          qsiexpq = qsiexpq * qsiexpq;
        }
        on_eb[i] += qsiexpq;
        if (newton_pair) on_eb[j] += qsiexpq;
      }
    }
  }

  // communicate the squared bonding energy between the various bins

  if (newton_pair) comm->reverse_comm(this);

  // Support Loop: take the square root of the bonding energy and
  // accumulate it in the energy accumulator if needed the store the
  // reciprocal in on_eb in order to not do it in the SECOND LOOP

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (i < nlocal) {
      eb_i = sqrt(on_eb[i]);
      if (eb_i != 0.0) {
        on_eb[i] = 1.0 / eb_i;
      } else {
        on_eb[i] = 0.0;
      }
      // if needed the bonding energy is accumulated:
      if (eflag_either) {
        if (eflag_atom) { eatom[i] -= eb_i; }
        if (eflag_global) { eng_vdwl -= eb_i; }
      }
    }
  }
  // this communication stores the denominators in the ghosts atoms,
  // this is needed because of how forces are calculated
  comm->forward_comm(this);

  // SECOND LOOP: given on_eb[i] calculates forces and energies
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      del[0] = xtmp - x[j][0];
      del[1] = ytmp - x[j][1];
      del[2] = ztmp - x[j][2];

      dijsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];
      if (dijsq < cutOffEnd2) {
        dij = sqrt(dijsq);
        if (dij < cutOffStart) {
          espo = 1.0 - dij / r0;
          aexpp = exp(p * espo) * A;
          Fr = (2.0 * aexpp) * (p / r0);
          qsiexpq = (QSI * QSI) * exp(2.0 * q * espo);
          Fb = -qsiexpq * q / r0;
        } else {
          polyval = dij - cutOffEnd;
          polyval2 = polyval * polyval;
          polyval3 = polyval2 * polyval;
          polyval4 = polyval3 * polyval;
          polyval5 = polyval4 * polyval;
          aexpp = a5 * polyval5 + a4 * polyval4 + a3 * polyval3;
          Fr = -2.0 * (5.0 * a5 * polyval4 + 4.0 * a4 * polyval3 + 3.0 * a3 * polyval2);
          qsiexpq = x5 * polyval5 + x4 * polyval4 + x3 * polyval3;
          Fb = ((5.0 * x5 * polyval4 + 4.0 * x4 * polyval3 + 3.0 * x3 * polyval2)) * qsiexpq;
        }

        // calculates the module of the pair energy between i and j
        fpair = (Fb * (on_eb[i] + on_eb[j]) + Fr) / dij;

        f[i][0] += del[0] * fpair;
        f[i][1] += del[1] * fpair;
        f[i][2] += del[2] * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= del[0] * fpair;
          f[j][1] -= del[1] * fpair;
          f[j][2] -= del[2] * fpair;
        }
        if (evflag) {
          ev_tally(i, j, nlocal, newton_pair, 2.0 * aexpp, 0.0, fpair, del[0], del[1], del[2]);
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSMATBSingle::settings(int narg, char **)
{
  if (narg > 0) error->all(FLERR, "Illegal pair_style command: smatb/single accepts no options");
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSMATBSingle::allocate()
{
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair_smatb:setflag");
  for (int i = 1; i < np1; i++) {
    for (int j = i; j < np1; j++) { setflag[i][j] = 0; }
  }

  memory->create(cutsq, np1, np1, "pair_smatb:cutsq");
  allocated = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSMATBSingle::coeff(int narg, char **arg)
{
  if (!allocated) { allocate(); }
  if (narg != 9) utils::missing_cmd_args(FLERR, "pair_style smatb/single", error);

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  r0 = utils::numeric(FLERR, arg[2], false, lmp);
  p = utils::numeric(FLERR, arg[3], false, lmp);
  q = utils::numeric(FLERR, arg[4], false, lmp);
  A = utils::numeric(FLERR, arg[5], false, lmp);
  QSI = utils::numeric(FLERR, arg[6], false, lmp);
  cutOffStart = utils::numeric(FLERR, arg[7], false, lmp);
  cutOffEnd = utils::numeric(FLERR, arg[8], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;

      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ------------------------------------------------------------------------ */

void PairSMATBSingle::init_style()
{
  if (force->newton_pair == 0) error->all(FLERR, "Pair style smatb/single requires newton pair on");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSMATBSingle::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  //calculating the polynomial linking to zero
  double es = cutOffEnd - cutOffStart;
  double es2 = es * es;
  double es3 = es2 * es;

  //variables for poly for p and A
  double expp = A * exp(p * (1.0 - cutOffStart / r0));
  double ap = -1.0 / es3;
  double bp = p / (r0 * es2);
  double cp = -(p * p) / (es * r0 * r0);

  a5 = expp * (12.0 * ap + 6.0 * bp + cp) / (2.0 * es2);
  a4 = expp * (15.0 * ap + 7.0 * bp + cp) / es;
  a3 = expp * (20.0 * ap + 8.0 * bp + cp) / 2.0;

  //variables for poly for q and qsi
  double expq = QSI * exp(q * (1.0 - cutOffStart / r0));
  double aq = -1 / es3;
  double bq = q / (es2 * r0);
  double cq = -(q * q) / (es * r0 * r0);

  x5 = expq * (12.0 * aq + 6.0 * bq + cq) / (2.0 * es2);
  x4 = expq * (15.0 * aq + 7.0 * bq + cq) / es;
  x3 = expq * (20.0 * aq + 8.0 * bq + cq) / 2.0;

  cutOffEnd2 = cutOffEnd * cutOffEnd;
  return cutOffEnd;
}

/* ---------------------------------------------------------------------- */

int PairSMATBSingle::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                       int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; ++i) {
    j = list[i];
    buf[m++] = on_eb[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; ++i) { on_eb[i] = buf[m++]; }
}

/* ---------------------------------------------------------------------- */

int PairSMATBSingle::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; ++i) { buf[m++] = on_eb[i]; }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    on_eb[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

//write binary data of this simulation:
void PairSMATBSingle::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&r0, sizeof(double), 1, fp);
        fwrite(&p, sizeof(double), 1, fp);
        fwrite(&q, sizeof(double), 1, fp);
        fwrite(&A, sizeof(double), 1, fp);
        fwrite(&QSI, sizeof(double), 1, fp);
        fwrite(&cutOffStart, sizeof(double), 1, fp);
        fwrite(&cutOffEnd, sizeof(double), 1, fp);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &r0, sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &p, sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &q, sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &A, sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &QSI, sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cutOffStart, sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cutOffEnd, sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&r0, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&p, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&q, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&A, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&QSI, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cutOffStart, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cutOffEnd, 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    fprintf(fp, "%d %g %g %g %g %g %g %g\n", i, r0, p, q, A, QSI, cutOffStart, cutOffEnd);
  }
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingle::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      fprintf(fp, "%d %d %g %g %g %g %g %g %g\n", i, j, r0, p, q, A, QSI, cutOffStart, cutOffEnd);
    }
  }
}
