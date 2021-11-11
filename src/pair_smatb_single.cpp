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
This pair style is written by Daniele Rapetti (iximiel@gmail.com)
------------------------------------------------------------------------- */
#include "pair_smatb_single.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;

#define MAXLINE 1024

PairSMATBSingle::PairSMATBSingle(LAMMPS *lmp) :
    Pair(lmp), nmax(0), on_eb(nullptr), r0(0), p(0), A(0), q(0), QSI(0), cutOffStart(0),
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

PairSMATBSingle::~PairSMATBSingle()
{
  if (copymode) { return; }
  memory->destroy(on_eb);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

void PairSMATBSingle::compute(int eflag, int vflag)
{    //workhorse routine that computes pairwise interactions
  //eflag means compute energy
  //vflag means compute virial
  int i, j, ii, jj, jnum;    //,itype,jtype;
  double xtmp, ytmp, ztmp, del[3], fpair;
  double dijsq, dij;
  double espo, aexpp, qsiexpq, eb_i, Fb, Fr;
  double polyval, polyval2, polyval3, polyval4, polyval5;
  //sets up the flags for energy caclulations
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
    eng_vdwl = 0;
  } else {
    evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
  }

  // grow on_eb array if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(on_eb, nmax, "pair_smatb:on_eb");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int newton_pair = force->newton_pair;

  // zero out on_eb
  if (newton_pair) {
    memset(on_eb, 0, nall * sizeof(on_eb[0]));
  } else {
    memset(on_eb, 0, nlocal * sizeof(on_eb[0]));
  }

  int inum = list->inum;
  int *ilist = list->ilist;
  int *jlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  //FIRST LOOP: CALCULATES the squared bounding energy and accumulate it in on_eb for each atom
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    //itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;
      //jtype = type[j];
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
        //if (newton_pair || j < nlocal) {
        on_eb[j] += qsiexpq;
        //}
      }
    }
  }

  //communicate the squared bounding energy between the various bins

  comm->reverse_comm_pair(this);

  //Support Loop: take the square root of the bounding energy and accumulate it in the energy accumulator if needed
  // the store the reciprocal in on_eb in order to not do it in the SECOND LOOP

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (i < nlocal) {
      eb_i = sqrt(on_eb[i]);
      if (eb_i != 0.0) {
        on_eb[i] = 1.0 / eb_i;
      } else {
        on_eb[i] = 0.0;
      }
      //if needed the bounding energy is accumulated:
      if (eflag_either) {
        if (eflag_atom) { eatom[i] -= eb_i; }
        if (eflag_global) { eng_vdwl -= eb_i; }
      }
    }
  }
  //this communication stores the denominators in the ghosts atoms, this is needed because of how forces are calculated
  comm->forward_comm_pair(this);

  //SECOND LOOP: given on_eb[i] calculates forces and energies
  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    //itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      //jtype = type[j];

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
        //if needed the repulsive energy is accumulated:
        if (eflag_either) {
          if (eflag_atom) {
            eatom[i] += aexpp;
            if (newton_pair || j < nlocal) { eatom[j] += aexpp; }
          }
          if (eflag_global) {
            if (newton_pair || j < nlocal) {
              eng_vdwl += 2.0 * (aexpp);
            } else {
              eng_vdwl += aexpp;
            }
          }
        }
        //calculates the module of the pair energy between i and j
        fpair = (Fb * (on_eb[i] + on_eb[j]) + Fr) / dij;

        f[i][0] += del[0] * fpair;
        f[i][1] += del[1] * fpair;
        f[i][2] += del[2] * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= del[0] * fpair;
          f[j][1] -= del[1] * fpair;
          f[j][2] -= del[2] * fpair;
        }
        if (vflag_atom) {
          ev_tally(i, j, nlocal, newton_pair, 0.0,
                   0.0,    //Energy is tally'd in the other parts of the potential
                   fpair, del[0], del[1], del[2]);
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
{    //reads the input script line with arguments you define
  if (narg > 0) error->all(FLERR, "Illegal pair_style command: smatb accepts no options");
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSMATBSingle::allocate()
{
  int n = atom->ntypes;
  int natoms = atom->natoms;

  memory->create(setflag, n + 1, n + 1, "pair_smatb:setflag");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) { setflag[i][j] = 0; }
  }

  memory->create(cutsq, n + 1, n + 1, "pair_smatb:cutsq");

  allocated = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSMATBSingle::coeff(int narg, char **arg)
{    //set coefficients for one i,j type pair
  if (!allocated) { allocate(); }
  if (narg != 9) {
    error->all(
        FLERR,
        "Incorrect args for pair coefficients:\n SMATB needs \"i j r0 p q A QSI CO_start CO_end\"");
  }
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  //reading parameters from input
  double myr0 = utils::numeric(FLERR, arg[2], false, lmp),
         myp = utils::numeric(FLERR, arg[3], false, lmp),
         myq = utils::numeric(FLERR, arg[4], false, lmp),
         myA = utils::numeric(FLERR, arg[5], false, lmp),
         myQSI = utils::numeric(FLERR, arg[6], false, lmp),
         mycutOffStart = utils::numeric(FLERR, arg[7], false, lmp),
         mycutOffEnd = utils::numeric(FLERR, arg[8], false, lmp);
  int count = 0;

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      r0 = myr0;
      p = myp;
      A = myA;
      q = myq;
      QSI = myQSI;
      cutOffStart = mycutOffStart;
      cutOffEnd = mycutOffEnd;
      setflag[i][j] = 1;

      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */
/*
  void PairSMATBSingle::init_style() {//initialization specific to this pair style
  neighbor->request(this,instance_me);
  }
*/

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSMATBSingle::init_one(int i, int j)
{    //perform initialization for one i,j type pair
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  //calculating the polynomial linking to zero
  double es = cutOffEnd - cutOffStart;
  double es2 = es * es;
  double es3 = es2 * es;

  //variables for poly for p and A
  double expp = A * exp(p * (1. - cutOffStart / r0));
  double ap = -1. / es3;
  double bp = p / (r0 * es2);
  double cp = -(p * p) / (es * r0 * r0);

  a5 = expp * (12. * ap + 6. * bp + cp) / (2. * es2);
  a4 = expp * (15. * ap + 7. * bp + cp) / es;
  a3 = expp * (20. * ap + 8. * bp + cp) / 2.;

  //variables for poly for q and qsi
  double expq = QSI * exp(q * (1. - cutOffStart / r0));
  double aq = -1 / es3;
  double bq = q / (es2 * r0);
  double cq = -(q * q) / (es * r0 * r0);

  x5 = expq * (12. * aq + 6. * bq + cq) / (2. * es2);
  x4 = expq * (15. * aq + 7. * bq + cq) / es;
  x3 = expq * (20. * aq + 8. * bq + cq) / 2.;

  cutOffEnd2 = cutOffEnd * cutOffEnd;
  if (i != j) {
    setflag[j][i] = 1;
    cutOffEnd2 = cutOffEnd2;

    r0 = r0;
    p = p;
    q = q;
    A = A;
    QSI = QSI;
    cutOffStart = cutOffStart;
    cutOffEnd = cutOffEnd;

    a3 = a3;
    a4 = a4;
    a5 = a5;
    x3 = x3;
    x4 = x4;
    x5 = x5;
  }
  //returns cutOffEnd to calculate cutforce and cutsq
  return cutOffEnd;
}

/* ---------------------------------------------------------------------- */

int PairSMATBSingle::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
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

//write binary data of this simulation:
void PairSMATBSingle::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

void PairSMATBSingle::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  size_t result;
  if (me == 0) {
    result = fread(&offset_flag, sizeof(int), 1, fp);
    result = fread(&mix_flag, sizeof(int), 1, fp);
    result = fread(&tail_flag, sizeof(int), 1, fp);
  }
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

void PairSMATBSingle::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  //"I J r0 p q A QSI CO_start CO_end"
  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
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
  //maybe we need to save also the values of the various polynomials
}

void PairSMATBSingle::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  size_t result;

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) { result = fread(&setflag[i][j], sizeof(int), 1, fp); }
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          result = fread(&r0, sizeof(double), 1, fp);
          result = fread(&p, sizeof(double), 1, fp);
          result = fread(&q, sizeof(double), 1, fp);
          result = fread(&A, sizeof(double), 1, fp);
          result = fread(&QSI, sizeof(double), 1, fp);
          result = fread(&cutOffStart, sizeof(double), 1, fp);
          result = fread(&cutOffEnd, sizeof(double), 1, fp);
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

void PairSMATBSingle::write_data(FILE *fp)
{
  //smatb needs I J r0 p q A QSI CO_start CO_end
  for (int i = 1; i <= atom->ntypes; i++) {
    fprintf(fp, "%d %g %g %g %g %g %g %g\n", i, r0, p, q, A, QSI, cutOffStart, cutOffEnd);
  }
}

void PairSMATBSingle::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      fprintf(fp, "%d %d %g %g %g %g %g %g %g\n", i, j, r0, p, q, A, QSI, cutOffStart, cutOffEnd);
    }
  }
}
