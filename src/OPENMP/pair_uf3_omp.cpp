/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_uf3_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>
#include <vector>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using MathConst::THIRD;

/* ---------------------------------------------------------------------- */

PairUF3OMP::PairUF3OMP(LAMMPS *lmp) : PairUF3(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairUF3OMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, cvatom, thr);

    // the base style requires newton pair on

    if (evflag) {
      if (eflag) {
        eval<1, 1>(ifrom, ito, thr);
      } else {
        eval<1, 0>(ifrom, ito, thr);
      }
    } else {
      eval<0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ----------------------------------------------------------------------
   same as PairUF3::compute() for the neighbor list entries iifrom to
   iito-1, with a thread-local short neighbor list for the 3-body terms
------------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairUF3OMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, k, ii, jj, kk, jnum, itype, jtype, ktype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fx, fy, fz;
  double del_rji[3], del_rki[3], del_rkj[3];
  double fij[3], fik[3], fjk[3];
  double fji[3], fki[3], fkj[3];
  double Fi[3], Fj[3], Fk[3];
  double rsq, rij, rik, rjk;
  double rij_sq, rik_sq, rjk_sq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  const double *const *const x = atom->x;
  double *const *const f = thr->get_f();
  double *const *const cvatom_thr = thr->get_cvatom_pair();
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  std::vector<int> neighshort;

  // loop over neighbors of my atoms
  for (ii = iifrom; ii < iito; ii++) {
    evdwl = 0;
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;
    if (jnum > (int) neighshort.size()) neighshort.resize(jnum);
    for (jj = 0; jj < jnum; jj++) {
      fx = 0;
      fy = 0;
      fz = 0;
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (rsq < cutsq[itype][jtype]) {
        rij = sqrt(rsq);

        if (pot_3b) {
          if (rij <= cut_3b_list[itype][jtype]) {
            neighshort[numshort] = j;
            numshort = numshort + 1;
          }
        }

        int knot_start_index = (this->*get_starting_index_2b)(itype, jtype, rij);

        double force_2b = cached_constants_2b_deri[itype][jtype][knot_start_index - 1][0];
        force_2b += rij * cached_constants_2b_deri[itype][jtype][knot_start_index - 1][1];
        force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 1][2];
        force_2b += cached_constants_2b_deri[itype][jtype][knot_start_index - 2][3];
        force_2b += rij * cached_constants_2b_deri[itype][jtype][knot_start_index - 2][4];
        force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 2][5];
        force_2b += cached_constants_2b_deri[itype][jtype][knot_start_index - 3][6];
        force_2b += rij * cached_constants_2b_deri[itype][jtype][knot_start_index - 3][7];
        force_2b += rsq * cached_constants_2b_deri[itype][jtype][knot_start_index - 3][8];

        fpair = -1 * force_2b / rij;

        fx = delx * fpair;
        fy = dely * fpair;
        fz = delz * fpair;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;

        if (EFLAG) {
          double rth = rsq * rij;
          evdwl = cached_constants_2b[itype][jtype][knot_start_index][0];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index][1];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index][2];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index][3];
          evdwl += cached_constants_2b[itype][jtype][knot_start_index - 1][4];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index - 1][5];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index - 1][6];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index - 1][7];
          evdwl += cached_constants_2b[itype][jtype][knot_start_index - 2][8];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index - 2][9];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index - 2][10];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index - 2][11];
          evdwl += cached_constants_2b[itype][jtype][knot_start_index - 3][12];
          evdwl += rij * cached_constants_2b[itype][jtype][knot_start_index - 3][13];
          evdwl += rsq * cached_constants_2b[itype][jtype][knot_start_index - 3][14];
          evdwl += rth * cached_constants_2b[itype][jtype][knot_start_index - 3][15];
        };

        if (EVFLAG) {
          ev_tally_xyz_thr(this, i, j, nlocal, /* newton_pair */ 1, evdwl, 0.0, fx, fy, fz, delx,
                           dely, delz, thr);

          // Centroid Stress
          if (vflag_either && cvflag_atom) {
            double v[6];

            v[0] = delx * fx;
            v[1] = dely * fy;
            v[2] = delz * fz;
            v[3] = delx * fy;
            v[4] = delx * fz;
            v[5] = dely * fz;

            cvatom_thr[i][0] += 0.5 * v[0];
            cvatom_thr[i][1] += 0.5 * v[1];
            cvatom_thr[i][2] += 0.5 * v[2];
            cvatom_thr[i][3] += 0.5 * v[3];
            cvatom_thr[i][4] += 0.5 * v[4];
            cvatom_thr[i][5] += 0.5 * v[5];
            cvatom_thr[i][6] += 0.5 * v[3];
            cvatom_thr[i][7] += 0.5 * v[4];
            cvatom_thr[i][8] += 0.5 * v[5];

            cvatom_thr[j][0] += 0.5 * v[0];
            cvatom_thr[j][1] += 0.5 * v[1];
            cvatom_thr[j][2] += 0.5 * v[2];
            cvatom_thr[j][3] += 0.5 * v[3];
            cvatom_thr[j][4] += 0.5 * v[4];
            cvatom_thr[j][5] += 0.5 * v[5];
            cvatom_thr[j][6] += 0.5 * v[3];
            cvatom_thr[j][7] += 0.5 * v[4];
            cvatom_thr[j][8] += 0.5 * v[5];
          }
        }
      }
    }

    // 3-body interaction
    // jth atom
    jnum = numshort - 1;
    for (jj = 0; jj < jnum; jj++) {
      fij[0] = fji[0] = 0;
      fij[1] = fji[1] = 0;
      fij[2] = fji[2] = 0;
      j = neighshort[jj];
      jtype = type[j];
      del_rji[0] = x[j][0] - xtmp;
      del_rji[1] = x[j][1] - ytmp;
      del_rji[2] = x[j][2] - ztmp;
      rij_sq = (del_rji[0] * del_rji[0]) + (del_rji[1] * del_rji[1]) + (del_rji[2] * del_rji[2]);
      rij = sqrt(rij_sq);

      // kth atom
      for (kk = jj + 1; kk < numshort; kk++) {

        fik[0] = fki[0] = 0;
        fik[1] = fki[1] = 0;
        fik[2] = fki[2] = 0;

        fjk[0] = fkj[0] = 0;
        fjk[1] = fkj[1] = 0;
        fjk[2] = fkj[2] = 0;

        k = neighshort[kk];
        ktype = type[k];
        del_rki[0] = x[k][0] - xtmp;
        del_rki[1] = x[k][1] - ytmp;
        del_rki[2] = x[k][2] - ztmp;
        rik_sq = (del_rki[0] * del_rki[0]) + (del_rki[1] * del_rki[1]) + (del_rki[2] * del_rki[2]);
        rik = sqrt(rik_sq);

        if ((rij <= cut_3b[itype][jtype][ktype]) && (rik <= cut_3b[itype][ktype][jtype]) &&
            (rij >= min_cut_3b[itype][jtype][ktype][2]) &&
            (rik >= min_cut_3b[itype][jtype][ktype][1])) {

          del_rkj[0] = x[k][0] - x[j][0];
          del_rkj[1] = x[k][1] - x[j][1];
          del_rkj[2] = x[k][2] - x[j][2];

          rjk_sq =
              (del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2]);
          rjk = sqrt(rjk_sq);

          if (rjk >= min_cut_3b[itype][jtype][ktype][0]) {
            double rij_th = rij * rij_sq;
            double rik_th = rik * rik_sq;
            double rjk_th = rjk * rjk_sq;

            int map_to = map_3b[itype][jtype][ktype];
            int knot_start_index_ij = (this->*get_starting_index_3b)(itype, jtype, ktype, rij, 2);
            int knot_start_index_ik = (this->*get_starting_index_3b)(itype, jtype, ktype, rik, 1);
            int knot_start_index_jk = (this->*get_starting_index_3b)(itype, jtype, ktype, rjk, 0);
            double basis_ij[4];
            double basis_ik[4];
            double basis_jk[4];
            double basis_ij_der[3];
            double basis_ik_der[3];
            double basis_jk_der[3];

            //--------------basis_ij
            basis_ij[0] = cached_constants_3b[map_to][0][knot_start_index_ij - 3][12];
            basis_ij[0] += rij * cached_constants_3b[map_to][0][knot_start_index_ij - 3][13];
            basis_ij[0] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij - 3][14];
            basis_ij[0] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij - 3][15];

            basis_ij[1] = cached_constants_3b[map_to][0][knot_start_index_ij - 2][8];
            basis_ij[1] += rij * cached_constants_3b[map_to][0][knot_start_index_ij - 2][9];
            basis_ij[1] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij - 2][10];
            basis_ij[1] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij - 2][11];

            basis_ij[2] = cached_constants_3b[map_to][0][knot_start_index_ij - 1][4];
            basis_ij[2] += rij * cached_constants_3b[map_to][0][knot_start_index_ij - 1][5];
            basis_ij[2] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij - 1][6];
            basis_ij[2] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij - 1][7];

            basis_ij[3] = cached_constants_3b[map_to][0][knot_start_index_ij][0];
            basis_ij[3] += rij * cached_constants_3b[map_to][0][knot_start_index_ij][1];
            basis_ij[3] += rij_sq * cached_constants_3b[map_to][0][knot_start_index_ij][2];
            basis_ij[3] += rij_th * cached_constants_3b[map_to][0][knot_start_index_ij][3];

            //--------------basis_ik
            basis_ik[0] = cached_constants_3b[map_to][1][knot_start_index_ik - 3][12];
            basis_ik[0] += rik * cached_constants_3b[map_to][1][knot_start_index_ik - 3][13];
            basis_ik[0] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik - 3][14];
            basis_ik[0] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik - 3][15];

            basis_ik[1] = cached_constants_3b[map_to][1][knot_start_index_ik - 2][8];
            basis_ik[1] += rik * cached_constants_3b[map_to][1][knot_start_index_ik - 2][9];
            basis_ik[1] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik - 2][10];
            basis_ik[1] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik - 2][11];

            basis_ik[2] = cached_constants_3b[map_to][1][knot_start_index_ik - 1][4];
            basis_ik[2] += rik * cached_constants_3b[map_to][1][knot_start_index_ik - 1][5];
            basis_ik[2] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik - 1][6];
            basis_ik[2] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik - 1][7];

            basis_ik[3] = cached_constants_3b[map_to][1][knot_start_index_ik][0];
            basis_ik[3] += rik * cached_constants_3b[map_to][1][knot_start_index_ik][1];
            basis_ik[3] += rik_sq * cached_constants_3b[map_to][1][knot_start_index_ik][2];
            basis_ik[3] += rik_th * cached_constants_3b[map_to][1][knot_start_index_ik][3];

            //--------------basis_jk
            basis_jk[0] = cached_constants_3b[map_to][2][knot_start_index_jk - 3][12];
            basis_jk[0] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk - 3][13];
            basis_jk[0] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk - 3][14];
            basis_jk[0] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk - 3][15];

            basis_jk[1] = cached_constants_3b[map_to][2][knot_start_index_jk - 2][8];
            basis_jk[1] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk - 2][9];
            basis_jk[1] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk - 2][10];
            basis_jk[1] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk - 2][11];

            basis_jk[2] = cached_constants_3b[map_to][2][knot_start_index_jk - 1][4];
            basis_jk[2] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk - 1][5];
            basis_jk[2] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk - 1][6];
            basis_jk[2] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk - 1][7];

            basis_jk[3] = cached_constants_3b[map_to][2][knot_start_index_jk][0];
            basis_jk[3] += rjk * cached_constants_3b[map_to][2][knot_start_index_jk][1];
            basis_jk[3] += rjk_sq * cached_constants_3b[map_to][2][knot_start_index_jk][2];
            basis_jk[3] += rjk_th * cached_constants_3b[map_to][2][knot_start_index_jk][3];

            //----------------basis_ij_der
            basis_ij_der[0] = cached_constants_3b_deri[map_to][0][knot_start_index_ij - 3][6];
            basis_ij_der[0] +=
                rij * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 3][7];
            basis_ij_der[0] +=
                rij_sq * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 3][8];

            basis_ij_der[1] = cached_constants_3b_deri[map_to][0][knot_start_index_ij - 2][3];
            basis_ij_der[1] +=
                rij * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 2][4];
            basis_ij_der[1] +=
                rij_sq * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 2][5];

            basis_ij_der[2] = cached_constants_3b_deri[map_to][0][knot_start_index_ij - 1][0];
            basis_ij_der[2] +=
                rij * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 1][1];
            basis_ij_der[2] +=
                rij_sq * cached_constants_3b_deri[map_to][0][knot_start_index_ij - 1][2];

            //----------------basis_ik_der
            basis_ik_der[0] = cached_constants_3b_deri[map_to][1][knot_start_index_ik - 3][6];
            basis_ik_der[0] +=
                rik * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 3][7];
            basis_ik_der[0] +=
                rik_sq * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 3][8];

            basis_ik_der[1] = cached_constants_3b_deri[map_to][1][knot_start_index_ik - 2][3];
            basis_ik_der[1] +=
                rik * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 2][4];
            basis_ik_der[1] +=
                rik_sq * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 2][5];

            basis_ik_der[2] = cached_constants_3b_deri[map_to][1][knot_start_index_ik - 1][0];
            basis_ik_der[2] +=
                rik * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 1][1];
            basis_ik_der[2] +=
                rik_sq * cached_constants_3b_deri[map_to][1][knot_start_index_ik - 1][2];

            //----------------basis_jk_der
            basis_jk_der[0] = cached_constants_3b_deri[map_to][2][knot_start_index_jk - 3][6];
            basis_jk_der[0] +=
                rjk * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 3][7];
            basis_jk_der[0] +=
                rjk_sq * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 3][8];

            basis_jk_der[1] = cached_constants_3b_deri[map_to][2][knot_start_index_jk - 2][3];
            basis_jk_der[1] +=
                rjk * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 2][4];
            basis_jk_der[1] +=
                rjk_sq * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 2][5];

            basis_jk_der[2] = cached_constants_3b_deri[map_to][2][knot_start_index_jk - 1][0];
            basis_jk_der[2] +=
                rjk * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 1][1];
            basis_jk_der[2] +=
                rjk_sq * cached_constants_3b_deri[map_to][2][knot_start_index_jk - 1][2];

            double triangle_eval[4] = {0, 0, 0, 0};

            int iknot_ij = knot_start_index_ij - 3;
            int iknot_ik = knot_start_index_ik - 3;
            int iknot_jk = knot_start_index_jk - 3;

            for (int l = 0; l < 3; l++) {
              const double basis_ij_der_i = basis_ij_der[l];
              for (int m = 0; m < 4; m++) {
                const double factor = basis_ij_der_i * basis_ik[m];
                const double *slice =
                    &coeff_for_der_ij[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                double tmp[4];
                tmp[0] = slice[0] * basis_jk[0];
                tmp[1] = slice[1] * basis_jk[1];
                tmp[2] = slice[2] * basis_jk[2];
                tmp[3] = slice[3] * basis_jk[3];
                double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
                triangle_eval[1] += factor * sum;
              }
            }

            for (int l = 0; l < 4; l++) {
              const double basis_ij_i = basis_ij[l];
              for (int m = 0; m < 3; m++) {
                const double factor = basis_ij_i * basis_ik_der[m];
                const double *slice =
                    &coeff_for_der_ik[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                double tmp[4];
                tmp[0] = slice[0] * basis_jk[0];
                tmp[1] = slice[1] * basis_jk[1];
                tmp[2] = slice[2] * basis_jk[2];
                tmp[3] = slice[3] * basis_jk[3];
                double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
                triangle_eval[2] += factor * sum;
              }
            }

            for (int l = 0; l < 4; l++) {
              const double basis_ij_i = basis_ij[l];
              for (int m = 0; m < 4; m++) {
                const double factor = basis_ij_i * basis_ik[m];
                const double *slice =
                    &coeff_for_der_jk[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                double tmp[3];
                tmp[0] = slice[0] * basis_jk_der[0];
                tmp[1] = slice[1] * basis_jk_der[1];
                tmp[2] = slice[2] * basis_jk_der[2];
                double sum = tmp[0] + tmp[1] + tmp[2];
                triangle_eval[3] += factor * sum;
              }
            }

            fij[0] = *(triangle_eval + 1) * (del_rji[0] / rij);
            fji[0] = -fij[0];
            fik[0] = *(triangle_eval + 2) * (del_rki[0] / rik);
            fki[0] = -fik[0];
            fjk[0] = *(triangle_eval + 3) * (del_rkj[0] / rjk);
            fkj[0] = -fjk[0];

            fij[1] = *(triangle_eval + 1) * (del_rji[1] / rij);
            fji[1] = -fij[1];
            fik[1] = *(triangle_eval + 2) * (del_rki[1] / rik);
            fki[1] = -fik[1];
            fjk[1] = *(triangle_eval + 3) * (del_rkj[1] / rjk);
            fkj[1] = -fjk[1];

            fij[2] = *(triangle_eval + 1) * (del_rji[2] / rij);
            fji[2] = -fij[2];
            fik[2] = *(triangle_eval + 2) * (del_rki[2] / rik);
            fki[2] = -fik[2];
            fjk[2] = *(triangle_eval + 3) * (del_rkj[2] / rjk);
            fkj[2] = -fjk[2];

            Fi[0] = fij[0] + fik[0];
            Fi[1] = fij[1] + fik[1];
            Fi[2] = fij[2] + fik[2];
            f[i][0] += Fi[0];
            f[i][1] += Fi[1];
            f[i][2] += Fi[2];

            Fj[0] = fji[0] + fjk[0];
            Fj[1] = fji[1] + fjk[1];
            Fj[2] = fji[2] + fjk[2];
            f[j][0] += Fj[0];
            f[j][1] += Fj[1];
            f[j][2] += Fj[2];

            Fk[0] = fki[0] + fkj[0];
            Fk[1] = fki[1] + fkj[1];
            Fk[2] = fki[2] + fkj[2];
            f[k][0] += Fk[0];
            f[k][1] += Fk[1];
            f[k][2] += Fk[2];

            if (EFLAG) {
              for (int l = 0; l < 4; l++) {
                const double basis_ij_i = basis_ij[l];
                for (int m = 0; m < 4; m++) {
                  const double factor = basis_ij_i * basis_ik[m];
                  const double *slice =
                      &n3b_coeff_array[map_to][iknot_ij + l][iknot_ik + m][iknot_jk];
                  double tmp[4];
                  tmp[0] = slice[0] * basis_jk[0];
                  tmp[1] = slice[1] * basis_jk[1];
                  tmp[2] = slice[2] * basis_jk[2];
                  tmp[3] = slice[3] * basis_jk[3];
                  double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
                  triangle_eval[0] += factor * sum;
                }
              }
              evdwl = *triangle_eval;
            }

            if (EVFLAG) {
              ev_tally3_thr(this, i, j, k, evdwl, 0, Fj, Fk, del_rji, del_rki, thr);
              // Centroid stress 3-body term
              if (vflag_either && cvflag_atom) {
                double ric[3];
                ric[0] = THIRD * (-del_rji[0] - del_rki[0]);
                ric[1] = THIRD * (-del_rji[1] - del_rki[1]);
                ric[2] = THIRD * (-del_rji[2] - del_rki[2]);

                cvatom_thr[i][0] += ric[0] * Fi[0];
                cvatom_thr[i][1] += ric[1] * Fi[1];
                cvatom_thr[i][2] += ric[2] * Fi[2];
                cvatom_thr[i][3] += ric[0] * Fi[1];
                cvatom_thr[i][4] += ric[0] * Fi[2];
                cvatom_thr[i][5] += ric[1] * Fi[2];
                cvatom_thr[i][6] += ric[1] * Fi[0];
                cvatom_thr[i][7] += ric[2] * Fi[0];
                cvatom_thr[i][8] += ric[2] * Fi[1];

                double rjc[3];
                rjc[0] = THIRD * (del_rji[0] - del_rkj[0]);
                rjc[1] = THIRD * (del_rji[1] - del_rkj[1]);
                rjc[2] = THIRD * (del_rji[2] - del_rkj[2]);

                cvatom_thr[j][0] += rjc[0] * Fj[0];
                cvatom_thr[j][1] += rjc[1] * Fj[1];
                cvatom_thr[j][2] += rjc[2] * Fj[2];
                cvatom_thr[j][3] += rjc[0] * Fj[1];
                cvatom_thr[j][4] += rjc[0] * Fj[2];
                cvatom_thr[j][5] += rjc[1] * Fj[2];
                cvatom_thr[j][6] += rjc[1] * Fj[0];
                cvatom_thr[j][7] += rjc[2] * Fj[0];
                cvatom_thr[j][8] += rjc[2] * Fj[1];

                double rkc[3];
                rkc[0] = THIRD * (del_rki[0] + del_rkj[0]);
                rkc[1] = THIRD * (del_rki[1] + del_rkj[1]);
                rkc[2] = THIRD * (del_rki[2] + del_rkj[2]);

                cvatom_thr[k][0] += rkc[0] * Fk[0];
                cvatom_thr[k][1] += rkc[1] * Fk[1];
                cvatom_thr[k][2] += rkc[2] * Fk[2];
                cvatom_thr[k][3] += rkc[0] * Fk[1];
                cvatom_thr[k][4] += rkc[0] * Fk[2];
                cvatom_thr[k][5] += rkc[1] * Fk[2];
                cvatom_thr[k][6] += rkc[1] * Fk[0];
                cvatom_thr[k][7] += rkc[2] * Fk[0];
                cvatom_thr[k][8] += rkc[2] * Fk[1];
              }
            }
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairUF3OMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairUF3::memory_usage();

  return bytes;
}
