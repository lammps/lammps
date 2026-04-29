/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_e3b_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;

#define DIM 3
#define NUMH 2
#define NUMO 2

/* ---------------------------------------------------------------------- */

PairE3BOMP::PairE3BOMP(LAMMPS *lmp) : PairE3B(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairE3BOMP::compute(int eflag, int vflag)
{
  if (natoms != atom->natoms) error->all(FLERR, "pair E3B requires a fixed number of atoms");

  ev_init(eflag, vflag);
  memset(sumExp, 0, sizeof(double) * maxID);

  pvector[0] = pvector[1] = pvector[2] = pvector[3] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int npair = 0;

  // SERIAL: first pass - build pair lists, accumulate sumExp, compute 2-body

  double fixx, fiyy, fizz, fxtmp, fytmp, fztmp;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair, rsq, tmpexp;
  double delxh, delyh, delzh, rsqh, tmpr;
  double scFact1, scFact2, scEng, scDer;
  bool addedH;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (type[i] != typeO) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fixx = fiyy = fizz = 0.0;

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      if (type[j] != typeO) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < rc2sq) {
        tmpr = sqrt(rsq);
        tmpexp = e2 * exp(-k2 * tmpr);
        fpair = k2 * tmpexp / tmpr;

        fxtmp = delx * fpair;
        fytmp = dely * fpair;
        fztmp = delz * fpair;
        fixx += fxtmp;
        fiyy += fytmp;
        fizz += fztmp;
        f[j][0] -= fxtmp;
        f[j][1] -= fytmp;
        f[j][2] -= fztmp;

        if (evflag) {
          ev_tally(i, j, nlocal, newton_pair, tmpexp, 0.0, fpair, delx, dely, delz);
          pvector[0] += tmpexp;
        }
      }

      if (rsq < rc3deltaSq) {
        pairO[npair][0] = i;
        pairO[npair][1] = j;
        addedH = false;

        for (int kk = 0; kk < NUMO; kk++) {
          int k = pairO[npair][kk];
          int otherO = pairO[npair][(kk + 1) % 2];
          for (int hh = 0; hh < NUMH; hh++) {
            int h = atom->map(tag[otherO] + hh + 1);
            h = domain->closest_image(otherO, h);
            pairH[npair][kk][hh] = h;

            delxh = x[k][0] - x[h][0];
            delyh = x[k][1] - x[h][1];
            delzh = x[k][2] - x[h][2];
            rsqh = delxh * delxh + delyh * delyh + delzh * delzh;

            if (rsqh < rc3sq) {
              tmpr = sqrt(rsqh);
              tmpexp = exp(-k3 * tmpr);
              if (tmpr > rs) {
                scFact1 = rc3 - tmpr;
                scFact2 = sc_num + 2 * tmpr;
                scEng = scFact1 * scFact1 * scFact2 * sc_denom;
                scDer = k3 * scEng - 6 * scFact1 * (rs - tmpr) * sc_denom;
              } else {
                scDer = k3;
                scEng = 1.0;
              }

              fpair3[npair][kk][hh] = scDer * tmpexp / tmpr;
              tmpexp *= scEng;
              exps[npair][kk][hh] = tmpexp;
              del3[npair][kk][hh][0] = delxh;
              del3[npair][kk][hh][1] = delyh;
              del3[npair][kk][hh][2] = delzh;

              sumExp[tag[k] - 1] += tmpexp;
              sumExp[tag[h] - 1] += tmpexp;

              addedH = true;
            } else {
              exps[npair][kk][hh] = 0.0;
              fpair3[npair][kk][hh] = 0.0;
            }
          }
        }
        if (addedH) {
          npair++;
          if (npair >= pairmax) error->one(FLERR, "neigh is too small");
        }
      }
    }

    f[i][0] += fixx;
    f[i][1] += fiyy;
    f[i][2] += fizz;
  }

  // SERIAL: communicate sumExp
  MPI_Allreduce(MPI_IN_PLACE, sumExp, maxID, MPI_DOUBLE, MPI_SUM, world);

  // PARALLEL: second pass - 3-body forces

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag, npair)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, npair, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    eval(ifrom, ito, thr, npair);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

void PairE3BOMP::eval(int iifrom, int iito, ThrData *const thr, int /*npair*/)
{
  double **f_thr = thr->get_f();
  tagint *tag = atom->tag;
  const int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;

  double fxtmp, fytmp, fztmp;
  double evdwl, fpair, partA, partB, partC;

  for (int ii = iifrom; ii < iito; ++ii) {
    for (int kk = 0; kk < NUMO; kk++) {
      const int i = pairO[ii][kk];
      const int otherO = (kk + 1) % 2;
      partB = eb *
          (sumExp[tag[pairO[ii][otherO]] - 1] + sumExp[tag[pairH[ii][otherO][0]] - 1] +
           sumExp[tag[pairH[ii][otherO][1]] - 1] -
           2 * (exps[ii][otherO][0] + exps[ii][otherO][1]));
      partC = ec * (sumExp[tag[i] - 1] - exps[ii][kk][0] - exps[ii][kk][1]);

      for (int hh = 0; hh < NUMH; hh++) {
        const int j = pairH[ii][kk][hh];
        const int otherH = (hh + 1) % 2;
        const int j2 = pairH[ii][kk][otherH];

        // type A
        partA = ea * (sumExp[tag[j2] - 1] - exps[ii][kk][otherH]);
        fpair = partA * fpair3[ii][kk][hh];
        fxtmp = fpair * del3[ii][kk][hh][0];
        fytmp = fpair * del3[ii][kk][hh][1];
        fztmp = fpair * del3[ii][kk][hh][2];

        f_thr[i][0] += fxtmp;
        f_thr[i][1] += fytmp;
        f_thr[i][2] += fztmp;
        f_thr[j][0] -= fxtmp;
        f_thr[j][1] -= fytmp;
        f_thr[j][2] -= fztmp;

        if (evflag) {
          evdwl = partA * exps[ii][kk][hh] * 0.5;
          ev_tally_thr(this, i, j, nlocal, newton_pair, evdwl, 0.0, fpair, del3[ii][kk][hh][0],
                       del3[ii][kk][hh][1], del3[ii][kk][hh][2], thr);
#pragma omp atomic
          pvector[1] += evdwl;
        }

        // type B
        fpair = partB * fpair3[ii][kk][hh];
        fxtmp = fpair * del3[ii][kk][hh][0];
        fytmp = fpair * del3[ii][kk][hh][1];
        fztmp = fpair * del3[ii][kk][hh][2];

        f_thr[i][0] += fxtmp;
        f_thr[i][1] += fytmp;
        f_thr[i][2] += fztmp;
        f_thr[j][0] -= fxtmp;
        f_thr[j][1] -= fytmp;
        f_thr[j][2] -= fztmp;

        if (evflag) {
          evdwl = partB * exps[ii][kk][hh] * 0.5;
          ev_tally_thr(this, i, j, nlocal, newton_pair, evdwl, 0.0, fpair, del3[ii][kk][hh][0],
                       del3[ii][kk][hh][1], del3[ii][kk][hh][2], thr);
#pragma omp atomic
          pvector[2] += evdwl;
        }

        // type C
        fpair = partC * fpair3[ii][kk][hh];
        fxtmp = fpair * del3[ii][kk][hh][0];
        fytmp = fpair * del3[ii][kk][hh][1];
        fztmp = fpair * del3[ii][kk][hh][2];

        f_thr[i][0] += fxtmp;
        f_thr[i][1] += fytmp;
        f_thr[i][2] += fztmp;
        f_thr[j][0] -= fxtmp;
        f_thr[j][1] -= fytmp;
        f_thr[j][2] -= fztmp;

        if (evflag) {
          evdwl = partC * exps[ii][kk][hh] * 0.5;
          ev_tally_thr(this, i, j, nlocal, newton_pair, evdwl, 0.0, fpair, del3[ii][kk][hh][0],
                       del3[ii][kk][hh][1], del3[ii][kk][hh][2], thr);
#pragma omp atomic
          pvector[3] += evdwl;
        }
      }    // end for hh in NUMH
    }      // end for kk in NUMO
  }        // end for ii in npair
}

/* ---------------------------------------------------------------------- */

double PairE3BOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairE3B::memory_usage();
  return bytes;
}
