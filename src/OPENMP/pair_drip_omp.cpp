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

#include "pair_drip_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

static constexpr double HALF = 0.5;

/* ---------------------------------------------------------------------- */

PairDRIPOMP::PairDRIPOMP(LAMMPS *lmp) : PairDRIP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairDRIPOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // find_nearest3neigh must be called before the OMP parallel region
  find_nearest3neigh();

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
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (vflag_either)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (vflag_either)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else
      eval<0, 0, 0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int VFLAG_EITHER>
void PairDRIPOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  double **x = atom->x;
  double **f_thr = thr->get_f();
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double ni[3];
  double dni_dri[3][3], dni_drnb1[3][3];
  double dni_drnb2[3][3], dni_drnb3[3][3];

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    if (nearest3neigh[i][0] == -1) continue;
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = map[type[i]];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    // normal and its derivatives w.r.t. atom i and its 3 nearest neighbors
    calc_normal(i, ni, dni_dri, dni_drnb1, dni_drnb2, dni_drnb3);

    double fi[3] = {0., 0., 0.};

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (nearest3neigh[j][0] == -1) continue;
      const int jtype = map[type[j]];

      const double delx = x[j][0] - xtmp;
      const double dely = x[j][1] - ytmp;
      const double delz = x[j][2] - ztmp;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int iparam_ij = elem2param[itype][jtype];
      Param &p = params[iparam_ij];

      // only include the interaction between different layers
      if (rsq < p.rcutsq && atom->molecule[i] != atom->molecule[j]) {

        double fj[3] = {0., 0., 0.};
        const double rvec[3] = {delx, dely, delz};

        const double phi_attr = calc_attractive(p, rsq, rvec, fi, fj);

        const double phi_repul = calc_repulsive_thr(i, j, p, rsq, rvec, ni, dni_dri, dni_drnb1,
                                                    dni_drnb2, dni_drnb3, fi, fj, f_thr, thr);

        double evdwl = 0.0;
        if (EFLAG) evdwl = HALF * (phi_repul + phi_attr);
        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1, evdwl, 0.0, 0, 0, 0, 0, thr);

        f_thr[j][0] += fj[0];
        f_thr[j][1] += fj[1];
        f_thr[j][2] += fj[2];
        if (VFLAG_EITHER) v_tally2_newton_thr(this, j, fj, x[j], thr);
      }
    }    // loop over jj

    f_thr[i][0] += fi[0];
    f_thr[i][1] += fi[1];
    f_thr[i][2] += fi[2];
    if (VFLAG_EITHER) v_tally2_newton_thr(this, i, fi, x[i], thr);
  }    // loop over ii
}

/* ---------------------------------------------------------------------- */

double PairDRIPOMP::calc_repulsive_thr(int const i, int const j, Param &p, double const rsq,
                                       double const *rvec, double const *ni, V3 const *dni_dri,
                                       V3 const *dni_drnb1, V3 const *dni_drnb2,
                                       V3 const *dni_drnb3, double *const fi, double *const fj,
                                       double **f_thr, ThrData *const thr)
{
  double **x = atom->x;

  const double C0 = p.C0;
  const double C2 = p.C2;
  const double C4 = p.C4;
  const double C = p.C;
  const double delta = p.delta;
  const double lambda = p.lambda;
  const double z0 = p.z0;
  const double cutoff = p.rcut;

  // nearest 3 neighbors of atoms i and j
  const int nbi1 = nearest3neigh[i][0];
  const int nbi2 = nearest3neigh[i][1];
  const int nbi3 = nearest3neigh[i][2];
  const int nbj1 = nearest3neigh[j][0];
  const int nbj2 = nearest3neigh[j][1];
  const int nbj3 = nearest3neigh[j][2];

  double fnbi1[3], fnbi2[3], fnbi3[3];
  double fnbj1[3], fnbj2[3], fnbj3[3];
  V3 dgij_dri, dgij_drj;
  V3 dgij_drk1, dgij_drk2, dgij_drk3;
  V3 dgij_drl1, dgij_drl2, dgij_drl3;
  V3 drhosqij_dri, drhosqij_drj;
  V3 drhosqij_drnb1, drhosqij_drnb2, drhosqij_drnb3;

  const double r = sqrt(rsq);

  get_drhosqij(rvec, ni, dni_dri, dni_drnb1, dni_drnb2, dni_drnb3, drhosqij_dri, drhosqij_drj,
               drhosqij_drnb1, drhosqij_drnb2, drhosqij_drnb3);

  double rhosqij, dtdij;
  const double tdij = td(C0, C2, C4, delta, rvec, r, ni, rhosqij, dtdij);

  double dgij_drhosq;
  const double gij = dihedral(i, j, p, rhosqij, dgij_drhosq, dgij_dri, dgij_drj, dgij_drk1,
                               dgij_drk2, dgij_drk3, dgij_drl1, dgij_drl2, dgij_drl3);

  const double V2 = C + tdij + gij;

  double dtp;
  const double tp = tap(r, cutoff, dtp);

  const double V1 = exp(-lambda * (r - z0));
  const double dV1 = -V1 * lambda;

  const double phi = tp * V1 * V2;

  for (int k = 0; k < 3; k++) {
    const double tmp = HALF * (dtp * V1 + tp * dV1) * V2 * rvec[k] / r;
    fi[k] += tmp;
    fj[k] -= tmp;

    fi[k] -= HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_dri[k] + dgij_dri[k]);
    fj[k] -= HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drj[k] + dgij_drj[k]);
    fnbi1[k] = -HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drnb1[k] + dgij_drk1[k]);
    fnbi2[k] = -HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drnb2[k] + dgij_drk2[k]);
    fnbi3[k] = -HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drnb3[k] + dgij_drk3[k]);
    fnbj1[k] = -HALF * tp * V1 * dgij_drl1[k];
    fnbj2[k] = -HALF * tp * V1 * dgij_drl2[k];
    fnbj3[k] = -HALF * tp * V1 * dgij_drl3[k];
  }

  for (int k = 0; k < 3; k++) {
    f_thr[nbi1][k] += fnbi1[k];
    f_thr[nbi2][k] += fnbi2[k];
    f_thr[nbi3][k] += fnbi3[k];
    f_thr[nbj1][k] += fnbj1[k];
    f_thr[nbj2][k] += fnbj2[k];
    f_thr[nbj3][k] += fnbj3[k];
  }

  if (vflag_either) {
    v_tally2_newton_thr(this, nbi1, fnbi1, x[nbi1], thr);
    v_tally2_newton_thr(this, nbi2, fnbi2, x[nbi2], thr);
    v_tally2_newton_thr(this, nbi3, fnbi3, x[nbi3], thr);
    v_tally2_newton_thr(this, nbj1, fnbj1, x[nbj1], thr);
    v_tally2_newton_thr(this, nbj2, fnbj2, x[nbj2], thr);
    v_tally2_newton_thr(this, nbj3, fnbj3, x[nbj3], thr);
  }

  return phi;
}

/* ---------------------------------------------------------------------- */

double PairDRIPOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDRIP::memory_usage();
  return bytes;
}
