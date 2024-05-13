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
   Contributing authors: Aidan Thompson (SNL), Axel Kohlmeyer (Temple U)

   Tomas Oppelstrup (LLNL): Optimization which reduces the number
   of iterations in the L,m1,m2 loops (by a factor of up to 10), and
   avoids evaluation of Ylm functions of negative m
------------------------------------------------------------------------- */

#include "compute_orientorder_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <utility>

using namespace LAMMPS_NS;
using namespace MathConst;
using MathSpecial::factorial;

#ifdef DBL_EPSILON
static constexpr double MY_EPSILON = (10.0 * DBL_EPSILON);
#else
static constexpr double MY_EPSILON = (10.0 * 2.220446049250313e-16);
#endif

static constexpr double QEPSILON = 1.0e-6;

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtom::ComputeOrientOrderAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), qlist(nullptr), qnormfac(nullptr), qnormfac2(nullptr), distsq(nullptr),
    nearest(nullptr), rlist(nullptr), qnarray(nullptr), qnm_r(nullptr), qnm_i(nullptr),
    w3jlist(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal compute orientorder/atom command");

  // set default values for optional args

  nnn = 12;
  cutsq = 0.0;
  wlflag = 0;
  wlhatflag = 0;
  qlcompflag = 0;
  chunksize = 16384;

  // specify which orders to request

  nqlist = 5;
  memory->create(qlist, nqlist, "orientorder/atom:qlist");
  qlist[0] = 4;
  qlist[1] = 6;
  qlist[2] = 8;
  qlist[3] = 10;
  qlist[4] = 12;
  qmax = 12;

  // process optional args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "nnn") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      if (strcmp(arg[iarg + 1], "NULL") == 0) {
        nnn = 0;
      } else {
        nnn = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        if (nnn <= 0) error->all(FLERR, "Illegal compute orientorder/atom command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "degrees") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      nqlist = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (nqlist <= 0) error->all(FLERR, "Illegal compute orientorder/atom command");
      memory->destroy(qlist);
      memory->create(qlist, nqlist, "orientorder/atom:qlist");
      iarg += 2;
      if (iarg + nqlist > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      qmax = 0;
      for (int il = 0; il < nqlist; il++) {
        qlist[il] = utils::numeric(FLERR, arg[iarg + il], false, lmp);
        if (qlist[il] < 0) error->all(FLERR, "Illegal compute orientorder/atom command");
        if (qlist[il] > qmax) qmax = qlist[il];
      }
      iarg += nqlist;
    } else if (strcmp(arg[iarg], "wl") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      wlflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "wl/hat") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      wlhatflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "components") == 0) {
      qlcompflag = 1;
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      qlcomp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iqlcomp = -1;
      for (int il = 0; il < nqlist; il++)
        if (qlcomp == qlist[il]) {
          iqlcomp = il;
          break;
        }
      if (iqlcomp == -1) error->all(FLERR, "Illegal compute orientorder/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      double cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (cutoff <= 0.0) error->all(FLERR, "Illegal compute orientorder/atom command");
      cutsq = cutoff * cutoff;
      iarg += 2;
    } else if (strcmp(arg[iarg], "chunksize") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute orientorder/atom command");
      chunksize = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (chunksize <= 0) error->all(FLERR, "Illegal compute orientorder/atom command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute orientorder/atom command");
  }

  ncol = nqlist;
  if (wlflag) ncol += nqlist;
  if (wlhatflag) ncol += nqlist;
  if (qlcompflag) ncol += 2 * (2 * qlcomp + 1);

  peratom_flag = 1;
  size_peratom_cols = ncol;

  nmax = 0;
  maxneigh = 0;

  memory->create(qnormfac, nqlist, "orientorder/atom:qnormfac");
  memory->create(qnormfac2, nqlist, "orientorder/atom:qnormfac2");
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    qnormfac[il] = sqrt(MY_4PI / (2.0 * l + 1.0));
    qnormfac2[il] = sqrt(2.0 * l + 1.0);
  }
}

/* --------------------------------------------------------------------- */

ComputeOrientOrderAtom::~ComputeOrientOrderAtom()
{
  if (copymode) return;

  memory->destroy(qnarray);
  memory->destroy(distsq);
  memory->destroy(rlist);
  memory->destroy(nearest);
  memory->destroy(qlist);
  memory->destroy(qnormfac);
  memory->destroy(qnormfac2);
  memory->destroy(qnm_r);
  memory->destroy(qnm_i);
  memory->destroy(w3jlist);
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute orientorder/atom requires a pair style be defined");
  if (cutsq == 0.0)
    cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR, "Compute orientorder/atom cutoff is longer than pairwise cutoff");

  memory->destroy(qnm_r);
  memory->destroy(qnm_i);
  memory->create(qnm_r, nqlist, qmax + 1, "orientorder/atom:qnm_r");
  memory->create(qnm_i, nqlist, qmax + 1, "orientorder/atom:qnm_i");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if ((modify->get_compute_by_style("orientorder/atom").size() > 1) && (comm->me == 0))
    error->warning(FLERR, "More than one instance of compute orientorder/atom");

  if (wlflag || wlhatflag) init_wigner3j();
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::compute_peratom()
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(qnarray);
    nmax = atom->nmax;
    memory->create(qnarray, nmax, ncol, "orientorder/atom:qnarray");
    array_atom = qnarray;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;
  memset(&qnarray[0][0], 0, sizeof(double) * nmax * ncol);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double *qn = qnarray[i];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // ensure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsq);
        memory->destroy(rlist);
        memory->destroy(nearest);
        maxneigh = jnum;
        memory->create(distsq, maxneigh, "orientorder/atom:distsq");
        memory->create(rlist, maxneigh, 3, "orientorder/atom:rlist");
        memory->create(nearest, maxneigh, "orientorder/atom:nearest");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // rlist[] = distance vector to each
      // nearest[] = atom indices of neighbors

      int ncount = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          distsq[ncount] = rsq;
          rlist[ncount][0] = delx;
          rlist[ncount][1] = dely;
          rlist[ncount][2] = delz;
          nearest[ncount++] = j;
        }
      }

      // if not nnn neighbors, order parameter = 0;

      if ((ncount == 0) || (ncount < nnn)) {
        for (jj = 0; jj < ncol; jj++) qn[jj] = 0.0;
        continue;
      }

      // if nnn > 0, use only nearest nnn neighbors

      if (nnn > 0) {
        select3(nnn, ncount, distsq, nearest, rlist);
        ncount = nnn;
      }

      calc_boop(rlist, ncount, qn, qlist, nqlist);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::memory_usage()
{
  double bytes = (double) ncol * nmax * sizeof(double);
  bytes += (double) (qmax * (2 * qmax + 1) + maxneigh * 4) * sizeof(double);
  bytes += (double) (nqlist + maxneigh) * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   select3 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary arrays at same time
------------------------------------------------------------------------- */

// Use no-op do while to create single statement

#define SWAP3(a, b)            \
  do {                         \
    std::swap((a)[0], (b)[0]); \
    std::swap((a)[1], (b)[1]); \
    std::swap((a)[2], (b)[2]); \
  } while (0)

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::select3(int k, int n, double *arr, int *iarr, double **arr3)
{
  int i, ir, j, l, mid, ia;
  double a, a3[3];

  arr--;
  iarr--;
  arr3--;
  l = 1;
  ir = n;
  while (true) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
        SWAP3(arr3[l], arr3[ir]);
      }
      return;
    } else {
      mid = (l + ir) >> 1;
      std::swap(arr[mid], arr[l + 1]);
      std::swap(iarr[mid], iarr[l + 1]);
      SWAP3(arr3[mid], arr3[l + 1]);
      if (arr[l] > arr[ir]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
        SWAP3(arr3[l], arr3[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
        std::swap(arr[l + 1], arr[ir]);
        std::swap(iarr[l + 1], iarr[ir]);
        SWAP3(arr3[l + 1], arr3[ir]);
      }
      if (arr[l] > arr[l + 1]) {
        std::swap(arr[l], arr[l + 1]);
        std::swap(iarr[l], iarr[l + 1]);
        SWAP3(arr3[l], arr3[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      ia = iarr[l + 1];
      a3[0] = arr3[l + 1][0];
      a3[1] = arr3[l + 1][1];
      a3[2] = arr3[l + 1][2];
      while (true) {
        do i++;
        while (arr[i] < a);
        do j--;
        while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i], arr[j]);
        std::swap(iarr[i], iarr[j]);
        SWAP3(arr3[i], arr3[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      iarr[l + 1] = iarr[j];
      iarr[j] = ia;
      arr3[l + 1][0] = arr3[j][0];
      arr3[l + 1][1] = arr3[j][1];
      arr3[l + 1][2] = arr3[j][2];
      arr3[j][0] = a3[0];
      arr3[j][1] = a3[1];
      arr3[j][2] = a3[2];
      if (j >= k) ir = j - 1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate the bond orientational order parameters
------------------------------------------------------------------------- */

void ComputeOrientOrderAtom::calc_boop(double **rlist, int ncount, double qn[], int qlist[],
                                       int nqlist)
{

  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for (int m = 0; m < l + 1; m++) {
      qnm_r[il][m] = 0.0;
      qnm_i[il][m] = 0.0;
    }
  }

  for (int ineigh = 0; ineigh < ncount; ineigh++) {
    const double *const r = rlist[ineigh];
    double rmag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    if (rmag <= MY_EPSILON) { return; }

    double costheta = r[2] / rmag;
    double expphi_r = r[0];
    double expphi_i = r[1];
    double rxymag = sqrt(expphi_r * expphi_r + expphi_i * expphi_i);
    if (rxymag <= MY_EPSILON) {
      expphi_r = 1.0;
      expphi_i = 0.0;
    } else {
      double rxymaginv = 1.0 / rxymag;
      expphi_r *= rxymaginv;
      expphi_i *= rxymaginv;
    }

    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];

      // calculate spherical harmonics
      // Ylm, -l <= m <= l
      // sign convention: sign(Yll(0,0)) = (-1)^l

      qnm_r[il][0] += polar_prefactor(l, 0, costheta);
      double expphim_r = expphi_r;
      double expphim_i = expphi_i;
      for (int m = 1; m <= +l; m++) {

        double prefactor = polar_prefactor(l, m, costheta);
        double ylm_r = prefactor * expphim_r;
        double ylm_i = prefactor * expphim_i;
        qnm_r[il][m] += ylm_r;
        qnm_i[il][m] += ylm_i;
        // Skip calculation of qnm for m<0 due to symmetry
        double tmp_r = expphim_r * expphi_r - expphim_i * expphi_i;
        double tmp_i = expphim_r * expphi_i + expphim_i * expphi_r;
        expphim_r = tmp_r;
        expphim_i = tmp_i;
      }
    }
  }

  // convert sums to averages

  double facn = 1.0 / ncount;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for (int m = 0; m < l + 1; m++) {
      qnm_r[il][m] *= facn;
      qnm_i[il][m] *= facn;
    }
  }

  // calculate Q_l
  // NOTE: optional W_l_hat and components of Q_qlcomp use these stored Q_l values

  int jj = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    double qm_sum = qnm_r[il][0] * qnm_r[il][0];
    for (int m = 1; m < l + 1; m++)
      qm_sum += 2.0 * (qnm_r[il][m] * qnm_r[il][m] + qnm_i[il][m] * qnm_i[il][m]);
    qn[jj++] = qnormfac[il] * sqrt(qm_sum);
  }

  // calculate W_l

  int nterms = 0;
  int widx_count = 0;
  if (wlflag || wlhatflag) {
    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];
      double wlsum = 0.0;
      for (int m1 = -l; m1 <= 0; m1++) {
        const int sgn = 1 - 2 * (m1 & 1);    // sgn = (-1)^m1
        for (int m2 = 0; m2 <= ((-m1) >> 1); m2++) {
          const int m3 = -(m1 + m2);
          // Loop enforces -L <= m1 <= 0 <= m2 <= m3 <= L, and m1 + m2 + m3 = 0

          // For even L, W3j is invariant under permutation of
          // (m1, m2, m3) and (m1, m2, m3) -> (-m1, -m2, -m3). The loop
          // structure enforces visiting only one member of each
          // such symmetry (invariance) group.

          // m1 <= 0, and Qlm[-m] = (-1)^m * conjg(Qlm[m])
          const double Q1Q2_r =
              (qnm_r[il][-m1] * qnm_r[il][m2] + qnm_i[il][-m1] * qnm_i[il][m2]) * sgn;
          const double Q1Q2_i =
              (qnm_r[il][-m1] * qnm_i[il][m2] - qnm_i[il][-m1] * qnm_r[il][m2]) * sgn;
          const double Q1Q2Q3 = Q1Q2_r * qnm_r[il][m3] - Q1Q2_i * qnm_i[il][m3];
          const double c = w3jlist[widx_count++];
          wlsum += Q1Q2Q3 * c;
        }
      }
      qn[jj++] = wlsum / qnormfac2[il];
      nterms++;
    }
  }

  // calculate W_l_hat

  if (wlhatflag) {
    const int jptr = jj - nterms;
    if (!wlflag) jj = jptr;
    for (int il = 0; il < nqlist; il++) {
      if (qn[il] < QEPSILON)
        qn[jj++] = 0.0;
      else {
        double qnfac = qnormfac[il] / qn[il];
        qn[jj++] = qn[jptr + il] * (qnfac * qnfac * qnfac) * qnormfac2[il];
      }
    }
  }

  // Calculate components of Q_l/|Q_l|, for l=qlcomp

  if (qlcompflag) {
    int il = iqlcomp;
    int l = qlcomp;
    if (qn[il] < QEPSILON)
      for (int m = 0; m < 2 * l + 1; m++) {
        qn[jj++] = 0.0;
        qn[jj++] = 0.0;
      }
    else {
      double qnfac = qnormfac[il] / qn[il];
      for (int m = -l; m < 0; m++) {
        // Computed only qnm for m>=0.
        // qnm[-m] = (-1)^m * conjg(qnm[m])
        const int sgn = 1 - 2 * (m & 1);    // sgn = (-1)^m
        qn[jj++] = qnm_r[il][-m] * qnfac * sgn;
        qn[jj++] = -qnm_i[il][-m] * qnfac * sgn;
      }
      for (int m = 0; m < l + 1; m++) {
        qn[jj++] = qnm_r[il][m] * qnfac;
        qn[jj++] = qnm_i[il][m] * qnfac;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   polar prefactor for spherical harmonic Y_l^m, where
   Y_l^m (theta, phi) = prefactor(l, m, cos(theta)) * exp(i*m*phi)
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::polar_prefactor(int l, int m, double costheta)
{
  const int mabs = abs(m);

  double prefactor = 1.0;
  for (int i = l - mabs + 1; i < l + mabs + 1; ++i) prefactor *= static_cast<double>(i);

  prefactor = sqrt(static_cast<double>(2 * l + 1) / (MY_4PI * prefactor)) *
      associated_legendre(l, mabs, costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
   sign convention: P(l,l) = (2l-1)!!(-sqrt(1-x^2))^l
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::associated_legendre(int l, int m, double x)
{
  if (l < m) return 0.0;

  double p(1.0), pm1(0.0), pm2(0.0);

  if (m != 0) {
    const double msqx = -sqrt(1.0 - x * x);
    for (int i = 1; i < m + 1; ++i) p *= static_cast<double>(2 * i - 1) * msqx;
  }

  for (int i = m + 1; i < l + 1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = (static_cast<double>(2 * i - 1) * x * pm1 - static_cast<double>(i + m - 1) * pm2) /
        static_cast<double>(i - m);
  }

  return p;
}

/* ----------------------------------------------------------------------
  Initialize table of Wigner 3j symbols
------------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init_wigner3j()
{
  int widx_count = 0;

  for (int il = 0; il < nqlist; il++) {
    const int l = qlist[il];

    for (int m1 = -l; m1 <= 0; m1++) {
      for (int m2 = 0; m2 <= ((-m1) >> 1); m2++) { widx_count++; }
    }
  }
  widx_max = widx_count;
  memory->destroy(w3jlist);
  memory->create(w3jlist, widx_max, "computeorientorderatom:w3jlist");

  widx_count = 0;

  for (int il = 0; il < nqlist; il++) {
    const int l = qlist[il];

    for (int m1 = -l; m1 <= 0; m1++) {
      for (int m2 = 0; m2 <= ((-m1) >> 1); m2++) {
        const int m3 = -(m1 + m2);
        // Loop enforces -L<=m1<=0<=m2<=m3<=L, and m1+m2+m3=0

        // For even L, W3j is invariant under permutation of
        // (m1,m2,m3) and (m1,m2,m3)->(-m1,-m2,-m3). The loop
        // structure enforces visiting only one member of each
        // such symmetry (invariance) group.

        // Determine number of elements in symmetry group of (m1,m2,m3)
        // Concise determination exploiting (m1,m2,m3) loop structure.
        int pfac;
        if (m1 == 0)
          pfac = 1;    // m1 = m2 = m3 = 0
        else if (m2 == 0 || m2 == m3) {
          // reduced group when only 3 permutations, or sign inversion
          // is equivalent to permutation
          pfac = 6;
        } else
          pfac = 12;    // 6 permutations * 2 signs

        w3jlist[widx_count] = w3j(l, m1, m2, m3) * pfac;
        widx_count++;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double ComputeOrientOrderAtom::triangle_coeff(const int a, const int b, const int c)
{
  return factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c) /
      factorial(a + b + c + 1);
}

/* ---------------------------------------------------------------------- */

double ComputeOrientOrderAtom::w3j(const int lmax, const int j1, const int j2, const int j3)
{
  const int a = lmax, b = lmax, c = lmax;
  const int alpha = j1, beta = j2, gamma = j3;
  struct {
    double operator()(const int a, const int b, const int c, const int alpha, const int beta,
                      const int t)
    {
      return factorial(t) * factorial(c - b + t + alpha) * factorial(c - a + t - beta) *
          factorial(a + b - c - t) * factorial(a - t - alpha) * factorial(b - t + beta);
    }
  } x;
  const double sgn = 1 - 2 * ((a - b - gamma) & 1);
  const double g = sqrt(triangle_coeff(lmax, lmax, lmax)) *
      sqrt(factorial(a + alpha) * factorial(a - alpha) * factorial(b + beta) * factorial(b - beta) *
           factorial(c + gamma) * factorial(c - gamma));
  double s = 0;
  int t = 0;
  while (c - b + t + alpha < 0 || c - a + t - beta < 0) t++;
  //     ^^ t>=-j1       ^^ t>=j2
  while (true) {
    if (a + b - c - t < 0) break;    // t<=lmax
    if (a - t - alpha < 0) break;    // t<=lmax-j1
    if (b - t + beta < 0) break;     // t<=lmax+j2
    const int m1t = 1 - 2 * (t & 1);
    s += m1t / x(lmax, lmax, lmax, alpha, beta, t);
    t++;
  }
  return sgn * g * s;
}
