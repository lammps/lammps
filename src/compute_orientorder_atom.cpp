/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:  Aidan Thompson (SNL)
                         Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "compute_orientorder_atom.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define QEPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtom::ComputeOrientOrderAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  qlist(NULL), distsq(NULL), nearest(NULL), rlist(NULL),
  qnarray(NULL), qnm_r(NULL), qnm_i(NULL), cglist(NULL)
{
  if (narg < 3 ) error->all(FLERR,"Illegal compute orientorder/atom command");

  // set default values for optional args

  nnn = 12;
  cutsq = 0.0;
  wlflag = 0;
  wlhatflag = 0;
  qlcompflag = 0;
  chunksize = 16384;

  // specify which orders to request

  nqlist = 5;
  memory->create(qlist,nqlist,"orientorder/atom:qlist");
  qlist[0] = 4;
  qlist[1] = 6;
  qlist[2] = 8;
  qlist[3] = 10;
  qlist[4] = 12;
  qmax = 12;

  // process optional args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nnn") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      if (strcmp(arg[iarg+1],"NULL") == 0) {
        nnn = 0;
      } else {
        nnn = force->numeric(FLERR,arg[iarg+1]);
        if (nnn <= 0)
          error->all(FLERR,"Illegal compute orientorder/atom command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"degrees") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      nqlist = force->numeric(FLERR,arg[iarg+1]);
      if (nqlist <= 0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      memory->destroy(qlist);
      memory->create(qlist,nqlist,"orientorder/atom:qlist");
      iarg += 2;
      if (iarg+nqlist > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      qmax = 0;
      for (int il = 0; il < nqlist; il++) {
        qlist[il] = force->numeric(FLERR,arg[iarg+il]);
        if (qlist[il] < 0)
          error->all(FLERR,"Illegal compute orientorder/atom command");
        if (qlist[il] > qmax) qmax = qlist[il];
      }
      iarg += nqlist;
    } else if (strcmp(arg[iarg],"wl") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) wlflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wlflag = 0;
      else error->all(FLERR,"Illegal compute orientorder/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"wl/hat") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) wlhatflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wlhatflag = 0;
      else error->all(FLERR,"Illegal compute orientorder/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"components") == 0) {
      qlcompflag = 1;
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      qlcomp = force->numeric(FLERR,arg[iarg+1]);
      iqlcomp = -1;
      for (int il = 0; il < nqlist; il++)
        if (qlcomp == qlist[il]) {
          iqlcomp = il;
          break;
        }
      if (iqlcomp == -1)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      double cutoff = force->numeric(FLERR,arg[iarg+1]);
      if (cutoff <= 0.0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      cutsq = cutoff*cutoff;
      iarg += 2;
    } else if (strcmp(arg[iarg],"chunksize") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      chunksize = force->numeric(FLERR,arg[iarg+1]);
      if (chunksize <= 0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute orientorder/atom command");
  }

  ncol = nqlist;
  if (wlflag) ncol += nqlist;
  if (wlhatflag) ncol += nqlist;
  if (qlcompflag) ncol += 2*(2*qlcomp+1);

  peratom_flag = 1;
  size_peratom_cols = ncol;

  nmax = 0;
  maxneigh = 0;
}

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtom::~ComputeOrientOrderAtom()
{
  if (copymode) return;

  memory->destroy(qnarray);
  memory->destroy(distsq);
  memory->destroy(rlist);
  memory->destroy(nearest);
  memory->destroy(qlist);
  memory->destroy(qnm_r);
  memory->destroy(qnm_i);
  memory->destroy(cglist);
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute orientorder/atom requires a "
               "pair style be defined");
  if (cutsq == 0.0) cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute orientorder/atom cutoff is "
               "longer than pairwise cutoff");

  memory->create(qnm_r,nqlist,2*qmax+1,"orientorder/atom:qnm_r");
  memory->create(qnm_i,nqlist,2*qmax+1,"orientorder/atom:qnm_i");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"orientorder/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute orientorder/atom");

  if (wlflag || wlhatflag) init_clebsch_gordan();
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(qnarray);
    nmax = atom->nmax;
    memory->create(qnarray,nmax,ncol,"orientorder/atom:qnarray");
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

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double* qn = qnarray[i];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // insure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsq);
        memory->destroy(rlist);
        memory->destroy(nearest);
        maxneigh = jnum;
        memory->create(distsq,maxneigh,"orientorder/atom:distsq");
        memory->create(rlist,maxneigh,3,"orientorder/atom:rlist");
        memory->create(nearest,maxneigh,"orientorder/atom:nearest");
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
        rsq = delx*delx + dely*dely + delz*delz;
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
        for (int jj = 0; jj < ncol; jj++)
          qn[jj] = 0.0;
        continue;
      }

      // if nnn > 0, use only nearest nnn neighbors

      if (nnn > 0) {
        select3(nnn,ncount,distsq,nearest,rlist);
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
  double bytes = ncol*nmax * sizeof(double);
  bytes += (qmax*(2*qmax+1)+maxneigh*4) * sizeof(double);
  bytes += (nqlist+maxneigh) * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   select3 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary arrays at same time
------------------------------------------------------------------------- */

// Use no-op do while to create single statement

#define SWAP(a,b) do {       \
    tmp = a; a = b; b = tmp; \
  } while(0)

#define ISWAP(a,b) do {        \
    itmp = a; a = b; b = itmp; \
  } while(0)

#define SWAP3(a,b) do {                  \
    tmp = a[0]; a[0] = b[0]; b[0] = tmp; \
    tmp = a[1]; a[1] = b[1]; b[1] = tmp; \
    tmp = a[2]; a[2] = b[2]; b[2] = tmp; \
  } while(0)

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::select3(int k, int n, double *arr, int *iarr, double **arr3)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp,a3[3];

  arr--;
  iarr--;
  arr3--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir]);
        ISWAP(iarr[l],iarr[ir]);
        SWAP3(arr3[l],arr3[ir]);
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1]);
      ISWAP(iarr[mid],iarr[l+1]);
      SWAP3(arr3[mid],arr3[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
        ISWAP(iarr[l],iarr[ir]);
        SWAP3(arr3[l],arr3[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
        ISWAP(iarr[l+1],iarr[ir]);
        SWAP3(arr3[l+1],arr3[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
        ISWAP(iarr[l],iarr[l+1]);
        SWAP3(arr3[l],arr3[l+1]);
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      a3[0] = arr3[l+1][0];
      a3[1] = arr3[l+1][1];
      a3[2] = arr3[l+1][2];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j]);
        ISWAP(iarr[i],iarr[j]);
        SWAP3(arr3[i],arr3[j]);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      arr3[l+1][0] = arr3[j][0];
      arr3[l+1][1] = arr3[j][1];
      arr3[l+1][2] = arr3[j][2];
      arr3[j][0] = a3[0];
      arr3[j][1] = a3[1];
      arr3[j][2] = a3[2];
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate the bond orientational order parameters
------------------------------------------------------------------------- */

void ComputeOrientOrderAtom::calc_boop(double **rlist,
                                       int ncount, double qn[],
                                       int qlist[], int nqlist) {

  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m = 0; m < 2*l+1; m++) {
      qnm_r[il][m] = 0.0;
      qnm_i[il][m] = 0.0;
    }
  }

  for(int ineigh = 0; ineigh < ncount; ineigh++) {
    const double * const r = rlist[ineigh];
    double rmag = dist(r);
    if(rmag <= MY_EPSILON) {
      return;
    }

    double costheta = r[2] / rmag;
    double expphi_r = r[0];
    double expphi_i = r[1];
    double rxymag = sqrt(expphi_r*expphi_r+expphi_i*expphi_i);
    if(rxymag <= MY_EPSILON) {
      expphi_r = 1.0;
      expphi_i = 0.0;
    } else {
      double rxymaginv = 1.0/rxymag;
      expphi_r *= rxymaginv;
      expphi_i *= rxymaginv;
    }

    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];

      qnm_r[il][l] += polar_prefactor(l, 0, costheta);
      double expphim_r = expphi_r;
      double expphim_i = expphi_i;
      for(int m = 1; m <= +l; m++) {
        double prefactor = polar_prefactor(l, m, costheta);
        double c_r = prefactor * expphim_r;
        double c_i = prefactor * expphim_i;
        qnm_r[il][m+l] += c_r;
        qnm_i[il][m+l] += c_i;
        if(m & 1) {
          qnm_r[il][-m+l] -= c_r;
          qnm_i[il][-m+l] += c_i;
        } else {
          qnm_r[il][-m+l] += c_r;
          qnm_i[il][-m+l] -= c_i;
        }
        double tmp_r = expphim_r*expphi_r - expphim_i*expphi_i;
        double tmp_i = expphim_r*expphi_i + expphim_i*expphi_r;
        expphim_r = tmp_r;
        expphim_i = tmp_i;
      }

    }
  }

  // convert sums to averages

  double facn = 1.0 / ncount;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m = 0; m < 2*l+1; m++) {
      qnm_r[il][m] *= facn;
      qnm_i[il][m] *= facn;
    }
  }

  // calculate Q_l
  // NOTE: optional W_l_hat and components of Q_qlcomp use these stored Q_l values

  int jj = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    double qnormfac = sqrt(MY_4PI/(2*l+1));
    double qm_sum = 0.0;
    for(int m = 0; m < 2*l+1; m++)
      qm_sum += qnm_r[il][m]*qnm_r[il][m] + qnm_i[il][m]*qnm_i[il][m];
    qn[jj++] = qnormfac * sqrt(qm_sum);
  }

  // TODO:
  // 1. [done]Need to allocate extra memory in qnarray[] for this option
  // 2. [done]Need to add keyword option
  // 3. [done]Need to calculate Clebsch-Gordan/Wigner 3j coefficients
  //     (Can try getting them from boop.py first)
  // 5. [done]Compare to bcc values in /Users/athomps/netapp/codes/MatMiner/matminer/matminer/featurizers/boop.py
  // 6. [done]I get the right answer for W_l, but need to make sure that factor of 1/sqrt(l+1) is right for cglist
  // 7. Add documentation
  // 8. [done] run valgrind
  // 9. [done] Add Wlhat
  // 10. Update memory_usage()
  // 11. Add exact FCC values for W_4, W_4_hat

  // calculate W_l

  if (wlflag) {
    int idxcg_count = 0;
    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];
      double wlsum = 0.0;
      for(int m1 = 0; m1 < 2*l+1; m1++) {
        for(int m2 = MAX(0,l-m1); m2 < MIN(2*l+1,3*l-m1+1); m2++) {
          int m = m1 + m2 - l;
          double qm1qm2_r = qnm_r[il][m1]*qnm_r[il][m2] - qnm_i[il][m1]*qnm_i[il][m2];
          double qm1qm2_i = qnm_r[il][m1]*qnm_i[il][m2] + qnm_i[il][m1]*qnm_r[il][m2];
          wlsum += (qm1qm2_r*qnm_r[il][m] + qm1qm2_i*qnm_i[il][m])*cglist[idxcg_count];
          idxcg_count++;
        }
      }
      qn[jj++] = wlsum/sqrt(2*l+1);
    }
  }

  // calculate W_l_hat

  if (wlhatflag) {
    int idxcg_count = 0;
    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];
      double wlsum = 0.0;
      for(int m1 = 0; m1 < 2*l+1; m1++) {
        for(int m2 = MAX(0,l-m1); m2 < MIN(2*l+1,3*l-m1+1); m2++) {
          int m = m1 + m2 - l;
          double qm1qm2_r = qnm_r[il][m1]*qnm_r[il][m2] - qnm_i[il][m1]*qnm_i[il][m2];
          double qm1qm2_i = qnm_r[il][m1]*qnm_i[il][m2] + qnm_i[il][m1]*qnm_r[il][m2];
          wlsum += (qm1qm2_r*qnm_r[il][m] + qm1qm2_i*qnm_i[il][m])*cglist[idxcg_count];
          idxcg_count++;
        }
      }
  //      Whats = [w/(q/np.sqrt(np.pi * 4 / (2 * l + 1)))**3 if abs(q) > 1.0e-6 else 0.0 for l,q,w in zip(range(1,max_l+1),Qs,Ws)]
      if (qn[il] < QEPSILON)
        qn[jj++] = 0.0;
      else {
        double qnormfac = sqrt(MY_4PI/(2*l+1));
        double qnfac = qnormfac/qn[il];
        qn[jj++] = wlsum/sqrt(2*l+1)*(qnfac*qnfac*qnfac);
      }
    }
  }

  // Calculate components of Q_l, for l=qlcomp

  if (qlcompflag) {
    int il = iqlcomp;
    int l = qlcomp;
    if (qn[il] < QEPSILON)
      for(int m = 0; m < 2*l+1; m++) {
        qn[jj++] = 0.0;
        qn[jj++] = 0.0;
      }
    else {
      double qnormfac = sqrt(MY_4PI/(2*l+1));
      double qnfac = qnormfac/qn[il];
      for(int m = 0; m < 2*l+1; m++) {
        qn[jj++] = qnm_r[il][m] * qnfac;
        qn[jj++] = qnm_i[il][m] * qnfac;
      }
    }
  }

}

/* ----------------------------------------------------------------------
   calculate scalar distance
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::dist(const double r[])
{
  return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

/* ----------------------------------------------------------------------
   polar prefactor for spherical harmonic Y_l^m, where
   Y_l^m (theta, phi) = prefactor(l, m, cos(theta)) * exp(i*m*phi)
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::polar_prefactor(int l, int m, double costheta)
{
  const int mabs = abs(m);

  double prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= static_cast<double>(i);

  prefactor = sqrt(static_cast<double>(2*l+1)/(MY_4PI*prefactor))
    * associated_legendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::associated_legendre(int l, int m, double x)
{
  if (l < m) return 0.0;

  double p(1.0), pm1(0.0), pm2(0.0);

  if (m != 0) {
    const double sqx = sqrt(1.0-x*x);
    for (int i=1; i < m+1; ++i)
      p *= static_cast<double>(2*i-1) * sqx;
  }

  for (int i=m+1; i < l+1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = (static_cast<double>(2*i-1)*x*pm1
         - static_cast<double>(i+m-1)*pm2) / static_cast<double>(i-m);
  }

  return p;
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients
   using the quasi-binomial formula VMK 8.2.1(3)
   specialized for case j1=j2=j=l
------------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init_clebsch_gordan()
{
  double sum,dcg,sfaccg, sfac1, sfac2;
  int m, aa2, bb2, cc2;
  int ifac, idxcg_count;

  idxcg_count = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m1 = 0; m1 < 2*l+1; m1++)
      for(int m2 = MAX(0,l-m1); m2 < MIN(2*l+1,3*l-m1+1); m2++)
        idxcg_count++;
  }
  idxcg_max = idxcg_count;
  memory->create(cglist, idxcg_max, "computeorientorderatom:cglist");

  idxcg_count = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m1 = 0; m1 < 2*l+1; m1++) {
        aa2 = m1 - l;
        for(int m2 = MAX(0,l-m1); m2 < MIN(2*l+1,3*l-m1+1); m2++) {
          bb2 = m2 - l;
          m = aa2 + bb2 + l;

          sum = 0.0;
          for (int z = MAX(0, MAX(-aa2, bb2));
               z <= MIN(l, MIN(l - aa2, l + bb2)); z++) {
            ifac = z % 2 ? -1 : 1;
            sum += ifac /
              (factorial(z) *
               factorial(l - z) *
               factorial(l - aa2 - z) *
               factorial(l + bb2 - z) *
               factorial(aa2 + z) *
               factorial(-bb2 + z));
          }

          cc2 = m - l;
          sfaccg = sqrt(factorial(l + aa2) *
                        factorial(l - aa2) *
                        factorial(l + bb2) *
                        factorial(l - bb2) *
                        factorial(l + cc2) *
                        factorial(l - cc2) *
                        (2*l + 1));

          sfac1 = factorial(3*l + 1);
          sfac2 = factorial(l);
          dcg = sqrt(sfac2*sfac2*sfac2 / sfac1);

          cglist[idxcg_count] = sum * dcg * sfaccg;
          idxcg_count++;
        }
      }
  }
}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

double ComputeOrientOrderAtom::factorial(int n)
{
  if (n < 0 || n > nmaxfactorial) {
    char str[128];
    sprintf(str, "Invalid argument to factorial %d", n);
    error->all(FLERR, str);
  }

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   factorial n table, size SNA::nmaxfactorial+1
------------------------------------------------------------------------- */

const double ComputeOrientOrderAtom::nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // nmaxfactorial = 167
};

