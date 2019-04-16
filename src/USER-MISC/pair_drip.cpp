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
   Contributing author: Mingjian Wen (University of Minnesota)
   e-mail: wenxx151@umn.edu
   based on "pair_style kolmogorov/crespi/full" by Wengen Ouyang

   This implements the DRIP model as described in
   M. Wen, S. Carr, S. Fang, E. Kaxiras, and E. B. Tadmor, Phys. Rev. B, 98,
   235404 (2018).
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include "pair_drip.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1
#define DIM 3
#define HALF 0.5

/* ---------------------------------------------------------------------- */

PairDRIP::PairDRIP(LAMMPS *lmp) : Pair(lmp)
{
  // initialize element to parameter maps
  single_enable = 0;
  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  cutmax = 0.0;
  nmax = 0;
  maxlocal = 0;
  DRIP_numneigh = NULL;
  DRIP_firstneigh = NULL;
  ipage = NULL;
  pgsize = oneatom = 0;

  normal = NULL;
  dnormal = NULL;
  dnormdri = NULL;




  nearest3neigh = NULL;



  // set comm size needed by this Pair
  comm_forward = 39;
  tap_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairDRIP::~PairDRIP()
{
  memory->destroy(DRIP_numneigh);
  memory->sfree(DRIP_firstneigh);
  delete [] ipage;
  memory->destroy(normal);
  memory->destroy(dnormal);
  memory->destroy(dnormdri);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);
  if (allocated) delete [] map;


  memory->destroy(nearest3neigh);

}


/* ---------------------------------------------------------------------- */

void PairDRIP::compute(int eflag, int vflag)
{

  int i,j,ii,jj,inum,jnum,itype,jtype,k,l,kk,ll;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,fpair1,fpair2, r, rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  int nbi1, nbi2, nbi3;
  double ni[DIM];
  double dni_dri[DIM][DIM], dni_drnb1[DIM][DIM];
  double dni_drnb2[DIM][DIM], dni_drnb3[DIM][DIM];


  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;





  double prodnorm1,prodnorm2,fkcx,fkcy,fkcz;
  double rhosq1,rhosq2,exp0,exp1,exp2,r2inv,r6inv,r8inv,Tap,dTap,Vkc;
  double frho1,frho2,sumC1,sumC2,sumC11,sumC22,sumCff,fsum,rdsq1,rdsq2;
  int *DRIP_neighs_i,*DRIP_neighs_j;


  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double dprodnorm2[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fp2[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double fprod2[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fl[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double delli[3] = {0.0, 0.0, 0.0};




  // find nearest 3 neighbors of each atom
  nearest3neigh();

  //TODO what does this comm do?
  // communicate the normal vector and its derivatives
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];


    // normal and its derivatives w.r.t. atom i and its 3 nearest neighs
    calc_normal(i, nbi1, nbi2, nbi3, ni, dni_dri,dni_drnb1, dni_drnb2, dni_drnb3);


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      jtag = tag[j];

//      // two-body interactions from full neighbor list, skip half of them
//      if (itag > jtag) {
//        if ((itag+jtag) % 2 == 0) continue;
//      } else if (itag < jtag) {
//        if ((itag+jtag) % 2 == 1) continue;
//      } else {
//        if (x[j][2] < ztmp) continue;
//        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
//        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
//      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      int iparam_ij = elem2param[itype][jtype];
      Param& p = params[iparam_ij];
      double rcutsq = p.rcutsq;


      // only include the interation between different layers
      if (rsq>1e-20 && rsq < rcutsq && atom->molecule[i] != atom->molecule[j]) {

        double phi_attr = calc_attractive(vflag, i,j,p, delx, dely, delz, rsq);

        double phi_repul = calc_repulsive(vflag);



      }
    }  //loop over jj
  }  // loop over ii

}


/* ---------------------------------------------------------------------- */

double PairDRIP::calc_attractive(int const i, int const j, Param& p,
    double const rsq, double const * rvec)
{

  double **f = atom->f;

  double const z0 = p.z0;
  double const A = p.A;
  double const cutoff = p.rcut;
  double const r = sqrt(rsq);

  double roz0_sq = rsq / (z0 * z0);
  double dtp;
  double tp = tap(r, cutoff, dtp);
  double r6 = A / (roz0_sq * roz0_sq * roz0_sq);
  double dr6 = -6 * r6 / r;
  double phi = -r6 * tp;

  double fpair = HALF * (r6 * dtp + dr6 * tp);
  f[i][0] += rvec[0] * fpair / r;
  f[i][1] += rvec[1] * fpair / r;
  f[i][2] += rvec[2] * fpair / r;
  f[j][0] -= rvec[0] * fpair / r;
  f[j][1] -= rvec[1] * fpair / r;
  f[j][2] -= rvec[2] * fpair / r;

  return phi;
}


/* ---------------------------------------------------------------------- */

double PairDRIP::calc_repulsive(int const i, int const j, Param& p,
    double const rsq, double const * rvec,
    int const nbi1, int const nbi2, int const nbi3, double const * ni,
    double const * dni_dr[DIM], double const * dni_drnb1[DIM],
    double const * dni_drnb2[DIM], double const * dni_drnb3[DIM])
{
  double **f = atom->f;
  double r = sqrt(rsq);

  // params
  double C0 = p.C0;
  double C2 = p.C2;
  double C4 = p.C4;
  double C = p.C;
  double delta = p.delta;
  double lambda = p.lambda;
  double z0 = p.z0;
  double cutoff = p.rcut;

  // nearest 3 neighbors of atom j
  int nbj1 = nearest3neigh[j][0];
  int nbj2 = nearest3neigh[j][1];
  int nbj3 = nearest3neigh[j][2];

  double[DIM] dgij_dri;
  double[DIM] dgij_drj;
  double[DIM] dgij_drk1;
  double[DIM] dgij_drk2;
  double[DIM] dgij_drk3;
  double[DIM] dgij_drl1;
  double[DIM] dgij_drl2;
  double[DIM] dgij_drl3;
  double[DIM] drhosqij_dri;
  double[DIM] drhosqij_drj;
  double[DIM] drhosqij_drnb1;
  double[DIM] drhosqij_drnb2;
  double[DIM] drhosqij_drnb3;


  // derivative of rhosq w.r.t coordinates of atoms i, j, and the nearests 3
  // neighs of i
  get_drhosqij(rvec, ni, dni_dri, dni_drnb1, dni_drnb2, dni_drnb3, drhosqij_dri,
      drhosqij_drj, drhosqij_drnb1, drhosqij_drnb2, drhosqij_drnb3);

  // transverse decay function f(rho) and its derivative w.r.t. rhosq
  double rhosqij;
  double dtdij;
  double tdij = td(C0, C2, C4, delta, rvec, r, ni, rhosqij, dtdij);

  // dihedral angle function and its derivateives
  double dgij_drhosq;
  double gij = dihedral(i, j, p, rhosqij, dgij_drhosq, dgij_dri, dgij_drj,
      dgij_drk1, dgij_drk2, dgij_drk3, dgij_drl1, dgij_drl2, dgij_drl3);

  double V2 = C + tdij + gij;

  // tap part
  double dtp;
  double tp = tap(r, cutoff, dtp);

  /* exponential part */
  double V1 = exp(-lambda * (r - z0));
  double dV1 = -V1 * lambda;

  double phi = tp * V1 * V2;

  for (int k = 0; k < DIM; k++) {

    // forces due to derivatives of tap and V1
    double tmp = -HALF * (dtp * V1 + tp * dV1) * V2 * rvec[k] / r;
    f[i][k] += tmp;
    f[j][k] -= tmp;

    // the following incldue the transverse decay part tdij and the dihedral part gij
    // derivative of V2 contribute to atoms i, j
    f[i][k] += HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_dri[k] + dgij_dri[k]);
    f[j][k] += HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drj[k] + dgij_drj[k]);

    // derivative of V2 contribute to neighs of atom i
    f[nbi1][k] += HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drnb1[k] + dgij_drk1[k]);
    f[nbi2][k] += HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drnb2[k] + dgij_drk2[k]);
    f[nbi3][k] += HALF * tp * V1 * ((dtdij + dgij_drhosq) * drhosqij_drnb3[k] + dgij_drk3[k]);

    // derivative of V2 contribute to neighs of atom j
    f[nbj1][k] += HALF * tp * V1 * dgij_drl1[k];
    f[nbj2][k] += HALF * tp * V1 * dgij_drl2[k];
    f[nbj3][k] += HALF * tp * V1 * dgij_drl3[k];
  }

  return phi;
}



/* ---------------------------------------------------------------------- */

void PairDRIP::find_nearest_3_neigh()
{

  int i,j,ii,jj,n,allnum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;


  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  memory->destroy(nearest3neigh);
  memory->create(nearest3neigh, allnum, 3, "DRIP:nearest3neigh");

  // store all DRIP neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];


    // init nb1 to be the 1st nearest neigh, nb3 the 3rd nearest
    int nb1 = -1;
    int nb2 = -1;
    int nb3 = -1;
    double nb1_rsq = 1.1e10;
    double nb2_rsq = 2.0e10;
    double nb3_rsq = 3.0e10;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      int iparam_ij = elem2param[itype][jtype];
      double rcutsq = params[iparam_ij].rcutsq;

      if (rsq > 1e-20 && rsq < rcutsq && atom->molecule[i] == atom->molecule[j]) {

        // find the 3 nearest neigh
        if (rsq < nb1_rsq) {
          nb3 = nb2;
          nb2 = nb1;
          nb1 = j;
          nb3_rsq = nb2_rsq;
          nb2_rsq = nb1_rsq;
          nb1_rsq = rsq;
        }
        else if (rsq < nb2_rsq) {
          nb3 = nb2;
          nb2 = j;
          nb3_rsq = nb2_rsq;
          nb2_rsq = rsq;
        }
        else if (rsq < nb3_rsq) {
          nb3 = j;
          nb3_rsq = rsq;
        }

      }
    }  // loop over jj

    // store neighbors to be used later to compute normal
    if (nb1_rsq >= 1.0e10 || nb2_rsq >= 1.0e10 || nb3_rsq >= 1.0e10) {
      error->one(FLERR,"No enough neighbors to construct normal.");
    } else{
      nearest3neigh[i][0] = nb1;
      nearest3neigh[i][1] = nb2;
      nearest3neigh[i][2] = nb3;
    }

  } // loop over ii
}


/* ---------------------------------------------------------------------- */

void PairDRIP::calc_normal(int const i, double ** const nearest3neigh,
    int& k1, int& k2, int& k3, double * const normal,
    double ** const dn_dri, double ** const dn_drk1,
    double ** const dn_drk2, double ** const dn_drk3)
{

  k1 = nearest3neigh[i][0];
  k2 = nearest3neigh[i][1];
  k3 = nearest3neigh[i][2];

  // normal does not depend on i, setting to zero
  for (int j = 0; j < DIM; j++) {
    for (int k = 0; k < DIM; k++) {
      dn_dri[j][k] = 0.0;
    }
  }

  // get normal and derives of normal w.r.t to its 3 nearest neighbors
  double **x = atom->x;
  deriv_cross(x[k1], x[k2], x[k3], normal, dn_drk1, dn_drk2, dn_drk3);
}


/* ---------------------------------------------------------------------- */
void PairDRIP::get_drhosqij(
    double const* const rij,
    double const* const ni,
    VectorOfSizeDIM const* const dni_dri,
    VectorOfSizeDIM const* const dni_drn1,
    VectorOfSizeDIM const* const dni_drn2,
    VectorOfSizeDIM const* const dni_drn3,
    double* const drhosq_dri,
    double* const drhosq_drj,
    double* const drhosq_drn1,
    double* const drhosq_drn2,
    double* const drhosq_drn3) const
{
  int k;
  double ni_dot_rij = 0;
  double dni_dri_dot_rij[DIM];
  double dni_drn1_dot_rij[DIM];
  double dni_drn2_dot_rij[DIM];
  double dni_drn3_dot_rij[DIM];

  ni_dot_rij = dot(ni, rij);
  mat_dot_vec(dni_dri, rij, dni_dri_dot_rij);
  mat_dot_vec(dni_drn1, rij, dni_drn1_dot_rij);
  mat_dot_vec(dni_drn2, rij, dni_drn2_dot_rij);
  mat_dot_vec(dni_drn3, rij, dni_drn3_dot_rij);

  for (k = 0; k < DIM; k++) {
    drhosq_dri[k] = -2 * rij[k] - 2 * ni_dot_rij * (-ni[k] + dni_dri_dot_rij[k]);
    drhosq_drj[k] = 2 * rij[k] - 2 * ni_dot_rij * ni[k];
    drhosq_drn1[k] = -2 * ni_dot_rij * dni_drn1_dot_rij[k];
    drhosq_drn2[k] = -2 * ni_dot_rij * dni_drn2_dot_rij[k];
    drhosq_drn3[k] = -2 * ni_dot_rij * dni_drn3_dot_rij[k];
  }
}


/* ---------------------------------------------------------------------- */
// Compute the normalized cross product of two vector rkl, rkm, and the
// derivates w.r.t rk, rl, rm.
// NOTE, the dcross_drk, dcross_drl, and dcross_drm is actually the transpose
// of the actual one.

void PairDRIP::deriv_cross( double const* rk, double const* rl, double const* rm,
    double* const cross, double ** const dcross_drk,
    double ** const dcross_drl, double ** const dcross_drm)
{
  double x[DIM];
  double y[DIM];
  double p[DIM];
  double q;
  double q_cubic;
  double d_invq_d_x0;
  double d_invq_d_x1;
  double d_invq_d_x2;
  double d_invq_d_y0;
  double d_invq_d_y1;
  double d_invq_d_y2;

  int i, j;


  // get x = rkl and y = rkm
  for (i = 0; i < DIM; i++) {
    x[i] = rl[i] - rk[i];
    y[i] = rm[i] - rk[i];
  }

  // cross product
  p[0] = x[1] * y[2] - x[2] * y[1];
  p[1] = x[2] * y[0] - x[0] * y[2];
  p[2] = x[0] * y[1] - x[1] * y[0];

  q = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

  // normalized cross
  cross[0] = p[0] / q;
  cross[1] = p[1] / q;
  cross[2] = p[2] / q;

  // compute derivatives
  // derivative of inverse q (i.e. 1/q) w.r.t x and y
  q_cubic = q * q * q;
  d_invq_d_x0 = (+p[1] * y[2] - p[2] * y[1]) / q_cubic;
  d_invq_d_x1 = (-p[0] * y[2] + p[2] * y[0]) / q_cubic;
  d_invq_d_x2 = (p[0] * y[1] - p[1] * y[0]) / q_cubic;
  d_invq_d_y0 = (-p[1] * x[2] + p[2] * x[1]) / q_cubic;
  d_invq_d_y1 = (p[0] * x[2] - p[2] * x[0]) / q_cubic;
  d_invq_d_y2 = (-p[0] * x[1] + p[1] * x[0]) / q_cubic;

  // dcross/drl transposed
  dcross_drl[0][0] = p[0] * d_invq_d_x0;
  dcross_drl[0][1] = -y[2] / q + p[1] * d_invq_d_x0;
  dcross_drl[0][2] = y[1] / q + p[2] * d_invq_d_x0;

  dcross_drl[1][0] = y[2] / q + p[0] * d_invq_d_x1;
  dcross_drl[1][1] = p[1] * d_invq_d_x1;
  dcross_drl[1][2] = -y[0] / q + p[2] * d_invq_d_x1;

  dcross_drl[2][0] = -y[1] / q + p[0] * d_invq_d_x2;
  dcross_drl[2][1] = y[0] / q + p[1] * d_invq_d_x2;
  dcross_drl[2][2] = p[2] * d_invq_d_x2;

  // dcross/drm transposed
  dcross_drm[0][0] = p[0] * d_invq_d_y0;
  dcross_drm[0][1] = x[2] / q + p[1] * d_invq_d_y0;
  dcross_drm[0][2] = -x[1] / q + p[2] * d_invq_d_y0;

  dcross_drm[1][0] = -x[2] / q + p[0] * d_invq_d_y1;
  dcross_drm[1][1] = p[1] * d_invq_d_y1;
  dcross_drm[1][2] = x[0] / q + p[2] * d_invq_d_y1;

  dcross_drm[2][0] = x[1] / q + p[0] * d_invq_d_y2;
  dcross_drm[2][1] = -x[0] / q + p[1] * d_invq_d_y2;
  dcross_drm[2][2] = p[2] * d_invq_d_y2;

  // dcross/drk transposed
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      dcross_drk[i][j] = -(dcross_drl[i][j] + dcross_drm[i][j]);
    }
  }

}

/* ---------------------------------------------------------------------- */


// derivartive of transverse decay function f(rho) w.r.t rho
double PairDRIP::td(double C0, double C2, double C4, double delta,
    double const* const rvec, double r,
    const double* const n,
    double& rho_sq, double& dtd) const
{
  double n_dot_r = dot(n, rvec);

  rho_sq = r * r - n_dot_r * n_dot_r;

  if (rho_sq < 0) {   // in case n is [0, 0, 1] and rho_sq is negative due to numerical error
    rho_sq = 0;
  }

  double del_sq = delta * delta;
  double rod_sq = rho_sq / del_sq;
  double td = exp(-rod_sq) * (C0 + rod_sq * (C2 + rod_sq * C4));
  dtd = -td / del_sq + exp(-rod_sq) * (C2 + 2 * C4 * rod_sq) / del_sq;

  return td;
}


/* ---------------------------------------------------------------------- */
// derivartive of dihedral angle func gij w.r.t rho, and atom positions
double PairDRIP::dihedral(
    const int i, const int j, Param& p, double const rhosq, double& d_drhosq,
    double* const d_dri, double* const d_drj,
    double* const d_drk1, double* const d_drk2, double* const d_drk3,
    double* const d_drl1, double* const d_drl2, double* const d_drl3)
{
  // get parameter
  double B = p.B;
  double eta = p.eta;
  double cut_rhosq = p.rhocutsq;

  // local vars
  double cos_kl[3][3];          // cos_omega_k1ijl1, cos_omega_k1ijl2 ...
  double d_dcos_kl[3][3];       // deriv of dihedral w.r.t to cos_omega_kijl
  double dcos_kl[3][3][4][DIM]; // 4 indicates k, i, j, l, e.g. dcoskl[0][1][0] means
                                // dcos_omega_k1ijl2 / drk


  // if larger than cutoff of rho, return 0
  if (rhosq >= cut_rhosq) {
    d_drhosq = 0;
    for (int dim = 0; dim < DIM; dim++) {
      d_dri[dim] = 0;
      d_drj[dim] = 0;
      d_drk1[dim] = 0;
      d_drk2[dim] = 0;
      d_drk3[dim] = 0;
      d_drl1[dim] = 0;
      d_drl2[dim] = 0;
      d_drl3[dim] = 0;
    }
    double dihe = 0.0;
    return dihe;
  }
  // 3 neighs of atoms i and j
  int k[3];
  int l[3];
  for (int m = 0; m < 3; m++) {
    k[m] = nearest3neigh[i][m];
    l[m] = nearest3neigh[j][m];
  }

  // cos_omega_kijl and the derivatives w.r.t coordinates
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      cos_kl[m][n] = deriv_cos_omega(
          coordinates[k[m]], coordinates[i], coordinates[j], coordinates[l[n]],
          dcos_kl[m][n][0], dcos_kl[m][n][1], dcos_kl[m][n][2], dcos_kl[m][n][3]);
    }
  }

  double epart1 = exp(-eta * cos_kl[0][0] * cos_kl[0][1] * cos_kl[0][2]);
  double epart2 = exp(-eta * cos_kl[1][0] * cos_kl[1][1] * cos_kl[1][2]);
  double epart3 = exp(-eta * cos_kl[2][0] * cos_kl[2][1] * cos_kl[2][2]);
  double D2 = epart1 + epart2 + epart3;

  // cutoff function
  double d_drhosq_tap;
  double D0 = B * tap_rho(rhosq, cut_rhosq, d_drhosq_tap);

  // dihedral energy
  double dihe = D0 * D2;

  // deriv of dihedral w.r.t rhosq
  d_drhosq = B * d_drhosq_tap * D2;

  // deriv of dihedral w.r.t cos_omega_kijl
  d_dcos_kl[0][0] = -D0 * epart1 * eta * cos_kl[0][1] * cos_kl[0][2];
  d_dcos_kl[0][1] = -D0 * epart1 * eta * cos_kl[0][0] * cos_kl[0][2];
  d_dcos_kl[0][2] = -D0 * epart1 * eta * cos_kl[0][0] * cos_kl[0][1];
  d_dcos_kl[1][0] = -D0 * epart2 * eta * cos_kl[1][1] * cos_kl[1][2];
  d_dcos_kl[1][1] = -D0 * epart2 * eta * cos_kl[1][0] * cos_kl[1][2];
  d_dcos_kl[1][2] = -D0 * epart2 * eta * cos_kl[1][0] * cos_kl[1][1];
  d_dcos_kl[2][0] = -D0 * epart3 * eta * cos_kl[2][1] * cos_kl[2][2];
  d_dcos_kl[2][1] = -D0 * epart3 * eta * cos_kl[2][0] * cos_kl[2][2];
  d_dcos_kl[2][2] = -D0 * epart3 * eta * cos_kl[2][0] * cos_kl[2][1];

  // initialization to be zero and later add values
  for (int dim = 0; dim < DIM; dim++) {
    d_drk1[dim] = 0.;
    d_drk2[dim] = 0.;
    d_drk3[dim] = 0.;
    d_dri[dim] = 0.;
    d_drj[dim] = 0.;
    d_drl1[dim] = 0.;
    d_drl2[dim] = 0.;
    d_drl3[dim] = 0.;
  }

  for (int m = 0; m < 3; m++) {
    for (int dim = 0; dim < 3; dim++) {
      d_drk1[dim] += d_dcos_kl[0][m] * dcos_kl[0][m][0][dim];
      d_drk2[dim] += d_dcos_kl[1][m] * dcos_kl[1][m][0][dim];
      d_drk3[dim] += d_dcos_kl[2][m] * dcos_kl[2][m][0][dim];
      d_drl1[dim] += d_dcos_kl[m][0] * dcos_kl[m][0][3][dim];
      d_drl2[dim] += d_dcos_kl[m][1] * dcos_kl[m][1][3][dim];
      d_drl3[dim] += d_dcos_kl[m][2] * dcos_kl[m][2][3][dim];
    }
    for (int n = 0; n < 3; n++) {
      for (int dim = 0; dim < 3; dim++) {
        d_dri[dim] += d_dcos_kl[m][n] * dcos_kl[m][n][1][dim];
        d_drj[dim] += d_dcos_kl[m][n] * dcos_kl[m][n][2][dim];
      }
    }
  }

  return dihe;
}
























/* ---------------------------------------------------------------------- */

void PairDRIP::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,k,l,kk,ll;
  tagint itag,jtag;
  double prodnorm1,prodnorm2,fkcx,fkcy,fkcz;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,fpair1,fpair2;
  double rsq,r,rhosq1,rhosq2,exp0,exp1,exp2,r2inv,r6inv,r8inv,Tap,dTap,Vkc;
  double frho1,frho2,sumC1,sumC2,sumC11,sumC22,sumCff,fsum,rdsq1,rdsq2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *DRIP_neighs_i,*DRIP_neighs_j;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double dprodnorm2[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fp2[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double fprod2[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fl[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double delli[3] = {0.0, 0.0, 0.0};

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // Build full neighbor list
  DRIP_neigh();
  // Calculate the normals
  calc_normal();

  // communicate the normal vector and its derivatives
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      jtag = tag[j];

      // two-body interactions from full neighbor list, skip half of them
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // only include the interation between different layers
      if (rsq < cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param& p = params[iparam_ij];

        r = sqrt(rsq);
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        r8inv = r2inv*r6inv;
	// turn on/off taper function
	if (tap_flag) {
	  Tap = calc_Tap(r,sqrt(cutsq[itype][jtype]));
	  dTap = calc_dTap(r,sqrt(cutsq[itype][jtype]));
	} else {Tap = 1.0; dTap = 0.0;}

        // Calculate the transverse distance
        // note that rho_ij does not equal to rho_ji except when normals are all along z
        prodnorm1 = normal[i][0]*delx + normal[i][1]*dely + normal[i][2]*delz;
        prodnorm2 = normal[j][0]*delx + normal[j][1]*dely + normal[j][2]*delz;
        rhosq1 = rsq - prodnorm1*prodnorm1;  // rho_ij
        rhosq2 = rsq - prodnorm2*prodnorm2;  // rho_ji
        rdsq1 = rhosq1*p.delta2inv; // (rho_ij/delta)^2
        rdsq2 = rhosq2*p.delta2inv; // (rho_ji/delta)^2

        // store exponents
        exp0 = exp(-p.lambda*(r-p.z0));
        exp1 = exp(-rdsq1);
        exp2 = exp(-rdsq2);

        sumC1 = p.C0 + p.C2*rdsq1 + p.C4*rdsq1*rdsq1;
        sumC2 = p.C0 + p.C2*rdsq2 + p.C4*rdsq2*rdsq2;
        sumC11 = (p.C2 + 2.0*p.C4*rdsq1)*p.delta2inv;
        sumC22 = (p.C2 + 2.0*p.C4*rdsq2)*p.delta2inv;
        frho1 = exp1*sumC1;
        frho2 = exp2*sumC2;
        sumCff = p.C + frho1 + frho2;
	Vkc = -p.A*p.z06*r6inv + exp0*sumCff;

        // derivatives
        fpair = -6.0*p.A*p.z06*r8inv + p.lambda*exp0/r*sumCff;
        fpair1 = 2.0*exp0*exp1*(p.delta2inv*sumC1 - sumC11);
        fpair2 = 2.0*exp0*exp2*(p.delta2inv*sumC2 - sumC22);
        fsum = fpair + fpair1 + fpair2;
        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0][i]*delx + dnormdri[1][0][i]*dely + dnormdri[2][0][i]*delz;
        dprodnorm1[1] = dnormdri[0][1][i]*delx + dnormdri[1][1][i]*dely + dnormdri[2][1][i]*delz;
        dprodnorm1[2] = dnormdri[0][2][i]*delx + dnormdri[1][2][i]*dely + dnormdri[2][2][i]*delz;
        // derivatives of the product of rji and nj, the result is a vector
        dprodnorm2[0] = dnormdri[0][0][j]*delx + dnormdri[1][0][j]*dely + dnormdri[2][0][j]*delz;
        dprodnorm2[1] = dnormdri[0][1][j]*delx + dnormdri[1][1][j]*dely + dnormdri[2][1][j]*delz;
        dprodnorm2[2] = dnormdri[0][2][j]*delx + dnormdri[1][2][j]*dely + dnormdri[2][2][j]*delz;
        fp1[0] = prodnorm1*normal[i][0]*fpair1;
        fp1[1] = prodnorm1*normal[i][1]*fpair1;
        fp1[2] = prodnorm1*normal[i][2]*fpair1;
        fp2[0] = prodnorm2*normal[j][0]*fpair2;
        fp2[1] = prodnorm2*normal[j][1]*fpair2;
        fp2[2] = prodnorm2*normal[j][2]*fpair2;
        fprod1[0] = prodnorm1*dprodnorm1[0]*fpair1;
        fprod1[1] = prodnorm1*dprodnorm1[1]*fpair1;
        fprod1[2] = prodnorm1*dprodnorm1[2]*fpair1;
        fprod2[0] = prodnorm2*dprodnorm2[0]*fpair2;
        fprod2[1] = prodnorm2*dprodnorm2[1]*fpair2;
        fprod2[2] = prodnorm2*dprodnorm2[2]*fpair2;
        fkcx = (delx*fsum - fp1[0] - fp2[0])*Tap - Vkc*dTap*delx/r;
        fkcy = (dely*fsum - fp1[1] - fp2[1])*Tap - Vkc*dTap*dely/r;
        fkcz = (delz*fsum - fp1[2] - fp2[2])*Tap - Vkc*dTap*delz/r;

        f[i][0] += fkcx - fprod1[0]*Tap;
        f[i][1] += fkcy - fprod1[1]*Tap;
        f[i][2] += fkcz - fprod1[2]*Tap;
        f[j][0] -= fkcx + fprod2[0]*Tap;
        f[j][1] -= fkcy + fprod2[1]*Tap;
        f[j][2] -= fkcz + fprod2[2]*Tap;

	// calculate the forces acted on the neighbors of atom i from atom j
	DRIP_neighs_i = DRIP_firstneigh[i];
  	for (kk = 0; kk < DRIP_numneigh[i]; kk++) {
	  k = DRIP_neighs_i[kk];
          if (k == i) continue;
          // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
          dprodnorm1[0] = dnormal[0][0][kk][i]*delx + dnormal[1][0][kk][i]*dely + dnormal[2][0][kk][i]*delz;
          dprodnorm1[1] = dnormal[0][1][kk][i]*delx + dnormal[1][1][kk][i]*dely + dnormal[2][1][kk][i]*delz;
          dprodnorm1[2] = dnormal[0][2][kk][i]*delx + dnormal[1][2][kk][i]*dely + dnormal[2][2][kk][i]*delz;
          fk[0] = (-prodnorm1*dprodnorm1[0]*fpair1)*Tap;
          fk[1] = (-prodnorm1*dprodnorm1[1]*fpair1)*Tap;
          fk[2] = (-prodnorm1*dprodnorm1[2]*fpair1)*Tap;
          f[k][0] += fk[0];
          f[k][1] += fk[1];
          f[k][2] += fk[2];
          delkj[0] = x[k][0] - x[j][0];
          delkj[1] = x[k][1] - x[j][1];
          delkj[2] = x[k][2] - x[j][2];
          if (evflag) ev_tally_xyz(k,j,nlocal,newton_pair,0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
	}

	// calculate the forces acted on the neighbors of atom j from atom i
	DRIP_neighs_j = DRIP_firstneigh[j];
  	for (ll = 0; ll < DRIP_numneigh[j]; ll++) {
	  l = DRIP_neighs_j[ll];
          if (l == j) continue;
          // derivatives of the product of rji and nj respect to rl, l=0,1,2, where atom l is the neighbors of atom j
          dprodnorm2[0] = dnormal[0][0][ll][j]*delx + dnormal[1][0][ll][j]*dely + dnormal[2][0][ll][j]*delz;
          dprodnorm2[1] = dnormal[0][1][ll][j]*delx + dnormal[1][1][ll][j]*dely + dnormal[2][1][ll][j]*delz;
          dprodnorm2[2] = dnormal[0][2][ll][j]*delx + dnormal[1][2][ll][j]*dely + dnormal[2][2][ll][j]*delz;
          fl[0] = (-prodnorm2*dprodnorm2[0]*fpair2)*Tap;
          fl[1] = (-prodnorm2*dprodnorm2[1]*fpair2)*Tap;
          fl[2] = (-prodnorm2*dprodnorm2[2]*fpair2)*Tap;
          f[l][0] += fl[0];
          f[l][1] += fl[1];
          f[l][2] += fl[2];
          delli[0] = x[l][0] - x[i][0];
          delli[1] = x[l][1] - x[i][1];
          delli[2] = x[l][2] - x[i][2];
          if (evflag) ev_tally_xyz(l,i,nlocal,newton_pair,0.0,0.0,fl[0],fl[1],fl[2],delli[0],delli[1],delli[2]);
	}

        if (eflag) {
          if (tap_flag) evdwl = Tap*Vkc;
          else  evdwl = Vkc;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0,fkcx,fkcy,fkcz,delx,dely,delz);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   Calculate the normals for each atom
------------------------------------------------------------------------- */
void PairDRIP::calc_normal()
{
  int i,j,ii,jj,inum,jnum;
  int cont,id,ip,m;
  double nn,xtp,ytp,ztp,delx,dely,delz,nn2;
  int *ilist,*jlist;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];

  double **x = atom->x;

  // grow normal array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(normal);
    memory->destroy(dnormal);
    memory->destroy(dnormdri);
    nmax = atom->nmax;
    memory->create(normal,nmax,3,"DRIP:normal");
    memory->create(dnormdri,3,3,nmax,"DRIP:dnormdri");
    memory->create(dnormal,3,3,3,nmax,"DRIP:dnormal");
  }

  inum = list->inum;
  ilist = list->ilist;
  //Calculate normals
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtp = x[i][0];
    ytp = x[i][1];
    ztp = x[i][2];

    //   Initialize the arrays
    for (id = 0; id < 3; id++){
      pv12[id] = 0.0;
      pv31[id] = 0.0;
      pv23[id] = 0.0;
      n1[id] = 0.0;
      dni[id] = 0.0;
      normal[i][id] = 0.0;
      for (ip = 0; ip < 3; ip++){
        vet[ip][id] = 0.0;
        dnn[ip][id] = 0.0;
        dpvdri[ip][id] = 0.0;
        dnormdri[ip][id][i] = 0.0;
        for (m = 0; m < 3; m++){
          dpv12[ip][id][m] = 0.0;
          dpv31[ip][id][m] = 0.0;
          dpv23[ip][id][m] = 0.0;
          dn1[ip][id][m] = 0.0;
          dnormal[ip][id][m][i] = 0.0;
        }
      }
    }

    cont = 0;
    jlist = DRIP_firstneigh[i];
    jnum = DRIP_numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = x[j][0] - xtp;
      dely = x[j][1] - ytp;
      delz = x[j][2] - ztp;
      vet[cont][0] = delx;
      vet[cont][1] = dely;
      vet[cont][2] = delz;
      cont++;
    }

    if (cont <= 1) {
      normal[i][0] = 0.0;
      normal[i][1] = 0.0;
      normal[i][2] = 1.0;
      // derivatives of normal vector is zero
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dnormdri[id][ip][i] = 0.0;
          for (m = 0; m < 3; m++){
            dnormal[id][ip][m][i] = 0.0;
          }
        }
      }
    }
    else if (cont == 2) {
      // for the atoms at the edge who has only two neighbor atoms
      pv12[0] = vet[0][1]*vet[1][2] - vet[1][1]*vet[0][2];
      pv12[1] = vet[0][2]*vet[1][0] - vet[1][2]*vet[0][0];
      pv12[2] = vet[0][0]*vet[1][1] - vet[1][0]*vet[0][1];
      dpvdri[0][0] = 0.0;
      dpvdri[0][1] = vet[0][2]-vet[1][2];
      dpvdri[0][2] = vet[1][1]-vet[0][1];
      dpvdri[1][0] = vet[1][2]-vet[0][2];
      dpvdri[1][1] = 0.0;
      dpvdri[1][2] = vet[0][0]-vet[1][0];
      dpvdri[2][0] = vet[0][1]-vet[1][1];
      dpvdri[2][1] = vet[1][0]-vet[0][0];
      dpvdri[2][2] = 0.0;

      // derivatives respect to the first neighbor, atom k
      dpv12[0][0][0] =  0.0;
      dpv12[0][1][0] =  vet[1][2];
      dpv12[0][2][0] = -vet[1][1];
      dpv12[1][0][0] = -vet[1][2];
      dpv12[1][1][0] =  0.0;
      dpv12[1][2][0] =  vet[1][0];
      dpv12[2][0][0] =  vet[1][1];
      dpv12[2][1][0] = -vet[1][0];
      dpv12[2][2][0] =  0.0;

      // derivatives respect to the second neighbor, atom l
      dpv12[0][0][1] =  0.0;
      dpv12[0][1][1] = -vet[0][2];
      dpv12[0][2][1] =  vet[0][1];
      dpv12[1][0][1] =  vet[0][2];
      dpv12[1][1][1] =  0.0;
      dpv12[1][2][1] = -vet[0][0];
      dpv12[2][0][1] = -vet[0][1];
      dpv12[2][1][1] =  vet[0][0];
      dpv12[2][2][1] =  0.0;

      // derivatives respect to the third neighbor, atom n
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv12[id][ip][2] = 0.0;
        }
      }

      n1[0] = pv12[0];
      n1[1] = pv12[1];
      n1[2] = pv12[2];
      // the magnitude of the normal vector
      nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
      nn = sqrt(nn2);
      if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[i][0] = n1[0]/nn;
      normal[i][1] = n1[1]/nn;
      normal[i][2] = n1[2]/nn;
      // derivatives of nn, dnn:3x1 vector
      dni[0] = (n1[0]*dpvdri[0][0] + n1[1]*dpvdri[1][0] + n1[2]*dpvdri[2][0])/nn;
      dni[1] = (n1[0]*dpvdri[0][1] + n1[1]*dpvdri[1][1] + n1[2]*dpvdri[2][1])/nn;
      dni[2] = (n1[0]*dpvdri[0][2] + n1[1]*dpvdri[1][2] + n1[2]*dpvdri[2][2])/nn;
      // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dnormdri[id][ip][i] = dpvdri[id][ip]/nn - n1[id]*dni[ip]/nn2;
        }
      }

      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          for (m = 0; m < 3; m++){
            dn1[id][ip][m] = dpv12[id][ip][m];
          }
        }
      }
      // derivatives of nn, dnn:3x3 vector
      // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
      // r[id][m]: the id's component of atom m
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          dnn[id][m] = (n1[0]*dn1[0][id][m] + n1[1]*dn1[1][id][m] + n1[2]*dn1[2][id][m])/nn;
        }
      }
      // dnormal[id][ip][m][i]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
      // for atom m, which is a neighbor atom of atom i, m=0,jnum-1
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          for (ip = 0; ip < 3; ip++){
            dnormal[id][ip][m][i] = dn1[id][ip][m]/nn - n1[id]*dnn[ip][m]/nn2;
          }
        }
      }
    }

    else if(cont == 3) {
      // for the atoms at the edge who has only two neighbor atoms
      pv12[0] = vet[0][1]*vet[1][2] - vet[1][1]*vet[0][2];
      pv12[1] = vet[0][2]*vet[1][0] - vet[1][2]*vet[0][0];
      pv12[2] = vet[0][0]*vet[1][1] - vet[1][0]*vet[0][1];
      // derivatives respect to the first neighbor, atom k
      dpv12[0][0][0] =  0.0;
      dpv12[0][1][0] =  vet[1][2];
      dpv12[0][2][0] = -vet[1][1];
      dpv12[1][0][0] = -vet[1][2];
      dpv12[1][1][0] =  0.0;
      dpv12[1][2][0] =  vet[1][0];
      dpv12[2][0][0] =  vet[1][1];
      dpv12[2][1][0] = -vet[1][0];
      dpv12[2][2][0] =  0.0;
      // derivatives respect to the second neighbor, atom l
      dpv12[0][0][1] =  0.0;
      dpv12[0][1][1] = -vet[0][2];
      dpv12[0][2][1] =  vet[0][1];
      dpv12[1][0][1] =  vet[0][2];
      dpv12[1][1][1] =  0.0;
      dpv12[1][2][1] = -vet[0][0];
      dpv12[2][0][1] = -vet[0][1];
      dpv12[2][1][1] =  vet[0][0];
      dpv12[2][2][1] =  0.0;

      // derivatives respect to the third neighbor, atom n
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv12[id][ip][2] = 0.0;
        }
      }

      pv31[0] = vet[2][1]*vet[0][2] - vet[0][1]*vet[2][2];
      pv31[1] = vet[2][2]*vet[0][0] - vet[0][2]*vet[2][0];
      pv31[2] = vet[2][0]*vet[0][1] - vet[0][0]*vet[2][1];
      // derivatives respect to the first neighbor, atom k
      dpv31[0][0][0] =  0.0;
      dpv31[0][1][0] = -vet[2][2];
      dpv31[0][2][0] =  vet[2][1];
      dpv31[1][0][0] =  vet[2][2];
      dpv31[1][1][0] =  0.0;
      dpv31[1][2][0] = -vet[2][0];
      dpv31[2][0][0] = -vet[2][1];
      dpv31[2][1][0] =  vet[2][0];
      dpv31[2][2][0] =  0.0;
      // derivatives respect to the third neighbor, atom n
      dpv31[0][0][2] =  0.0;
      dpv31[0][1][2] =  vet[0][2];
      dpv31[0][2][2] = -vet[0][1];
      // derivatives of pv13[1] to rn
      dpv31[1][0][2] = -vet[0][2];
      dpv31[1][1][2] =  0.0;
      dpv31[1][2][2] =  vet[0][0];
      // derivatives of pv13[2] to rn
      dpv31[2][0][2] =  vet[0][1];
      dpv31[2][1][2] = -vet[0][0];
      dpv31[2][2][2] =  0.0;

      // derivatives respect to the second neighbor, atom l
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv31[id][ip][1] = 0.0;
        }
      }

      pv23[0] = vet[1][1]*vet[2][2] - vet[2][1]*vet[1][2];
      pv23[1] = vet[1][2]*vet[2][0] - vet[2][2]*vet[1][0];
      pv23[2] = vet[1][0]*vet[2][1] - vet[2][0]*vet[1][1];
      // derivatives respect to the second neighbor, atom k
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv23[id][ip][0] = 0.0;
        }
      }
      // derivatives respect to the second neighbor, atom l
      dpv23[0][0][1] =  0.0;
      dpv23[0][1][1] =  vet[2][2];
      dpv23[0][2][1] = -vet[2][1];
      dpv23[1][0][1] = -vet[2][2];
      dpv23[1][1][1] =  0.0;
      dpv23[1][2][1] =  vet[2][0];
      dpv23[2][0][1] =  vet[2][1];
      dpv23[2][1][1] = -vet[2][0];
      dpv23[2][2][1] =  0.0;
      // derivatives respect to the third neighbor, atom n
      dpv23[0][0][2] =  0.0;
      dpv23[0][1][2] = -vet[1][2];
      dpv23[0][2][2] =  vet[1][1];
      dpv23[1][0][2] =  vet[1][2];
      dpv23[1][1][2] =  0.0;
      dpv23[1][2][2] = -vet[1][0];
      dpv23[2][0][2] = -vet[1][1];
      dpv23[2][1][2] =  vet[1][0];
      dpv23[2][2][2] =  0.0;

//############################################################################################
      // average the normal vectors by using the 3 neighboring planes
      n1[0] = (pv12[0] + pv31[0] + pv23[0])/cont;
      n1[1] = (pv12[1] + pv31[1] + pv23[1])/cont;
      n1[2] = (pv12[2] + pv31[2] + pv23[2])/cont;
      // the magnitude of the normal vector
      nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
      nn = sqrt(nn2);
      if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[i][0] = n1[0]/nn;
      normal[i][1] = n1[1]/nn;
      normal[i][2] = n1[2]/nn;

      // for the central atoms, dnormdri is always zero
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dnormdri[id][ip][i] = 0.0;
        }
      } // end of derivatives of normals respect to atom i

      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          for (m = 0; m < 3; m++){
            dn1[id][ip][m] = (dpv12[id][ip][m] + dpv23[id][ip][m] + dpv31[id][ip][m])/cont;
          }
        }
      }
      // derivatives of nn, dnn:3x3 vector
      // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
      // r[id][m]: the id's component of atom m
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          dnn[id][m] = (n1[0]*dn1[0][id][m] + n1[1]*dn1[1][id][m] + n1[2]*dn1[2][id][m])/nn;
        }
      }
      // dnormal[id][ip][m][i]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
      // for atom m, which is a neighbor atom of atom i, m=0,jnum-1
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          for (ip = 0; ip < 3; ip++){
            dnormal[id][ip][m][i] = dn1[id][ip][m]/nn - n1[id]*dnn[ip][m]/nn2;
          }
        }
      }
    }
    else {
      error->one(FLERR,"There are too many neighbors for calculating normals");
    }
  }

}


/* ----------------------------------------------------------------------
 create neighbor list from main neighbor list for calculating the normals
------------------------------------------------------------------------- */

void PairDRIP::DRIP_neigh()
{
  int i,j,ii,jj,n,allnum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(DRIP_numneigh);
    memory->sfree(DRIP_firstneigh);
    memory->create(DRIP_numneigh,maxlocal,"DRIP:numneigh");
    DRIP_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                             "DRIP:firstneigh");
  }

  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all DRIP neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;


      int iparam_ij = elem2param[itype][jtype];
      double rcutsq = params[iparam_ij].rcutsq;

      if (rsq != 0 && rsq < rcutsq && atom->molecule[i] == atom->molecule[j]) {
        neighptr[n++] = j;
      }
    }

    DRIP_firstneigh[i] = neighptr;
    DRIP_numneigh[i] = n;
    if (n > 3) error->one(FLERR,"There are too many neighbors for some atoms, please check your configuration");
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDRIP::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style drip requires newton pair on");
  if (!atom->molecule_flag)
    error->all(FLERR,"Pair style drip requires atom attribute molecule");

  // need a full neighbor list, including neighbors of ghosts

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  // local DRIP neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage= comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom,pgsize,PGDELTA);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDRIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  // MOVE init of setflag ot other places; se pair_sw
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[atom->ntypes+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDRIP::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
  if (strcmp(force->pair_style,"hybrid/overlay")!=0)
    error->all(FLERR,"ERROR: requires hybrid/overlay pair_style");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDRIP::coeff(int narg, char **arg)
{
  int i,j,n;

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }


  read_file(arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDRIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ----------------------------------------------------------------------
   read DRIP file
------------------------------------------------------------------------- */

void PairDRIP::read_file(char *filename)
{
  int params_per_line = 14;
  char **words = new char*[params_per_line+1];
  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open DRIP potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int i,j,n,m,nwords,ielement,jelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Insufficient format in DRIP potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement = 1st args
    // if these 2 args are in element list, then parse this line
    // else skip to next line (continue)

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].C0       = atof(words[2]);
    params[nparams].C2       = atof(words[3]);
    params[nparams].C4       = atof(words[4]);
    params[nparams].C        = atof(words[5]);
    params[nparams].delta    = atof(words[6]);
    params[nparams].lambda   = atof(words[7]);
    params[nparams].A        = atof(words[8]);
    params[nparams].z0       = atof(words[9]);
    params[nparams].B        = atof(words[10]);
    params[nparams].eta      = atof(words[11]);
    params[nparams].rhocut   = atof(words[12]);
    params[nparams].rcut     = atof(words[13]);

    // convenient precomputations
    params[nparams].rhocutsq = params[nparams].rhocut * params[nparams].rhocut;
    params[nparams].rcutsq   = params[nparams].rcut * params[nparams].rcut;

    // set max cutoff
    if(params[nparams].rcut > cutmax) cutmax = params[nparams].rcut;


    nparams++;
    //if(nparams >= pow(atom->ntypes,3)) break;
  }

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      n = -1;
      for (m = 0; m < nparams; m++) {
        if (i == params[m].ielement && j == params[m].jelement) {
          if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
          n = m;
        }
      }
      if (n < 0) error->all(FLERR,"Potential file is missing an entry");
      elem2param[i][j] = n;
    }
  }
  delete [] words;
}

/* ---------------------------------------------------------------------- */

int PairDRIP::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m,l,ip,id;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = normal[j][0];
    buf[m++] = normal[j][1];
    buf[m++] = normal[j][2];
    buf[m++] = dnormdri[0][0][j];
    buf[m++] = dnormdri[0][1][j];
    buf[m++] = dnormdri[0][2][j];
    buf[m++] = dnormdri[1][0][j];
    buf[m++] = dnormdri[1][1][j];
    buf[m++] = dnormdri[1][2][j];
    buf[m++] = dnormdri[2][0][j];
    buf[m++] = dnormdri[2][1][j];
    buf[m++] = dnormdri[2][2][j];
    for (l = 0; l < 3; l++){
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
	  buf[m++] = dnormal[id][ip][l][j];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairDRIP::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last,l,ip,id;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    normal[i][0] = buf[m++];
    normal[i][1] = buf[m++];
    normal[i][2] = buf[m++];
    dnormdri[0][0][i] = buf[m++];
    dnormdri[0][1][i] = buf[m++];
    dnormdri[0][2][i] = buf[m++];
    dnormdri[1][0][i] = buf[m++];
    dnormdri[1][1][i] = buf[m++];
    dnormdri[1][2][i] = buf[m++];
    dnormdri[2][0][i] = buf[m++];
    dnormdri[2][1][i] = buf[m++];
    dnormdri[2][2][i] = buf[m++];
    for (l = 0; l < 3; l++){
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
	  dnormal[id][ip][l][i] = buf[m++];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
