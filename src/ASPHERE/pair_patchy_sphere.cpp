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
   Contributing author: Stefan Paquay (Brandeis University)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_patchy_sphere.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

constexpr const bool nine_pairs = false;

/* ---------------------------------------------------------------------- */

PairPatchySphere::PairPatchySphere(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  // writedata = 1; // Add later.
  allocated = 0;
  setflag = NULL;

}

/* ---------------------------------------------------------------------- */


PairPatchySphere::~PairPatchySphere()
{
  if (allocated) {
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(inner_cut);

    memory->destroy(offset);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(epsilon_b);
    memory->destroy(phi_m);
    memory->destroy(theta_m);

    memory->destroy(bond_pairs);
    memory->destroy(bond_vectors);

    memory->destroy(setflag);

  }
}

/* ---------------------------------------------------------------------- */


void PairPatchySphere::allocate()
{
  int n = atom->ntypes;

  memory->create(cutsq, n+1,n+1, "pair:cutsq");
  memory->create(cut, n+1,n+1, "pair:cut");
  memory->create(inner_cut, n+1,n+1, "pair:inner_cut");
  memory->create(offset, n+1,n+1, "pair:offset");
  memory->create(epsilon, n+1,n+1, "pair:epsilon");
  memory->create(sigma, n+1,n+1, "pair:sigma");
  memory->create(epsilon_b, n+1,n+1, "pair:epsilon_b");
  memory->create(phi_m, n+1, n+1, "pair:phi_m");
  memory->create(theta_m, n+1, n+1, "pair:theta_m");

  // For now the complementary pairs are fixed.
  n_bonds = 3;
  if (nine_pairs){
    n_bond_pairs = 9;
  } else {
    n_bond_pairs = 3;
  }
  memory->create(bond_pairs,   n_bond_pairs, 6, "pair:bond_pairs");
  memory->create(bond_vectors, n_bonds,      3, "pair:bond_vectors");

  if( comm->me == 0 ) fprintf(screen, "n is %d\n", n);
  memory->create<int>(setflag, n+1, n+1, "pair:setflag");


  allocated = 1;
}


/* ---------------------------------------------------------------------- */


void PairPatchySphere::settings(int narg, char **arg)
{
  // params: phi_m, theta_m, cut-off
  if (narg != 1 && narg != 2) error->all(FLERR,"Illegal pair_style command");

  if (debug_output && comm->me == 0) {
    fprintf(screen, "Got %d args:", narg);
    for (int i = 0; i < narg; ++i) {
      fprintf(screen, " %s", arg[i]);
    }
    fprintf(screen, "\n");
  }
  cut_global = force->numeric(FLERR, arg[0]);

  if (narg == 2)
    debug_output = force->inumeric(FLERR, arg[1]);

  if (allocated){
    for (int i = 1; i < atom->ntypes; ++i){
      for (int j = 1; j <= atom->ntypes; ++j){
        cut[i][j] = cut_global;
        cutsq[i][j] = cut[i][j]*cut[i][j];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairPatchySphere::coeff(int narg, char **arg)
{
  // Coeffs are epsilon, sigma, epsilon_b, phi_m, theta_m
  if (narg < 7 || narg > 8)
    error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();

  double cut_one = cut_global;
  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one   = force->numeric(FLERR, arg[2]);
  double sigma_one     = force->numeric(FLERR, arg[3]);
  double epsilon_b_one = force->numeric(FLERR, arg[4]);
  double phi_m_one     = force->numeric(FLERR, arg[5]);
  double theta_m_one   = force->numeric(FLERR, arg[6]);

  double rc_one = sigma_one * pow(2.0,1.0/6.0);

  if (narg == 8) cut_one = force->numeric(FLERR, arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; ++i){
    for (int j = jlo; j <= jhi; ++j ){
      epsilon[i][j]   = epsilon_one;
      sigma[i][j]     = sigma_one;
      epsilon_b[i][j] = epsilon_b_one;
      phi_m[i][j]     = phi_m_one;
      theta_m[i][j]   = theta_m_one;
      cut[i][j]       = cut_one;
      cutsq[i][j]     = cut_one*cut_one;
      inner_cut[i][j] = rc_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

void PairPatchySphere::init_style()
{
  avec = static_cast<AtomVecEllipsoid*>(atom->style_match("ellipsoid"));
  if (!avec) error->all(FLERR, "Pair patchy_sphere requires atom style ellipsoid");

  neighbor->request(this,instance_me);

  // Set up the bond vectors here:
  bond_pairs[0][0] = 0;
  bond_pairs[0][1] = 1;
  bond_pairs[0][2] = 1;
  bond_pairs[0][3] = 0;
  bond_pairs[0][4] = 2;
  bond_pairs[0][5] = 2;

  bond_pairs[1][0] = 1;
  bond_pairs[1][1] = 0;
  bond_pairs[1][2] = 0;
  bond_pairs[1][3] = 1;
  bond_pairs[1][4] = 2;
  bond_pairs[1][5] = 2;

  bond_pairs[2][0] = 2;
  bond_pairs[2][1] = 2;
  bond_pairs[2][2] = 0;
  bond_pairs[2][3] = 1;
  bond_pairs[2][4] = 1;
  bond_pairs[2][5] = 0;

  if (nine_pairs) {
    bond_pairs[3][0] = 0;
    bond_pairs[3][1] = 2;
    bond_pairs[3][2] = 2;
    bond_pairs[3][3] = 0;
    bond_pairs[3][4] = 1;
    bond_pairs[3][5] = 1;

    bond_pairs[4][0] = 1;
    bond_pairs[4][1] = 1;
    bond_pairs[4][2] = 2;
    bond_pairs[4][3] = 0;
    bond_pairs[4][4] = 0;
    bond_pairs[4][5] = 2;

    bond_pairs[5][0] = 2;
    bond_pairs[5][1] = 0;
    bond_pairs[5][2] = 0;
    bond_pairs[5][3] = 2;
    bond_pairs[5][4] = 1;
    bond_pairs[5][5] = 1;

    bond_pairs[6][0] = 0;
    bond_pairs[6][1] = 0;
    bond_pairs[6][2] = 1;
    bond_pairs[6][3] = 2;
    bond_pairs[6][4] = 2;
    bond_pairs[6][5] = 1;

    bond_pairs[7][0] = 1;
    bond_pairs[7][1] = 2;
    bond_pairs[7][2] = 0;
    bond_pairs[7][3] = 0;
    bond_pairs[7][4] = 2;
    bond_pairs[7][5] = 1;

    bond_pairs[8][0] = 2;
    bond_pairs[8][1] = 1;
    bond_pairs[8][2] = 0;
    bond_pairs[8][3] = 0;
    bond_pairs[8][4] = 1;
    bond_pairs[8][5] = 2;
  }

  bond_vectors[2][0] = 1.0;
  bond_vectors[2][1] = 0.0;
  bond_vectors[2][2] = 0.0;

  bond_vectors[1][0] = -0.309017;
  bond_vectors[1][1] = 0.951057;
  bond_vectors[1][2] = 0.0;

  bond_vectors[0][0] = -0.309017;
  bond_vectors[0][1] = -0.425325;
  bond_vectors[0][2] =  0.850651;

  for (int i = 0; i < n_bonds; ++i) {
    double bx  = bond_vectors[i][0];
    double by  = bond_vectors[i][1];
    double bz  = bond_vectors[i][2];

    double nb2 = bx*bx + by*by + bz*bz;
    if ( fabs(nb2 - 1.0) > 1e-4 ){
      char msg[2048];
      snprintf(msg, 2048, "Bond vector %d is not of unit length! norm^2 - 1.0 = %g",
               i+1, nb2 - 1.0);
      error->all(FLERR,msg);
    }
  }

}

/* ---------------------------------------------------------------------- */

double PairPatchySphere::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, "Not all patchy sphere coefficients set!");

  epsilon[j][i]   = epsilon[i][j];
  sigma[j][i]     = sigma[i][j];
  epsilon_b[j][i] = epsilon_b[i][j];
  phi_m[j][i]     = phi_m[i][j];
  theta_m[j][i]   = theta_m[i][j];
  cut[j][i]       = cut[i][j];
  setflag[j][i]   = setflag[i][j];
  inner_cut[j][i] = inner_cut[i][j];


  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairPatchySphere::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double evdwl, one_eng, rsq, r2inv, r6inv, forcelj, factor_lj;
  double fforce[3], ttor[3], rtor[3], r12[3];
  double a1[3][3], a2[3][3];
  int *ilist, *jlist, *numneigh, **firstneigh;
  double *iquat, *jquat;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  // For now always assume you are dealing with patchy <--> patchy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    iquat = bonus[ellipsoid[i]].quat;
    MathExtra::quat_to_mat_trans(iquat,a1); // a1 is rotation matrix
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {

        jquat = bonus[ellipsoid[j]].quat;
        MathExtra::quat_to_mat_trans(jquat,a2);
        one_eng = bond_vec_analytic(i,j,a1,a2,r12,rsq,
                                    fforce,ttor,rtor);


        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
        tor[i][1] += ttor[1];
        tor[i][2] += ttor[2];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
          tor[j][1] += rtor[1];
          tor[j][2] += rtor[2];
        }

        if (eflag) evdwl = factor_lj*one_eng;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                 -r12[0],-r12[1],-r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ------------------------------------------------------------ */

double PairPatchySphere::bond_vec_analytic (const int i, const int j,
                                            double a1[3][3], double a2[3][3],
                                            double *r12, const double rsq,
                                            double *fforce, double *ttor,
                                            double *rtor)
{
  // At this point, rsq and r12 contain the distance vector and distance
  // squared, a1 and a2 contain the rotation matrices.
  // We need to calculate fforce, ttor and rtor.

  // For each bond vector pair:
  int *type = atom->type;
  int itype = type[i];
  int jtype = type[j];
  double s = sigma[itype][jtype];
  double sigma2 = s*s;
  double inner_rc = inner_cut[itype][jtype];
  double inner_rc2 = inner_rc*inner_rc;
  double energy = 0.0;

  ttor[0] = ttor[1] = ttor[2] = 0.0;
  rtor[0] = rtor[1] = rtor[2] = 0.0;
  fforce[0] = fforce[1] = fforce[2] = 0.0;

  // Distance criterion for repulsive LJ part:
  if (rsq < inner_rc2) {
    double rinv2 = 1.0 / rsq;
    double s6 = sigma2*sigma2*sigma2;
    double r6 = rinv2*rinv2*rinv2;
    s6 *= r6;

    double force_lj = 48.0*epsilon[itype][jtype]*s6*(s6-0.5);
    force_lj *= -rinv2;

    // Add force to the force vector:
    fforce[0] += r12[0] * force_lj;
    fforce[1] += r12[1] * force_lj;
    fforce[2] += r12[2] * force_lj;

    energy += 4.0*epsilon[itype][jtype]*( s6*(s6-1.0) + 0.25 );
  }

  // Calculate the rotated vectors here:
  double b_rot[3*6];

  MathExtra::matvec(a1, bond_vectors[0], b_rot);
  MathExtra::matvec(a1, bond_vectors[1], b_rot+3);
  MathExtra::matvec(a1, bond_vectors[2], b_rot+6);

  MathExtra::matvec(a2, bond_vectors[0], b_rot+9);
  MathExtra::matvec(a2, bond_vectors[1], b_rot+12);
  MathExtra::matvec(a2, bond_vectors[2], b_rot+15);

  // Contributions of the bond vectors:
  for (int pair_i = 0; pair_i < n_bond_pairs; ++pair_i){
    bool bail = false;
    // Calculate for each pair the attractive patch part:
    int alpha_i   = bond_pairs[pair_i][0];
    int beta_i    = bond_pairs[pair_i][1];
    int gamma_i   = bond_pairs[pair_i][2];
    int epsilon_i = bond_pairs[pair_i][3];
    int eta_i     = bond_pairs[pair_i][4];
    int nu_i      = bond_pairs[pair_i][5];

    /*
    if (debug_output && comm->me == 0 && i == 1) {
      fprintf( screen, "At i=1: pairs: %d %d %d %d %d %d\n", alpha_i, beta_i,
               gamma_i, epsilon_i, eta_i, nu_i );
    }
    */

    // Alpha, gamma and eta belong to particle i,
    // beta, epsilon and nu belong to particle j.
    double *b_alp = b_rot + 3*alpha_i;
    double *b_bet = b_rot + 3*beta_i + 9;
    double *b_gam = b_rot + 3*gamma_i;
    double *b_eps = b_rot + 3*epsilon_i + 9;
    double *b_eta = b_rot + 3*eta_i;
    double *b_nu  = b_rot + 3*nu_i + 9;

    // Distance between bond vectors, attractive force:
    double patch_dist[3];
    const double patch_b = 0.56123102415468651;

    patch_dist[0] = r12[0] + patch_b*(b_alp[0] - b_bet[0]);
    patch_dist[1] = r12[1] + patch_b*(b_alp[1] - b_bet[1]);
    patch_dist[2] = r12[2] + patch_b*(b_alp[2] - b_bet[2]);
    double patch_r2 = patch_dist[0]*patch_dist[0] + patch_dist[1]*patch_dist[1]
	    + patch_dist[2]*patch_dist[2];
    double criterion = cut[itype][jtype] - s*1.122462048309373;
    if (debug_output) {
      fprintf(screen,"Distance criterion is %g\n", criterion);
    }
    double crit2 = criterion*criterion;

    if (patch_r2 >= crit2) {
      // The attraction is out of range. Ignore this contribution.
      bail = true;
    }

    double patch_r = sqrt(patch_r2);
    if (debug_output) {
      fprintf(screen, "r_ij is %g (vector: %g %g %g)\n",
              sqrt(rsq), r12[0], r12[1], r12[2]);
      fprintf(screen, "b_alpha = (%g %g %g)\n",
              b_alp[0], b_alp[1], b_alp[2]);
      fprintf(screen, "b_beta  = (%g %g %g)\n",
              b_bet[0], b_bet[1], b_bet[2]);
      fprintf(screen, "Patch distance is %g (vector: %g %g %g)\n",
              patch_r, patch_dist[0], patch_dist[1], patch_dist[2]);
    }


    bool calc_theta = true;

    // The total force, torque and reverse torque from this bond vector pair
    // will be tallied in these:
    double acc_force[3] = { 0.0, 0.0, 0.0 };
    double acc_ttor[3]  = { 0.0, 0.0, 0.0 };
    double acc_rtor[3]  = { 0.0, 0.0, 0.0 };

    double u_theta = 2.0;
    double theta_force[3] = {0.0,0.0,0.0};
    double theta_ttor[3] = {0.0,0.0,0.0};
    double theta_rtor[3] = {0.0,0.0,0.0};

    if (calc_theta && !bail) {
      // Switching function: Theta
      u_theta = bond_vec_theta_part(i, j, itype, jtype, r12, b_alp, b_bet,
                                    theta_force, theta_ttor, theta_rtor);
      if (debug_output) {
        fprintf(screen, "u_theta = %g\n", u_theta);
      }
      if (u_theta == 0) {
        // No contribution from this bond vector set.
        bail = true;
        theta_force[0] = theta_force[1] = theta_force[2] = 0.0;
        theta_ttor[0] = theta_ttor[1] = theta_ttor[2] = 0.0;
        theta_rtor[0] = theta_rtor[1] = theta_rtor[2] = 0.0;
      }
    }

    if (debug_output)
      fprintf(screen, "u_theta = %g\n", u_theta);

    // Tally force and torque:
    bool calc_phi1 = true;
    bool calc_phi2 = true;

    // phi_1:
    double u_phi1 = 2.0;
    double phi1_force[3] = {0.0,0.0,0.0};
    double phi1_ttor[3]  = {0.0,0.0,0.0};
    double phi1_rtor[3]  = {0.0,0.0,0.0};

    if (!bail && calc_phi1) {
      u_phi1 = bond_vec_phi_part(i, j, itype, jtype, r12, rsq, b_alp, b_bet,
                                 b_gam, b_eps, phi1_force, phi1_ttor, phi1_rtor);
      if (u_phi1 == 0) {
        // No contribution from this bond vector pair.
        bail = true;
      }
    }else{
      phi1_force[0] = phi1_force[1] = phi1_force[2] = 0.0;
      phi1_ttor[0] = phi1_ttor[1] = phi1_ttor[2] = 0.0;
      phi1_rtor[0] = phi1_rtor[1] = phi1_rtor[2] = 0.0;
      // If do not calc, ignore its contributino to forces and torques
      // but put u_phi1 to 2.0 so that the phi factor does not change anything.
      u_phi1 = 2.0;
    }

    double phi2_force[3] = {0.0,0.0,0.0};
    double phi2_ttor[3]  = {0.0,0.0,0.0};
    double phi2_rtor[3]  = {0.0,0.0,0.0};

    double u_phi2 = 2.0;

    if (!bail && calc_phi2) {
      u_phi2 = bond_vec_phi_part(i, j, itype, jtype, r12, rsq, b_alp, b_bet,
                                 b_eta, b_nu, phi2_force, phi2_ttor, phi2_rtor);
      if (u_phi2 == 0) {
        // No contribution from this bond vector pair.
        bail = true;
      }
    }else{
      phi2_force[0] = phi2_force[1] = phi2_force[2] = 0.0;
      phi2_ttor[0] = phi2_ttor[1] = phi2_ttor[2] = 0.0;
      phi2_rtor[0] = phi2_rtor[1] = phi2_rtor[2] = 0.0;
      // If do not calc, ignore its contributino to forces and torques
      // but put u_phi1 to 2.0 so that the phi factor does not change anything.
      u_phi2 = 2.0;
    }

    // attractive LJ force. Has weird shape of 2^{1/6} + patch_dist.
    double rr  = patch_r + inner_cut[itype][jtype];
    double sr2 = sigma2 / (rr*rr);
    double sr6 = sr2*sr2*sr2;

    // 0.5 instead of 4 because we omit the factors 0.5 for switching functions.
    double bond_eps = 0.5*epsilon_b[itype][jtype];
    // Add the offset.
    double src2 = sigma2 / cutsq[itype][jtype];
    double src6 = src2*src2*src2;
    double u_LJ = bond_eps * ( sr6 * ( sr6 - 1.0 ) - src6 * ( src6 - 1.0 ) );

    if (debug_output) {
      fprintf(screen, "patch distance was %g, u_LJ = %g\n",
              patch_r, u_LJ);
    }

    // LJ force:
    double lj_force[3] = {0.0, 0.0, 0.0};

    if (!bail) {
      // bond_eps is already a factor 8 smaller so this is regular LJ:
      double pair_f = 48.0*bond_eps*sr6*(sr6 - 0.5)/patch_r2;
      lj_force[0] = pair_f*patch_dist[0];
      lj_force[1] = pair_f*patch_dist[1];
      lj_force[2] = pair_f*patch_dist[2];

      const double c1 = u_LJ*u_phi1*u_phi2;
      const double c2 = u_theta*u_phi1*u_phi2;
      const double c3 = u_LJ*u_theta*u_phi2;
      const double c4 = u_LJ*u_theta*u_phi1;

      // Theta part produces no force.
      acc_force[0] = lj_force[0]*c2 + phi1_force[0]*c3 + phi2_force[0]*c4;
      acc_force[1] = lj_force[1]*c2 + phi1_force[1]*c3 + phi2_force[1]*c4;
      acc_force[2] = lj_force[2]*c2 + phi1_force[2]*c3 + phi2_force[2]*c4;

      // There is another torque coming from the lj force.
      double force_ttor[3], force_rtor[3];
      MathExtra::cross3(b_alp, lj_force, force_ttor);
      MathExtra::cross3(b_bet, lj_force, force_rtor);

      acc_ttor[0] = force_ttor[0]*c2 + theta_ttor[0]*c1 + phi1_ttor[0]*c3 + phi2_ttor[0]*c4;
      acc_ttor[1] = force_ttor[1]*c2 + theta_ttor[1]*c1 + phi1_ttor[1]*c3 + phi2_ttor[1]*c4;
      acc_ttor[2] = force_ttor[2]*c2 + theta_ttor[2]*c1 + phi1_ttor[2]*c3 + phi2_ttor[2]*c4;

      acc_rtor[0] = force_rtor[0]*c2 + theta_rtor[0]*c1 + phi1_rtor[0]*c3 + phi2_rtor[0]*c4;
      acc_rtor[1] = force_rtor[1]*c2 + theta_rtor[1]*c1 + phi1_rtor[1]*c3 + phi2_rtor[1]*c4;
      acc_rtor[2] = force_rtor[2]*c2 + theta_rtor[2]*c1 + phi1_rtor[2]*c3 + phi2_rtor[2]*c4;

      double dE = u_LJ*u_phi1*u_phi2*u_theta;
      if (debug_output)
        fprintf(screen, "  Energy contribution of this pair is %g\n", dE);
      energy += dE;

      fforce[0] += acc_force[0];
      fforce[1] += acc_force[1];
      fforce[2] += acc_force[2];

      ttor[0] += acc_ttor[0];
      ttor[1] += acc_ttor[1];
      ttor[2] += acc_ttor[2];

      rtor[0] += acc_rtor[0];
      rtor[1] += acc_rtor[1];
      rtor[2] += acc_rtor[2];
    }

    if (debug_output && comm->me == 0 && !bail) {
      fprintf(screen, "Energy contributions: %g %g %g %g\n",
              u_LJ,u_theta, u_phi1, u_phi2);
      fprintf(screen, "Forces: (%g %g %g) (%g %g %g)\n",
              acc_force[0], acc_force[1], acc_force[2],
              fforce[0], fforce[1], fforce[2]);
      fprintf(screen, "Torque (t): (%g %g %g) (%g %g %g)\n",
              acc_ttor[0], acc_ttor[1], acc_ttor[2],
              ttor[0], ttor[1], ttor[2]);
      fprintf(screen, "Torque (r): (%g %g %g) (%g %g %g)\n\n",
              acc_rtor[0], acc_rtor[1], acc_rtor[2],
              rtor[0], rtor[1], rtor[2]);
      fprintf(screen, "pair %d, bond vectors: b1.b2 = %g, b1.b3 = %g, b2.b3 = %g\n\n\n",
              pair_i, MathExtra::dot3(b_alp, b_bet), MathExtra::dot3(b_alp,b_gam),
              MathExtra::dot3(b_bet, b_gam));
    }
  }

  return energy;
}



double PairPatchySphere::bond_vec_theta_part(const int i, const int j,
                                             const int itype, const int jtype,
                                             const double *r12, const double *b_alp,
                                             const double *b_bet, double *fforce,
                                             double *ttor, double *rtor)
{
  // Calculate the force and torque acquired due to the gradient in theta_12
  double cos_theta = -MathExtra::dot3(b_alp, b_bet);
  const double theta_param = MathConst::MY_PI / theta_m[itype][jtype];

  const double cos_theta_tol = cos(theta_m[itype][jtype]);
  if (cos_theta <= cos_theta_tol) {
    // There is no theta-contribution to the potential energy.
    // This means the switching function is 0 and hence the entire potential
    // is 0 too.
    fforce[0] = fforce[1] = fforce[2] = 0.0;
    ttor[0] = ttor[1] = ttor[2] = 0.0;
    rtor[0] = rtor[1] = rtor[2] = 0.0;
    return 0;
  }

  double theta = (cos_theta >= 1.0 ? 0.0 :acos(cos_theta));

  // Determine the forces and torques coming from theta:

  const double patch_b = 0.56123102415468651;
  if (theta == 0.0) {
    ttor[0] = ttor[1] = ttor[2] = 0.0;
    rtor[0] = rtor[1] = rtor[2] = 0.0;
    fforce[0] = fforce[1] = fforce[2] = 0.0;
    return 2.0; // At theta == 0 the switching function has its maximum, 2.
  }else{
    double c[3];
    MathExtra::cross3(b_alp, b_bet, c);
    double factor = patch_b*patch_b*theta_param * sin(theta_param*theta)/(2.0*sin(theta));
    ttor[0] = factor*c[0];
    ttor[1] = factor*c[1];
    ttor[2] = factor*c[2];

    rtor[0] = -ttor[0];
    rtor[1] = -ttor[1];
    rtor[2] = -ttor[2];

    fforce[0] = fforce[1] = fforce[2] = 0.0;
    return cos(theta * theta_param) + 1.0;
  }
}


/**
   Calculates the phi part to the force and torque.
   It returns the value of the switching function (in range of [0,2]).
   The force and torque needs to be multiplied by the other switching
   functions and potentials u_LJ, u_phi (the other one) and u_theta.

   If this return 0, you know there is no torque or force for this bond pair.
   If it return 2, the potential is maximal and there is no torque or force from
   this phi but there might be from others.
*/
double PairPatchySphere::bond_vec_phi_part( const int i, const int j,
                                            const int itype, const int jtype,
                                            const double *r12, const double rsq,
                                            const double *b_alp, const double *b_bet,
                                            const double *b_gam, const double *b_eps,
                                            double *fforce, double *ttor, double *rtor)
{

  // Calculate the force and torque acquired due to the gradient in phi

  const double patch_b  = 0.56123102415468651;
  const double patch_b2 = patch_b*patch_b;

  const double gamma_dot_rij = patch_b*MathExtra::dot3(b_gam, r12);
  const double eps_dot_rij   = patch_b*MathExtra::dot3(b_eps, r12);
  const double gamma_dot_eps = patch_b2*MathExtra::dot3(b_gam, b_eps);

  const double r_eps_rij = rsq - eps_dot_rij*eps_dot_rij;
  const double r_gam_rij = rsq - gamma_dot_rij*gamma_dot_rij;

  constexpr const double tol = 1e-4;
  if (r_eps_rij < tol || r_gam_rij < tol) {
    // This means phi is close to 0. In that case, the switching function
    // Is close to 2.0 and there is no torque or force coming from this.
    fforce[0] = fforce[1] = fforce[2] = 0.0;
    ttor[0] = ttor[1] = ttor[2] = 0.0;
    rtor[0] = rtor[1] = rtor[2] = 0.0;
    return 2.0;
  }

  double inv_sq = 1.0 / sqrt(r_gam_rij * r_eps_rij);
  double cos_phi = gamma_dot_eps * rsq - gamma_dot_rij * eps_dot_rij;
  cos_phi *= inv_sq;

  double phi = 0.0;
  if (cos_phi >= 1.0) {
    phi = 0.0;
  }else if (cos_phi < -1.0) {
    phi = MathConst::MY_PI;
  }else{
    phi = acos(cos_phi);
  }

  const double phi_tol = phi_m[itype][jtype];
  if (phi >= phi_tol) {
    fforce[0] = fforce[1] = fforce[2] = 0.0;
    ttor[0] = ttor[1] = ttor[2] = 0.0;
    rtor[0] = rtor[1] = rtor[2];
    // In this case the switching function is 0.
    return 0.0;
  }

  double psi = gamma_dot_rij * eps_dot_rij - gamma_dot_eps * rsq;
  double grad_cos_phi[3];
  double grad_cos_phi_gam[3];
  double fbj[3];

  /*
  grad_cos_phi = b_eps * gamma_dot_rij + b_gam * eps_dot_rij
	  - 2.0*r12*gamma_dot_eps - (r12 - b_eps * eps_dot_rij)*(psi/r_eps_rij)
	  - (r12 - b_gam * gamma_dot_rij)*(psi/r_gam_rij);
  */
  const double p_eps_rij = psi/r_eps_rij;
  const double p_gam_rij = psi/r_gam_rij;

  const double eps_p_eps = eps_dot_rij * p_eps_rij;
  const double gam_p_gam = gamma_dot_rij * p_gam_rij;

  double c1 = gamma_dot_rij + eps_p_eps;
  double c2 = eps_dot_rij + gam_p_gam;
  double c3 = 2.0 * gamma_dot_eps + p_eps_rij + p_gam_rij;

  c1 *= inv_sq;
  c2 *= inv_sq;
  c3 *= inv_sq;

  grad_cos_phi[0] = c1*b_eps[0] + c2*b_gam[0] - c3*r12[0];
  grad_cos_phi[1] = c1*b_eps[1] + c2*b_gam[1] - c3*r12[1];
  grad_cos_phi[2] = c1*b_eps[2] + c2*b_gam[2] - c3*r12[2];

  c1 = eps_dot_rij + gamma_dot_rij * p_gam_rij;
  c2 = rsq * p_gam_rij;
  c1 *= inv_sq;
  c2 *= inv_sq;

  grad_cos_phi_gam[0] = r12[0] * c1 - b_eps[0] * rsq - b_gam[0] * c2;
  grad_cos_phi_gam[1] = r12[1] * c1 - b_eps[1] * rsq - b_gam[1] * c2;
  grad_cos_phi_gam[2] = r12[2] * c1 - b_eps[2] * rsq - b_gam[2] * c2;

  c1 = gamma_dot_rij + eps_dot_rij * p_eps_rij;
  c2 = rsq * p_eps_rij;
  c1 *= inv_sq;
  c2 *= inv_sq;

  fbj[0] = r12[0] * c1 - b_gam[0] * rsq - b_eps[0] * c2;
  fbj[1] = r12[1] * c1 - b_gam[1] * rsq - b_eps[1] * c2;
  fbj[2] = r12[2] * c1 - b_gam[2] * rsq - b_eps[2] * c2;

  const double phi_param = MathConst::MY_PI / phi_m[itype][jtype];
  double sinc = phi > 0.001 ? sin(phi*phi_param)/sin(phi) : phi_param;

  // The torque is
  double chain = 0.5*phi_param * sinc;

  fforce[0] = grad_cos_phi[0] * chain;
  fforce[1] = grad_cos_phi[1] * chain;
  fforce[2] = grad_cos_phi[2] * chain;

  double c[3];
  MathExtra::cross3(b_gam, grad_cos_phi_gam, c);

  ttor[0] = c[0] * chain;
  ttor[1] = c[1] * chain;
  ttor[2] = c[2] * chain;

  double k[3];
  MathExtra::cross3(b_eps, fbj, k);
  rtor[0] = k[0] * chain;
  rtor[1] = k[1] * chain;
  rtor[2] = k[2] * chain;

  fforce[0] = grad_cos_phi[0] * chain;
  fforce[1] = grad_cos_phi[1] * chain;
  fforce[2] = grad_cos_phi[2] * chain;

  return cos(phi*phi_param) + 1.0;
}
