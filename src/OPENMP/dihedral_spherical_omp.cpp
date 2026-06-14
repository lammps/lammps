// clang-format off
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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "dihedral_spherical_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "neighbor.h"

#include <cmath>

#include "omp_compat.h"
#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;

// --------------------------------------------
// ------- Calculate the dihedral angle -------
// --------------------------------------------
static const int g_dim = 3;

static void norm3safe(double *v)
{
  double inv_scale = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  double scale = 1.0;
  if (inv_scale > 0.0) scale = 1.0 / inv_scale;
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

static double Phi(double const *x1,    //array holding x,y,z coords atom 1
                  double const *x2,    // :       :      :      :        2
                  double const *x3,    // :       :      :      :        3
                  double const *x4,    // :       :      :      :        4
                  Domain *domain,      //<-periodic boundary information
                  double *vb12,        //<-preallocated vector will store x2-x1
                  double *vb23,        //<-preallocated vector will store x3-x2
                  double *vb34,        //<-preallocated vector will store x4-x3
                  double *n123,        //<-will store normal to plane x1,x2,x3
                  double *n234)        //<-will store normal to plane x2,x3,x4
{
  for (int d = 0; d < g_dim; ++d) {
    vb12[d] = x2[d] - x1[d];    // 1st bond
    vb23[d] = x3[d] - x2[d];    // 2nd bond
    vb34[d] = x4[d] - x3[d];    // 3rd bond
  }

  //Consider periodic boundary conditions:
  domain->minimum_image(FLERR, vb12[0],vb12[1],vb12[2]);
  domain->minimum_image(FLERR, vb23[0],vb23[1],vb23[2]);
  domain->minimum_image(FLERR, vb34[0],vb34[1],vb34[2]);

  //--- Compute the normal to the planes formed by atoms 1,2,3 and 2,3,4 ---

  cross3(vb23, vb12, n123);    // <- n123=vb23 x vb12
  cross3(vb23, vb34, n234);    // <- n234=vb23 x vb34

  norm3safe(n123);
  norm3safe(n234);

  double cos_phi = -dot3(n123, n234);

  if (cos_phi > 1.0)
    cos_phi = 1.0;
  else if (cos_phi < -1.0)
    cos_phi = -1.0;

  double phi = acos(cos_phi);

  if (dot3(n123, vb34) > 0.0) {
    phi = -phi;       //(Note: Negative dihedral angles are possible only in 3-D.)
    phi += MY_2PI;    //<- This ensure phi is always in the range 0 to 2*PI
  }
  return phi;
}    // Phi()

/* ---------------------------------------------------------------------- */

DihedralSphericalOMP::DihedralSphericalOMP(class LAMMPS *lmp)
  : DihedralSpherical(lmp), ThrOMP(lmp,THR_DIHEDRAL)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void DihedralSphericalOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->ndihedrallist;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, cvatom, thr);

    if (inum > 0) {
      if (evflag) {
        if (eflag) {
          if (force->newton_bond) eval<1,1,1>(ifrom, ito, thr);
          else eval<1,1,0>(ifrom, ito, thr);
        } else {
          if (force->newton_bond) eval<1,0,1>(ifrom, ito, thr);
          else eval<1,0,0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_bond) eval<0,0,1>(ifrom, ito, thr);
        else eval<0,0,0>(ifrom, ito, thr);
      }
    }
    thr->timer(Timer::BOND);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void DihedralSphericalOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,n,type;
  double edihedral,f1[3],f2[3],f3[3],f4[3];

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const * const dihedrallist = neighbor->dihedrallist;
  const int nlocal = atom->nlocal;

  double vb12[g_dim];    // displacement vector from atom i1 towards atom i2
  double vb23[g_dim];    // displacement vector from atom i2 towards atom i3
  double vb34[g_dim];    // displacement vector from atom i3 towards atom i4

  double n123[g_dim];    //n123=vb23 x vb12 / |vb23 x vb12|
  double n234[g_dim];    //n234=vb23 x vb34 / |vb23 x vb34|

  double proj12on23[g_dim];
  double proj34on23[g_dim];
  double perp12on23[g_dim];
  double perp34on23[g_dim];

  edihedral = 0.0;

  for (n = nfrom; n < nto; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // ------ Step 1: Compute the dihedral angle "phi" ------

    double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain,
                     vb12, vb23, vb34, n123, n234);

    // Step 2: Compute the gradients of phi, theta1, theta2 with atom position

    // ===================== Step2a) phi dependence: ========================

    double dphi_dx1[g_dim];
    double dphi_dx2[g_dim];
    double dphi_dx3[g_dim];
    double dphi_dx4[g_dim];

    double dot123 = dot3(vb12, vb23);
    double dot234 = dot3(vb23, vb34);

    double L23sqr = dot3(vb23, vb23);
    double L23 = sqrt(L23sqr);

    double inv_L23sqr = 0.0;
    double inv_L23 = 0.0;
    if (L23sqr != 0.0) {
      inv_L23sqr = 1.0 / L23sqr;
      inv_L23 = 1.0 / L23;
    }

    double neg_inv_L23 = -inv_L23;
    double dot123_over_L23sqr = dot123 * inv_L23sqr;
    double dot234_over_L23sqr = dot234 * inv_L23sqr;

    for (int d = 0; d < g_dim; ++d) {
      proj12on23[d] = vb23[d] * dot123_over_L23sqr;
      proj34on23[d] = vb23[d] * dot234_over_L23sqr;
      perp12on23[d] = vb12[d] - proj12on23[d];
      perp34on23[d] = vb34[d] - proj34on23[d];
    }

    double perp12on23_len = sqrt(dot3(perp12on23, perp12on23));
    double perp34on23_len = sqrt(dot3(perp34on23, perp34on23));

    double inv_perp12on23 = 0.0;
    if (perp12on23_len != 0.0) inv_perp12on23 = 1.0 / perp12on23_len;
    double inv_perp34on23 = 0.0;
    if (perp34on23_len != 0.0) inv_perp34on23 = 1.0 / perp34on23_len;

    for (int d = 0; d < g_dim; ++d) {
      dphi_dx1[d] = n123[d] * inv_perp12on23;
      dphi_dx4[d] = n234[d] * inv_perp34on23;
    }

    double proj12on23_len = dot123 * inv_L23;
    double proj34on23_len = dot234 * inv_L23;

    double dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
    double dphi234_dx2_coef = inv_L23 * proj34on23_len;

    double dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
    double dphi123_dx3_coef = inv_L23 * proj12on23_len;

    for (int d = 0; d < g_dim; ++d) {
      dphi_dx2[d] = dphi123_dx2_coef*dphi_dx1[d] + dphi234_dx2_coef*dphi_dx4[d];
      dphi_dx3[d] = dphi123_dx3_coef*dphi_dx1[d] + dphi234_dx3_coef*dphi_dx4[d];
    }

    // ============= Step2b) theta1 and theta2 dependence: =============

    double dth1_dx1[g_dim];
    double dth1_dx2[g_dim];
    double dth1_dx3[g_dim];

    double dth2_dx2[g_dim];
    double dth2_dx3[g_dim];
    double dth2_dx4[g_dim];

    double L12sqr = dot3(vb12, vb12);
    double L12 = sqrt(L12sqr);
    double L34sqr = dot3(vb34, vb34);
    double L34 = sqrt(L34sqr);
    double inv_L12sqr = 0.0;
    double inv_L12 = 0.0;
    double inv_L34sqr = 0.0;
    double inv_L34 = 0.0;
    if (L12sqr != 0.0) {
      inv_L12sqr = 1.0 / L12sqr;
      inv_L12 = 1.0 / L12;
    }
    if (L34sqr != 0.0) {
      inv_L34sqr = 1.0 / L34sqr;
      inv_L34 = 1.0 / L34;
    }

    double proj23on12[g_dim];
    double perp23on12[g_dim];
    double proj23on34[g_dim];
    double perp23on34[g_dim];

    double dot123_over_L12sqr = dot123 * inv_L12sqr;
    double dot234_over_L34sqr = dot234 * inv_L34sqr;

    for (int d = 0; d < g_dim; ++d) {
      proj23on12[d] = vb12[d] * dot123_over_L12sqr;
      proj23on34[d] = vb34[d] * dot234_over_L34sqr;
      perp23on12[d] = vb23[d] - proj23on12[d];
      perp23on34[d] = vb23[d] - proj23on34[d];
    }

    double perp23on12_len = sqrt(dot3(perp23on12, perp23on12));
    double perp23on34_len = sqrt(dot3(perp23on34, perp23on34));

    double inv_perp23on12 = 0.0;
    if (perp23on12_len != 0.0) inv_perp23on12 = 1.0 / perp23on12_len;
    double inv_perp23on34 = 0.0;
    if (perp23on34_len != 0.0) inv_perp23on34 = 1.0 / perp23on34_len;

    double coeff_dth1_dx1 = -inv_perp23on12 * inv_L12;
    double coeff_dth1_dx3 = inv_perp12on23 * inv_L23;
    double coeff_dth2_dx2 = -inv_perp34on23 * inv_L23;
    double coeff_dth2_dx4 = inv_perp23on34 * inv_L34;

    for (int d = 0; d < g_dim; ++d) {
      dth1_dx1[d] = perp23on12[d] * coeff_dth1_dx1;
      dth1_dx3[d] = perp12on23[d] * coeff_dth1_dx3;
      dth1_dx2[d] = -(dth1_dx1[d] + dth1_dx3[d]);

      dth2_dx2[d] = perp34on23[d] * coeff_dth2_dx2;
      dth2_dx4[d] = perp23on34[d] * coeff_dth2_dx4;
      dth2_dx3[d] = -(dth2_dx2[d] + dth2_dx4[d]);
    }

    double ct1 = -dot123 * inv_L12 * inv_L23;
    if (ct1 < -1.0) ct1 = -1.0;
    else if (ct1 > 1.0) ct1 = 1.0;
    double theta1 = acos(ct1);

    double ct2 = -dot234 * inv_L23 * inv_L34;
    if (ct2 < -1.0) ct2 = -1.0;
    else if (ct2 > 1.0) ct2 = 1.0;
    double theta2 = acos(ct2);

    // - Step 3: Calculate the energy and force in the phi & theta1/2 directions

    double u = 0.0;
    double m_du_dth1 = 0.0;
    double m_du_dth2 = 0.0;
    double m_du_dphi = 0.0;

    u = CalcGeneralizedForces(type, phi, theta1, theta2, &m_du_dth1, &m_du_dth2, &m_du_dphi);

    if (EFLAG) edihedral = u;

    // ----- Step 4: Calculate the force direction in real space -----

    for (int d = 0; d < g_dim; ++d) {
      f1[d] = m_du_dphi*dphi_dx1[d] + m_du_dth1*dth1_dx1[d];
      f2[d] = m_du_dphi*dphi_dx2[d] + m_du_dth1*dth1_dx2[d] + m_du_dth2*dth2_dx2[d];
      f3[d] = m_du_dphi*dphi_dx3[d] + m_du_dth1*dth1_dx3[d] + m_du_dth2*dth2_dx3[d];
      f4[d] = m_du_dphi*dphi_dx4[d] + m_du_dth2*dth2_dx4[d];
    }

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (EVFLAG)
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,edihedral,f1,f3,f4,
                   vb12[0],vb12[1],vb12[2],vb23[0],vb23[1],vb23[2],
                   vb34[0],vb34[1],vb34[2],thr);
  }
}
