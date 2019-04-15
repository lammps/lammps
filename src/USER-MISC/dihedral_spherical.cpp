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
   Contributing author: Andrew Jewett (Caltech)
   [ using code borrowed from Loukas D. Peristeras (Scienomics SARL)
     and Paul Crozier (SNL) ]
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include "dihedral_spherical.h"

using namespace std;
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

DihedralSpherical::DihedralSpherical(LAMMPS *lmp) : Dihedral(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralSpherical::~DihedralSpherical()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(nterms);

    for (int i=1; i<= atom->ndihedraltypes; i++) {
      if ( Ccoeff[i] ) delete [] Ccoeff[i];
      if ( phi_mult[i] ) delete [] phi_mult[i];
      if ( phi_shift[i] ) delete [] phi_shift[i];
      if ( phi_offset[i] ) delete [] phi_offset[i];
      if ( theta1_mult[i] ) delete [] theta1_mult[i];
      if ( theta1_shift[i] ) delete [] theta1_shift[i];
      if ( theta1_offset[i] ) delete [] theta1_offset[i];
      if ( theta2_mult[i] ) delete [] theta2_mult[i];
      if ( theta2_shift[i] ) delete [] theta2_shift[i];
      if ( theta2_offset[i] ) delete [] theta2_offset[i];
    }
    delete [] Ccoeff;
    delete [] phi_mult;
    delete [] phi_shift;
    delete [] phi_offset;
    delete [] theta1_mult;
    delete [] theta1_shift;
    delete [] theta1_offset;
    delete [] theta2_mult;
    delete [] theta2_shift;
    delete [] theta2_offset;
  }
}


static void norm3safe(double *v) {
  double inv_scale = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  double scale = 1.0;
  if (inv_scale > 0.0)
    scale = 1.0 / inv_scale;
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}


// --------------------------------------------
// ------- Calculate the dihedral angle -------
// --------------------------------------------
static const int g_dim=3;

static double Phi(double const *x1, //array holding x,y,z coords atom 1
                  double const *x2, // :       :      :      :        2
                  double const *x3, // :       :      :      :        3
                  double const *x4, // :       :      :      :        4
                  Domain *domain, //<-periodic boundary information
                  double *vb12,   //<-preallocated vector will store x2-x1
                  double *vb23,   //<-preallocated vector will store x3-x2
                  double *vb34,   //<-preallocated vector will store x4-x3
                  double *n123,   //<-will store normal to plane x1,x2,x3
                  double *n234)   //<-will store normal to plane x2,x3,x4
{

  for (int d=0; d < g_dim; ++d) {
    vb12[d] = x2[d] - x1[d]; // 1st bond
    vb23[d] = x3[d] - x2[d]; // 2nd bond
    vb34[d] = x4[d] - x3[d]; // 3rd bond
  }

  //Consider periodic boundary conditions:
  domain->minimum_image(vb12[0],vb12[1],vb12[2]);
  domain->minimum_image(vb23[0],vb23[1],vb23[2]);
  domain->minimum_image(vb34[0],vb34[1],vb34[2]);

  //--- Compute the normal to the planes formed by atoms 1,2,3 and 2,3,4 ---

  cross3(vb23, vb12, n123);        // <- n123=vb23 x vb12
  cross3(vb23, vb34, n234);        // <- n234=vb23 x vb34

  norm3safe(n123);
  norm3safe(n234);

  double cos_phi = -dot3(n123, n234);

  if (cos_phi > 1.0)
    cos_phi = 1.0;
  else if (cos_phi < -1.0)
    cos_phi = -1.0;

  double phi = acos(cos_phi);

  if (dot3(n123, vb34) > 0.0) {
    phi = -phi;   //(Note: Negative dihedral angles are possible only in 3-D.)
    phi += MY_2PI; //<- This insure phi is always in the range 0 to 2*PI
  }
  return phi;
} // DihedralSpherical::Phi()



/* ---------------------------------------------------------------------- */

void DihedralSpherical::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double edihedral,f1[3],f2[3],f3[3],f4[3];

  double **x = atom->x;
  double **f = atom->f;

  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  // The dihedral angle "phi" is the angle between n123 and n234
  // the planes defined by atoms i1,i2,i3, and i2,i3,i4.
  //
  // Definitions of vectors: vb12, vb23, vb34, perp12on23
  //                         proj12on23, perp43on32, proj43on32
  //
  //  Note: The positions of the 4 atoms are labeled x[i1], x[i2], x[i3], x[i4]
  //        (which are also vectors)
  //
  //             proj12on23                          proj34on23
  //             --------->                         ----------->
  //
  //
  //
  //                  x[i2]                         x[i3]
  //    .                __@----------vb23-------->@                   .
  //   /|\                /|                        \                  |
  //    |                /                           \                 |
  //    |               /                             \                |
  // perp12on23        /                               \               |
  //    |             /                                 \          perp34on23
  //    |          vb12                                  \             |
  //    |           /                                   vb34           |
  //    |          /                                       \           |
  //    |         /                                         \          |
  //    |        /                                           \         |
  //            @                                             \        |
  //                                                          _\|     \|/
  //         x[i1]                                              @
  //
  //                                                           x[i4]
  //

  double vb12[g_dim]; // displacement vector from atom i1 towards atom i2
  //     vb12[d]       = x[i2][d] - x[i1][d]      (for d=0,1,2)
  double vb23[g_dim]; // displacement vector from atom i2 towards atom i3
  //     vb23[d]       = x[i3][d] - x[i2][d]      (for d=0,1,2)
  double vb34[g_dim]; // displacement vector from atom i3 towards atom i4
  //     vb34[d]       = x[i4][d] - x[i3][d]      (for d=0,1,2)

  //  n123 & n234: These two unit vectors are normal to the planes
  //               defined by atoms 1,2,3 and 2,3,4.
  double n123[g_dim]; //n123=vb23 x vb12 / |vb23 x vb12|  ("x" is cross product)
  double n234[g_dim]; //n234=vb23 x vb34 / |vb23 x vb34|  ("x" is cross product)

  // The next 4 vectors are needed to calculate  dphi_dx  = d phi / dx
  double proj12on23[g_dim];
  //    proj12on23[d] = (vb23[d]/|vb23|) * dot3(vb12,vb23)/|vb12|*|vb23|
  double proj34on23[g_dim];
  //    proj34on23[d] = (vb34[d]/|vb23|) * dot3(vb34,vb23)/|vb34|*|vb23|
  double perp12on23[g_dim];
  //    perp12on23[d] = v12[d] - proj12on23[d]
  double perp34on23[g_dim];
  //    perp34on23[d] = v34[d] - proj34on23[d]

  edihedral = 0.0;
  ev_init(eflag,vflag);


  for (n = 0; n < ndihedrallist; n++) {

    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // ------ Step 1: Compute the dihedral angle "phi" ------
    //

    // Phi() calculates the dihedral angle.
    // This function also calculates the vectors:
    // vb12, vb23, vb34, n123, and n234, which we will need later.


    double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain,
                     vb12, vb23, vb34, n123, n234);

    // Step 2: Compute the gradients of phi, theta1, theta2 with atom position:


    // ===================== Step2a) phi dependence: ========================
    //
    // Gradient variables:
    //
    // dphi_dx1, dphi_dx2, dphi_dx3, dphi_dx4 are the gradients of phi with
    // respect to the atomic positions of atoms i1, i2, i3, i4, respectively.
    // As an example, consider dphi_dx1.  The d'th element is:
    double dphi_dx1[g_dim]; //                 d phi
    double dphi_dx2[g_dim]; // dphi_dx1[d] = ----------    (partial derivatives)
    double dphi_dx3[g_dim]; //               d x[i1][d]
    double dphi_dx4[g_dim]; //where d=0,1,2 corresponds to x,y,z    (g_dim==3)

    double dot123             = dot3(vb12, vb23);
    double dot234             = dot3(vb23, vb34);

    double L23sqr             = dot3(vb23, vb23);
    double L23                = sqrt(L23sqr);     // (central bond length)

    double inv_L23sqr = 0.0;
    double inv_L23    = 0.0;
    if (L23sqr != 0.0) {
      inv_L23sqr = 1.0 / L23sqr;
      inv_L23 = 1.0 / L23;
    }

    double neg_inv_L23        = -inv_L23;
    double dot123_over_L23sqr = dot123 * inv_L23sqr;
    double dot234_over_L23sqr = dot234 * inv_L23sqr;

    for (int d=0; d < g_dim; ++d) {
      // See figure above for a visual definitions of these vectors:
      proj12on23[d] = vb23[d] * dot123_over_L23sqr;
      proj34on23[d] = vb23[d] * dot234_over_L23sqr;
      perp12on23[d] = vb12[d] - proj12on23[d];
      perp34on23[d] = vb34[d] - proj34on23[d];
    }

    // --- Compute the gradient vectors dphi/dx1 and dphi/dx4: ---

    // These two gradients point in the direction of n123 and n234,
    // and are scaled by the distances of atoms 1 and 4 from the central axis.
    // Distance of atom 1 to central axis:
    double perp12on23_len = sqrt(dot3(perp12on23, perp12on23));
    // Distance of atom 4 to central axis:
    double perp34on23_len = sqrt(dot3(perp34on23, perp34on23));

    double inv_perp12on23 = 0.0;
    if (perp12on23_len != 0.0) inv_perp12on23 = 1.0 / perp12on23_len;
    double inv_perp34on23 = 0.0;
    if (perp34on23_len != 0.0) inv_perp34on23 = 1.0 / perp34on23_len;

    for (int d=0; d < g_dim; ++d) {
      dphi_dx1[d] = n123[d] * inv_perp12on23;
      dphi_dx4[d] = n234[d] * inv_perp34on23;
    }

    // --- Compute the gradient vectors dphi/dx2 and dphi/dx3: ---
    //
    // This is more tricky because atoms 2 and 3 are shared by both planes
    // 123 and 234 (the angle between which defines "phi").  Moving either
    // one of these atoms effects both the 123 and 234 planes
    // Both the 123 and 234 planes intersect with the plane perpendicular to the
    // central bond axis (vb23).  The two lines where these intersections occur
    // will shift when you move either atom 2 or atom 3.  The angle between
    // these lines is the dihedral angle, phi.  We can define four quantities:
    // dphi123_dx2 is the change in "phi" due to the movement of the 123 plane
    //             ...as a result of moving atom 2.
    // dphi234_dx2 is the change in "phi" due to the movement of the 234 plane
    //             ...as a result of moving atom 2.
    // dphi123_dx3 is the change in "phi" due to the movement of the 123 plane
    //             ...as a result of moving atom 3.
    // dphi234_dx3 is the change in "phi" due to the movement of the 234 plane
    //             ...as a result of moving atom 3.

    double proj12on23_len = dot123 * inv_L23;
    double proj34on23_len = dot234 * inv_L23;
    // Interpretation:
    //The magnitude of "proj12on23_len" is the length of the proj12on23 vector.
    //The sign is positive if it points in the same direction as the central
    //bond (vb23).  Otherwise it is negative.  The same goes for "proj34on23".
    //(In the example figure in the comment above, both variables are positive.)

    // The following 8 lines of code are used to calculate the gradient of phi
    // with respect to the two "middle" atom positions (x[i2] and x[i3]).
    // For an explanation of the formula used below, download the file
    // "dihedral_table_2011-8-02.tar.gz" at the bottom of this post:
    //    http://lammps.sandia.gov/threads/msg22233.html
    // Unpack it and go to this subdirectory:
    //    "supporting_information/doc/gradient_formula_explanation/"
    double dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
    double dphi234_dx2_coef = inv_L23 * proj34on23_len;

    double dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
    double dphi123_dx3_coef = inv_L23 * proj12on23_len;

    for (int d=0; d < g_dim; ++d) {
      // Recall that the n123 and n234 plane normal vectors are proportional to
      // the dphi/dx1 and dphi/dx2 gradients vectors
      // It turns out we can save slightly more CPU cycles by expressing
      // dphi/dx2 and dphi/dx3 as linear combinations of dphi/dx1 and dphi/dx2
      // which we computed already (instead of n123 & n234).
      dphi_dx2[d] = dphi123_dx2_coef*dphi_dx1[d] + dphi234_dx2_coef*dphi_dx4[d];
      dphi_dx3[d] = dphi123_dx3_coef*dphi_dx1[d] + dphi234_dx3_coef*dphi_dx4[d];
    }


    // ============= Step2b) theta1 and theta2 dependence: =============

    // --- Compute the gradient vectors dtheta1/dx1 and dtheta2/dx4: ---

    // These two gradients point in the direction of n123 and n234,
    // and are scaled by the distances of atoms 1 and 4 from the central axis.
    // Distance of atom 1 to central axis:
    double dth1_dx1[g_dim]; //                d theta1      (partial
    double dth1_dx2[g_dim]; // dth1_dx1[d] = ----------     derivative)
    double dth1_dx3[g_dim]; //               d x[i1][d]
    //Note dth1_dx4 = 0

    //Note dth2_dx1 = 0
    double dth2_dx2[g_dim]; //                   d theta2      (partial
    double dth2_dx3[g_dim]; // dth2_dx1[d] = ----------     derivative)
    double dth2_dx4[g_dim]; //                  d x[i1][d]
                            //where d=0,1,2 corresponds to x,y,z  (g_dim==3)

    double L12sqr     = dot3(vb12, vb12);
    double L12        = sqrt(L12sqr);
    double L34sqr     = dot3(vb34, vb34);
    double L34        = sqrt(L34sqr);
    double inv_L12sqr = 0.0;
    double inv_L12    = 0.0;
    double inv_L34sqr = 0.0;
    double inv_L34    = 0.0;
    if (L12sqr != 0.0) {
      inv_L12sqr = 1.0 / L12sqr;
      inv_L12 = 1.0 / L12;
    }
    if (L34sqr != 0.0) {
      inv_L34sqr = 1.0 / L34sqr;
      inv_L34 = 1.0 / L34;
    }

    // The next 2 vectors are needed for calculating dth1_dx = d theta1 / d x
    double proj23on12[g_dim];
    //    proj23on12[d] = (vb12[d]/|vb12|) * dot3(vb23,vb12)/|vb23|*|vb12|
    double perp23on12[g_dim];
    //    perp23on12[d] = v23[d] - proj23on12[d]

    // The next 2 vectors are needed for calculating dth2_dx = d theta2 / d x
    double proj23on34[g_dim];
    //    proj23on34[d] = (vb23[d]/|vb34|) * dot3(vb23,vb34)/|vb23|*|vb34|
    double perp23on34[g_dim];
    //    perp23on34[d] = v23[d] - proj23on34[d]

    double dot123_over_L12sqr = dot123 * inv_L12sqr;
    double dot234_over_L34sqr = dot234 * inv_L34sqr;

    /*                           __            .
     *               proj23on12   .\            .
     *                           .               .  proj23on34
     *                          .                 .
     *                         .                   .
     *                  x[i2] .                    _./   x[i3]
     *                     __@----------vb23-------->@
     *                      /|    /              \    \
     *                     /   theta1          theta2  \
     *                    /  <-'                    `-> \
     *                   /                               \
     *                  /                                 \
     *               vb12                                  \
     *                /                                   vb34
     *               /                                       \
     *              /                                         \
     *             /                                           \
     *            @                                             \
     *                                                          _\|
     *         x[i1]                                              @
     *
     *                                                           x[i4]
     */

    for (int d=0; d < g_dim; ++d) {
      // See figure above for a visual definitions of these vectors:
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
    double coeff_dth1_dx3 =  inv_perp12on23 * inv_L23;
    double coeff_dth2_dx2 = -inv_perp34on23 * inv_L23;
    double coeff_dth2_dx4 =  inv_perp23on34 * inv_L34;

    for (int d=0; d < g_dim; ++d) {
      dth1_dx1[d] = perp23on12[d] * coeff_dth1_dx1;
      dth1_dx3[d] = perp12on23[d] * coeff_dth1_dx3;
      dth1_dx2[d] = -(dth1_dx1[d] + dth1_dx3[d]);
      //dtheta1_dx4 = 0

      //dtheta2_dx1 = 0
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

    double u=0.0;            // u = energy
    double m_du_dth1 = 0.0;  // m_du_dth1 = -du / d theta1
    double m_du_dth2 = 0.0;  // m_du_dth2 = -du / d theta2
    double m_du_dphi = 0.0;  // m_du_dphi = -du / d phi

    u = CalcGeneralizedForces(type,
                              phi, theta1, theta2,
                              &m_du_dth1, &m_du_dth2, &m_du_dphi);

    if (eflag) edihedral = u;

    // ----- Step 4: Calculate the force direction in real space -----

    // chain rule:
    //          d U      d U    d phi    d U    d theta1     d U    d theta2
    // -f  =   -----  = ----- * ----- + -------*-------  + --------*--------
    //          d x     d phi    d x    d theta1   d X     d theta2   d X
    for(int d=0; d < g_dim; ++d) {
      f1[d] = m_du_dphi*dphi_dx1[d]+m_du_dth1*dth1_dx1[d];
                                                           //note: dth2_dx1[d]=0
      f2[d] = m_du_dphi*dphi_dx2[d]+m_du_dth1*dth1_dx2[d]+m_du_dth2*dth2_dx2[d];
      f3[d] = m_du_dphi*dphi_dx3[d]+m_du_dth1*dth1_dx3[d]+m_du_dth2*dth2_dx3[d];
      f4[d] = m_du_dphi*dphi_dx4[d]          +            m_du_dth2*dth2_dx4[d];
                                      //note: dth1_dx4[d] = 0
    }

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,
               nlocal,newton_bond,edihedral,
               f1,f3,f4,
               vb12[0],vb12[1],vb12[2],
               vb23[0],vb23[1],vb23[2],
               vb34[0],vb34[1],vb34[2]);
  }
} // void DihedralSpherical::compute()






// ---           CalcGeneralizedForces()                                    ---
// --- Calculate the energy as a function of theta1, theta2, and phi        ---
// --- as well as its derivatives (with respect to theta1, theta2, and phi) ---



// The code above above is sufficiently general that it can work with any
// any function of the angles theta1, theta2, and phi.  However the
// function below calculates the energy and force according to this specific
// formula:
//
// E(\theta_1,\theta_2,\phi) =
//    \sum_{i=1}^N C_i \Theta_{1i}(\theta_1) \Theta_{2i}(\theta_2) \Phi_i(\phi)
// where:
// \Theta_{1i}(\theta_1)     =  cos((\theta_1-a_i)K_i) + u_i
// \Theta_{2i}(\theta_2)     =  cos((\theta_2-b_i)L_i) + v_i
// \Phi_i(\phi)              =  cos((\phi  -  c_i)M_i) + w_i




double DihedralSpherical::
CalcGeneralizedForces(int type,
                      double phi,
                      double theta1,
                      double theta2,
                      double *m_du_dth1,
                      double *m_du_dth2,
                      double *m_du_dphi)
{
  double energy = 0.0;
  assert(m_du_dphi && m_du_dphi && m_du_dphi);
  *m_du_dphi = 0.0;
  *m_du_dth1 = 0.0;
  *m_du_dth2 = 0.0;

  int i = type;
  for (int j = 0; j < nterms[i]; j++) {

    // (It's common that some terms in an expansion have phi_multi[i][j]=0.
    //  When this happens, perhaps it will speed up the calculation to avoid
    //  unnecessary calls to the cos() and sin() functions. Check this below)
    // I also check whether theta1_mult[i][j] and theta2_mult[i][j] are 0.
    double cp = 1.0;
    double sp = 0.0;
    if (phi_mult[i][j] != 0.0) {
      double p   = phi_mult[i][j]  * (phi    - phi_shift[i][j]);
      cp = cos(p);
      sp = sin(p);
    }

    double ct1 = 1.0;
    double st1 = 0.0;
    if (theta1_mult[i][j] != 0.0) {
      double t1  = theta1_mult[i][j]*(theta1 - theta1_shift[i][j]);
      ct1 = cos(t1);
      st1 = sin(t1);
    }

    double ct2 = 1.0;
    double st2 = 0.0;
    if (theta2_mult[i][j] != 0.0) {
      double t2  = theta2_mult[i][j]*(theta2 - theta2_shift[i][j]);
      ct2 = cos(t2);
      st2 = sin(t2);
    }

    energy     +=  Ccoeff[i][j] * (phi_offset[i][j]    - cp) *
                                  (theta1_offset[i][j] - ct1) *
                                  (theta2_offset[i][j] - ct2);

    // Forces:
    *m_du_dphi += -Ccoeff[i][j] *  sp * phi_mult[i][j] *
                                  (theta1_offset[i][j] - ct1) *
                                  (theta2_offset[i][j] - ct2);

    *m_du_dth1 += -Ccoeff[i][j] * (phi_offset[i][j]    - cp) *
                                   st1 * theta1_mult[i][j] *
                                  (theta2_offset[i][j] - ct2);

    *m_du_dth2 += -Ccoeff[i][j] * (phi_offset[i][j]    - cp) *
                                  (theta1_offset[i][j] - ct1) *
                                   st2 * theta2_mult[i][j];


    // Things to consider later:
    // To speed up the computation, one could try to simplify the expansion:
    //  IE by factoring out common terms, and precomputing trig functions once:
    //     cos(K*theta1), sin(K*theta1),
    //     cos(L*theta2), sin(L*theta2), and
    //     cos(M*phi), sin(M*phi)
    // Also: For integer K,L,M, the trig functions cos(M*phi) and sin(M*phi)
    //       can be calculated more efficiently using polynomials of
    //       cos(phi) and sin(phi)

  } //for (int j = 0; j < nterms[i]; j++) {

  return energy;

} //CalcGeneralizedForces()







void DihedralSpherical::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(nterms,n+1,"dihedral:nterms");

  Ccoeff = new double * [n+1];
  phi_mult = new double * [n+1];
  phi_shift = new double * [n+1];
  phi_offset = new double * [n+1];
  theta1_mult = new double * [n+1];
  theta1_shift = new double * [n+1];
  theta1_offset = new double * [n+1];
  theta2_mult = new double * [n+1];
  theta2_shift = new double * [n+1];
  theta2_offset = new double * [n+1];
  for (int i = 1; i <= n; i++) {
    Ccoeff[i] = NULL;
    phi_mult[i] = NULL;
    phi_shift[i] = NULL;
    phi_offset[i] = NULL;
    theta1_mult[i] = NULL;
    theta1_shift[i] = NULL;
    theta1_offset[i] = NULL;
    theta2_mult[i] = NULL;
    theta2_shift[i] = NULL;
    theta2_offset[i] = NULL;
  }

  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralSpherical::coeff(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  int nterms_one = force->inumeric(FLERR,arg[1]);

  if (nterms_one < 1)
    error->all(FLERR,"Incorrect number of terms arg for dihedral coefficients");

  if (2+10*nterms_one < narg)
    error->all(FLERR,"Incorrect number of arguments for dihedral coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    nterms[i] = nterms_one;
    Ccoeff[i] = new double [nterms_one];
    phi_mult[i] = new double [nterms_one];
    phi_shift[i] = new double [nterms_one];
    phi_offset[i] = new double [nterms_one];
    theta1_mult[i] = new double [nterms_one];
    theta1_shift[i] = new double [nterms_one];
    theta1_offset[i] = new double [nterms_one];
    theta2_mult[i] = new double [nterms_one];
    theta2_shift[i] = new double [nterms_one];
    theta2_offset[i] = new double [nterms_one];
    for (int j = 0; j < nterms_one; j++) {
      int offset = 1+10*j;
      Ccoeff[i][j] = force->numeric(FLERR,arg[offset+1]);
      phi_mult[i][j] = force->numeric(FLERR,arg[offset+2]);
      phi_shift[i][j] = force->numeric(FLERR,arg[offset+3]) * MY_PI/180.0;
      phi_offset[i][j] = force->numeric(FLERR,arg[offset+4]);
      theta1_mult[i][j] = force->numeric(FLERR,arg[offset+5]);
      theta1_shift[i][j] = force->numeric(FLERR,arg[offset+6]) * MY_PI/180.0;
      theta1_offset[i][j] = force->numeric(FLERR,arg[offset+7]);
      theta2_mult[i][j] = force->numeric(FLERR,arg[offset+8]);
      theta2_shift[i][j] = force->numeric(FLERR,arg[offset+9]) * MY_PI/180.0;
      theta2_offset[i][j] = force->numeric(FLERR,arg[offset+10]);
    }
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralSpherical::write_restart(FILE *fp)
{

  fwrite(&nterms[1],sizeof(int),atom->ndihedraltypes,fp);
  for(int i = 1; i <= atom->ndihedraltypes; i++) {
    fwrite(Ccoeff[i],sizeof(double),nterms[i],fp);
    fwrite(phi_mult[i],sizeof(double),nterms[i],fp);
    fwrite(phi_shift[i],sizeof(double),nterms[i],fp);
    fwrite(phi_offset[i],sizeof(double),nterms[i],fp);
    fwrite(theta1_mult[i],sizeof(double),nterms[i],fp);
    fwrite(theta1_shift[i],sizeof(double),nterms[i],fp);
    fwrite(theta1_offset[i],sizeof(double),nterms[i],fp);
    fwrite(theta2_mult[i],sizeof(double),nterms[i],fp);
    fwrite(theta2_shift[i],sizeof(double),nterms[i],fp);
    fwrite(theta2_offset[i],sizeof(double),nterms[i],fp);
  }

}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralSpherical::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0)
    fread(&nterms[1],sizeof(int),atom->ndihedraltypes,fp);

  MPI_Bcast(&nterms[1],atom->ndihedraltypes,MPI_INT,0,world);

  // allocate
  for (int i=1; i<=atom->ndihedraltypes; i++) {
    Ccoeff[i] = new double [nterms[i]];
    phi_mult[i] = new double [nterms[i]];
    phi_shift[i] = new double [nterms[i]];
    phi_offset[i] = new double [nterms[i]];
    theta1_mult[i] = new double [nterms[i]];
    theta1_shift[i] = new double [nterms[i]];
    theta1_offset[i] = new double [nterms[i]];
    theta2_mult[i] = new double [nterms[i]];
    theta2_shift[i] = new double [nterms[i]];
    theta2_offset[i] = new double [nterms[i]];
  }

  if (comm->me == 0) {
    for (int i=1; i<=atom->ndihedraltypes; i++) {
      fread(Ccoeff[i],sizeof(double),nterms[i],fp);
      fread(phi_mult[i],sizeof(double),nterms[i],fp);
      fread(phi_shift[i],sizeof(double),nterms[i],fp);
      fread(phi_offset[i],sizeof(double),nterms[i],fp);
      fread(theta1_mult[i],sizeof(double),nterms[i],fp);
      fread(theta1_shift[i],sizeof(double),nterms[i],fp);
      fread(theta1_offset[i],sizeof(double),nterms[i],fp);
      fread(theta2_mult[i],sizeof(double),nterms[i],fp);
      fread(theta2_shift[i],sizeof(double),nterms[i],fp);
      fread(theta2_offset[i],sizeof(double),nterms[i],fp);
    }
  }

  for (int i=1; i<=atom->ndihedraltypes; i++) {
    MPI_Bcast(Ccoeff[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(phi_mult[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(phi_shift[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(phi_offset[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta1_mult[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta1_shift[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta1_offset[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta2_mult[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta2_shift[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta2_offset[i],nterms[i],MPI_DOUBLE,0,world);
  }

  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}




/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralSpherical::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    fprintf(fp,"%d %d ", i , nterms[i]);
    for (int j = 0; j < nterms[i]; j++) {
      fprintf(fp, "%g %g %g %g %g %g %g %g %g %g ", Ccoeff[i][j],
              phi_mult[i][j], phi_shift[i][j]*180.0/MY_PI, phi_offset[i][j],
              theta1_mult[i][j], theta1_shift[i][j]*180.0/MY_PI,
              theta1_offset[i][j], theta2_mult[i][j],
              theta2_shift[i][j]*180.0/MY_PI, theta2_offset[i][j]);
    }
    fprintf(fp,"\n");
  }
}






// Not needed?
// single() calculates the dihedral-angle energy of atoms i1, i2, i3, i4.
//double DihedralSpherical::single(int type, int i1, int i2, int i3, int i4)
//{
//  //variables we will need
//  double vb12[g_dim];
//  double vb23[g_dim];
//  double vb34[g_dim];
//
//  // Some functions calculate numbers we don't care about. Store in variables:
//  double n123[g_dim]; // (will be ignored)
//  double n234[g_dim]; // (will be ignored)
//  double m_du_dth1;   // (will be ignored)
//  double m_du_dth2;   // (will be ignored)
//  double m_du_dphi;   // (will be ignored)
//
//  double **x = atom->x;
//
//  // Calculate the 4-body angle: phi
//  double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain,
//                   vb12, vb23, vb34, n123, n234);
//
//  // Calculate the 3-body angles: theta1 and theta2
//  double L12 = sqrt(dot3(vb12, vb12));
//  double L23 = sqrt(dot3(vb23, vb23));
//  double L34 = sqrt(dot3(vb34, vb34));
//
//  double ct1 = -dot3(vb12, vb23) / (L12 * L23);
//  if (ct1 < -1.0) ct1 = -1.0;
//  else if (ct1 > 1.0) ct1 = 1.0;
//  double theta1 = acos(ct1);
//
//  double ct2 = -dot3(vb23, vb34) / (L23 * L34);
//  if (ct2 < -1.0) ct2 = -1.0;
//  else if (ct2 > 1.0) ct2 = 1.0;
//  double theta2 = acos(ct2);
//
//  double u = CalcGeneralizedForces(type,
//                                   phi, theta1, theta2,
//                                   &m_du_dth1, &m_du_dth2, &m_du_dphi);
//  return u;
//}

