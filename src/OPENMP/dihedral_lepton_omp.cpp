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

#include "dihedral_lepton_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

#include "Lepton.h"
#include "lepton_utils.h"
#include "omp_compat.h"
using namespace LAMMPS_NS;
using MathExtra::dot3;

static constexpr int g_dim = 3;

/* ---------------------------------------------------------------------- */

DihedralLeptonOMP::DihedralLeptonOMP(class LAMMPS *lmp) :
    DihedralLepton(lmp), ThrOMP(lmp, THR_DIHEDRAL)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void DihedralLeptonOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->ndihedrallist;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (inum > 0) {
      if (evflag) {
        if (eflag) {
          if (force->newton_bond)
            eval<1, 1, 1>(ifrom, ito, thr);
          else
            eval<1, 1, 0>(ifrom, ito, thr);
        } else {
          if (force->newton_bond)
            eval<1, 0, 1>(ifrom, ito, thr);
          else
            eval<1, 0, 0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_bond)
          eval<0, 0, 1>(ifrom, ito, thr);
        else
          eval<0, 0, 0>(ifrom, ito, thr);
      }
    }
    thr->timer(Timer::BOND);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void DihedralLeptonOMP::eval(int nfrom, int nto, ThrData *const thr)
{
  std::vector<Lepton::CompiledExpression> dihedralforce;
  std::vector<Lepton::CompiledExpression> dihedralpot;
  std::vector<bool> has_ref;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, Pointers::lmp));
      dihedralforce.emplace_back(parsed.differentiate("phi").createCompiledExpression());
      has_ref.push_back(true);
      try {
        dihedralforce.back().getVariableReference("phi");
      } catch (Lepton::Exception &) {
        has_ref.back() = false;
      }
      if (EFLAG) dihedralpot.emplace_back(parsed.createCompiledExpression());
    }
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  const double *const *const x = atom->x;
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *const *const dihedrallist = neighbor->dihedrallist;
  const int nlocal = atom->nlocal;

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
  //                           .
  //                          .
  //                         .
  //                  x[i2] .                       x[i3]
  //    .                __@----------vb23-------->@ . . . .           .
  //   /|\                /|                        \                  |
  //    |                /                           \                 |
  //    |               /                             \                |
  // perp12vs23        /                               \               |
  //    |             /                                 \          perp34vs23
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

  double vb12[g_dim];    // displacement vector from atom i1 towards atom i2
  //     vb12[d]       = x[i2][d] - x[i1][d]      (for d=0,1,2)
  double vb23[g_dim];    // displacement vector from atom i2 towards atom i3
  //     vb23[d]       = x[i3][d] - x[i2][d]      (for d=0,1,2)
  double vb34[g_dim];    // displacement vector from atom i3 towards atom i4
  //     vb34[d]       = x[i4][d] - x[i3][d]      (for d=0,1,2)

  //  n123 & n234: These two unit vectors are normal to the planes
  //               defined by atoms 1,2,3 and 2,3,4.
  double n123[g_dim];    //n123=vb23 x vb12 / |vb23 x vb12|  ("x" is cross product)
  double n234[g_dim];    //n234=vb23 x vb34 / |vb23 x vb34|  ("x" is cross product)

  double proj12on23[g_dim];
  //    proj12on23[d] = (vb23[d]/|vb23|) * dot3(vb12,vb23)/|vb12|*|vb23|
  double proj34on23[g_dim];
  //    proj34on23[d] = (vb34[d]/|vb23|) * dot3(vb34,vb23)/|vb34|*|vb23|
  double perp12on23[g_dim];
  //    perp12on23[d] = v12[d] - proj12on23[d]
  double perp34on23[g_dim];
  //    perp34on23[d] = v34[d] - proj34on23[d]

  double f1[3], f2[3], f3[3], f4[3];

  for (int n = nfrom; n < nto; n++) {
    const int i1 = dihedrallist[n][0];
    const int i2 = dihedrallist[n][1];
    const int i3 = dihedrallist[n][2];
    const int i4 = dihedrallist[n][3];
    const int type = dihedrallist[n][4];

    // ------ Step 1: Compute the dihedral angle "phi" ------
    //

    // get_phi() calculates the dihedral angle.
    // This function also calculates the vectors:
    // vb12, vb23, vb34, n123, and n234, which we will need later.

    const double phi = get_phi(x[i1], x[i2], x[i3], x[i4], domain, vb12, vb23, vb34, n123, n234);

    // ------ Step 2: Compute the gradient of phi with atomic position: ------
    //
    // Gradient variables:
    //
    // dphi_dx1, dphi_dx2, dphi_dx3, dphi_dx4 are the gradients of phi with
    // respect to the atomic positions of atoms i1, i2, i3, i4, respectively.
    // As an example, consider dphi_dx1.  The d'th element is:
    double dphi_dx1[g_dim];    //                 d phi
    double dphi_dx2[g_dim];    // dphi_dx1[d] = ----------    (partial derivatives)
    double dphi_dx3[g_dim];    //               d x[i1][d]
    double dphi_dx4[g_dim];    //where d=0,1,2 corresponds to x,y,z  (if g_dim==3)

    double dot123 = dot3(vb12, vb23);
    double dot234 = dot3(vb23, vb34);
    double L23sqr = dot3(vb23, vb23);
    double L23 = sqrt(L23sqr);    // (central bond length)
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

    for (int d = 0; d < g_dim; ++d) {
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

    // The forumula used in the 8 lines below explained here:
    //   "supporting_information/doc/gradient_formula_explanation/"
    double dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
    double dphi234_dx2_coef = inv_L23 * proj34on23_len;

    double dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
    double dphi123_dx3_coef = inv_L23 * proj12on23_len;

    for (int d = 0; d < g_dim; ++d) {
      // Recall that the n123 and n234 plane normal vectors are proportional to
      // the dphi/dx1 and dphi/dx2 gradients vectors
      // It turns out we can save slightly more CPU cycles by expressing
      // dphi/dx2 and dphi/dx3 as linear combinations of dphi/dx1 and dphi/dx2
      // which we computed already (instead of n123 & n234).
      dphi_dx2[d] = dphi123_dx2_coef * dphi_dx1[d] + dphi234_dx2_coef * dphi_dx4[d];
      dphi_dx3[d] = dphi123_dx3_coef * dphi_dx1[d] + dphi234_dx3_coef * dphi_dx4[d];
    }

    const int idx = type2expression[type];
    if (has_ref[idx]) dihedralforce[idx].getVariableReference("phi") = phi;
    double m_du_dphi = -dihedralforce[idx].evaluate();

    // ----- Step 4: Calculate the force direction in real space -----

    // chain rule:
    //          d U          d U      d phi
    // -f  =   -----   =    -----  *  -----
    //          d x         d phi      d x
    for (int d = 0; d < g_dim; ++d) {
      f1[d] = m_du_dphi * dphi_dx1[d];
      f2[d] = m_du_dphi * dphi_dx2[d];
      f3[d] = m_du_dphi * dphi_dx3[d];
      f4[d] = m_du_dphi * dphi_dx4[d];
    }

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += f1[0];
      f[i1].y += f1[1];
      f[i1].z += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x += f2[0];
      f[i2].y += f2[1];
      f[i2].z += f2[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3].x += f3[0];
      f[i3].y += f3[1];
      f[i3].z += f3[2];
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4].x += f4[0];
      f[i4].y += f4[1];
      f[i4].z += f4[2];
    }

    double edihedral = 0.0;
    if (EFLAG) {
      try {
        dihedralpot[idx].getVariableReference("phi") = phi;
      } catch (Lepton::Exception &) {
        ;    // ignore -> constant potential
      }
      edihedral = dihedralpot[idx].evaluate();
    }
    if (EVFLAG)
      ev_tally_thr(this, i1, i2, i3, i4, nlocal, NEWTON_BOND, edihedral, f1, f3, f4, -vb12[0],
                   -vb12[1], -vb12[2], vb23[0], vb23[1], vb23[2], vb34[0], vb34[1], vb34[2], thr);
  }
}
