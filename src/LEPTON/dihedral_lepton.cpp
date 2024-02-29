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
   Using parts of dihedral_table.cpp by Andrew Jewett (jewett.aij at gmail)
------------------------------------------------------------------------- */

#include "dihedral_lepton.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>

#include "Lepton.h"
#include "lepton_utils.h"

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::MY_2PI;
using MathConst::RAD2DEG;
using MathExtra::cross3;
using MathExtra::dot3;
using MathExtra::norm3;

static constexpr int g_dim = 3;

/* ---------------------------------------------------------------------- */

DihedralLepton::DihedralLepton(LAMMPS *_lmp) : Dihedral(_lmp), type2expression(nullptr)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralLepton::~DihedralLepton()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(type2expression);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralLepton::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  ev_init(eflag, vflag);
  if (evflag) {
    if (eflag) {
      if (force->newton_bond)
        eval<1, 1, 1>();
      else
        eval<1, 1, 0>();
    } else {
      if (force->newton_bond)
        eval<1, 0, 1>();
      else
        eval<1, 0, 0>();
    }
  } else {
    if (force->newton_bond)
      eval<0, 0, 1>();
    else
      eval<0, 0, 0>();
  }
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_BOND> void DihedralLepton::eval()
{
  std::vector<Lepton::CompiledExpression> dihedralforce;
  std::vector<Lepton::CompiledExpression> dihedralpot;
  std::vector<bool> has_ref;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
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
  double *const *const f = atom->f;
  const int *const *const dihedrallist = neighbor->dihedrallist;
  const int ndihedrallist = neighbor->ndihedrallist;
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

  for (int n = 0; n < ndihedrallist; n++) {
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
      ev_tally(i1, i2, i3, i4, nlocal, NEWTON_BOND, edihedral, f1, f3, f4, -vb12[0], -vb12[1],
               -vb12[2], vb23[0], vb23[1], vb23[2], vb34[0], vb34[1], vb34[2]);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralLepton::allocate()
{
  allocated = 1;
  const int np1 = atom->ndihedraltypes + 1;

  memory->create(type2expression, np1, "dihedral:type2expression");
  memory->create(setflag, np1, "dihedral:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void DihedralLepton::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect number of args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  // remove whitespace and quotes from expression string and then
  // check if the expression can be parsed and evaluated without error
  std::string exp_one = LeptonUtils::condense(arg[1]);
  try {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp));
    auto dihedralpot = parsed.createCompiledExpression();
    auto dihedralforce = parsed.differentiate("phi").createCompiledExpression();
    try {
      dihedralpot.getVariableReference("phi") = 0.0;
    } catch (Lepton::Exception &) {
      if (comm->me == 0)
        error->warning(FLERR, "Lepton potential expression {} does not depend on 'phi'", exp_one);
    }
    try {
      dihedralforce.getVariableReference("phi") = 0.0;
    } catch (Lepton::Exception &) {
      if (comm->me == 0)
        error->warning(FLERR, "Force from Lepton expression {} does not depend on 'phi'", exp_one);
    }
    dihedralforce.evaluate();
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  std::size_t idx = 0;
  for (const auto &exp : expressions) {
    if (exp == exp_one) break;
    ++idx;
  }

  // if not found, add to list
  if ((expressions.size() == 0) || (idx == expressions.size())) expressions.push_back(exp_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    type2expression[i] = idx;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralLepton::write_restart(FILE *fp)
{
  fwrite(&type2expression[1], sizeof(int), atom->ndihedraltypes, fp);

  int num = expressions.size();
  int maxlen = 0;
  for (const auto &exp : expressions) maxlen = MAX(maxlen, (int) exp.size());
  ++maxlen;

  fwrite(&num, sizeof(int), 1, fp);
  fwrite(&maxlen, sizeof(int), 1, fp);
  for (const auto &exp : expressions) {
    int n = exp.size() + 1;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(exp.c_str(), sizeof(char), n, fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralLepton::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &type2expression[1], sizeof(int), atom->ndihedraltypes, fp, nullptr,
                  error);
  }
  MPI_Bcast(&type2expression[1], atom->ndihedraltypes, MPI_INT, 0, world);
  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;

  int num, maxlen, len;
  if (comm->me == 0) {
    utils::sfread(FLERR, &num, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &maxlen, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&num, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxlen, 1, MPI_INT, 0, world);

  char *buf = new char[maxlen];

  for (int i = 0; i < num; ++i) {
    if (comm->me == 0) {
      utils::sfread(FLERR, &len, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, buf, sizeof(char), len, fp, nullptr, error);
    }
    MPI_Bcast(buf, maxlen, MPI_CHAR, 0, world);
    expressions.emplace_back(buf);
  }

  delete[] buf;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralLepton::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp, "%d %s\n", i, expressions[type2expression[i]].c_str());
}

// --------------------------------------------
// ------- Calculate the dihedral angle -------
// --------------------------------------------

double DihedralLepton::get_phi(double const *x1,    //array holding x,y,z coords atom 1
                               double const *x2,    // :       :      :      :        2
                               double const *x3,    // :       :      :      :        3
                               double const *x4,    // :       :      :      :        4
                               Domain *domain,      //<-periodic boundary information
                               // The following arrays are of doubles with g_dim elements.
                               // (g_dim is a constant known at compile time, usually 3).
                               // Their contents is calculated by this function.
                               // Space for these vectors must be allocated in advance.
                               // (This is not hidden internally because these vectors
                               //  may be needed outside the function, later on.)
                               double *vb12,    // will store x2-x1
                               double *vb23,    // will store x3-x2
                               double *vb34,    // will store x4-x3
                               double *n123,    // will store normal to plane x1,x2,x3
                               double *n234)    // will store normal to plane x2,x3,x4
    const
{
  for (int d = 0; d < g_dim; ++d) {
    vb12[d] = x2[d] - x1[d];    // 1st bond
    vb23[d] = x3[d] - x2[d];    // 2nd bond
    vb34[d] = x4[d] - x3[d];    // 3rd bond
  }

  //Consider periodic boundary conditions:
  domain->minimum_image(vb12[0], vb12[1], vb12[2]);
  domain->minimum_image(vb23[0], vb23[1], vb23[2]);
  domain->minimum_image(vb34[0], vb34[1], vb34[2]);

  //--- Compute the normal to the planes formed by atoms 1,2,3 and 2,3,4 ---

  cross3(vb23, vb12, n123);    // <- n123=vb23 x vb12
  cross3(vb23, vb34, n234);    // <- n234=vb23 x vb34

  norm3(n123);
  norm3(n234);

  double cos_phi = -dot3(n123, n234);

  if (cos_phi > 1.0)
    cos_phi = 1.0;
  else if (cos_phi < -1.0)
    cos_phi = -1.0;

  double phi = acos(cos_phi);

  if (dot3(n123, vb34) > 0.0) {
    phi = -phi;       //(Note: Negative dihedral angles are possible only in 3-D.)
    phi += MY_2PI;    //<- This ensures phi is always in the range 0 to 2*PI
  }
  return phi;
}
