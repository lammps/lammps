// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <cstdlib>
#include <cstring>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarparse.h"

#ifdef COLVARS_LAMMPS
// Use open-source Jacobi implementation
#include "math_eigen_impl.h"
#else
// Fall back to NR routine
#include "nr_jacobi.h"
#endif


bool      colvarmodule::rotation::monitor_crossings = false;
cvm::real colvarmodule::rotation::crossing_threshold = 1.0E-02;


std::string cvm::rvector::to_simple_string() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(cvm::cv_prec);
  os << x << " " << y << " " << z;
  return os.str();
}


int cvm::rvector::from_simple_string(std::string const &s)
{
  std::stringstream stream(s);
  if ( !(stream >> x) ||
       !(stream >> y) ||
       !(stream >> z) ) {
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


std::ostream & operator << (std::ostream &os, colvarmodule::rvector const &v)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  os.width(2);
  os << "( ";
  os.width(w); os.precision(p);
  os << v.x << " , ";
  os.width(w); os.precision(p);
  os << v.y << " , ";
  os.width(w); os.precision(p);
  os << v.z << " )";
  return os;
}


std::istream & operator >> (std::istream &is, colvarmodule::rvector &v)
{
  std::streampos const start_pos = is.tellg();
  char sep;
  if ( !(is >> sep) || !(sep == '(') ||
       !(is >> v.x) || !(is >> sep)  || !(sep == ',') ||
       !(is >> v.y) || !(is >> sep)  || !(sep == ',') ||
       !(is >> v.z) || !(is >> sep)  || !(sep == ')') ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }
  return is;
}

std::string cvm::quaternion::to_simple_string() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(cvm::cv_prec);
  os << q0 << " " << q1 << " " << q2 << " " << q3;
  return os.str();
}

int cvm::quaternion::from_simple_string(std::string const &s)
{
  std::stringstream stream(s);
  if ( !(stream >> q0) ||
       !(stream >> q1) ||
       !(stream >> q2) ||
       !(stream >> q3) ) {
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

std::ostream & operator << (std::ostream &os, colvarmodule::quaternion const &q)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  os.width(2);
  os << "( ";
  os.width(w); os.precision(p);
  os << q.q0 << " , ";
  os.width(w); os.precision(p);
  os << q.q1 << " , ";
  os.width(w); os.precision(p);
  os << q.q2 << " , ";
  os.width(w); os.precision(p);
  os << q.q3 << " )";
  return os;
}


std::istream & operator >> (std::istream &is, colvarmodule::quaternion &q)
{
  std::streampos const start_pos = is.tellg();
  char sep;
  if ( !(is >> sep)  || !(sep == '(') ||
       !(is >> q.q0) || !(is >> sep)  || !(sep == ',') ||
       !(is >> q.q1) || !(is >> sep)  || !(sep == ',') ||
       !(is >> q.q2) || !(is >> sep)  || !(sep == ',') ||
       !(is >> q.q3) || !(is >> sep)  || !(sep == ')') ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
  }
  return is;
}


cvm::quaternion
cvm::quaternion::position_derivative_inner(cvm::rvector const &pos,
                                            cvm::rvector const &vec) const
{
  cvm::quaternion result(0.0, 0.0, 0.0, 0.0);


  result.q0 =   2.0 * pos.x * q0 * vec.x
               +2.0 * pos.y * q0 * vec.y
               +2.0 * pos.z * q0 * vec.z

               -2.0 * pos.y * q3 * vec.x
               +2.0 * pos.z * q2 * vec.x

               +2.0 * pos.x * q3 * vec.y
               -2.0 * pos.z * q1 * vec.y

               -2.0 * pos.x * q2 * vec.z
               +2.0 * pos.y * q1 * vec.z;


  result.q1 =  +2.0 * pos.x * q1 * vec.x
               -2.0 * pos.y * q1 * vec.y
               -2.0 * pos.z * q1 * vec.z

               +2.0 * pos.y * q2 * vec.x
               +2.0 * pos.z * q3 * vec.x

               +2.0 * pos.x * q2 * vec.y
               -2.0 * pos.z * q0 * vec.y

               +2.0 * pos.x * q3 * vec.z
               +2.0 * pos.y * q0 * vec.z;


  result.q2 =  -2.0 * pos.x * q2 * vec.x
               +2.0 * pos.y * q2 * vec.y
               -2.0 * pos.z * q2 * vec.z

               +2.0 * pos.y * q1 * vec.x
               +2.0 * pos.z * q0 * vec.x

               +2.0 * pos.x * q1 * vec.y
               +2.0 * pos.z * q3 * vec.y

               -2.0 * pos.x * q0 * vec.z
               +2.0 * pos.y * q3 * vec.z;


  result.q3 =  -2.0 * pos.x * q3 * vec.x
               -2.0 * pos.y * q3 * vec.y
               +2.0 * pos.z * q3 * vec.z

               -2.0 * pos.y * q0 * vec.x
               +2.0 * pos.z * q1 * vec.x

               +2.0 * pos.x * q0 * vec.y
               +2.0 * pos.z * q2 * vec.y

               +2.0 * pos.x * q1 * vec.z
               +2.0 * pos.y * q2 * vec.z;

  return result;
}




// Calculate the optimal rotation between two groups, and implement it
// as a quaternion.  Uses the method documented in: Coutsias EA,
// Seok C, Dill KA.  Using quaternions to calculate RMSD.  J Comput
// Chem. 25(15):1849-57 (2004) DOI: 10.1002/jcc.20110 PubMed: 15376254

#ifdef COLVARS_LAMMPS
namespace {
  inline void *new_Jacobi_solver(int size) {
    return reinterpret_cast<void *>(new MathEigen::Jacobi<cvm::real,
                                    cvm::vector1d<cvm::real> &,
                                    cvm::matrix2d<cvm::real> &>(4));
  }
}
#endif


int colvarmodule::rotation::init()
{
  b_debug_gradients = false;
  lambda = 0.0;
  cvm::main()->cite_feature("Optimal rotation via flexible fitting");
  return COLVARS_OK;
}


colvarmodule::rotation::rotation()
{
  init();
#ifdef COLVARS_LAMMPS
  jacobi = new_Jacobi_solver(4);
#else
  jacobi = NULL;
#endif
}


colvarmodule::rotation::rotation(cvm::quaternion const &qi)
  : q(qi)
{
  init();
#ifdef COLVARS_LAMMPS
  jacobi = new_Jacobi_solver(4);
#else
  jacobi = NULL;
#endif
}


colvarmodule::rotation::rotation(cvm::real angle, cvm::rvector const &axis)
{
  init();
  cvm::rvector const axis_n = axis.unit();
  cvm::real const sina = cvm::sin(angle/2.0);
  q = cvm::quaternion(cvm::cos(angle/2.0),
                      sina * axis_n.x, sina * axis_n.y, sina * axis_n.z);
#ifdef COLVARS_LAMMPS
  jacobi = new_Jacobi_solver(4);
#else
  jacobi = NULL;
#endif
}


colvarmodule::rotation::~rotation()
{
#ifdef COLVARS_LAMMPS
  delete reinterpret_cast<
    MathEigen::Jacobi<cvm::real,
                      cvm::vector1d<cvm::real> &,
                      cvm::matrix2d<cvm::real> &> *>(jacobi);
#endif
}


void colvarmodule::rotation::build_correlation_matrix(
                                        std::vector<cvm::atom_pos> const &pos1,
                                        std::vector<cvm::atom_pos> const &pos2)
{
  // build the correlation matrix
  size_t i;
  for (i = 0; i < pos1.size(); i++) {
    C.xx += pos1[i].x * pos2[i].x;
    C.xy += pos1[i].x * pos2[i].y;
    C.xz += pos1[i].x * pos2[i].z;
    C.yx += pos1[i].y * pos2[i].x;
    C.yy += pos1[i].y * pos2[i].y;
    C.yz += pos1[i].y * pos2[i].z;
    C.zx += pos1[i].z * pos2[i].x;
    C.zy += pos1[i].z * pos2[i].y;
    C.zz += pos1[i].z * pos2[i].z;
  }
}


void colvarmodule::rotation::compute_overlap_matrix()
{
  // build the "overlap" matrix, whose eigenvectors are stationary
  // points of the RMSD in the space of rotations
  S[0][0] =    C.xx + C.yy + C.zz;
  S[1][0] =    C.yz - C.zy;
  S[0][1] = S[1][0];
  S[2][0] =  - C.xz + C.zx ;
  S[0][2] = S[2][0];
  S[3][0] =    C.xy - C.yx;
  S[0][3] = S[3][0];
  S[1][1] =    C.xx - C.yy - C.zz;
  S[2][1] =    C.xy + C.yx;
  S[1][2] = S[2][1];
  S[3][1] =    C.xz + C.zx;
  S[1][3] = S[3][1];
  S[2][2] = - C.xx + C.yy - C.zz;
  S[3][2] =   C.yz + C.zy;
  S[2][3] = S[3][2];
  S[3][3] = - C.xx - C.yy + C.zz;
}


#ifndef COLVARS_LAMMPS
namespace {

void diagonalize_matrix(cvm::matrix2d<cvm::real> &m,
                        cvm::vector1d<cvm::real> &eigval,
                        cvm::matrix2d<cvm::real> &eigvec)
{
  eigval.resize(4);
  eigval.reset();
  eigvec.resize(4, 4);
  eigvec.reset();

  // diagonalize
  int jac_nrot = 0;
  if (NR_Jacobi::jacobi(m.c_array(), eigval.c_array(), eigvec.c_array(), &jac_nrot) !=
      COLVARS_OK) {
    cvm::error("Too many iterations in jacobi diagonalization.\n"
               "This is usually the result of an ill-defined set of atoms for "
               "rotational alignment (RMSD, rotateReference, etc).\n");
  }
  NR_Jacobi::eigsrt(eigval.c_array(), eigvec.c_array());
  // jacobi saves eigenvectors by columns
  NR_Jacobi::transpose(eigvec.c_array());

  // normalize eigenvectors
  for (size_t ie = 0; ie < 4; ie++) {
    cvm::real norm2 = 0.0;
    size_t i;
    for (i = 0; i < 4; i++) {
      norm2 += eigvec[ie][i] * eigvec[ie][i];
    }
    cvm::real const norm = cvm::sqrt(norm2);
    for (i = 0; i < 4; i++) {
      eigvec[ie][i] /= norm;
    }
  }
}

}
#endif


// Calculate the rotation, plus its derivatives

void colvarmodule::rotation::calc_optimal_rotation(
                                        std::vector<cvm::atom_pos> const &pos1,
                                        std::vector<cvm::atom_pos> const &pos2)
{
  C.reset();
  build_correlation_matrix(pos1, pos2);

  S.resize(4, 4);
  S.reset();
  compute_overlap_matrix();

  S_backup.resize(4, 4);
  S_backup = S;

  if (b_debug_gradients) {
    cvm::log("S     = "+cvm::to_str(S_backup, cvm::cv_width, cvm::cv_prec)+"\n");
  }

  S_eigval.resize(4);
  S_eigvec.resize(4, 4);

#ifdef COLVARS_LAMMPS
  MathEigen::Jacobi<cvm::real,
                    cvm::vector1d<cvm::real> &,
                    cvm::matrix2d<cvm::real> &> *ecalc =
    reinterpret_cast<MathEigen::Jacobi<cvm::real,
                                       cvm::vector1d<cvm::real> &,
                                       cvm::matrix2d<cvm::real> &> *>(jacobi);

  int ierror = ecalc->Diagonalize(S, S_eigval, S_eigvec);
  if (ierror) {
    cvm::error("Too many iterations in jacobi diagonalization.\n"
               "This is usually the result of an ill-defined set of atoms for "
               "rotational alignment (RMSD, rotateReference, etc).\n");
  }
#else
  diagonalize_matrix(S, S_eigval, S_eigvec);
#endif


  // eigenvalues and eigenvectors
  cvm::real const L0 = S_eigval[0];
  cvm::real const L1 = S_eigval[1];
  cvm::real const L2 = S_eigval[2];
  cvm::real const L3 = S_eigval[3];
  cvm::quaternion const Q0(S_eigvec[0]);
  cvm::quaternion const Q1(S_eigvec[1]);
  cvm::quaternion const Q2(S_eigvec[2]);
  cvm::quaternion const Q3(S_eigvec[3]);

  lambda = L0;
  q = Q0;

  if (cvm::rotation::monitor_crossings) {
    if (q_old.norm2() > 0.0) {
      q.match(q_old);
      if (q_old.inner(q) < (1.0 - crossing_threshold)) {
        cvm::log("Warning: one molecular orientation has changed by more than "+
                 cvm::to_str(crossing_threshold)+": discontinuous rotation ?\n");
      }
    }
    q_old = q;
  }

  if (b_debug_gradients) {
    cvm::log("L0 = "+cvm::to_str(L0, cvm::cv_width, cvm::cv_prec)+
             ", Q0 = "+cvm::to_str(Q0, cvm::cv_width, cvm::cv_prec)+
             ", Q0*Q0 = "+cvm::to_str(Q0.inner(Q0), cvm::cv_width, cvm::cv_prec)+
             "\n");
    cvm::log("L1 = "+cvm::to_str(L1, cvm::cv_width, cvm::cv_prec)+
             ", Q1 = "+cvm::to_str(Q1, cvm::cv_width, cvm::cv_prec)+
             ", Q0*Q1 = "+cvm::to_str(Q0.inner(Q1), cvm::cv_width, cvm::cv_prec)+
             "\n");
    cvm::log("L2 = "+cvm::to_str(L2, cvm::cv_width, cvm::cv_prec)+
             ", Q2 = "+cvm::to_str(Q2, cvm::cv_width, cvm::cv_prec)+
             ", Q0*Q2 = "+cvm::to_str(Q0.inner(Q2), cvm::cv_width, cvm::cv_prec)+
             "\n");
    cvm::log("L3 = "+cvm::to_str(L3, cvm::cv_width, cvm::cv_prec)+
             ", Q3 = "+cvm::to_str(Q3, cvm::cv_width, cvm::cv_prec)+
             ", Q0*Q3 = "+cvm::to_str(Q0.inner(Q3), cvm::cv_width, cvm::cv_prec)+
             "\n");
  }

  // calculate derivatives of L0 and Q0 with respect to each atom in
  // either group; note: if dS_1 is a null vector, nothing will be
  // calculated
  size_t ia;
  for (ia = 0; ia < dS_1.size(); ia++) {

    cvm::real const &a2x = pos2[ia].x;
    cvm::real const &a2y = pos2[ia].y;
    cvm::real const &a2z = pos2[ia].z;

    cvm::matrix2d<cvm::rvector> &ds_1 = dS_1[ia];

    // derivative of the S matrix
    ds_1.reset();
    ds_1[0][0].set( a2x,  a2y,  a2z);
    ds_1[1][0].set( 0.0,  a2z, -a2y);
    ds_1[0][1] = ds_1[1][0];
    ds_1[2][0].set(-a2z,  0.0,  a2x);
    ds_1[0][2] = ds_1[2][0];
    ds_1[3][0].set( a2y, -a2x,  0.0);
    ds_1[0][3] = ds_1[3][0];
    ds_1[1][1].set( a2x, -a2y, -a2z);
    ds_1[2][1].set( a2y,  a2x,  0.0);
    ds_1[1][2] = ds_1[2][1];
    ds_1[3][1].set( a2z,  0.0,  a2x);
    ds_1[1][3] = ds_1[3][1];
    ds_1[2][2].set(-a2x,  a2y, -a2z);
    ds_1[3][2].set( 0.0,  a2z,  a2y);
    ds_1[2][3] = ds_1[3][2];
    ds_1[3][3].set(-a2x, -a2y,  a2z);

    cvm::rvector                &dl0_1 = dL0_1[ia];
    cvm::vector1d<cvm::rvector> &dq0_1 = dQ0_1[ia];

    // matrix multiplications; derivatives of L_0 and Q_0 are
    // calculated using Hellmann-Feynman theorem (i.e. exploiting the
    // fact that the eigenvectors Q_i form an orthonormal basis)

    dl0_1.reset();
    for (size_t i = 0; i < 4; i++) {
      for (size_t j = 0; j < 4; j++) {
        dl0_1 += Q0[i] * ds_1[i][j] * Q0[j];
      }
    }

    dq0_1.reset();
    for (size_t p = 0; p < 4; p++) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          dq0_1[p] +=
            (Q1[i] * ds_1[i][j] * Q0[j]) / (L0-L1) * Q1[p] +
            (Q2[i] * ds_1[i][j] * Q0[j]) / (L0-L2) * Q2[p] +
            (Q3[i] * ds_1[i][j] * Q0[j]) / (L0-L3) * Q3[p];
        }
      }
    }
  }

  // do the same for the second group
  for (ia = 0; ia < dS_2.size(); ia++) {

    cvm::real const &a1x = pos1[ia].x;
    cvm::real const &a1y = pos1[ia].y;
    cvm::real const &a1z = pos1[ia].z;

    cvm::matrix2d<cvm::rvector> &ds_2 = dS_2[ia];

    ds_2.reset();
    ds_2[0][0].set( a1x,  a1y,  a1z);
    ds_2[1][0].set( 0.0, -a1z,  a1y);
    ds_2[0][1] = ds_2[1][0];
    ds_2[2][0].set( a1z,  0.0, -a1x);
    ds_2[0][2] = ds_2[2][0];
    ds_2[3][0].set(-a1y,  a1x,  0.0);
    ds_2[0][3] = ds_2[3][0];
    ds_2[1][1].set( a1x, -a1y, -a1z);
    ds_2[2][1].set( a1y,  a1x,  0.0);
    ds_2[1][2] = ds_2[2][1];
    ds_2[3][1].set( a1z,  0.0,  a1x);
    ds_2[1][3] = ds_2[3][1];
    ds_2[2][2].set(-a1x,  a1y, -a1z);
    ds_2[3][2].set( 0.0,  a1z,  a1y);
    ds_2[2][3] = ds_2[3][2];
    ds_2[3][3].set(-a1x, -a1y,  a1z);

    cvm::rvector                &dl0_2 = dL0_2[ia];
    cvm::vector1d<cvm::rvector> &dq0_2 = dQ0_2[ia];

    dl0_2.reset();
    for (size_t i = 0; i < 4; i++) {
      for (size_t j = 0; j < 4; j++) {
        dl0_2 += Q0[i] * ds_2[i][j] * Q0[j];
      }
    }

    dq0_2.reset();
    for (size_t p = 0; p < 4; p++) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          dq0_2[p] +=
            (Q1[i] * ds_2[i][j] * Q0[j]) / (L0-L1) * Q1[p] +
            (Q2[i] * ds_2[i][j] * Q0[j]) / (L0-L2) * Q2[p] +
            (Q3[i] * ds_2[i][j] * Q0[j]) / (L0-L3) * Q3[p];
        }
      }
    }

    if (b_debug_gradients) {

      cvm::matrix2d<cvm::real> S_new(4, 4);
      cvm::vector1d<cvm::real> S_new_eigval(4);
      cvm::matrix2d<cvm::real> S_new_eigvec(4, 4);

      // make an infitesimal move along each cartesian coordinate of
      // this atom, and solve again the eigenvector problem
      for (size_t comp = 0; comp < 3; comp++) {

        S_new = S_backup;
        // diagonalize the new overlap matrix
        for (size_t i = 0; i < 4; i++) {
          for (size_t j = 0; j < 4; j++) {
            S_new[i][j] +=
              colvarmodule::debug_gradients_step_size * ds_2[i][j][comp];
          }
        }

        //           cvm::log("S_new = "+cvm::to_str(cvm::to_str (S_new), cvm::cv_width, cvm::cv_prec)+"\n");

#ifdef COLVARS_LAMMPS
        ecalc->Diagonalize(S_new, S_new_eigval, S_new_eigvec);
#else
        diagonalize_matrix(S_new, S_new_eigval, S_new_eigvec);
#endif

        cvm::real const &L0_new = S_new_eigval[0];
        cvm::quaternion const Q0_new(S_new_eigvec[0]);

        cvm::real const DL0 = (dl0_2[comp]) * colvarmodule::debug_gradients_step_size;
        cvm::quaternion const DQ0(dq0_2[0][comp] * colvarmodule::debug_gradients_step_size,
                                  dq0_2[1][comp] * colvarmodule::debug_gradients_step_size,
                                  dq0_2[2][comp] * colvarmodule::debug_gradients_step_size,
                                  dq0_2[3][comp] * colvarmodule::debug_gradients_step_size);

        cvm::log(  "|(l_0+dl_0) - l_0^new|/l_0 = "+
                   cvm::to_str(cvm::fabs(L0+DL0 - L0_new)/L0, cvm::cv_width, cvm::cv_prec)+
                   ", |(q_0+dq_0) - q_0^new| = "+
                   cvm::to_str((Q0+DQ0 - Q0_new).norm(), cvm::cv_width, cvm::cv_prec)+
                   "\n");
      }
    }
  }
}



