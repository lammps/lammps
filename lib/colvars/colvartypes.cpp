// -*- c++ -*-

#include <stdlib.h>
#include <string.h>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarparse.h"


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
  size_t const start_pos = is.tellg();
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
  size_t const start_pos = is.tellg();

  std::string euler("");

  if ( (is >> euler) && (colvarparse::to_lower_cppstr(euler) ==
                         std::string("euler")) ) {

    // parse the Euler angles

    char sep;
    cvm::real phi, theta, psi;
    if ( !(is >> sep)   || !(sep == '(') ||
         !(is >> phi)   || !(is >> sep)  || !(sep == ',') ||
         !(is >> theta) || !(is >> sep)  || !(sep == ',') ||
         !(is >> psi)   || !(is >> sep)  || !(sep == ')') ) {
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    }

    q = colvarmodule::quaternion(phi, theta, psi);

  } else {

    // parse the quaternion components

    is.seekg(start_pos, std::ios::beg);
    char sep;
    if ( !(is >> sep)  || !(sep == '(') ||
         !(is >> q.q0) || !(is >> sep)  || !(sep == ',') ||
         !(is >> q.q1) || !(is >> sep)  || !(sep == ',') ||
         !(is >> q.q2) || !(is >> sep)  || !(sep == ',') ||
         !(is >> q.q3) || !(is >> sep)  || !(sep == ')') ) {
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    }
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

void colvarmodule::rotation::build_matrix(std::vector<cvm::atom_pos> const &pos1,
                                          std::vector<cvm::atom_pos> const &pos2,
                                          cvm::matrix2d<cvm::real>         &S)
{
  // build the correlation matrix
  C.resize(3, 3);
  C.reset();
  size_t i;
  for (i = 0; i < pos1.size(); i++) {
    C.xx() += pos1[i].x * pos2[i].x;
    C.xy() += pos1[i].x * pos2[i].y;
    C.xz() += pos1[i].x * pos2[i].z;
    C.yx() += pos1[i].y * pos2[i].x;
    C.yy() += pos1[i].y * pos2[i].y;
    C.yz() += pos1[i].y * pos2[i].z;
    C.zx() += pos1[i].z * pos2[i].x;
    C.zy() += pos1[i].z * pos2[i].y;
    C.zz() += pos1[i].z * pos2[i].z;
  }

  // build the "overlap" matrix, whose eigenvectors are stationary
  // points of the RMSD in the space of rotations
  S[0][0] =    C.xx() + C.yy() + C.zz();
  S[1][0] =    C.yz() - C.zy();
  S[0][1] = S[1][0];
  S[2][0] =  - C.xz() + C.zx() ;
  S[0][2] = S[2][0];
  S[3][0] =    C.xy() - C.yx();
  S[0][3] = S[3][0];
  S[1][1] =    C.xx() - C.yy() - C.zz();
  S[2][1] =    C.xy() + C.yx();
  S[1][2] = S[2][1];
  S[3][1] =    C.xz() + C.zx();
  S[1][3] = S[3][1];
  S[2][2] = - C.xx() + C.yy() - C.zz();
  S[3][2] =   C.yz() + C.zy();
  S[2][3] = S[3][2];
  S[3][3] = - C.xx() - C.yy() + C.zz();
}


void colvarmodule::rotation::diagonalize_matrix(cvm::matrix2d<cvm::real> &S,
                                                cvm::vector1d<cvm::real> &S_eigval,
                                                cvm::matrix2d<cvm::real> &S_eigvec)
{
  S_eigval.resize(4);
  S_eigval.reset();
  S_eigvec.resize(4,4);
  S_eigvec.reset();

  // diagonalize
  int jac_nrot = 0;
  jacobi(S.c_array(), S_eigval.c_array(), S_eigvec.c_array(), &jac_nrot);
  eigsrt(S_eigval.c_array(), S_eigvec.c_array());
  // jacobi saves eigenvectors by columns
  transpose(S_eigvec.c_array());

  // normalize eigenvectors
  for (size_t ie = 0; ie < 4; ie++) {
    cvm::real norm2 = 0.0;
    size_t i;
    for (i = 0; i < 4; i++) {
      norm2 += std::pow(S_eigvec[ie][i], int(2));
    }
    cvm::real const norm = std::sqrt(norm2);
    for (i = 0; i < 4; i++) {
      S_eigvec[ie][i] /= norm;
    }
  }
}


// Calculate the rotation, plus its derivatives

void colvarmodule::rotation::calc_optimal_rotation(std::vector<cvm::atom_pos> const &pos1,
                                                   std::vector<cvm::atom_pos> const &pos2)
{
  S.resize(4,4);
  S.reset();

  build_matrix(pos1, pos2, S);

  S_backup.resize(4,4);
  S_backup = S;

  if (cvm::debug()) {
    if (b_debug_gradients) {
      cvm::log("S     = "+cvm::to_str(cvm::to_str(S_backup), cvm::cv_width, cvm::cv_prec)+"\n");
    }
  }

  diagonalize_matrix(S, S_eigval, S_eigvec);

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

  if (q_old.norm2() > 0.0) {
    q.match(q_old);
    if (q_old.inner(q) < (1.0 - crossing_threshold)) {
      cvm::log("Warning: one molecular orientation has changed by more than "+
               cvm::to_str(crossing_threshold)+": discontinuous rotation ?\n");
    }
  }
  q_old = q;

  if (cvm::debug()) {
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

    if (cvm::debug()) {

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

          diagonalize_matrix(S_new, S_new_eigval, S_new_eigvec);

          cvm::real const &L0_new = S_new_eigval[0];
          cvm::quaternion const Q0_new(S_new_eigvec[0]);

          cvm::real const DL0 = (dl0_2[comp]) * colvarmodule::debug_gradients_step_size;
          cvm::quaternion const q0(Q0);
          cvm::quaternion const DQ0(dq0_2[0][comp] * colvarmodule::debug_gradients_step_size,
                                    dq0_2[1][comp] * colvarmodule::debug_gradients_step_size,
                                    dq0_2[2][comp] * colvarmodule::debug_gradients_step_size,
                                    dq0_2[3][comp] * colvarmodule::debug_gradients_step_size);

          cvm::log(  "|(l_0+dl_0) - l_0^new|/l_0 = "+
                     cvm::to_str(std::fabs(L0+DL0 - L0_new)/L0, cvm::cv_width, cvm::cv_prec)+
                     ", |(q_0+dq_0) - q_0^new| = "+
                     cvm::to_str((Q0+DQ0 - Q0_new).norm(), cvm::cv_width, cvm::cv_prec)+
                     "\n");
        }
      }
    }
  }
}



// Numerical Recipes routine for diagonalization

#define ROTATE(a,i,j,k,l) g=a[i][j]; \
  h=a[k][l];                         \
  a[i][j]=g-s*(h+g*tau);             \
  a[k][l]=h+s*(g-h*tau);

#define n 4

void jacobi(cvm::real **a, cvm::real *d, cvm::real **v, int *nrot)
{
  int j,iq,ip,i;
  cvm::real tresh,theta,tau,t,sm,s,h,g,c;

  cvm::vector1d<cvm::real> b(n);
  cvm::vector1d<cvm::real> z(n);

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) {
      v[ip][iq]=0.0;
    }
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=0;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
        sm += std::fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
        g=100.0*std::fabs(a[ip][iq]);
        if (i > 4 && (cvm::real)(std::fabs(d[ip])+g) == (cvm::real)std::fabs(d[ip])
            && (cvm::real)(std::fabs(d[iq])+g) == (cvm::real)std::fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (std::fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((cvm::real)(std::fabs(h)+g) == (cvm::real)std::fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(std::fabs(theta)+std::sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/std::sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0;j<=ip-1;j++) {
            ROTATE(a,j,ip,j,iq)
              }
          for (j=ip+1;j<=iq-1;j++) {
            ROTATE(a,ip,j,j,iq)
              }
          for (j=iq+1;j<n;j++) {
            ROTATE(a,ip,j,iq,j)
              }
          for (j=0;j<n;j++) {
            ROTATE(v,j,ip,j,iq)
              }
          ++(*nrot);
        }
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  cvm::error("Too many iterations in routine jacobi.\n");
}

void eigsrt(cvm::real *d, cvm::real **v)
{
  int k,j,i;
  cvm::real p;

  for (i=0;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=0;j<n;j++) {
        p=v[j][i];
        v[j][i]=v[j][k];
        v[j][k]=p;
      }
    }
  }
}

void transpose(cvm::real **v)
{
  cvm::real p;
  int i,j;
  for (i=0;i<n;i++) {
    for (j=i+1;j<n;j++) {
      p=v[i][j];
      v[i][j]=v[j][i];
      v[j][i]=p;
    }
  }
}

#undef n
#undef ROTATE
