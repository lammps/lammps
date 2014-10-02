/// -*- c++ -*-

#ifndef COLVARTYPES_H
#define COLVARTYPES_H

#include <cmath>

#ifndef PI
#define PI 3.14159265358979323846
#endif

// ----------------------------------------------------------------------
/// Linear algebra functions and data types used in the collective
/// variables implemented so far
// ----------------------------------------------------------------------


/// 1-dimensional vector of real numbers with three components
class colvarmodule::rvector {

public:

  cvm::real x, y, z;

  inline rvector()
    : x (0.0), y (0.0), z (0.0)
  {}

  inline rvector (cvm::real const &x_i,
                  cvm::real const &y_i,
                  cvm::real const &z_i)
    : x (x_i), y (y_i), z (z_i)
  {}

  inline rvector (cvm::real v)
    : x (v), y (v), z (v)
  {}

  /// \brief Set all components to a scalar value
  inline void set (cvm::real const &value = 0.0) {
    x = y = z = value;
  }

  /// \brief Set all components to zero
  inline void reset() {
    x = y = z = 0.0;
  }

  /// \brief Access cartesian components by index
  inline cvm::real & operator [] (int const &i) {
    return (i == 0) ? x : (i == 1) ? y : (i == 2) ? z : x;
  }

  /// \brief Access cartesian components by index
  inline cvm::real const & operator [] (int const &i) const {
    return (i == 0) ? x : (i == 1) ? y : (i == 2) ? z : x;
  }


  inline cvm::rvector & operator = (cvm::real const &v)
  {
    x = v;
    y = v;
    z = v;
    return *this;
  }

  inline void operator += (cvm::rvector const &v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
  }

  inline void operator -= (cvm::rvector const &v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }

  inline void operator *= (cvm::real const &v)
  {
    x *= v;
    y *= v;
    z *= v;
  }

  inline void operator /= (cvm::real const& v)
  {
    x /= v;
    y /= v;
    z /= v;
  }

  inline cvm::real norm2() const
  {
    return (x*x + y*y + z*z);
  }

  inline cvm::real norm() const
  {
    return std::sqrt (this->norm2());
  }

  inline cvm::rvector unit() const
  {
    const cvm::real n = this->norm();
    return (n > 0. ? cvm::rvector (x, y, z)/n : cvm::rvector (1., 0., 0.));
  }

  static inline size_t output_width (size_t const &real_width)
  {
    return 3*real_width + 10;
  }


  static inline cvm::rvector outer (cvm::rvector const &v1, cvm::rvector const &v2)
  {
    return cvm::rvector ( v1.y*v2.z - v2.y*v1.z,
                         -v1.x*v2.z + v2.x*v1.z,
                          v1.x*v2.y - v2.x*v1.y);
  }

  friend inline cvm::rvector operator - (cvm::rvector const &v)
  {
    return cvm::rvector (-v.x, -v.y, -v.z);
  }

  friend inline int operator == (cvm::rvector const &v1, cvm::rvector const &v2)
  {
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
  }

  friend inline int operator != (cvm::rvector const &v1, cvm::rvector const &v2)
  {
    return (v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z);
  }

  friend inline cvm::rvector operator + (cvm::rvector const &v1, cvm::rvector const &v2)
  {
    return cvm::rvector (v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
  }
  friend inline cvm::rvector operator - (cvm::rvector const &v1, cvm::rvector const &v2)
  {
    return cvm::rvector (v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
  }

  friend inline cvm::real operator * (cvm::rvector const &v1, cvm::rvector const &v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  friend inline cvm::rvector operator * (cvm::real const &a, cvm::rvector const &v)
  {
    return cvm::rvector (a*v.x, a*v.y, a*v.z);
  }

  friend inline cvm::rvector operator * (cvm::rvector const &v, cvm::real const &a)
  {
    return cvm::rvector (a*v.x, a*v.y, a*v.z);
  }

  friend inline cvm::rvector operator / (cvm::rvector const &v, cvm::real const &a)
  {
    return cvm::rvector (v.x/a, v.y/a, v.z/a);
  }


};



/// \brief Arbitrary size array (one dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T, size_t const length> class colvarmodule::vector1d
{
protected:

  /// Underlying C-array
  T *array;

public:

  /// Length of the array
  inline size_t size()
  {
    return length;
  }

  /// Default constructor
  inline vector1d (T const &t = T())
  {
    array = new T[length];
    reset();
  }

  /// Set all elements to zero
  inline void reset()
  {
    for (size_t i = 0; i < length; i++) {
      array[i] = T (0.0);
    }
  }

  /// Constructor from a 1-d C array
  inline vector1d (T const *v)
  {
    array = new T[length];
    for (size_t i = 0; i < length; i++) {
      array[i] = v[i];
    }
  }

  /// Copy constructor
  inline vector1d (vector1d<T, length> const &v)
  {
    array = new T[length];
    for (size_t i = 0; i < length; i++) {
      array[i] = v.array[i];
    }
  }

  /// Assignment
  inline vector1d<T, length> & operator = (vector1d<T, length> const &v)
  {
    for (size_t i = 0; i < length; i++) {
      this->array[i] = v.array[i];
    }
    return *this;
  }

  /// Destructor
  inline ~vector1d() {
    delete [] array;
  }

  /// Return the 1-d C array
  inline T *c_array() { return array; }

  /// Return the 1-d C array
  inline operator T *() { return array; }

  /// Inner product
  inline friend T operator * (vector1d<T, length> const &v1,
                              vector1d<T, length> const &v2)
  {
    T prod (0.0);
    for (size_t i = 0; i < length; i++) {
      prod += v1.array[i] * v2.array[i];
    }
    return prod;
  }

  /// Formatted output
  friend std::ostream & operator << (std::ostream &os,
                                     vector1d<T, length> const &v)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    os << "( ";
    for (size_t i = 0; i < length-1; i++) {
      os.width (w); os.precision (p);
      os << v.array[i] << " , ";
    }
    os.width (w); os.precision (p);
    os << v.array[length-1] << " )";
    return os;
  }

};



/// \brief Arbitrary size array (two dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T,
          size_t const outer_length,
          size_t const inner_length> class colvarmodule::matrix2d
{
protected:

  /// Underlying C array
  T **array;

public:

  /// Allocation routine, used by all constructors
  inline void alloc() {
    array = new T * [outer_length];
    for (size_t i = 0; i < outer_length; i++) {
      array[i] = new T [inner_length];
    }
  }

  /// Set all elements to zero
  inline void reset()
  {
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++) {
        array[i][j] = T (0.0);
      }
    }
  }

  /// Default constructor
  inline matrix2d()
  {
    this->alloc();
    reset();
  }

  /// Constructor from a 2-d C array
  inline matrix2d (T const **m)
  {
    this->alloc();
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++)
        array[i][j] = m[i][j];
    }
  }

  /// Copy constructor
  inline matrix2d (matrix2d<T, outer_length, inner_length> const &m)
  {
    this->alloc();
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++)
        this->array[i][j] = m.array[i][j];
    }
  }

  /// Assignment
  inline matrix2d<T, outer_length, inner_length> &
  operator = (matrix2d<T, outer_length, inner_length> const &m)
  {
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++)
        this->array[i][j] = m.array[i][j];
    }
    return *this;
  }

  /// Destructor
  inline ~matrix2d() {
    for (size_t i = 0; i < outer_length; i++) {
      delete [] array[i];
    }
    delete [] array;
  }

  /// Return the 2-d C array
  inline T **c_array() { return array; }

  /// Return the 2-d C array
  inline operator T **() { return array; }

//   /// Matrix addi(c)tion
//   inline friend matrix2d<T, outer_length, inner_length>
//   operator + (matrix2d<T, outer_length, inner_length> const &mat1,
//               matrix2d<T, outer_length, inner_length> const &mat2) {
//     matrix2d<T, outer_length, inner_length> sum;
//     for (size_t i = 0; i < outer_length; i++) {
//       for (size_t j = 0; j < inner_length; j++) {
//         sum[i][j] = mat1[i][j] + mat2[i][j];
//       }
//     }
//   }

//   /// Matrix subtraction
//   inline friend matrix2d<T, outer_length, inner_length>
//   operator - (matrix2d<T, outer_length, inner_length> const &mat1,
//               matrix2d<T, outer_length, inner_length> const &mat2) {
//     matrix2d<T, outer_length, inner_length> sum;
//     for (size_t i = 0; i < outer_length; i++) {
//       for (size_t j = 0; j < inner_length; j++) {
//         sum[i][j] = mat1[i][j] - mat2[i][j];
//       }
//     }
//   }

  /// Formatted output
  friend std::ostream & operator << (std::ostream &os,
                                     matrix2d<T, outer_length, inner_length> const &m)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    os << "(";
    for (size_t i = 0; i < outer_length; i++) {
      os << " ( ";
      for (size_t j = 0; j < inner_length-1; j++) {
        os.width (w);
        os.precision (p);
        os << m.array[i][j] << " , ";
      }
      os.width (w);
      os.precision (p);
      os << m.array[i][inner_length-1] << " )";
    }

    os << " )";
    return os;
  }

};


/// \brief 2-dimensional array of real numbers with three components
/// along each dimension (works with colvarmodule::rvector)
class colvarmodule::rmatrix
  : public colvarmodule::matrix2d<colvarmodule::real, 3, 3> {
private:

public:

  /// Return the xx element
  inline cvm::real & xx() { return array[0][0]; }
  /// Return the xy element
  inline cvm::real & xy() { return array[0][1]; }
  /// Return the xz element
  inline cvm::real & xz() { return array[0][2]; }
  /// Return the yx element
  inline cvm::real & yx() { return array[1][0]; }
  /// Return the yy element
  inline cvm::real & yy() { return array[1][1]; }
  /// Return the yz element
  inline cvm::real & yz() { return array[1][2]; }
  /// Return the zx element
  inline cvm::real & zx() { return array[2][0]; }
  /// Return the zy element
  inline cvm::real & zy() { return array[2][1]; }
  /// Return the zz element
  inline cvm::real & zz() { return array[2][2]; }

  /// Return the xx element
  inline cvm::real xx() const { return array[0][0]; }
  /// Return the xy element
  inline cvm::real xy() const { return array[0][1]; }
  /// Return the xz element
  inline cvm::real xz() const { return array[0][2]; }
  /// Return the yx element
  inline cvm::real yx() const { return array[1][0]; }
  /// Return the yy element
  inline cvm::real yy() const { return array[1][1]; }
  /// Return the yz element
  inline cvm::real yz() const { return array[1][2]; }
  /// Return the zx element
  inline cvm::real zx() const { return array[2][0]; }
  /// Return the zy element
  inline cvm::real zy() const { return array[2][1]; }
  /// Return the zz element
  inline cvm::real zz() const { return array[2][2]; }

  /// Constructor from a 2-d C array
  inline rmatrix (cvm::real const **m)
    : cvm::matrix2d<cvm::real, 3, 3> (m)
  {}

  /// Default constructor
  inline rmatrix()
    : cvm::matrix2d<cvm::real, 3, 3>()
  {}

  /// Constructor component by component
  inline rmatrix (cvm::real const &xxi,
                  cvm::real const &xyi,
                  cvm::real const &xzi,
                  cvm::real const &yxi,
                  cvm::real const &yyi,
                  cvm::real const &yzi,
                  cvm::real const &zxi,
                  cvm::real const &zyi,
                  cvm::real const &zzi)
    : cvm::matrix2d<cvm::real, 3, 3>()
  {
    this->xx() = xxi;
    this->xy() = xyi;
    this->xz() = xzi;
    this->yx() = yxi;
    this->yy() = yyi;
    this->yz() = yzi;
    this->zx() = zxi;
    this->zy() = zyi;
    this->zz() = zzi;
  }

  /// Destructor
  inline ~rmatrix()
  {}

  /// Return the determinant
  inline cvm::real determinant() const
  {
    return
      (  xx() * (yy()*zz() - zy()*yz()))
      - (yx() * (xy()*zz() - zy()*xz()))
      + (zx() * (xy()*yz() - yy()*xz()));
  }

  inline cvm::rmatrix transpose() const
  {
    return cvm::rmatrix (this->xx(),
                         this->yx(),
                         this->zx(),
                         this->xy(),
                         this->yy(),
                         this->zy(),
                         this->xz(),
                         this->yz(),
                         this->zz());
  }

  friend cvm::rvector operator * (cvm::rmatrix const &m, cvm::rvector const &r);

  //   friend inline cvm::rmatrix const operator * (cvm::rmatrix const &m1, cvm::rmatrix const &m2) {
  //     return cvm::rmatrix (m1.xx()*m2.xx() + m1.xy()*m2.yx() + m1.xz()*m2.yz(),
  //                     m1.xx()*m2.xy() + m1.xy()*m2.yy() + m1.xz()*m2.zy(),
  //                     m1.xx()*m2.xz() + m1.xy()*m2.yz() + m1.xz()*m2.zz(),
  //                     m1.yx()*m2.xx() + m1.yy()*m2.yx() + m1.yz()*m2.yz(),
  //                     m1.yx()*m2.xy() + m1.yy()*m2.yy() + m1.yz()*m2.yy(),
  //                     m1.yx()*m2.xz() + m1.yy()*m2.yz() + m1.yz()*m2.yz(),
  //                     m1.zx()*m2.xx() + m1.zy()*m2.yx() + m1.zz()*m2.yz(),
  //                     m1.zx()*m2.xy() + m1.zy()*m2.yy() + m1.zz()*m2.yy(),
  //                     m1.zx()*m2.xz() + m1.zy()*m2.yz() + m1.zz()*m2.yz());
  //   }

};


inline cvm::rvector operator * (cvm::rmatrix const &m,
                                cvm::rvector const &r)
{
  return cvm::rvector (m.xx()*r.x + m.xy()*r.y + m.xz()*r.z,
                       m.yx()*r.x + m.yy()*r.y + m.yz()*r.z,
                       m.zx()*r.x + m.zy()*r.y + m.zz()*r.z);
}


/// Numerical recipes diagonalization
void jacobi (cvm::real **a, cvm::real d[], cvm::real **v, int *nrot);

/// Eigenvector sort
void eigsrt (cvm::real d[], cvm::real **v);

/// Transpose the matrix
void transpose (cvm::real **v);




/// \brief 1-dimensional vector of real numbers with four components and
/// a quaternion algebra
class colvarmodule::quaternion {

public:

  cvm::real q0, q1, q2, q3;

  /// Constructor from a 3-d vector
  inline quaternion (cvm::real const &x, cvm::real const &y, cvm::real const &z)
    : q0 (0.0), q1 (x), q2 (y), q3 (z)
  {}

  /// Constructor component by component
  inline quaternion (cvm::real const qv[4])
    : q0 (qv[0]), q1 (qv[1]), q2 (qv[2]), q3 (qv[3])
  {}

  /// Constructor component by component
  inline quaternion (cvm::real const &q0i,
                     cvm::real const &q1i,
                     cvm::real const &q2i,
                     cvm::real const &q3i)
    : q0 (q0i), q1 (q1i), q2 (q2i), q3 (q3i)
  {}

  /// "Constructor" after Euler angles (in radians)
  ///
  /// http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
  inline void set_from_euler_angles (cvm::real const &phi_in,
                                     cvm::real const &theta_in,
                                     cvm::real const &psi_in)
  {
    q0 = ( (std::cos (phi_in/2.0)) * (std::cos (theta_in/2.0)) * (std::cos (psi_in/2.0)) +
           (std::sin (phi_in/2.0)) * (std::sin (theta_in/2.0)) * (std::sin (psi_in/2.0)) );

    q1 = ( (std::sin (phi_in/2.0)) * (std::cos (theta_in/2.0)) * (std::cos (psi_in/2.0)) -
           (std::cos (phi_in/2.0)) * (std::sin (theta_in/2.0)) * (std::sin (psi_in/2.0)) );

    q2 = ( (std::cos (phi_in/2.0)) * (std::sin (theta_in/2.0)) * (std::cos (psi_in/2.0)) +
           (std::sin (phi_in/2.0)) * (std::cos (theta_in/2.0)) * (std::sin (psi_in/2.0)) );

    q3 = ( (std::cos (phi_in/2.0)) * (std::cos (theta_in/2.0)) * (std::sin (psi_in/2.0)) -
           (std::sin (phi_in/2.0)) * (std::sin (theta_in/2.0)) * (std::cos (psi_in/2.0)) );
  }

  /// \brief Default constructor
  inline quaternion()
  {
    reset();
  }

  /// \brief Set all components to a scalar
  inline void set (cvm::real const &value = 0.0)
  {
    q0 = q1 = q2 = q3 = value;
  }

  /// \brief Set all components to zero (null quaternion)
  inline void reset()
  {
    q0 = q1 = q2 = q3 = 0.0;
  }

  /// \brief Set the q0 component to 1 and the others to 0 (quaternion
  /// representing no rotation)
  inline void reset_rotation()
  {
    q0 = 1.0;
    q1 = q2 = q3 = 0.0;
  }

  /// Tell the number of characters required to print a quaternion, given that of a real number
  static inline size_t output_width (size_t const &real_width)
  {
    return 4*real_width + 13;
  }

  /// \brief Formatted output operator
  friend std::ostream & operator << (std::ostream &os, cvm::quaternion const &q);
  /// \brief Formatted input operator
  friend std::istream & operator >> (std::istream &is, cvm::quaternion &q);

  /// Access the quaternion as a 4-d array (return a reference)
  inline cvm::real & operator [] (int const &i) {
    switch (i) {
    case 0:
      return this->q0;
    case 1:
      return this->q1;
    case 2:
      return this->q2;
    case 3:
      return this->q3;
    default:
      cvm::error ("Error: incorrect quaternion component.\n");
      return q0;
    }
  }

  /// Access the quaternion as a 4-d array (return a value)
  inline cvm::real operator [] (int const &i) const {
    switch (i) {
    case 0:
      return this->q0;
    case 1:
      return this->q1;
    case 2:
      return this->q2;
    case 3:
      return this->q3;
    default:
      cvm::error ("Error: trying to access a quaternion "
                 "component which is not between 0 and 3.\n");
      return 0.0;
    }
  }

  /// Square norm of the quaternion
  inline cvm::real norm2() const
  {
    return q0*q0 + q1*q1 + q2*q2 + q3*q3;
  }

  /// Norm of the quaternion
  inline cvm::real norm() const
  {
    return std::sqrt (this->norm2());
  }

  /// Return the conjugate quaternion
  inline cvm::quaternion conjugate() const
  {
    return cvm::quaternion (q0, -q1, -q2, -q3);
  }

  inline void operator *= (cvm::real const &a)
  {
    q0 *= a; q1 *= a; q2 *= a; q3 *= a;
  }

  inline void operator /= (cvm::real const &a)
  {
    q0 /= a; q1 /= a; q2 /= a; q3 /= a;
  }

  inline void set_positive()
  {
    if (q0 > 0.0) return;
    q0 = -q0;
    q1 = -q1;
    q2 = -q2;
    q3 = -q3;
  }

  inline void operator += (cvm::quaternion const &h)
  {
    q0+=h.q0; q1+=h.q1; q2+=h.q2; q3+=h.q3;
  }
  inline void operator -= (cvm::quaternion const &h)
  {
    q0-=h.q0; q1-=h.q1; q2-=h.q2; q3-=h.q3;
  }

  /// Promote a 3-vector to a quaternion
  static inline cvm::quaternion promote (cvm::rvector const &v)
  {
    return cvm::quaternion (0.0, v.x, v.y, v.z);
  }
  /// Return the vector component
  inline cvm::rvector get_vector() const
  {
    return cvm::rvector (q1, q2, q3);
  }


  friend inline cvm::quaternion operator + (cvm::quaternion const &h, cvm::quaternion const &q)
  {
    return cvm::quaternion (h.q0+q.q0, h.q1+q.q1, h.q2+q.q2, h.q3+q.q3);
  }

  friend inline cvm::quaternion operator - (cvm::quaternion const &h, cvm::quaternion const &q)
  {
    return cvm::quaternion (h.q0-q.q0, h.q1-q.q1, h.q2-q.q2, h.q3-q.q3);
  }

  /// \brief Provides the quaternion product.  \b NOTE: for the inner
  /// product use: \code h.inner (q); \endcode
  friend inline cvm::quaternion operator * (cvm::quaternion const &h, cvm::quaternion const &q)
  {
    return cvm::quaternion (h.q0*q.q0 - h.q1*q.q1 - h.q2*q.q2 - h.q3*q.q3,
                            h.q0*q.q1 + h.q1*q.q0 + h.q2*q.q3 - h.q3*q.q2,
                            h.q0*q.q2 + h.q2*q.q0 + h.q3*q.q1 - h.q1*q.q3,
                            h.q0*q.q3 + h.q3*q.q0 + h.q1*q.q2 - h.q2*q.q1);
  }

  friend inline cvm::quaternion operator * (cvm::real const &c,
                                            cvm::quaternion const &q)
  {
    return cvm::quaternion (c*q.q0, c*q.q1, c*q.q2, c*q.q3);
  }
  friend inline cvm::quaternion operator * (cvm::quaternion const &q,
                                            cvm::real const &c)
  {
    return cvm::quaternion (q.q0*c, q.q1*c, q.q2*c, q.q3*c);
  }
  friend inline cvm::quaternion operator / (cvm::quaternion const &q,
                                            cvm::real const &c)
  {
    return cvm::quaternion (q.q0/c, q.q1/c, q.q2/c, q.q3/c);
  }


  /// \brief Rotate v through this quaternion (put it in the rotated
  /// reference frame)
  inline cvm::rvector rotate (cvm::rvector const &v) const
  {
    return ((*this) * promote (v) * ((*this).conjugate())).get_vector();
  }

  /// \brief Rotate Q2 through this quaternion (put it in the rotated
  /// reference frame)
  inline cvm::quaternion rotate (cvm::quaternion const &Q2) const
  {
    cvm::rvector const vq_rot = this->rotate (Q2.get_vector());
    return cvm::quaternion (Q2.q0, vq_rot.x, vq_rot.y, vq_rot.z);
  }

  /// Return the 3x3 matrix associated to this quaternion
  inline cvm::rmatrix rotation_matrix() const
  {
    cvm::rmatrix R;

    R.xx() = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R.yy() = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    R.zz() = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    R.xy() = 2.0 * (q1*q2 - q0*q3);
    R.xz() = 2.0 * (q0*q2 + q1*q3);

    R.yx() = 2.0 * (q0*q3 + q1*q2);
    R.yz() = 2.0 * (q2*q3 - q0*q1);

    R.zx() = 2.0 * (q1*q3 - q0*q2);
    R.zy() = 2.0 * (q0*q1 + q2*q3);

    return R;
  }


  /// \brief Multiply the given vector by the derivative of the given
  /// (rotated) position with respect to the quaternion
  cvm::quaternion position_derivative_inner (cvm::rvector const &pos,
                                             cvm::rvector const &vec) const;


  /// \brief Return the cosine between the orientation frame
  /// associated to this quaternion and another
  inline cvm::real cosine (cvm::quaternion const &q) const
  {
    cvm::real const iprod = this->inner (q);
    return 2.0*iprod*iprod - 1.0;
  }

  /// \brief Square distance from another quaternion on the
  /// 4-dimensional unit sphere: returns the square of the angle along
  /// the shorter of the two geodesics
  inline cvm::real dist2 (cvm::quaternion const &Q2) const
  {
    cvm::real const cos_omega = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;

    cvm::real const omega = std::acos ( (cos_omega > 1.0) ? 1.0 :
                                     ( (cos_omega < -1.0) ? -1.0 : cos_omega) );

    // get the minimum distance: x and -x are the same quaternion
    if (cos_omega > 0.0)
      return omega * omega;
    else
      return (PI-omega) * (PI-omega);
  }

  /// Gradient of the square distance: returns a 4-vector equivalent
  /// to that provided by slerp
  inline cvm::quaternion dist2_grad (cvm::quaternion const &Q2) const
  {
    cvm::real const cos_omega = this->q0*Q2.q0 + this->q1*Q2.q1 + this->q2*Q2.q2 + this->q3*Q2.q3;
    cvm::real const omega = std::acos ( (cos_omega > 1.0) ? 1.0 :
                                     ( (cos_omega < -1.0) ? -1.0 : cos_omega) );
    cvm::real const sin_omega = std::sin (omega);

    if (std::fabs (sin_omega) < 1.0E-14) {
      // return a null 4d vector
      return cvm::quaternion (0.0, 0.0, 0.0, 0.0);
    }

    cvm::quaternion const
      grad1 ((-1.0)*sin_omega*Q2.q0 + cos_omega*(this->q0-cos_omega*Q2.q0)/sin_omega,
             (-1.0)*sin_omega*Q2.q1 + cos_omega*(this->q1-cos_omega*Q2.q1)/sin_omega,
             (-1.0)*sin_omega*Q2.q2 + cos_omega*(this->q2-cos_omega*Q2.q2)/sin_omega,
             (-1.0)*sin_omega*Q2.q3 + cos_omega*(this->q3-cos_omega*Q2.q3)/sin_omega);

    if (cos_omega > 0.0) {
      return 2.0*omega*grad1;
    }
    else {
      return -2.0*(PI-omega)*grad1;
    }
  }

  /// \brief Choose the closest between Q2 and -Q2 and save it back.
  /// Not required for dist2() and dist2_grad()
  inline void match (cvm::quaternion &Q2) const
  {
    cvm::real const cos_omega = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;
    if (cos_omega < 0.0) Q2 *= -1.0;
  }

  /// \brief Inner product (as a 4-d vector) with Q2; requires match()
  /// if the largest overlap is looked for
  inline cvm::real inner (cvm::quaternion const &Q2) const
  {
    cvm::real const prod = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;
    return prod;
  }


};


/// \brief A rotation between two sets of coordinates (for the moment
/// a wrapper for colvarmodule::quaternion)
class colvarmodule::rotation
{
public:

  /// \brief The rotation itself (implemented as a quaternion)
  cvm::quaternion q;

  /// \brief Eigenvalue corresponding to the optimal rotation
  cvm::real lambda;

  /// \brief Perform gradient tests
  bool b_debug_gradients;

  /// \brief Positions to superimpose: the rotation should brings pos1
  /// into pos2
  std::vector< cvm::atom_pos > pos1, pos2;

  /// Derivatives of S
  std::vector< cvm::matrix2d<cvm::rvector, 4, 4> > dS_1,  dS_2;
  /// Derivatives of leading eigenvalue
  std::vector< cvm::rvector >                      dL0_1, dL0_2;
  /// Derivatives of leading eigenvector
  std::vector< cvm::vector1d<cvm::rvector, 4> >    dQ0_1, dQ0_2;

  /// Allocate space for the derivatives of the rotation
  inline void request_group1_gradients (size_t const &n)
  {
    dS_1.resize  (n, cvm::matrix2d<cvm::rvector, 4, 4>());
    dL0_1.resize (n, cvm::rvector (0.0, 0.0, 0.0));
    dQ0_1.resize (n, cvm::vector1d<cvm::rvector, 4>());
  }

  inline void request_group2_gradients (size_t const &n)
  {
    dS_2.resize  (n, cvm::matrix2d<cvm::rvector, 4, 4>());
    dL0_2.resize (n, cvm::rvector (0.0, 0.0, 0.0));
    dQ0_2.resize (n, cvm::vector1d<cvm::rvector, 4>());
  }

  /// \brief Calculate the optimal rotation and store the
  /// corresponding eigenvalue and eigenvector in the arguments l0 and
  /// q0; if the gradients have been previously requested, calculate
  /// them as well
  ///
  /// The method to derive the optimal rotation is defined in:
  /// Coutsias EA, Seok C, Dill KA.
  /// Using quaternions to calculate RMSD.
  /// J Comput Chem. 25(15):1849-57 (2004)
  /// DOI: 10.1002/jcc.20110  PubMed: 15376254
  void calc_optimal_rotation (std::vector<atom_pos> const &pos1,
                              std::vector<atom_pos> const &pos2);

  /// Default constructor
  inline rotation()
    : b_debug_gradients (false)
  {}

  /// Constructor after a quaternion
  inline rotation (cvm::quaternion const &qi)
    : q (qi),
      b_debug_gradients (false)
  {
  }

  /// Constructor after an axis of rotation and an angle (in radians)
  inline rotation (cvm::real const &angle, cvm::rvector const &axis)
    : b_debug_gradients (false)
  {
    cvm::rvector const axis_n = axis.unit();
    cvm::real const sina = std::sin (angle/2.0);
    q = cvm::quaternion (std::cos (angle/2.0),
                         sina * axis_n.x, sina * axis_n.y, sina * axis_n.z);
  }

  /// Destructor
  inline ~rotation()
  {}

  /// Return the rotated vector
  inline cvm::rvector rotate (cvm::rvector const &v) const
  {
    return q.rotate (v);
  }

  /// Return the inverse of this rotation
  inline cvm::rotation inverse() const
  {
    return cvm::rotation (this->q.conjugate());
  }

  /// Return the associated 3x3 matrix
  inline cvm::rmatrix matrix() const
  {
    return q.rotation_matrix();
  }


  /// \brief Return the spin angle (in degrees) with respect to the
  /// provided axis (which MUST be normalized)
  inline cvm::real spin_angle (cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real alpha = (180.0/PI) * 2.0 * std::atan2 (axis * q_vec, q.q0);
    while (alpha >  180.0) alpha -= 360;
    while (alpha < -180.0) alpha += 360;
    return alpha;
  }

  /// \brief Return the derivative of the spin angle with respect to
  /// the quaternion
  inline cvm::quaternion dspin_angle_dq (cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real const iprod = axis * q_vec;

    if (q.q0 != 0.0) {

      // cvm::real const x = iprod/q.q0;

      cvm::real const dspindx = (180.0/PI) * 2.0 * (1.0 / (1.0 + (iprod*iprod)/(q.q0*q.q0)));

      return
        cvm::quaternion ( dspindx * (iprod * (-1.0) / (q.q0*q.q0)),
                          dspindx * ((1.0/q.q0) * axis.x),
                          dspindx * ((1.0/q.q0) * axis.y),
                          dspindx * ((1.0/q.q0) * axis.z));
    } else {
      // (1/(1+x^2)) ~ (1/x)^2
      return
        cvm::quaternion ((180.0/PI) * 2.0 * ((-1.0)/iprod), 0.0, 0.0, 0.0);
      // XX TODO: What if iprod == 0? XX
    }
  }

  /// \brief Return the projection of the orientation vector onto a
  /// predefined axis
  inline cvm::real cos_theta (cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real const alpha =
      (180.0/PI) * 2.0 * std::atan2 (axis * q_vec, q.q0);

    cvm::real const cos_spin_2 = std::cos (alpha * (PI/180.0) * 0.5);
    cvm::real const cos_theta_2 = ( (cos_spin_2 != 0.0) ?
                                    (q.q0 / cos_spin_2) :
                                    (0.0) );
    // cos(2t) = 2*cos(t)^2 - 1
    return 2.0 * (cos_theta_2*cos_theta_2) - 1.0;
  }

   /// Return the derivative of the tilt wrt the quaternion
  inline cvm::quaternion dcos_theta_dq (cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real const iprod = axis * q_vec;

    cvm::real const cos_spin_2 = std::cos (std::atan2 (iprod, q.q0));

    if (q.q0 != 0.0)  {

      cvm::real const d_cos_theta_dq0 =
        (4.0 * q.q0 / (cos_spin_2*cos_spin_2)) *
        (1.0 - (iprod*iprod)/(q.q0*q.q0) / (1.0 + (iprod*iprod)/(q.q0*q.q0)));

      cvm::real const d_cos_theta_dqn =
        (4.0 * q.q0 / (cos_spin_2*cos_spin_2) *
         (iprod/q.q0) / (1.0 + (iprod*iprod)/(q.q0*q.q0)));

      return cvm::quaternion (d_cos_theta_dq0,
                              d_cos_theta_dqn * axis.x,
                              d_cos_theta_dqn * axis.y,
                              d_cos_theta_dqn * axis.z);
    } else {

      cvm::real const d_cos_theta_dqn =
        (4.0 / (cos_spin_2*cos_spin_2 * iprod));

      return cvm::quaternion (0.0,
                              d_cos_theta_dqn * axis.x,
                              d_cos_theta_dqn * axis.y,
                              d_cos_theta_dqn * axis.z);
    }
  }

  /// \brief Threshold for the eigenvalue crossing test
  static cvm::real crossing_threshold;

protected:

  /// \brief Previous value of the rotation (used to warn the user
  /// when the structure changes too much, and there may be an
  /// eigenvalue crossing)
  cvm::quaternion q_old;

  /// Build the overlap matrix S (used by calc_optimal_rotation())
  void build_matrix (std::vector<cvm::atom_pos> const &pos1,
                     std::vector<cvm::atom_pos> const &pos2,
                     cvm::matrix2d<real, 4, 4>        &S);

  /// Diagonalize the overlap matrix S (used by calc_optimal_rotation())
  void diagonalize_matrix (cvm::matrix2d<cvm::real, 4, 4> &S,
                           cvm::real                       S_eigval[4],
                           cvm::matrix2d<cvm::real, 4, 4> &S_eigvec);
};


#endif
