/// -*- c++ -*-

#ifndef COLVARVALUE_H
#define COLVARVALUE_H

#include "colvarmodule.h"


/// \brief Value of a collective variable: this is a metatype which
/// can be set at runtime.  By default it is set to be a scalar
/// number, and can be treated as such in all operations (this is
/// done by most \link cvc \endlink implementations).
///
/// \link colvarvalue \endlink allows \link colvar \endlink to be
/// treat different data types.  By default, a \link colvarvalue
/// \endlink variable is a scalar number.  To use it as
/// another type, declare and initialize it as
/// \code colvarvalue x(colvarvalue::type_xxx), use \link x.type
/// (colvarvalue::type_xxx) \endlink at a later stage, or if unset,
//  assign the type with \code x = y; \endcode, provided y is correctly set.
///
/// All operators (either unary or binary) on a \link
/// colvarvalue \endlink object performs one or more checks on the
/// \link Type \endlink, except when reading from a stream, when there is no way to
/// detect the \link Type \endlink.  To use  \code is >> x; \endcode x \b MUST
/// already have a type correcly set up for properly parsing the
/// stream. No problem of course with the output streams: \code os << x;
///
/// \em Note \em on \em performance: to avoid type checks in a long array of \link
/// colvarvalue \endlink objects, use one of the the "_opt" functions or implement a new one

class colvarvalue {

public:

  /// \brief Possible types of value
  ///
  /// These three cover most possibilities of data type one can
  /// devise.  If you need to implement a new colvar with a very
  /// complex data type, it's better to put an allocatable class here
  enum Type {
    /// Undefined type
    type_notset,
    /// Scalar number, implemented as \link colvarmodule::real \endlink (default)
    type_scalar,
    /// 3-dimensional vector, implemented as \link colvarmodule::rvector \endlink
    type_vector,
    /// 3-dimensional unit vector, implemented as \link colvarmodule::rvector \endlink
    type_unitvector,
    /// 3-dimensional vector that is a derivative of a unitvector
    type_unitvectorderiv,
    /// 4-dimensional unit vector representing a rotation, implemented as \link colvarmodule::quaternion \endlink
    type_quaternion,
    /// 4-dimensional vector that is a derivative of a quaternion
    type_quaternionderiv
  };

  /// Runtime description of value types
  inline std::string static const type_desc(Type t)
  {
    switch (t) {
    case colvarvalue::type_notset:
    default:
      return "not set"; break;
    case colvarvalue::type_scalar:
      return "scalar number"; break;
    case colvarvalue::type_vector:
      return "3-dimensional vector"; break;
    case colvarvalue::type_unitvector:
      return "3-dimensional unit vector"; break;
    case colvarvalue::type_unitvectorderiv:
      return "derivative of a 3-dimensional unit vector"; break;
    case colvarvalue::type_quaternion:
      return "4-dimensional unit quaternion"; break;
    case colvarvalue::type_quaternionderiv:
      return "4-dimensional tangent vector"; break;
    }
  }

  /// User keywords for specifying value types in the configuration
  inline std::string static const type_keyword(Type t)
  {
    switch (t) {
    case colvarvalue::type_notset:
    default:
      return "not_set"; break;
    case colvarvalue::type_scalar:
      return "scalar"; break;
    case colvarvalue::type_vector:
      return "vector"; break;
    case colvarvalue::type_unitvector:
      return "unit_vector"; break;
    case colvarvalue::type_unitvectorderiv:
      return ""; break;
    case colvarvalue::type_quaternion:
      return "unit_quaternion"; break;
    case colvarvalue::type_quaternionderiv:
      return ""; break;
    }
  }

  /// Number of degrees of freedom for each type
  inline size_t static num_df(Type t)
  {
    switch (t) {
    case colvarvalue::type_notset:
    default:
      return 0; break;
    case colvarvalue::type_scalar:
      return 1; break;
    case colvarvalue::type_vector:
      return 3; break;
    case colvarvalue::type_unitvector:
      return 2; break;
    case colvarvalue::type_unitvectorderiv:
      return 2; break;
    case colvarvalue::type_quaternion:
      return 3; break;
    case colvarvalue::type_quaternionderiv:
      return 3; break;
    }
  }

  /// \brief Real data member
  cvm::real       real_value;

  /// \brief Vector data member
  cvm::rvector    rvector_value;

  /// \brief Quaternion data member
  cvm::quaternion quaternion_value;

  /// Current type of this colvarvalue object
  Type  value_type;

  static inline bool type_checking()
  {
    return true;
  }

  /// \brief Default constructor: this class defaults to a scalar
  /// number and always behaves like it unless you change its type
  inline colvarvalue()
    : real_value(0.0), value_type(type_scalar)
  {}

  /// Constructor from a type specification
  inline colvarvalue(Type const &vti)
    : value_type(vti)
  {
    reset();
  }

  /// Copy constructor from real base type
  inline colvarvalue(cvm::real const &x)
    : real_value(x), value_type(type_scalar)
  {}

  /// \brief Copy constructor from rvector base type (Note: this sets
  /// automatically a type \link type_vector \endlink , if you want a
  /// \link type_unitvector \endlink you must set it explicitly)
  inline colvarvalue(cvm::rvector const &v)
    : rvector_value(v), value_type(type_vector)
  {}

  /// \brief Copy constructor from rvector base type (additional
  /// argument to make possible to choose a \link type_unitvector
  /// \endlink
  inline colvarvalue(cvm::rvector const &v, Type const &vti)
    : rvector_value(v), value_type(vti)
  {}

  /// \brief Copy constructor from quaternion base type
  inline colvarvalue(cvm::quaternion const &q)
    : quaternion_value(q), value_type(type_quaternion)
  {}

  /// Copy constructor from another \link colvarvalue \endlink
  inline colvarvalue(colvarvalue const &x)
    : value_type(x.value_type)
  {
    reset();

    switch (x.value_type) {
    case type_scalar:
      real_value = x.real_value;
      break;
    case type_vector:
    case type_unitvector:
      rvector_value = x.rvector_value;
      break;
    case type_quaternion:
      quaternion_value = x.quaternion_value;
      break;
    case type_notset:
    default:
      break;
    }
  }


  /// Set to the null value for the data type currently defined
  void reset();

  /// \brief If the variable has constraints (e.g. unitvector or
  /// quaternion), transform it to satisfy them; this function needs
  /// to be called only when the \link colvarvalue \endlink
  /// is calculated outside of \link cvc \endlink objects
  void apply_constraints();

  /// Get the current type
  inline Type type() const
  {
    return value_type;
  }

  /// Set the type explicitly
  inline void type(Type const &vti)
  {
    reset();
    value_type = vti;
    reset();
  }

  /// Set the type after another \link colvarvalue \endlink
  inline void type(colvarvalue const &x)
  {
    reset();
    value_type = x.value_type;
    reset();
  }

  /// Make the type a derivative of the original type
  /// (constraints do not apply on time derivatives of vector values)
  inline void is_derivative();

  /// Square norm of this colvarvalue
  cvm::real norm2() const;

  /// Norm of this colvarvalue
  inline cvm::real norm() const
  {
    return std::sqrt(this->norm2());
  }

  /// \brief Return the value whose scalar product with this value is
  /// 1
  inline colvarvalue inverse() const;

  /// Square distance between this \link colvarvalue \endlink and another
  cvm::real dist2(colvarvalue const &x2) const;

  /// Derivative with respect to this \link colvarvalue \endlink of the square distance
  colvarvalue dist2_grad(colvarvalue const &x2) const;

  /// Assignment operator (type of x is checked)
  colvarvalue & operator = (colvarvalue const &x);

  void operator += (colvarvalue const &x);
  void operator -= (colvarvalue const &x);
  void operator *= (cvm::real const &a);
  void operator /= (cvm::real const &a);


  // Casting operator
  inline operator cvm::real() const
  {
    if (value_type != type_scalar) {
      error_rside(type_scalar);
    }
    return real_value;
  }

  // Casting operator
  inline operator cvm::rvector() const
  {
    if ((value_type != type_vector) &&
        (value_type != type_unitvector) &&
        (value_type != type_unitvectorderiv)) {
      error_rside(type_vector);
    }
    return rvector_value;
  }

  // Casting operator
  inline operator cvm::quaternion() const
  {
    if ((value_type != type_quaternion) &&
        (value_type != type_quaternionderiv)) {
      error_rside(type_quaternion);
    }
    return quaternion_value;
  }

  /// Special case when the variable is a real number, and all operations are defined
  inline bool is_scalar() const
  {
    return (value_type == type_scalar);
  }

  /// Ensure that the two types are the same within a binary operator
  void static check_types(colvarvalue const &x1, colvarvalue const &x2);

  /// Undefined operation
  void undef_op() const;

  /// Trying to assign this \link colvarvalue \endlink object to
  /// another object set with a different type
  void error_lside(Type const &vt) const;

  /// Trying to assign another \link colvarvalue \endlink object set
  /// with a different type to this object
  void error_rside(Type const &vt) const;

  /// Give the number of characters required to output this
  /// colvarvalue, given the current type assigned and the number of
  /// characters for a real number
  size_t output_width(size_t const &real_width) const;

  /// Formats value as a script-friendly string (space separated list)
  std::string to_simple_string() const;

  /// Parses value from a script-friendly string (space separated list)
  int from_simple_string(std::string const &s);

  // optimized routines for operations with an array; xv and inner as
  // vectors are assumed to have the same number of elements (i.e. the
  // end for inner comes together with xv_end in the loop)

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  void static inner_opt(colvarvalue const                        &x,
                        std::vector<colvarvalue>::iterator       &xv,
                        std::vector<colvarvalue>::iterator const &xv_end,
                        std::vector<cvm::real>::iterator         &inner);

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  void static inner_opt(colvarvalue const                        &x,
                        std::list<colvarvalue>::iterator         &xv,
                        std::list<colvarvalue>::iterator const   &xv_end,
                        std::vector<cvm::real>::iterator         &inner);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial, (3cos^2(w)-1)/2, of one collective variable with an
  /// array
  void static p2leg_opt(colvarvalue const                        &x,
                        std::vector<colvarvalue>::iterator       &xv,
                        std::vector<colvarvalue>::iterator const &xv_end,
                        std::vector<cvm::real>::iterator         &inner);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial of one collective variable with an array
  void static p2leg_opt(colvarvalue const                        &x,
                        std::list<colvarvalue>::iterator         &xv,
                        std::list<colvarvalue>::iterator const   &xv_end,
                        std::vector<cvm::real>::iterator         &inner);

};






inline void colvarvalue::reset()
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value = cvm::real(0.0);
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    rvector_value = cvm::rvector(0.0, 0.0, 0.0);
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value = cvm::quaternion(0.0, 0.0, 0.0, 0.0);
    break;
  case colvarvalue::type_notset:
  default:
    break;
  }
}


inline void colvarvalue::apply_constraints()
{
  switch (value_type) {
  case colvarvalue::type_scalar:
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvectorderiv:
  case colvarvalue::type_quaternionderiv:
   break;
  case colvarvalue::type_unitvector:
    rvector_value /= std::sqrt(rvector_value.norm2());
    break;
  case colvarvalue::type_quaternion:
    quaternion_value /= std::sqrt(quaternion_value.norm2());
    break;
  case colvarvalue::type_notset:
  default:
    break;
  }
}


inline void colvarvalue::is_derivative()
{
  switch (value_type) {
  case colvarvalue::type_scalar:
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvectorderiv:
  case colvarvalue::type_quaternionderiv:
    break;
  case colvarvalue::type_unitvector:
    type(colvarvalue::type_unitvectorderiv);
    break;
  case colvarvalue::type_quaternion:
    type(colvarvalue::type_quaternionderiv);
    break;
  case colvarvalue::type_notset:
  default:
    break;
  }
}



inline cvm::real colvarvalue::norm2() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    return (this->real_value)*(this->real_value);
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return (this->rvector_value).norm2();
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return (this->quaternion_value).norm2();
  case colvarvalue::type_notset:
  default:
    return 0.0;
  }
}


inline colvarvalue colvarvalue::inverse() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(1.0/real_value);
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return colvarvalue(cvm::rvector(1.0/rvector_value.x,
                                    1.0/rvector_value.y,
                                    1.0/rvector_value.z));
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(cvm::quaternion(1.0/quaternion_value.q0,
                                       1.0/quaternion_value.q1,
                                       1.0/quaternion_value.q2,
                                       1.0/quaternion_value.q3));
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
  return colvarvalue();
}

inline colvarvalue & colvarvalue::operator = (colvarvalue const &x)
{
  if (this->value_type != type_notset)
    if (this->value_type != x.value_type)
      error_lside(x.value_type);

  this->value_type = x.value_type;

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    this->real_value = x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    this->rvector_value = x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    this->quaternion_value = x.quaternion_value;
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
    break;
  }
  return *this;
}

inline void colvarvalue::operator += (colvarvalue const &x)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(*this, x);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    this->real_value += x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    this->rvector_value += x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    this->quaternion_value += x.quaternion_value;
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}

inline void colvarvalue::operator -= (colvarvalue const &x)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(*this, x);

  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value -= x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    rvector_value -= x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value -= x.quaternion_value;
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}

inline void colvarvalue::operator *= (cvm::real const &a)
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value *= a;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvectorderiv:
    rvector_value *= a;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value *= a;
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}

inline void colvarvalue::operator /= (cvm::real const &a)
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value /= a; break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    rvector_value /= a; break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value /= a; break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}


// binary operations between two colvarvalues

inline colvarvalue operator + (colvarvalue const &x1,
                               colvarvalue const &x2)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(x1.real_value + x2.real_value);
  case colvarvalue::type_vector:
    return colvarvalue(x1.rvector_value + x2.rvector_value);
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return colvarvalue(x1.rvector_value + x2.rvector_value,
                       colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x1.quaternion_value + x2.quaternion_value);
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}

inline colvarvalue operator - (colvarvalue const &x1,
                               colvarvalue const &x2)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(x1.real_value - x2.real_value);
  case colvarvalue::type_vector:
    return colvarvalue(x1.rvector_value - x2.rvector_value);
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return colvarvalue(x1.rvector_value - x2.rvector_value,
                       colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x1.quaternion_value - x2.quaternion_value);
  default:
    x1.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}


// binary operations with real numbers

inline colvarvalue operator * (cvm::real const &a,
                               colvarvalue const &x)
{
  switch (x.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(a * x.real_value);
  case colvarvalue::type_vector:
    return colvarvalue(a * x.rvector_value);
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return colvarvalue(a * x.rvector_value,
                       colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(a * x.quaternion_value);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  }
}

inline colvarvalue operator * (colvarvalue const &x,
                               cvm::real const &a)
{
  switch (x.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(x.real_value * a);
  case colvarvalue::type_vector:
    return colvarvalue(x.rvector_value * a);
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return colvarvalue(x.rvector_value * a,
                       colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x.quaternion_value * a);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  }
}

inline colvarvalue operator / (colvarvalue const &x,
                               cvm::real const &a)
{
  switch (x.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(x.real_value / a);
  case colvarvalue::type_vector:
    return colvarvalue(x.rvector_value / a);
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return colvarvalue(x.rvector_value / a,
                       colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x.quaternion_value / a);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  }
}


// inner product between two colvarvalues

inline cvm::real operator * (colvarvalue const &x1,
                             colvarvalue const &x2)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return (x1.real_value * x2.real_value);
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return (x1.rvector_value * x2.rvector_value);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    // the "*" product is the quaternion product, here the inner
    // member function is used instead
    return (x1.quaternion_value.inner(x2.quaternion_value));
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return 0.0;
  };
}



inline cvm::real colvarvalue::dist2(colvarvalue const &x2) const
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(*this, x2);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    return (this->real_value - x2.real_value)*(this->real_value - x2.real_value);
  case colvarvalue::type_vector:
    return (this->rvector_value - x2.rvector_value).norm2();
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    // angle between (*this) and x2 is the distance
    return std::acos(this->rvector_value * x2.rvector_value) * std::acos(this->rvector_value * x2.rvector_value);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    // angle between (*this) and x2 is the distance, the quaternion
    // object has it implemented internally
    return this->quaternion_value.dist2(x2.quaternion_value);
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return 0.0;
  };
}


inline colvarvalue colvarvalue::dist2_grad(colvarvalue const &x2) const
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types(*this, x2);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    return 2.0 * (this->real_value - x2.real_value);
  case colvarvalue::type_vector:
    return 2.0 * (this->rvector_value - x2.rvector_value);
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    {
      cvm::rvector const &v1 = this->rvector_value;
      cvm::rvector const &v2 = x2.rvector_value;
      cvm::real const cos_t = v1 * v2;
      cvm::real const sin_t = std::sqrt(1.0 - cos_t*cos_t);
      return colvarvalue( 2.0 * sin_t *
                          cvm::rvector((-1.0) * sin_t * v2.x +
                                       cos_t/sin_t * (v1.x - cos_t*v2.x),
                                       (-1.0) * sin_t * v2.y +
                                       cos_t/sin_t * (v1.y - cos_t*v2.y),
                                       (-1.0) * sin_t * v2.z +
                                       cos_t/sin_t * (v1.z - cos_t*v2.z)
                                       ),
                          colvarvalue::type_unitvector );
    }
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return this->quaternion_value.dist2_grad(x2.quaternion_value);
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}


inline void colvarvalue::check_types(colvarvalue const &x1,
                                     colvarvalue const &x2)
{
  if (x1.value_type != x2.value_type) {
    if (((x1.value_type == type_unitvector) &&
         (x2.value_type == type_unitvectorderiv)) ||
        ((x2.value_type == type_unitvector) &&
         (x1.value_type == type_unitvectorderiv)) ||
        ((x1.value_type == type_quaternion) &&
         (x2.value_type == type_quaternionderiv)) ||
        ((x2.value_type == type_quaternion) &&
         (x1.value_type == type_quaternionderiv))) {
      return;
    }
    cvm::error("Performing an operation between two colvar "
               "values with different types, \""+
               colvarvalue::type_desc(x1.value_type)+
               "\" and \""+
               colvarvalue::type_desc(x2.value_type)+
               "\".\n");
  }
}

inline void colvarvalue::undef_op() const
{
  cvm::error("Error: Undefined operation on a colvar of type \""+
             type_desc(this->value_type)+"\".\n");
}

inline void colvarvalue::error_rside
(colvarvalue::Type const &vt) const
{
  cvm::error("Trying to assign a colvar value with type \""+
             type_desc(this->value_type)+"\" to one with type \""+
             type_desc(vt)+"\".\n");
}

inline void colvarvalue::error_lside(colvarvalue::Type const &vt) const
{
  cvm::error("Trying to use a colvar value with type \""+
             type_desc(vt)+"\" as one of type \""+
             type_desc(this->value_type)+"\".\n");
}



inline void colvarvalue::inner_opt(colvarvalue                        const &x,
                            std::vector<colvarvalue>::iterator       &xv,
                            std::vector<colvarvalue>::iterator const &xv_end,
                            std::vector<cvm::real>::iterator         &inner)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::vector<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator    &ii = inner;

  switch (x.value_type) {
  case colvarvalue::type_scalar:
    while (xvi != xv_end) {
      *(ii++) += (xvi++)->real_value * x.real_value;
    }
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    while (xvi != xv_end) {
      *(ii++) += (xvi++)->rvector_value * x.rvector_value;
    }
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    while (xvi != xv_end) {
      *(ii++) += ((xvi++)->quaternion_value).cosine(x.quaternion_value);
    }
    break;
  default:
    x.undef_op();
  };
}

inline void colvarvalue::inner_opt(colvarvalue const                      &x,
                            std::list<colvarvalue>::iterator       &xv,
                            std::list<colvarvalue>::iterator const &xv_end,
                            std::vector<cvm::real>::iterator       &inner)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::list<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator  &ii = inner;

  switch (x.value_type) {
  case colvarvalue::type_scalar:
    while (xvi != xv_end) {
      *(ii++) += (xvi++)->real_value * x.real_value;
    }
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    while (xvi != xv_end) {
      *(ii++) += (xvi++)->rvector_value * x.rvector_value;
    }
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    while (xvi != xv_end) {
      *(ii++) += ((xvi++)->quaternion_value).cosine(x.quaternion_value);
    }
    break;
  default:
    x.undef_op();
  };
}


inline void colvarvalue::p2leg_opt(colvarvalue const                        &x,
                            std::vector<colvarvalue>::iterator       &xv,
                            std::vector<colvarvalue>::iterator const &xv_end,
                            std::vector<cvm::real>::iterator         &inner)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::vector<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator    &ii = inner;

  switch (x.value_type) {
  case colvarvalue::type_scalar:
    cvm::error("Error: cannot calculate Legendre polynomials "
               "for scalar variables.\n");
    return;
    break;
  case colvarvalue::type_vector:
    while (xvi != xv_end) {
      cvm::real const cosine =
        ((xvi)->rvector_value * x.rvector_value) /
        ((xvi)->rvector_value.norm() * x.rvector_value.norm());
      xvi++;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    while (xvi != xv_end) {
      cvm::real const cosine = (xvi++)->rvector_value * x.rvector_value;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    while (xvi != xv_end) {
      cvm::real const cosine = (xvi++)->quaternion_value.cosine(x.quaternion_value);
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  default:
    x.undef_op();
  };
}

inline void colvarvalue::p2leg_opt(colvarvalue const                        &x,
                            std::list<colvarvalue>::iterator         &xv,
                            std::list<colvarvalue>::iterator const   &xv_end,
                            std::vector<cvm::real>::iterator         &inner)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::list<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator  &ii = inner;

  switch (x.value_type) {
  case colvarvalue::type_scalar:
    cvm::error("Error: cannot calculate Legendre polynomials "
               "for scalar variables.\n");
    break;
  case colvarvalue::type_vector:
    while (xvi != xv_end) {
      cvm::real const cosine =
        ((xvi)->rvector_value * x.rvector_value) /
        ((xvi)->rvector_value.norm() * x.rvector_value.norm());
      xvi++;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    while (xvi != xv_end) {
      cvm::real const cosine = (xvi++)->rvector_value * x.rvector_value;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    while (xvi != xv_end) {
      cvm::real const cosine = (xvi++)->quaternion_value.cosine(x.quaternion_value);
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  default:
    x.undef_op();
  };
}

inline std::string colvarvalue::to_simple_string() const
{
  switch (type()) {
  case colvarvalue::type_scalar:
    return cvm::to_str(real_value, 0, cvm::cv_prec);
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return rvector_value.to_simple_string();
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return quaternion_value.to_simple_string();
    break;
  case colvarvalue::type_notset:
    undef_op();
    break;
  }
  return std::string();
}

inline int colvarvalue::from_simple_string(std::string const &s)
{
  switch (type()) {
  case colvarvalue::type_scalar:
    return ((std::istringstream(s) >> real_value)
            ? COLVARS_OK : COLVARS_ERROR);
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return rvector_value.from_simple_string(s);
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return quaternion_value.from_simple_string(s);
    break;
  case colvarvalue::type_notset:
    break;
  default:
    undef_op();
  }
  return COLVARS_ERROR;
}

inline std::ostream & operator << (std::ostream &os, colvarvalue const &x)
{
  switch (x.type()) {
  case colvarvalue::type_scalar:
    os << x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    os << x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    os << x.quaternion_value;
    break;
  case colvarvalue::type_notset:
    os << "not set"; break;
  }
  return os;
}


inline std::ostream & operator << (std::ostream &os, std::vector<colvarvalue> const &v)
{
  for (size_t i = 0; i < v.size(); i++) {
    os << v[i];
  }
  return os;
}


inline std::istream & operator >> (std::istream &is, colvarvalue &x)
{
  if (x.type() == colvarvalue::type_notset) {
    cvm::error("Trying to read from a stream a colvarvalue, "
               "which has not yet been assigned a data type.\n");
    return is;
  }

  switch (x.type()) {
  case colvarvalue::type_scalar:
    is >> x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvectorderiv:
    is >> x.rvector_value;
    break;
  case colvarvalue::type_unitvector:
    is >> x.rvector_value;
    x.apply_constraints();
    break;
  case colvarvalue::type_quaternion:
    is >> x.quaternion_value;
    x.apply_constraints();
    break;
  case colvarvalue::type_quaternionderiv:
    is >> x.quaternion_value;
    break;
  default:
    x.undef_op();
  }
  return is;
}


inline size_t colvarvalue::output_width(size_t const &real_width) const
{
  switch (this->value_type) {
  case colvarvalue::type_scalar:
    return real_width;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
  case colvarvalue::type_unitvectorderiv:
    return cvm::rvector::output_width(real_width);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return cvm::quaternion::output_width(real_width);
  case colvarvalue::type_notset:
  default:
    return 0;
  }
}



#endif
