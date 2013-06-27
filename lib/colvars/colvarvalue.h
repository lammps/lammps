#ifndef COLVARVALUE_H
#define COLVARVALUE_H

#include "colvarmodule.h"


/// \brief Value of a collective variable: this is a metatype which
/// can be set at runtime.  By default it is set to be a scalar
/// number, and can be treated like that in all operations (this is
/// done by most \link cvc \endlink implementations).
///
/// \link colvarvalue \endlink allows \link colvar \endlink to be
/// treat different data types.  By default, a \link colvarvalue
/// \endlink variable is a scalar number.  If you want to use it as
/// another type, you should declare and initialize a variable as
/// \code colvarvalue x (colvarvalue::type_xxx); \endcode where
/// type_xxx is a value within the \link Type \endlink enum.
/// Alternatively, initialize x with \link x.type
/// (colvarvalue::type_xxx) \endlink at a later stage.
///
/// Given a colvarvalue variable x which is not yet assigned (and
/// thus has not yet a type) it is also possible to correctly assign
/// the type with \code x = y; \endcode if y is correctly set.
/// Otherwise, an error will be raised if the \link Type \endlink of x
/// is different from the \link Type \endlink of y.
///
/// Also, every operator (either unary or binary) on a \link
/// colvarvalue \endlink object performs one or more checks on the
/// \link Type \endlink to avoid errors such as initializing a
/// three-dimensional vector with a scalar number (legal otherwise).
///
/// \b Special \b case: when reading from a stream, there is no way to
/// detect the \link Type \endlink and safely handle errors at the
/// same time.  Hence, when using \code is >> x; \endcode x \b MUST
/// already have a type correcly set up for properly parsing the
/// stream.  An error will be raised otherwise.  Usually this is not
/// the problem, because \link colvarvalue \endlink objects are first
/// initialized in the configuration, and the stream operation will be
/// performed only when reading restart files.
///
/// No problem of course with the output streams: \code os << x;
/// \endcode will print a different output according to the value of
/// colvarvalue::value_type, and the width of such output is returned
/// by colvarvalue::output_width()
///
/// \em Note \em on \em performance: at every operation between two
/// \link colvarvalue \endlink objects, their two \link Type \endlink
/// flags will be checked for a match.  In a long array of \link
/// colvarvalue \endlink objects this is time consuming: a few static
/// functions are defined ("xxx_opt" functions) within the \link
/// colvarvalue \endlink class, which only check the matching once for
/// a large array, and execute different loops according to the type.
/// You should do the same for every time consuming loop involving
/// operations on \link colvarvalue \endlink objects if you want
/// e.g. to optimize your colvar bias.

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
    /// 4-dimensional unit vector representing a rotation, implemented as \link colvarmodule::quaternion \endlink
    type_quaternion
  };

  /// Runtime description of value types
  std::string static const type_desc[colvarvalue::type_quaternion+1];

  /// Number of degrees of freedom for each type
  size_t static const      dof_num[  colvarvalue::type_quaternion+1];

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
    : real_value (0.0), value_type (type_scalar)
  {}

  /// Constructor from a type specification
  inline colvarvalue (Type const &vti)
    : value_type (vti)
  {
    reset();
  }

  /// Copy constructor from real base type
  inline colvarvalue (cvm::real const &x)
    : real_value (x), value_type (type_scalar)
  {}

  /// \brief Copy constructor from rvector base type (Note: this sets
  /// automatically a type \link type_vector \endlink , if you want a
  /// \link type_unitvector \endlink you must set it explicitly)
  inline colvarvalue (cvm::rvector const &v)
    : rvector_value (v), value_type (type_vector)
  {}

  /// \brief Copy constructor from rvector base type (additional
  /// argument to make possible to choose a \link type_unitvector
  /// \endlink
  inline colvarvalue (cvm::rvector const &v, Type const &vti)
    : rvector_value (v), value_type (vti)
  {}

  /// \brief Copy constructor from quaternion base type
  inline colvarvalue (cvm::quaternion const &q)
    : quaternion_value (q), value_type (type_quaternion)
  {}

  /// Copy constructor from another \link colvarvalue \endlink
  inline colvarvalue (colvarvalue const &x)
    : value_type (x.value_type)
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
  /// quaternion), transform it to satisfy them; use it when the \link
  /// colvarvalue \endlink is not calculated from \link cvc
  /// \endlink objects, but manipulated by you
  void apply_constraints();


  /// Get the current type
  inline Type type() const
  {
    return value_type;
  }

  /// Set the type explicitly
  inline void type (Type const &vti)
  {
    reset();
    value_type = vti;
    reset();
  }

  /// Set the type after another \link colvarvalue \endlink
  inline void type (colvarvalue const &x)
  {
    reset();
    value_type = x.value_type;
    reset();
  }

  /// Square norm of this colvarvalue
  cvm::real norm2() const;

  /// Norm of this colvarvalue
  inline cvm::real norm() const
  {
    return std::sqrt (this->norm2());
  }

  /// \brief Return the value whose scalar product with this value is
  /// 1
  inline colvarvalue inverse() const;

  /// Square distance between this \link colvarvalue \endlink and another
  cvm::real dist2 (colvarvalue const &x2) const;

  /// Derivative with respect to this \link colvarvalue \endlink of the square distance
  colvarvalue dist2_grad (colvarvalue const &x2) const;

  /// Assignment operator (type of x is checked)
  colvarvalue & operator = (colvarvalue const &x);

  void operator += (colvarvalue const &x);
  void operator -= (colvarvalue const &x);
  void operator *= (cvm::real const &a);
  void operator /= (cvm::real const &a);


  // Casting operator
  inline operator cvm::real() const
  {
    if (value_type != type_scalar) error_rside (type_scalar);
    return real_value;
  }

  // Casting operator
  inline operator cvm::rvector() const
  {
    if ( (value_type != type_vector) && (value_type != type_unitvector))
      error_rside (type_vector);
    return rvector_value;
  }

  // Casting operator
  inline operator cvm::quaternion() const
  {
    if (value_type != type_quaternion) error_rside (type_quaternion);
    return quaternion_value;
  }

  /// Special case when the variable is a real number, and all operations are defined
  inline bool is_scalar()
  {
    return (value_type == type_scalar);
  }


  /// Ensure that the two types are the same within a binary operator
  void static check_types (colvarvalue const &x1, colvarvalue const &x2);

  /// Undefined operation
  void undef_op() const;

  /// Trying to assign this \link colvarvalue \endlink object to
  /// another object set with a different type
  void error_lside (Type const &vt) const;

  /// Trying to assign another \link colvarvalue \endlink object set
  /// with a different type to this object
  void error_rside (Type const &vt) const;

  /// Give the number of characters required to output this
  /// colvarvalue, given the current type assigned and the number of
  /// characters for a real number
  size_t output_width (size_t const &real_width) const;


  // optimized routines for operations with an array; xv and inner as
  // vectors are assumed to have the same number of elements (i.e. the
  // end for inner comes together with xv_end in the loop)

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  void static inner_opt (colvarvalue const                        &x,
                         std::vector<colvarvalue>::iterator       &xv,
                         std::vector<colvarvalue>::iterator const &xv_end,
                         std::vector<cvm::real>::iterator         &inner);

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  void static inner_opt (colvarvalue const                        &x,
                         std::list<colvarvalue>::iterator         &xv,
                         std::list<colvarvalue>::iterator const   &xv_end,
                         std::vector<cvm::real>::iterator         &inner);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial, (3cos^2(w)-1)/2, of one collective variable with an
  /// array
  void static p2leg_opt (colvarvalue const                        &x,
                         std::vector<colvarvalue>::iterator       &xv,
                         std::vector<colvarvalue>::iterator const &xv_end,
                         std::vector<cvm::real>::iterator         &inner);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial of one collective variable with an array
  void static p2leg_opt (colvarvalue const                        &x,
                         std::list<colvarvalue>::iterator         &xv,
                         std::list<colvarvalue>::iterator const   &xv_end,
                         std::vector<cvm::real>::iterator         &inner);

};



inline void colvarvalue::reset()
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value = cvm::real (0.0);
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    rvector_value = cvm::rvector (0.0, 0.0, 0.0);
    break;
  case colvarvalue::type_quaternion:
    quaternion_value = cvm::quaternion (0.0, 0.0, 0.0, 0.0);
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
    break;
  case colvarvalue::type_vector:
    break;
  case colvarvalue::type_unitvector:
    rvector_value /= std::sqrt (rvector_value.norm2());
    break;
  case colvarvalue::type_quaternion:
    quaternion_value /= std::sqrt (quaternion_value.norm2());
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
    return (this->rvector_value).norm2();
  case colvarvalue::type_quaternion:
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
    return colvarvalue (1.0/real_value);
  case colvarvalue::type_vector:
    return colvarvalue (cvm::rvector (1.0/rvector_value.x,
                                      1.0/rvector_value.y,
                                      1.0/rvector_value.z));
  case colvarvalue::type_quaternion:
    return colvarvalue (cvm::quaternion (1.0/quaternion_value.q0,
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
      error_lside (x.value_type);

  this->value_type = x.value_type;

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    this->real_value = x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    this->rvector_value = x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
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
    colvarvalue::check_types (*this, x);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    this->real_value += x.real_value;
    break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    this->rvector_value += x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
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
    colvarvalue::check_types (*this, x);

  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value -= x.real_value; break;
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    rvector_value -= x.rvector_value; break;
  case colvarvalue::type_quaternion:
    quaternion_value -= x.quaternion_value; break;
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
  case colvarvalue::type_unitvector:
    rvector_value *= a;
    break;
  case colvarvalue::type_quaternion:
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
    rvector_value /= a; break;
  case colvarvalue::type_quaternion:
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
    colvarvalue::check_types (x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue (x1.real_value + x2.real_value);
  case colvarvalue::type_vector:
    return colvarvalue (x1.rvector_value + x2.rvector_value);
  case colvarvalue::type_unitvector:
    return colvarvalue (x1.rvector_value + x2.rvector_value,
                        colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
    return colvarvalue (x1.quaternion_value + x2.quaternion_value);
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return colvarvalue (colvarvalue::type_notset);
  };
}

inline colvarvalue operator - (colvarvalue const &x1,
                               colvarvalue const &x2)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types (x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue (x1.real_value - x2.real_value);
  case colvarvalue::type_vector:
    return colvarvalue (x1.rvector_value - x2.rvector_value);
  case colvarvalue::type_unitvector:
    return colvarvalue (x1.rvector_value - x2.rvector_value,
                        colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
    return colvarvalue (x1.quaternion_value - x2.quaternion_value);
  default:
    x1.undef_op();
    return colvarvalue (colvarvalue::type_notset);
  };
}


// binary operations with real numbers

inline colvarvalue operator * (cvm::real const &a,
                               colvarvalue const &x)
{
  switch (x.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue (a * x.real_value);
  case colvarvalue::type_vector:
    return colvarvalue (a * x.rvector_value);
  case colvarvalue::type_unitvector:
    return colvarvalue (a * x.rvector_value,
                        colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
    return colvarvalue (a * x.quaternion_value);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue (colvarvalue::type_notset);
  }
}

inline colvarvalue operator * (colvarvalue const &x,
                               cvm::real const &a)
{
  switch (x.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue (x.real_value * a);
  case colvarvalue::type_vector:
    return colvarvalue (x.rvector_value * a);
  case colvarvalue::type_unitvector:
    return colvarvalue (x.rvector_value * a,
                        colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
    return colvarvalue (x.quaternion_value * a);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue (colvarvalue::type_notset);
  }
}

inline colvarvalue operator / (colvarvalue const &x,
                               cvm::real const &a)
{
  switch (x.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue (x.real_value / a);
  case colvarvalue::type_vector:
    return colvarvalue (x.rvector_value / a);
  case colvarvalue::type_unitvector:
    return colvarvalue (x.rvector_value / a,
                        colvarvalue::type_unitvector);
  case colvarvalue::type_quaternion:
    return colvarvalue (x.quaternion_value / a);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue (colvarvalue::type_notset);
  }
}


// inner product between two colvarvalues

inline cvm::real operator * (colvarvalue const &x1,
                             colvarvalue const &x2)
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types (x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return (x1.real_value * x2.real_value);
  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    return (x1.rvector_value * x2.rvector_value);
  case colvarvalue::type_quaternion:
    // the "*" product is the quaternion product, here the inner
    // member function is used instead
    return (x1.quaternion_value.inner (x2.quaternion_value));
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return 0.0;
  };
}



inline cvm::real colvarvalue::dist2 (colvarvalue const &x2) const
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types (*this, x2);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    return (this->real_value - x2.real_value)*(this->real_value - x2.real_value);
  case colvarvalue::type_vector:
    return (this->rvector_value - x2.rvector_value).norm2();
  case colvarvalue::type_unitvector:
    // angle between (*this) and x2 is the distance
    return std::acos (this->rvector_value * x2.rvector_value) * std::acos (this->rvector_value * x2.rvector_value);
  case colvarvalue::type_quaternion:
    // angle between (*this) and x2 is the distance, the quaternion
    // object has it implemented internally
    return this->quaternion_value.dist2 (x2.quaternion_value);
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return 0.0;
  };
}


inline colvarvalue colvarvalue::dist2_grad (colvarvalue const &x2) const
{
  if (colvarvalue::type_checking())
    colvarvalue::check_types (*this, x2);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    return 2.0 * (this->real_value - x2.real_value);
  case colvarvalue::type_vector:
    return 2.0 * (this->rvector_value - x2.rvector_value);
  case colvarvalue::type_unitvector:
    {
      cvm::rvector const &v1 = this->rvector_value;
      cvm::rvector const &v2 = x2.rvector_value;
      cvm::real const cos_t = v1 * v2;
      cvm::real const sin_t = std::sqrt (1.0 - cos_t*cos_t);
      return colvarvalue ( 2.0 * sin_t *
                           cvm::rvector ( (-1.0) * sin_t * v2.x +
                                          cos_t/sin_t * (v1.x - cos_t*v2.x),
                                          (-1.0) * sin_t * v2.y +
                                          cos_t/sin_t * (v1.y - cos_t*v2.y),
                                          (-1.0) * sin_t * v2.z +
                                          cos_t/sin_t * (v1.z - cos_t*v2.z)
                                          ),
                           colvarvalue::type_unitvector );
    }
  case colvarvalue::type_quaternion:
    return this->quaternion_value.dist2_grad (x2.quaternion_value);
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return colvarvalue (colvarvalue::type_notset);
  };
}


inline void colvarvalue::check_types (colvarvalue const &x1,
                                      colvarvalue const &x2)
{
  if (x1.value_type != x2.value_type) {
    cvm::log ("x1 type = "+cvm::to_str (x1.value_type)+
              ", x2 type = "+cvm::to_str (x2.value_type)+"\n");
    cvm::fatal_error ("Performing an operation between two colvar "
                      "values with different types, \""+
                      colvarvalue::type_desc[x1.value_type]+
                      "\" and \""+
                      colvarvalue::type_desc[x2.value_type]+
                      "\".\n");
  }
}




std::ostream & operator << (std::ostream &os, colvarvalue const &x);
std::ostream & operator << (std::ostream &os, std::vector<colvarvalue> const &v);

std::istream & operator >> (std::istream &is, colvarvalue &x);





#endif


// Emacs
// Local Variables:
// mode: C++
// End:
