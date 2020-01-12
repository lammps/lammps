// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARVALUE_H
#define COLVARVALUE_H

#include "colvarmodule.h"
#include "colvartypes.h"


/// \brief Value of a collective variable: this is a metatype which
/// can be set at runtime.  By default it is set to be a scalar
/// number, and can be treated as such in all operations (this is
/// done by most \link colvar::cvc \endlink implementations).
///
/// \link colvarvalue \endlink allows \link colvar \endlink to be
/// treat different data types.  By default, a \link colvarvalue
/// \endlink variable is a scalar number.  To use it as
/// another type, declare and initialize it as
/// `colvarvalue x(colvarvalue::type_xxx)`, use `x.type (colvarvalue::type_xxx)`
///  at a later stage, or if unset,
///  assign the type with `x = y;`, provided y is correctly set.
///
/// All operators (either unary or binary) on a \link
/// colvarvalue \endlink object performs one or more checks on the
/// \link Type \endlink, except when reading from a stream, when there is no way to
/// detect the \link Type \endlink.  To use  `is >> x;` x \b MUST
/// already have a type correcly set up for properly parsing the
/// stream. No problem of course with the output streams: `os << x;`
///
/// \em Note \em on \em performance: to avoid type checks in a long array of \link
/// colvarvalue \endlink objects, use one of the existing "_opt" functions or implement a new one


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
    type_3vector,
    /// 3-dimensional unit vector, implemented as \link colvarmodule::rvector \endlink
    type_unit3vector,
    /// 3-dimensional vector that is a derivative of a unitvector
    type_unit3vectorderiv,
    /// 4-dimensional unit vector representing a rotation, implemented as \link colvarmodule::quaternion \endlink
    type_quaternion,
    /// 4-dimensional vector that is a derivative of a quaternion
    type_quaternionderiv,
    /// vector (arbitrary dimension)
    type_vector,
    /// Needed to iterate through enum
    type_all
  };

  /// Current type of this colvarvalue object
  Type value_type;

  /// \brief Real data member
  cvm::real real_value;

  /// \brief 3-dimensional vector data member
  cvm::rvector rvector_value;

  /// \brief Quaternion data member
  cvm::quaternion quaternion_value;

  /// \brief Generic vector data member
  cvm::vector1d<cvm::real> vector1d_value;

  /// \brief If \link vector1d_value \endlink is a concatenation of colvarvalues,
  /// keep track of the individual types
  std::vector<Type> elem_types;

  /// \brief If \link vector1d_value \endlink is a concatenation of colvarvalues,
  /// these mark the initial components of each colvarvalue
  std::vector<int> elem_indices;

  /// \brief If \link vector1d_value \endlink is a concatenation of colvarvalues,
  /// these mark how many components for each colvarvalue
  std::vector<int> elem_sizes;

  /// \brief Whether or not the type check is enforced
  static inline bool type_checking()
  {
    return true;
  }

  /// Runtime description of value types
  static std::string const type_desc(Type t);

  /// User keywords for specifying value types in the configuration
  static std::string const type_keyword(Type t);

  /// Number of degrees of freedom for each supported type
  static size_t num_df(Type t);

  /// Number of dimensions for each supported type (used to allocate vector1d_value)
  static size_t num_dimensions(Type t);

  /// Number of dimensions of this variable
  size_t size() const;

  /// \brief Default constructor: this class defaults to a scalar
  /// number and always behaves like it unless you change its type
  inline colvarvalue()
    : value_type(type_scalar), real_value(0.0)
  {}

  /// Constructor from a type specification
  inline colvarvalue(Type const &vti)
    : value_type(vti)
  {
    reset();
  }

  /// Copy constructor from real base type
  inline colvarvalue(cvm::real const &x)
    : value_type(type_scalar), real_value(x)
  {}

  /// \brief Copy constructor from rvector base type (Note: this sets
  /// by default a type \link type_3vector \endlink , if you want a
  /// \link type_unit3vector \endlink you must set it explicitly)
  inline colvarvalue(cvm::rvector const &v, Type vti = type_3vector)
    : value_type(vti), rvector_value(v)
  {}

  /// \brief Copy constructor from quaternion base type
  inline colvarvalue(cvm::quaternion const &q, Type vti = type_quaternion)
    : value_type(vti), quaternion_value(q)
  {}

  /// Copy constructor from vector1d base type
  colvarvalue(cvm::vector1d<cvm::real> const &v, Type vti = type_vector);

  /// Copy constructor from another \link colvarvalue \endlink
  colvarvalue(colvarvalue const &x);


  /// Set to the null value for the data type currently defined
  void reset();

  /// \brief If the variable has constraints (e.g. unitvector or
  /// quaternion), transform it to satisfy them; this function needs
  /// to be called only when the \link colvarvalue \endlink
  /// is calculated outside of \link colvar::cvc \endlink objects
  void apply_constraints();

  /// Get the current type
  inline Type type() const
  {
    return value_type;
  }

  /// Set the type explicitly
  void type(Type const &vti);

  /// Set the type after another \link colvarvalue \endlink
  void type(colvarvalue const &x);

  /// Make the type a derivative of the original type
  /// (so that its constraints do not apply)
  void is_derivative();

  /// Square norm of this colvarvalue
  cvm::real norm2() const;

  /// Norm of this colvarvalue
  inline cvm::real norm() const
  {
    return cvm::sqrt(this->norm2());
  }

  /// Sum of the components of this colvarvalue (if more than one dimension)
  cvm::real sum() const;

  /// Return a colvarvalue object of the same type and all components set to 1
  colvarvalue ones() const;

  /// Square distance between this \link colvarvalue \endlink and another
  cvm::real dist2(colvarvalue const &x2) const;

  /// Derivative with respect to this \link colvarvalue \endlink of the square distance
  colvarvalue dist2_grad(colvarvalue const &x2) const;

  /// Return the midpoint between x1 and x2, optionally weighted by lambda
  /// (which must be between 0.0 and 1.0)
  static colvarvalue const interpolate(colvarvalue const &x1,
                                       colvarvalue const &x2,
                                       cvm::real const lambda = 0.5);

  /// Assignment operator (type of x is checked)
  colvarvalue & operator = (colvarvalue const &x);

  void operator += (colvarvalue const &x);
  void operator -= (colvarvalue const &x);
  void operator *= (cvm::real const &a);
  void operator /= (cvm::real const &a);

  // Binary operators (return values)
  friend colvarvalue operator + (colvarvalue const &x1, colvarvalue const &x2);
  friend colvarvalue operator - (colvarvalue const &x1, colvarvalue const &x2);
  friend colvarvalue operator * (colvarvalue const &x, cvm::real const &a);
  friend colvarvalue operator * (cvm::real const &a, colvarvalue const &x);
  friend colvarvalue operator / (colvarvalue const &x, cvm::real const &a);

  /// Inner product
  friend cvm::real operator * (colvarvalue const &x1, colvarvalue const &x2);

  // Cast to scalar
  inline operator cvm::real() const
  {
    if (value_type != type_scalar) {
      cvm::error("Error: trying to use a variable of type \""+
                 type_desc(value_type)+"\" as one of type \""+
                 type_desc(type_scalar)+"\".\n");
    }
    return real_value;
  }

  // Cast to 3-vector
  inline operator cvm::rvector() const
  {
    if ((value_type != type_3vector) &&
        (value_type != type_unit3vector) &&
        (value_type != type_unit3vectorderiv)) {
      cvm::error("Error: trying to use a variable of type \""+
                 type_desc(value_type)+"\" as one of type \""+
                 type_desc(type_3vector)+"\".\n");
    }
    return rvector_value;
  }

  // Cast to quaternion
  inline operator cvm::quaternion() const
  {
    if ((value_type != type_quaternion) &&
        (value_type != type_quaternionderiv)) {
      cvm::error("Error: trying to use a variable of type \""+
                 type_desc(value_type)+"\" as one of type \""+
                 type_desc(type_quaternion)+"\".\n");
    }
    return quaternion_value;
  }

  // Create a n-dimensional vector from one of the basic types, or return the existing vector
  cvm::vector1d<cvm::real> const as_vector() const;


  /// Whether this variable is a real number
  inline bool is_scalar() const
  {
    return (value_type == type_scalar);
  }


  /// Add an element to the vector (requires that type_vector is already set).
  /// This is only needed to use this object as a vector of "complex" colvar values.
  /// To use it instead as a plain n-dimensional vector, access vector1d_value directly.
  void add_elem(colvarvalue const &x);

  /// Get a single colvarvalue out of elements of the vector
  colvarvalue const get_elem(int const i_begin, int const i_end, Type const vt) const;

  /// Get a single colvarvalue out of elements of the vector
  colvarvalue const get_elem(int const icv) const;

  /// Set elements of the vector from a single colvarvalue (uses the rank of x
  /// to compute the length)
  void set_elem(int const icv, colvarvalue const &x);

  /// Set elements of the vector from a single colvarvalue
  void set_elem(int const i_begin, int const i_end, colvarvalue const &x);

  /// Make each element a random number in N(0,1)
  void set_random();

  /// Make each element equal to the given argument
  void set_ones(cvm::real assigned_value = 1.0);

  /// Get a scalar number out of an element of the vector
  cvm::real operator [] (int const i) const;

  /// Use an element of the vector as a scalar number
  cvm::real & operator [] (int const i);

  /// Ensure that the two types are the same within a binary operator
  static int check_types(colvarvalue const &x1, colvarvalue const &x2);

  /// Ensure that the two types are the same within an assignment, or that the left side is type_notset
  static int check_types_assign(Type const &vt1, Type const &vt2);

  /// Undefined operation
  void undef_op() const;


  /// \brief Formatted output operator
  friend std::ostream & operator << (std::ostream &os, colvarvalue const &q);

  /// \brief Formatted input operator
  friend std::istream & operator >> (std::istream &is, colvarvalue &q);

  /// Give the number of characters required to output this
  /// colvarvalue, given the current type assigned and the number of
  /// characters for a real number
  size_t output_width(size_t const &real_width) const;

  /// Formats value as a script-friendly string (space separated list)
  std::string to_simple_string() const;

  /// Parses value from a script-friendly string (space separated list)
  int from_simple_string(std::string const &s);


  // optimized routines for operations on arrays of colvar values;
  // xv and result are assumed to have the same number of elements

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  static void inner_opt(colvarvalue const                        &x,
                        std::vector<colvarvalue>::iterator       &xv,
                        std::vector<colvarvalue>::iterator const &xv_end,
                        std::vector<cvm::real>::iterator         &result);

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  static void inner_opt(colvarvalue const                        &x,
                        std::list<colvarvalue>::iterator         &xv,
                        std::list<colvarvalue>::iterator const   &xv_end,
                        std::vector<cvm::real>::iterator         &result);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial, (3cos^2(w)-1)/2, of one collective variable with an
  /// array
  static void p2leg_opt(colvarvalue const                        &x,
                        std::vector<colvarvalue>::iterator       &xv,
                        std::vector<colvarvalue>::iterator const &xv_end,
                        std::vector<cvm::real>::iterator         &result);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial of one collective variable with an array
  static void p2leg_opt(colvarvalue const                        &x,
                        std::list<colvarvalue>::iterator         &xv,
                        std::list<colvarvalue>::iterator const   &xv_end,
                        std::vector<cvm::real>::iterator         &result);

};


inline size_t colvarvalue::size() const
{
  switch (value_type) {
  case colvarvalue::type_notset:
  default:
    return 0; break;
  case colvarvalue::type_scalar:
    return 1; break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return 3; break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return 4; break;
  case colvarvalue::type_vector:
    return vector1d_value.size(); break;
  }
}


inline cvm::real colvarvalue::operator [] (int const i) const
{
  switch (value_type) {
  case colvarvalue::type_notset:
  default:
    cvm::error("Error: trying to access a colvar value "
               "that is not initialized.\n", BUG_ERROR);
    return 0.0; break;
  case colvarvalue::type_scalar:
    return real_value; break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return rvector_value[i]; break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return quaternion_value[i]; break;
  case colvarvalue::type_vector:
    return vector1d_value[i]; break;
  }
}


inline cvm::real & colvarvalue::operator [] (int const i)
{
  switch (value_type) {
  case colvarvalue::type_notset:
  default:
    cvm::error("Error: trying to access a colvar value "
               "that is not initialized.\n", BUG_ERROR);
    return real_value; break;
  case colvarvalue::type_scalar:
    return real_value; break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return rvector_value[i]; break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return quaternion_value[i]; break;
  case colvarvalue::type_vector:
    return vector1d_value[i]; break;
  }
}


inline int colvarvalue::check_types(colvarvalue const &x1,
                                    colvarvalue const &x2)
{
  if (!colvarvalue::type_checking()) {
    return COLVARS_OK;
  }

  if (x1.type() != x2.type()) {
    if (((x1.type() == type_unit3vector) &&
         (x2.type() == type_unit3vectorderiv)) ||
        ((x2.type() == type_unit3vector) &&
         (x1.type() == type_unit3vectorderiv)) ||
        ((x1.type() == type_quaternion) &&
         (x2.type() == type_quaternionderiv)) ||
        ((x2.type() == type_quaternion) &&
         (x1.type() == type_quaternionderiv))) {
      return COLVARS_OK;
    } else {
      cvm::error("Trying to perform an operation between two colvar "
                 "values with different types, \""+
                 colvarvalue::type_desc(x1.type())+
                 "\" and \""+
                 colvarvalue::type_desc(x2.type())+
                 "\".\n");
      return COLVARS_ERROR;
    }
  }

  if (x1.type() == type_vector) {
    if (x1.vector1d_value.size() != x2.vector1d_value.size()) {
      cvm::error("Trying to perform an operation between two vector colvar "
                 "values with different sizes, "+
                 cvm::to_str(x1.vector1d_value.size())+
                 " and "+
                 cvm::to_str(x2.vector1d_value.size())+
                 ".\n");
      return COLVARS_ERROR;
    }
  }
  return COLVARS_OK;
}


inline int colvarvalue::check_types_assign(colvarvalue::Type const &vt1,
                                           colvarvalue::Type const &vt2)
{
  if (!colvarvalue::type_checking()) {
    return COLVARS_OK;
  }

  if (vt1 != type_notset) {
    if (((vt1 == type_unit3vector) &&
         (vt2 == type_unit3vectorderiv)) ||
        ((vt2 == type_unit3vector) &&
         (vt1 == type_unit3vectorderiv)) ||
        ((vt1 == type_quaternion) &&
         (vt2 == type_quaternionderiv)) ||
        ((vt2 == type_quaternion) &&
         (vt1 == type_quaternionderiv))) {
      return COLVARS_OK;
    } else {
      if (vt1 != vt2) {
        cvm::error("Trying to assign a colvar value with type \""+
                   type_desc(vt2)+"\" to one with type \""+
                   type_desc(vt1)+"\".\n");
        return COLVARS_ERROR;
      }
    }
  }
  return COLVARS_OK;
}


inline colvarvalue & colvarvalue::operator = (colvarvalue const &x)
{
  check_types_assign(this->type(), x.type());
  value_type = x.type();

  switch (this->type()) {
  case colvarvalue::type_scalar:
    this->real_value = x.real_value;
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    this->rvector_value = x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    this->quaternion_value = x.quaternion_value;
    break;
  case colvarvalue::type_vector:
    vector1d_value = x.vector1d_value;
    elem_types = x.elem_types;
    elem_indices = x.elem_indices;
    elem_sizes = x.elem_sizes;
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
  colvarvalue::check_types(*this, x);

  switch (this->type()) {
  case colvarvalue::type_scalar:
    this->real_value += x.real_value;
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    this->rvector_value += x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    this->quaternion_value += x.quaternion_value;
    break;
  case colvarvalue::type_vector:
    this->vector1d_value += x.vector1d_value;
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}


inline void colvarvalue::operator -= (colvarvalue const &x)
{
  colvarvalue::check_types(*this, x);

  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value -= x.real_value;
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    rvector_value -= x.rvector_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value -= x.quaternion_value;
    break;
  case colvarvalue::type_vector:
    this->vector1d_value -= x.vector1d_value;
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vectorderiv:
    rvector_value *= a;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value *= a;
    break;
  case colvarvalue::type_vector:
    this->vector1d_value *= a;
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    rvector_value /= a; break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value /= a; break;
  case colvarvalue::type_vector:
    this->vector1d_value /= a;
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}


inline cvm::vector1d<cvm::real> const colvarvalue::as_vector() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    {
      cvm::vector1d<cvm::real> v(1);
      v[0] = real_value;
      return v;
    }
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return rvector_value.as_vector();
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return quaternion_value.as_vector();
  case colvarvalue::type_vector:
    return vector1d_value;
  case colvarvalue::type_notset:
  default:
    return cvm::vector1d<cvm::real>(0);
  }
}


inline cvm::real colvarvalue::norm2() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    return (this->real_value)*(this->real_value);
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return (this->rvector_value).norm2();
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return (this->quaternion_value).norm2();
  case colvarvalue::type_vector:
    if (elem_types.size() > 0) {
      // if we have information about non-scalar types, use it
      cvm::real result = 0.0;
      size_t i;
      for (i = 0; i < elem_types.size(); i++) {
        result += (this->get_elem(i)).norm2();
      }
      return result;
    } else {
      return vector1d_value.norm2();
    }
    break;
  case colvarvalue::type_notset:
  default:
    return 0.0;
  }
}


inline cvm::real colvarvalue::sum() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    return (this->real_value);
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return (this->rvector_value).x + (this->rvector_value).y +
      (this->rvector_value).z;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return (this->quaternion_value).q0 + (this->quaternion_value).q1 +
      (this->quaternion_value).q2 + (this->quaternion_value).q3;
  case colvarvalue::type_vector:
    return (this->vector1d_value).sum();
  case colvarvalue::type_notset:
  default:
    return 0.0;
  }
}


inline cvm::real colvarvalue::dist2(colvarvalue const &x2) const
{
  colvarvalue::check_types(*this, x2);

  switch (this->type()) {
  case colvarvalue::type_scalar:
    return (this->real_value - x2.real_value)*(this->real_value - x2.real_value);
  case colvarvalue::type_3vector:
    return (this->rvector_value - x2.rvector_value).norm2();
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    // angle between (*this) and x2 is the distance
    return cvm::acos(this->rvector_value * x2.rvector_value) * cvm::acos(this->rvector_value * x2.rvector_value);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    // angle between (*this) and x2 is the distance, the quaternion
    // object has it implemented internally
    return this->quaternion_value.dist2(x2.quaternion_value);
  case colvarvalue::type_vector:
    return (this->vector1d_value - x2.vector1d_value).norm2();
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return 0.0;
  };
}


#endif
