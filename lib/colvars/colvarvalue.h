/// -*- c++ -*-

#ifndef COLVARVALUE_H
#define COLVARVALUE_H

#include "colvarmodule.h"
#include "colvartypes.h"


/// \brief Value of a collective variable: this is a metatype which
/// can be set at runtime.  By default it is set to be a scalar
/// number, and can be treated as such in all operations (this is
/// done by most \link cvc \endlink implementations).
///
/// \link colvarvalue \endlink allows \link colvar \endlink to be
/// treat different data types.  By default, a \link colvarvalue
/// \endlink variable is a scalar number.  To use it as
/// another type, declare and initialize it as
/// \code colvarvalue x(colvarvalue::type_xxx)\endcode, use \link x.type
/// (colvarvalue::type_xxx) \endlink at a later stage, or if unset,
///  assign the type with \code x = y; \endcode, provided y is correctly set.
///
/// All operators (either unary or binary) on a \link
/// colvarvalue \endlink object performs one or more checks on the
/// \link Type \endlink, except when reading from a stream, when there is no way to
/// detect the \link Type \endlink.  To use  \code is >> x; \endcode x \b MUST
/// already have a type correcly set up for properly parsing the
/// stream. No problem of course with the output streams: \code os << x; \endcode
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
  /// automatically a type \link type_3vector \endlink , if you want a
  /// \link type_unit3vector \endlink you must set it explicitly)
  inline colvarvalue(cvm::rvector const &v)
    : value_type(type_3vector), rvector_value(v)
  {}

  /// \brief Copy constructor from rvector base type (additional
  /// argument to make possible to choose a \link type_unit3vector
  /// \endlink
  inline colvarvalue(cvm::rvector const &v, Type const &vti)
    : value_type(vti), rvector_value(v)
  {}

  /// \brief Copy constructor from quaternion base type
  inline colvarvalue(cvm::quaternion const &q)
    : value_type(type_quaternion), quaternion_value(q)
  {}

  /// Copy constructor from another \link colvarvalue \endlink
  colvarvalue(colvarvalue const &x);

  /// Copy constructor from vector1d base type
  colvarvalue(cvm::vector1d<cvm::real> const &v, Type const &vti);


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
    // reset the value based on the previous type
    reset();
    if ((value_type == type_vector) && (vti != type_vector)) {
      vector1d_value.resize(0);
    }
    value_type = vti;
  }

  /// Set the type after another \link colvarvalue \endlink
  inline void type(colvarvalue const &x)
  {
    // reset the value held from based on the previous type
    reset();
    if (x.type() == type_vector) {
      vector1d_value.resize(x.vector1d_value.size());
    } else {
      if (value_type == type_vector) {
        vector1d_value.resize(0);
      }
    }
    value_type = x.type();
  }

  /// Make the type a derivative of the original type
  /// (so that its constraints do not apply)
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

  /// Set elements of the vector from a single colvarvalue
  void set_elem(int const i_begin, int const i_end, colvarvalue const &x);

  /// Get a single colvarvalue out of elements of the vector
  colvarvalue const get_elem(int const icv) const;

  /// Set elements of the vector from a single colvarvalue
  void set_elem(int const icv, colvarvalue const &x);

  /// Get a scalar number out of an element of the vector
  inline cvm::real operator [] (int const i) const
  {
    if (vector1d_value.size() > 0) {
      return vector1d_value[i];
    } else {
      cvm::error("Error: trying to use as a vector a variable that is not initialized as such.\n");
      return 0.0;
    }
  }

  /// Use an element of the vector as a scalar number
  inline cvm::real & operator [] (int const i)
  {
    if (vector1d_value.size() > 0) {
      return vector1d_value[i];
    } else {
      cvm::error("Error: trying to use as a vector a variable that is not initialized as such.\n");
      real_value = 0.0;
      return real_value;
    }
  }


  /// Ensure that the two types are the same within a binary operator
  int static check_types(colvarvalue const &x1, colvarvalue const &x2);

  /// Ensure that the two types are the same within an assignment, or that the left side is type_notset
  int static check_types_assign(Type const &vt1, Type const &vt2);

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
  void static inner_opt(colvarvalue const                        &x,
                        std::vector<colvarvalue>::iterator       &xv,
                        std::vector<colvarvalue>::iterator const &xv_end,
                        std::vector<cvm::real>::iterator         &result);

  /// \brief Optimized routine for the inner product of one collective
  /// variable with an array
  void static inner_opt(colvarvalue const                        &x,
                        std::list<colvarvalue>::iterator         &xv,
                        std::list<colvarvalue>::iterator const   &xv_end,
                        std::vector<cvm::real>::iterator         &result);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial, (3cos^2(w)-1)/2, of one collective variable with an
  /// array
  void static p2leg_opt(colvarvalue const                        &x,
                        std::vector<colvarvalue>::iterator       &xv,
                        std::vector<colvarvalue>::iterator const &xv_end,
                        std::vector<cvm::real>::iterator         &result);

  /// \brief Optimized routine for the second order Legendre
  /// polynomial of one collective variable with an array
  void static p2leg_opt(colvarvalue const                        &x,
                        std::list<colvarvalue>::iterator         &xv,
                        std::list<colvarvalue>::iterator const   &xv_end,
                        std::vector<cvm::real>::iterator         &result);

};



inline std::string const colvarvalue::type_desc(Type t)
{
  switch (t) {
  case colvarvalue::type_scalar:
    return "scalar number"; break;
  case colvarvalue::type_3vector:
    return "3-dimensional vector"; break;
  case colvarvalue::type_unit3vector:
    return "3-dimensional unit vector"; break;
  case colvarvalue::type_unit3vectorderiv:
    return "derivative of a 3-dimensional unit vector"; break;
  case colvarvalue::type_quaternion:
    return "4-dimensional unit quaternion"; break;
  case colvarvalue::type_quaternionderiv:
    return "4-dimensional tangent vector"; break;
  case colvarvalue::type_vector:
    return "n-dimensional vector"; break;
  case colvarvalue::type_notset:
    // fallthrough
  default:
    return "not set"; break;
  }
}


inline std::string const colvarvalue::type_keyword(Type t)
{
  switch (t) {
  case colvarvalue::type_notset:
  default:
    return "not_set"; break;
  case colvarvalue::type_scalar:
    return "scalar"; break;
  case colvarvalue::type_3vector:
    return "vector3"; break;
  case colvarvalue::type_unit3vector:
    return "unit_vector3"; break;
  case colvarvalue::type_unit3vectorderiv:
    return ""; break;
  case colvarvalue::type_quaternion:
    return "unit_quaternion"; break;
  case colvarvalue::type_quaternionderiv:
    return ""; break;
  case colvarvalue::type_vector:
    return "vector"; break;
  }
}


inline size_t colvarvalue::num_df(Type t)
{
  switch (t) {
  case colvarvalue::type_notset:
  default:
    return 0; break;
  case colvarvalue::type_scalar:
    return 1; break;
  case colvarvalue::type_3vector:
    return 3; break;
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return 2; break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return 3; break;
  case colvarvalue::type_vector:
    // the size of a vector is unknown without its object
    return 0; break;
  }
}


inline size_t colvarvalue::num_dimensions(Type t)
{
  switch (t) {
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
    // the size of a vector is unknown without its object
    return 0; break;
  }
}


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


inline colvarvalue::colvarvalue(colvarvalue const &x)
  : value_type(x.type())
{
  switch (x.type()) {
  case type_scalar:
    real_value = x.real_value;
    break;
  case type_3vector:
  case type_unit3vector:
  case type_unit3vectorderiv:
    rvector_value = x.rvector_value;
    break;
  case type_quaternion:
  case type_quaternionderiv:
    quaternion_value = x.quaternion_value;
    break;
  case type_vector:
    vector1d_value = x.vector1d_value;
    elem_types = x.elem_types;
    elem_indices = x.elem_indices;
    elem_sizes = x.elem_sizes;
  case type_notset:
  default:
    break;
  }
}

inline colvarvalue::colvarvalue(cvm::vector1d<cvm::real> const &v, Type const &vti)
{
  if ((vti != type_vector) && (v.size() != num_dimensions(vti))) {
    cvm::error("Error: trying to initialize a variable of type \""+type_desc(vti)+
               "\" using a vector of size "+cvm::to_str(v.size())+
               ".\n");
    value_type = type_notset;
  } else {
    value_type = vti;
    switch (vti) {
    case type_scalar:
      real_value = v[0];
      break;
    case type_3vector:
    case type_unit3vector:
    case type_unit3vectorderiv:
      rvector_value = cvm::rvector(v);
      break;
    case type_quaternion:
    case type_quaternionderiv:
      quaternion_value = cvm::quaternion(v);
      break;
    case type_vector:
      vector1d_value = v;
      break;
    case type_notset:
    default:
      break;
    }
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
    if (vt1 != vt2) {
      cvm::error("Trying to assign a colvar value with type \""+
                 type_desc(vt2)+"\" to one with type \""+
                 type_desc(vt1)+"\".\n");
      return COLVARS_ERROR;
    }
  }
  return COLVARS_OK;
}


inline void colvarvalue::undef_op() const
{
  cvm::error("Error: Undefined operation on a colvar of type \""+
             type_desc(this->type())+"\".\n");
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


inline void colvarvalue::reset()
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    real_value = 0.0;
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    rvector_value.reset();
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    quaternion_value.reset();
    break;
  case colvarvalue::type_vector:
    vector1d_value.reset();
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vectorderiv:
  case colvarvalue::type_quaternionderiv:
    break;
  case colvarvalue::type_unit3vector:
    rvector_value /= std::sqrt(rvector_value.norm2());
    break;
  case colvarvalue::type_quaternion:
    quaternion_value /= std::sqrt(quaternion_value.norm2());
    break;
  case colvarvalue::type_vector:
    if (elem_types.size() > 0) {
      // if we have information about non-scalar types, use it
      size_t i;
      for (i = 0; i < elem_types.size(); i++) {
        if (elem_sizes[i] == 1) continue; // TODO this can be optimized further
        colvarvalue cvtmp(vector1d_value.slice(elem_indices[i],
                                               elem_indices[i] + elem_sizes[i]), elem_types[i]);
        cvtmp.apply_constraints();
        set_elem(i, cvtmp);
      }
    }
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vectorderiv:
  case colvarvalue::type_quaternionderiv:
    break;
  case colvarvalue::type_unit3vector:
    type(colvarvalue::type_unit3vectorderiv);
    break;
  case colvarvalue::type_quaternion:
    type(colvarvalue::type_quaternionderiv);
    break;
  case colvarvalue::type_vector:
    // TODO
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
    return std::acos(this->rvector_value * x2.rvector_value) * std::acos(this->rvector_value * x2.rvector_value);
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
