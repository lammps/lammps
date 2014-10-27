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
    type_vector
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
  /// automatically a type \link type_3vector \endlink , if you want a
  /// \link type_unit3vector \endlink you must set it explicitly)
  inline colvarvalue(cvm::rvector const &v)
    : rvector_value(v), value_type(type_3vector)
  {}

  /// \brief Copy constructor from rvector base type (additional
  /// argument to make possible to choose a \link type_unit3vector
  /// \endlink
  inline colvarvalue(cvm::rvector const &v, Type const &vti)
    : rvector_value(v), value_type(vti)
  {}

  /// \brief Copy constructor from quaternion base type
  inline colvarvalue(cvm::quaternion const &q)
    : quaternion_value(q), value_type(type_quaternion)
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

  // Cast to n-dimensional vector
  cvm::vector1d<cvm::real> const as_vector() const;


  /// Whether this variable is a real number
  inline bool is_scalar() const
  {
    return (value_type == type_scalar);
  }

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
  void static check_types(colvarvalue const &x1, colvarvalue const &x2);

  /// Ensure that the two types are the same within an assignment, or that the left side is type_notset
  void static check_types_assign(Type const &vt1, Type const &vt2);

  /// Undefined operation
  void undef_op() const;

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



inline std::string const colvarvalue::type_desc(Type t)
{
  switch (t) {
  case colvarvalue::type_notset:
  default:
    return "not set"; break;
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
  : value_type(x.value_type)
{
  reset();

  switch (x.value_type) {
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
  reset();
  if (v.size() != num_dimensions(value_type)) {
    cvm::error("Error: trying to initialize a variable of type \""+type_desc(vti)+
               "\" using a vector of size \""+cvm::to_str(v.size())+
               "\".\n");
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


inline colvarvalue const colvarvalue::get_elem(int const i_begin, int const i_end, Type const vt) const
{
  if (vector1d_value.size() > 0) {
    cvm::vector1d<cvm::real> const v(vector1d_value.slice(i_begin, i_end));
    return colvarvalue(v, vt);
  } else {
    cvm::error("Error: trying to get an element from a variable that is not a vector.\n");
    return colvarvalue(type_notset);
  }
}

inline void colvarvalue::set_elem(int const i_begin, int const i_end, colvarvalue const &x)
{
  if (vector1d_value.size() > 0) {
    vector1d_value.sliceassign(i_begin, i_end, x.as_vector());
  } else {
    cvm::error("Error: trying to set an element for a variable that is not a vector.\n");
  }
}


inline colvarvalue const colvarvalue::get_elem(int const icv) const
{
  if (elem_types.size() > 0) {
    return get_elem(elem_indices[icv], elem_indices[icv] + elem_sizes[icv],
                    elem_types[icv]);
  } else {
    cvm::error("Error: trying to get a colvarvalue element from a vector colvarvalue that was initialized as a plain array.\n");
    return colvarvalue(type_notset);
  }
}

inline void colvarvalue::set_elem(int const icv, colvarvalue const &x)
{
  if (elem_types.size() > 0) {
    check_types_assign(elem_types[icv], x.value_type);
    set_elem(elem_indices[icv], elem_indices[icv] + elem_sizes[icv], x);
  } else {
    cvm::error("Error: trying to set a colvarvalue element for a colvarvalue that was initialized as a plain array.\n");
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
      for (size_t i = 0; i < elem_types.size(); i++) {
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
      for (size_t i = 0; i < elem_types.size(); i++) {
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


inline colvarvalue colvarvalue::inverse() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(1.0/real_value);
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return colvarvalue(cvm::rvector(1.0/rvector_value.x,
                                    1.0/rvector_value.y,
                                    1.0/rvector_value.z));
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(cvm::quaternion(1.0/quaternion_value.q0,
                                       1.0/quaternion_value.q1,
                                       1.0/quaternion_value.q2,
                                       1.0/quaternion_value.q3));
    break;
  case colvarvalue::type_vector:
    {
      cvm::vector1d<cvm::real> result(vector1d_value);
      if (elem_types.size() > 0) {
        // if we have information about non-scalar types, use it
        for (size_t i = 0; i < elem_types.size(); i++) {
          result.sliceassign(elem_indices[i], elem_indices[i]+elem_sizes[i],
                             cvm::vector1d<cvm::real>((this->get_elem(i)).inverse()));
        }
      } else {
        for (size_t i = 0; i < result.size(); i++) {
          if (result[i] != 0.0) {
            result = 1.0/result[i];
          }
        }
      }
      return colvarvalue(result, type_vector);
    }
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
    break;
  }
  return colvarvalue();
}


inline colvarvalue & colvarvalue::operator = (colvarvalue const &x)
{
  check_types_assign(this->value_type, x.value_type);
  value_type = x.value_type;

  switch (this->value_type) {
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

  switch (this->value_type) {
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
    // TODO move this to vector1d
    for (size_t i = 0; i < vector1d_value.size(); i++) {
      vector1d_value[i] += x.vector1d_value[i];
    }
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
    // TODO move this to vector1d
    for (size_t i = 0; i < vector1d_value.size(); i++) {
      vector1d_value[i] -= x.vector1d_value[i];
    }
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
    // TODO move this to vector1d
    for (size_t i = 0; i < vector1d_value.size(); i++) {
      vector1d_value[i] *= a;
    }
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
    // TODO move this to vector1d
    for (size_t i = 0; i < vector1d_value.size(); i++) {
      vector1d_value[i] /= a;
    }
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
  }
}


// binary operations between two colvarvalues

inline colvarvalue operator + (colvarvalue const &x1,
                               colvarvalue const &x2)
{
  colvarvalue::check_types(x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(x1.real_value + x2.real_value);
  case colvarvalue::type_3vector:
    return colvarvalue(x1.rvector_value + x2.rvector_value);
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return colvarvalue(x1.rvector_value + x2.rvector_value,
                       colvarvalue::type_unit3vector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x1.quaternion_value + x2.quaternion_value);
  case colvarvalue::type_vector:
    // TODO
    break;
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}

inline colvarvalue operator - (colvarvalue const &x1,
                               colvarvalue const &x2)
{
  colvarvalue::check_types(x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return colvarvalue(x1.real_value - x2.real_value);
  case colvarvalue::type_3vector:
    return colvarvalue(x1.rvector_value - x2.rvector_value);
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return colvarvalue(x1.rvector_value - x2.rvector_value,
                       colvarvalue::type_unit3vector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x1.quaternion_value - x2.quaternion_value);
  case colvarvalue::type_vector:
    // TODO
    break;
  case colvarvalue::type_notset:
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
  case colvarvalue::type_3vector:
    return colvarvalue(a * x.rvector_value);
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return colvarvalue(a * x.rvector_value,
                       colvarvalue::type_unit3vector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(a * x.quaternion_value);
  case colvarvalue::type_vector:
    // TODO
    break;
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
  case colvarvalue::type_3vector:
    return colvarvalue(x.rvector_value * a);
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return colvarvalue(x.rvector_value * a,
                       colvarvalue::type_unit3vector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x.quaternion_value * a);
  case colvarvalue::type_vector:
    // TODO
    break;
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
  case colvarvalue::type_3vector:
    return colvarvalue(x.rvector_value / a);
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return colvarvalue(x.rvector_value / a,
                       colvarvalue::type_unit3vector);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return colvarvalue(x.quaternion_value / a);
  case colvarvalue::type_vector:
    // TODO
    break;
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
  colvarvalue::check_types(x1, x2);

  switch (x1.value_type) {
  case colvarvalue::type_scalar:
    return (x1.real_value * x2.real_value);
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return (x1.rvector_value * x2.rvector_value);
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    // the "*" product is the quaternion product, here the inner
    // member function is used instead
    return (x1.quaternion_value.inner(x2.quaternion_value));
  case colvarvalue::type_vector:
    // TODO
    break;
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return 0.0;
  };
}



inline cvm::real colvarvalue::dist2(colvarvalue const &x2) const
{
  colvarvalue::check_types(*this, x2);

  switch (this->value_type) {
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
    // TODO
    break;
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return 0.0;
  };
}


inline colvarvalue colvarvalue::dist2_grad(colvarvalue const &x2) const
{
  colvarvalue::check_types(*this, x2);

  switch (this->value_type) {
  case colvarvalue::type_scalar:
    return 2.0 * (this->real_value - x2.real_value);
  case colvarvalue::type_3vector:
    return 2.0 * (this->rvector_value - x2.rvector_value);
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
                          colvarvalue::type_unit3vector );
    }
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return this->quaternion_value.dist2_grad(x2.quaternion_value);
  case colvarvalue::type_vector:
    // TODO
    break;
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}



inline cvm::vector1d<cvm::real> const colvarvalue::as_vector() const
{
  switch (value_type) {
  case colvarvalue::type_scalar:
    return cvm::vector1d<cvm::real>(1, real_value);
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    return rvector_value.as_vector();
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return quaternion_value.as_vector();
    break;
  case colvarvalue::type_vector:
    if (vector1d_value.size() == 0) {
      cvm::error("Error: trying to use as vector a variable that is not yet initialized as such\".\n");
    }
    return vector1d_value;
    break;
  case colvarvalue::type_notset:
  default:
    break;
  }
}


inline void colvarvalue::check_types(colvarvalue const &x1,
                                     colvarvalue const &x2)
{
  if (!colvarvalue::type_checking()) return;

  if (x1.value_type != x2.value_type) {
    if (((x1.value_type == type_unit3vector) &&
         (x2.value_type == type_unit3vectorderiv)) ||
        ((x2.value_type == type_unit3vector) &&
         (x1.value_type == type_unit3vectorderiv)) ||
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

  if (x1.value_type == type_vector) {
    if (x1.vector1d_value.size() != x2.vector1d_value.size()) {
      cvm::error("Performing an operation between two vector colvar "
                 "values with different sizes, "+
                 cvm::to_str(x1.vector1d_value.size())+
                 " and "+
                 cvm::to_str(x2.vector1d_value.size())+
                 ".\n");
    }
  }
}

inline void colvarvalue::check_types_assign(colvarvalue::Type const &vt1,
                                            colvarvalue::Type const &vt2)
{
  if (vt1 != type_notset) {
    if (vt1 != vt2) {
      cvm::error("Trying to assign a colvar value with type \""+
                 type_desc(vt2)+"\" to one with type \""+
                 type_desc(vt1)+"\".\n");
    }
  }
}

inline void colvarvalue::undef_op() const
{
  cvm::error("Error: Undefined operation on a colvar of type \""+
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
    while (xvi != xv_end) {
      cvm::real const cosine =
        ((xvi)->rvector_value * x.rvector_value) /
        ((xvi)->rvector_value.norm() * x.rvector_value.norm());
      xvi++;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
    while (xvi != xv_end) {
      cvm::real const cosine =
        ((xvi)->rvector_value * x.rvector_value) /
        ((xvi)->rvector_value.norm() * x.rvector_value.norm());
      xvi++;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vectorderiv:
    is >> x.rvector_value;
    break;
  case colvarvalue::type_unit3vector:
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
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
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
