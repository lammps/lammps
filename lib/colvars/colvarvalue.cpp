// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <vector>
#include <sstream>
#include <iostream>

#include "colvarmodule.h"
#include "colvarvalue.h"



colvarvalue::colvarvalue()
  : value_type(type_scalar), real_value(0.0)
{}


colvarvalue::colvarvalue(Type const &vti)
  : value_type(vti), real_value(0.0)
{
  reset();
}

colvarvalue::colvarvalue(cvm::real const &x)
  : value_type(type_scalar), real_value(x)
{}


colvarvalue::colvarvalue(cvm::rvector const &v, colvarvalue::Type vti)
  : value_type(vti), real_value(0.0), rvector_value(v)
{}


colvarvalue::colvarvalue(cvm::quaternion const &q, colvarvalue::Type vti)
  : value_type(vti), real_value(0.0), quaternion_value(q)
{}


colvarvalue::colvarvalue(colvarvalue const &x)
  : value_type(x.type()), real_value(0.0)
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


colvarvalue::colvarvalue(cvm::vector1d<cvm::real> const &v,
                         colvarvalue::Type vti)
  : real_value(0.0)
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


std::string const colvarvalue::type_desc(Type t)
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


std::string const colvarvalue::type_keyword(Type t)
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


size_t colvarvalue::num_df(Type t)
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


size_t colvarvalue::num_dimensions(Type t)
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


void colvarvalue::reset()
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


void colvarvalue::apply_constraints()
{
  switch (value_type) {
  case colvarvalue::type_scalar:
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vectorderiv:
  case colvarvalue::type_quaternionderiv:
    break;
  case colvarvalue::type_unit3vector:
    rvector_value /= cvm::sqrt(rvector_value.norm2());
    break;
  case colvarvalue::type_quaternion:
    quaternion_value /= cvm::sqrt(quaternion_value.norm2());
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


void colvarvalue::type(Type const &vti)
{
  if (vti != value_type) {
    // reset the value based on the previous type
    reset();
    if ((value_type == type_vector) && (vti != type_vector)) {
      vector1d_value.clear();
    }
    value_type = vti;
  }
}


void colvarvalue::type(colvarvalue const &x)
{
  if (x.type() != value_type) {
    // reset the value based on the previous type
    reset();
    if (value_type == type_vector) {
      vector1d_value.clear();
    }
    value_type = x.type();
  }

  if (x.type() == type_vector) {
    vector1d_value.resize(x.vector1d_value.size());
  }
}


void colvarvalue::is_derivative()
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


void colvarvalue::add_elem(colvarvalue const &x)
{
  if (this->value_type != type_vector) {
    cvm::error("Error: trying to set an element for a variable that is not set to be a vector.\n");
    return;
  }
  size_t const n = vector1d_value.size();
  size_t const nd = num_dimensions(x.value_type);
  elem_types.push_back(x.value_type);
  elem_indices.push_back(n);
  elem_sizes.push_back(nd);
  vector1d_value.resize(n + nd);
  set_elem(n, x);
}


colvarvalue const colvarvalue::get_elem(int const i_begin, int const i_end, Type const vt) const
{
  if (vector1d_value.size() > 0) {
    cvm::vector1d<cvm::real> const v(vector1d_value.slice(i_begin, i_end));
    return colvarvalue(v, vt);
  } else {
    cvm::error("Error: trying to get an element from a variable that is not a vector.\n");
    return colvarvalue(type_notset);
  }
}


void colvarvalue::set_elem(int const i_begin, int const i_end, colvarvalue const &x)
{
  if (vector1d_value.size() > 0) {
    vector1d_value.sliceassign(i_begin, i_end, x.as_vector());
  } else {
    cvm::error("Error: trying to set an element for a variable that is not a vector.\n");
  }
}


colvarvalue const colvarvalue::get_elem(int const icv) const
{
  if (elem_types.size() > 0) {
    return get_elem(elem_indices[icv], elem_indices[icv] + elem_sizes[icv],
                    elem_types[icv]);
  } else {
    cvm::error("Error: trying to get a colvarvalue element from a vector colvarvalue that was initialized as a plain array.\n");
    return colvarvalue(type_notset);
  }
}


void colvarvalue::set_elem(int const icv, colvarvalue const &x)
{
  if (elem_types.size() > 0) {
    check_types_assign(elem_types[icv], x.value_type);
    set_elem(elem_indices[icv], elem_indices[icv] + elem_sizes[icv], x);
  } else {
    cvm::error("Error: trying to set a colvarvalue element for a colvarvalue that was initialized as a plain array.\n");
  }
}


void colvarvalue::set_random()
{
  size_t ic;
  switch (this->type()) {
  case colvarvalue::type_scalar:
    this->real_value = cvm::rand_gaussian();
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    this->rvector_value.x = cvm::rand_gaussian();
    this->rvector_value.y = cvm::rand_gaussian();
    this->rvector_value.z = cvm::rand_gaussian();
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    this->quaternion_value.q0 = cvm::rand_gaussian();
    this->quaternion_value.q1 = cvm::rand_gaussian();
    this->quaternion_value.q2 = cvm::rand_gaussian();
    this->quaternion_value.q3 = cvm::rand_gaussian();
    break;
  case colvarvalue::type_vector:
    for (ic = 0; ic < this->vector1d_value.size(); ic++) {
      this->vector1d_value[ic] = cvm::rand_gaussian();
    }
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
    break;
  }
}


void colvarvalue::set_ones(cvm::real assigned_value)
{
  size_t ic;
  switch (this->type()) {
  case colvarvalue::type_scalar:
    this->real_value = assigned_value;
    break;
  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    this->rvector_value.x = assigned_value;
    this->rvector_value.y = assigned_value;
    this->rvector_value.z = assigned_value;
    break;
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    this->quaternion_value.q0 = assigned_value;
    this->quaternion_value.q1 = assigned_value;
    this->quaternion_value.q2 = assigned_value;
    this->quaternion_value.q3 = assigned_value;
    break;
  case colvarvalue::type_vector:
    for (ic = 0; ic < this->vector1d_value.size(); ic++) {
      this->vector1d_value[ic] = assigned_value;
    }
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
    break;
  }
}


void colvarvalue::undef_op() const
{
  cvm::error("Error: Undefined operation on a colvar of type \""+
             type_desc(this->type())+"\".\n");
}


// binary operations between two colvarvalues

colvarvalue operator + (colvarvalue const &x1,
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
    return colvarvalue(x1.vector1d_value + x2.vector1d_value, colvarvalue::type_vector);
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}


colvarvalue operator - (colvarvalue const &x1,
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
    return colvarvalue(x1.vector1d_value - x2.vector1d_value, colvarvalue::type_vector);
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}


// binary operations with real numbers

colvarvalue operator * (cvm::real const &a,
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
    return colvarvalue(x.vector1d_value * a, colvarvalue::type_vector);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  }
}


colvarvalue operator * (colvarvalue const &x,
                        cvm::real const &a)
{
  return a * x;
}


colvarvalue operator / (colvarvalue const &x,
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
    return colvarvalue(x.vector1d_value / a, colvarvalue::type_vector);
  case colvarvalue::type_notset:
  default:
    x.undef_op();
    return colvarvalue(colvarvalue::type_notset);
  }
}


// inner product between two colvarvalues

cvm::real operator * (colvarvalue const &x1,
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
    return (x1.vector1d_value * x2.vector1d_value);
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    return 0.0;
  };
}


colvarvalue colvarvalue::dist2_grad(colvarvalue const &x2) const
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
      return colvarvalue(2.0 * (cos_t * v1 - v2), colvarvalue::type_unit3vectorderiv);
    }
  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    return this->quaternion_value.dist2_grad(x2.quaternion_value);
  case colvarvalue::type_vector:
    return colvarvalue(2.0 * (this->vector1d_value - x2.vector1d_value), colvarvalue::type_vector);
    break;
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
}


/// Return the midpoint between x1 and x2, optionally weighted by lambda
/// (which must be between 0.0 and 1.0)
colvarvalue const colvarvalue::interpolate(colvarvalue const &x1,
                                           colvarvalue const &x2,
                                           cvm::real const lambda)
{
  colvarvalue::check_types(x1, x2);

  if ((lambda < 0.0) || (lambda > 1.0)) {
    cvm::error("Error: trying to interpolate between two colvarvalues with a "
               "lamdba outside [0:1].\n", COLVARS_BUG_ERROR);
  }

  colvarvalue interp = ((1.0-lambda)*x1 + lambda*x2);
  cvm::real const d2 = x1.dist2(x2);

  switch (x1.type()) {
  case colvarvalue::type_scalar:
  case colvarvalue::type_3vector:
  case colvarvalue::type_vector:
  case colvarvalue::type_unit3vectorderiv:
  case colvarvalue::type_quaternionderiv:
    return interp;
    break;
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_quaternion:
    if (interp.norm()/cvm::sqrt(d2) < 1.0e-6) {
      cvm::error("Error: interpolation between "+cvm::to_str(x1)+" and "+
                 cvm::to_str(x2)+" with lambda = "+cvm::to_str(lambda)+
                 " is undefined: result = "+cvm::to_str(interp)+"\n",
                 COLVARS_INPUT_ERROR);
    }
    interp.apply_constraints();
    return interp;
    break;
  case colvarvalue::type_notset:
  default:
    x1.undef_op();
    break;
  }
  return colvarvalue(colvarvalue::type_notset);
}


std::string colvarvalue::to_simple_string() const
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
  case colvarvalue::type_vector:
    return vector1d_value.to_simple_string();
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
    break;
  }
  return std::string();
}


int colvarvalue::from_simple_string(std::string const &s)
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
  case colvarvalue::type_vector:
    return vector1d_value.from_simple_string(s);
    break;
  case colvarvalue::type_notset:
  default:
    undef_op();
    break;
  }
  return COLVARS_ERROR;
}

std::ostream & operator << (std::ostream &os, colvarvalue const &x)
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
  case colvarvalue::type_vector:
    os << x.vector1d_value;
    break;
  case colvarvalue::type_notset:
  default:
    os << "not set";
    break;
  }
  return os;
}


std::ostream & operator << (std::ostream &os, std::vector<colvarvalue> const &v)
{
  size_t i;
  for (i = 0; i < v.size(); i++) {
    os << v[i];
  }
  return os;
}


std::istream & operator >> (std::istream &is, colvarvalue &x)
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
  case colvarvalue::type_vector:
    is >> x.vector1d_value;
    break;
  case colvarvalue::type_notset:
  default:
    x.undef_op();
  }
  return is;
}


size_t colvarvalue::output_width(size_t const &real_width) const
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
  case colvarvalue::type_vector:
    // note how this depends on the vector's size
    return vector1d_value.output_width(real_width);
  case colvarvalue::type_notset:
  default:
    return 0;
  }
}


void colvarvalue::inner_opt(colvarvalue                        const &x,
                            std::vector<colvarvalue>::iterator       &xv,
                            std::vector<colvarvalue>::iterator const &xv_end,
                            std::vector<cvm::real>::iterator         &result)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::vector<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator    &ii = result;

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
  case colvarvalue::type_vector:
    while (xvi != xv_end) {
      *(ii++) += (xvi++)->vector1d_value * x.vector1d_value;
    }
    break;
  default:
    x.undef_op();
  };
}


void colvarvalue::inner_opt(colvarvalue const                      &x,
                            std::list<colvarvalue>::iterator       &xv,
                            std::list<colvarvalue>::iterator const &xv_end,
                            std::vector<cvm::real>::iterator       &result)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::list<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator  &ii = result;

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
  case colvarvalue::type_vector:
    while (xvi != xv_end) {
      *(ii++) += (xvi++)->vector1d_value * x.vector1d_value;
    }
    break;
  default:
    x.undef_op();
  };
}


void colvarvalue::p2leg_opt(colvarvalue const                        &x,
                            std::vector<colvarvalue>::iterator       &xv,
                            std::vector<colvarvalue>::iterator const &xv_end,
                            std::vector<cvm::real>::iterator         &result)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::vector<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator    &ii = result;

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
  case colvarvalue::type_vector:
    while (xvi != xv_end) {
      cvm::real const cosine =
        ((xvi)->vector1d_value * x.vector1d_value) /
        ((xvi)->vector1d_value.norm() * x.rvector_value.norm());
      xvi++;
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  default:
    x.undef_op();
  };
}


void colvarvalue::p2leg_opt(colvarvalue const                        &x,
                            std::list<colvarvalue>::iterator         &xv,
                            std::list<colvarvalue>::iterator const   &xv_end,
                            std::vector<cvm::real>::iterator         &result)
{
  // doing type check only once, here
  colvarvalue::check_types(x, *xv);

  std::list<colvarvalue>::iterator &xvi = xv;
  std::vector<cvm::real>::iterator  &ii = result;

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



