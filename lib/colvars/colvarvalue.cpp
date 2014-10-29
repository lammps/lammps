/// -*- c++ -*-

#include <vector>
#include <sstream>
#include <iostream>

#include "colvarmodule.h"
#include "colvarvalue.h"



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


colvarvalue colvarvalue::inverse() const
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
        size_t i;
        for (i = 0; i < elem_types.size(); i++) {
          result.sliceassign(elem_indices[i], elem_indices[i]+elem_sizes[i],
                             cvm::vector1d<cvm::real>((this->get_elem(i)).inverse()));
        }
      } else {
        size_t i;
        for (i = 0; i < result.size(); i++) {
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
    return colvarvalue(2.0 * (this->vector1d_value - x2.vector1d_value), colvarvalue::type_vector);
    break;
  case colvarvalue::type_notset:
  default:
    this->undef_op();
    return colvarvalue(colvarvalue::type_notset);
  };
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



