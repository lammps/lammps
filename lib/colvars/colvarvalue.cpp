/// -*- c++ -*-

#include <vector>

#include "colvarmodule.h"
#include "colvarvalue.h"



std::string const colvarvalue::type_desc[colvarvalue::type_all+1] =
  { "not_set",
    "scalar number",
    "3-dimensional vector",
    "3-dimensional unit vector",
    "3-dimensional tangent vector",
    "4-dimensional unit quaternion",
    "4-dimensional tangent vector",
  };

std::string const colvarvalue::type_keyword[colvarvalue::type_all+1] =
  { "not_set",
    "scalar",
    "vector",
    "unit_vector",
    "",
    "unit_quaternion",
    "",
  };

size_t const      colvarvalue::dof_num[  colvarvalue::type_all+1] =
  { 0, 1, 3, 2, 2, 3, 3 };


void colvarvalue::undef_op() const
{
  cvm::error ("Error: Undefined operation on a colvar of type \""+
              colvarvalue::type_desc[this->value_type]+"\".\n");
}

void colvarvalue::error_rside
(colvarvalue::Type const &vt) const
{
  cvm::error("Trying to assign a colvar value with type \""+
             type_desc[this->value_type]+"\" to one with type \""+
             type_desc[vt]+"\".\n");
}

void colvarvalue::error_lside(colvarvalue::Type const &vt) const
{
  cvm::error("Trying to use a colvar value with type \""+
             type_desc[vt]+"\" as one of type \""+
             type_desc[this->value_type]+"\".\n");
}



void colvarvalue::inner_opt(colvarvalue                        const &x,
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
      *(ii++) += ((xvi++)->quaternion_value).cosine (x.quaternion_value);
    }
    break;
  default:
    x.undef_op();
  };
}

void colvarvalue::inner_opt(colvarvalue const                      &x,
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
      *(ii++) += ((xvi++)->quaternion_value).cosine (x.quaternion_value);
    }
    break;
  default:
    x.undef_op();
  };
}


void colvarvalue::p2leg_opt(colvarvalue const                        &x,
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
      cvm::real const cosine = (xvi++)->quaternion_value.cosine (x.quaternion_value);
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
      cvm::real const cosine = (xvi++)->quaternion_value.cosine (x.quaternion_value);
      *(ii++) += 1.5*cosine*cosine - 0.5;
    }
    break;
  default:
    x.undef_op();
  };
}


std::ostream & operator << (std::ostream &os, colvarvalue const &x)
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


std::ostream & operator << (std::ostream &os, std::vector<colvarvalue> const &v)
{
  for (size_t i = 0; i < v.size(); i++) {
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


size_t colvarvalue::output_width(size_t const &real_width) const
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


