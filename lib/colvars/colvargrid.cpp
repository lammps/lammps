// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <ctime>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvargrid.h"
#include "colvargrid_def.h"



colvar_grid_count::colvar_grid_count()
  : colvar_grid<size_t>()
{
  mult = 1;
}

colvar_grid_count::colvar_grid_count(std::vector<colvar *>  &colvars,
                                     std::string config)
  : colvar_grid<size_t>(colvars, 0, 1, false, nullptr, config)
{}

colvar_grid_count::colvar_grid_count(std::vector<colvar *>  &colvars,
                                     std::shared_ptr<const colvar_grid_params> params)
  : colvar_grid<size_t>(colvars, 0, 1, false, params)
{}

std::string colvar_grid_count::get_state_params() const
{
  return colvar_grid<size_t>::get_state_params();
}

int colvar_grid_count::parse_params(std::string const &conf,
                                    colvarparse::Parse_Mode const parse_mode)
{
  return colvar_grid<size_t>::parse_params(conf, parse_mode);
}

std::istream & colvar_grid_count::read_restart(std::istream &is)
{
  return colvar_grid<size_t>::read_restart(is);
}

cvm::memory_stream & colvar_grid_count::read_restart(cvm::memory_stream &is)
{
  return colvar_grid<size_t>::read_restart(is);
}

std::ostream & colvar_grid_count::write_restart(std::ostream &os)
{
  return colvar_grid<size_t>::write_restart(os);
}

cvm::memory_stream & colvar_grid_count::write_restart(cvm::memory_stream &os)
{
  return colvar_grid<size_t>::write_restart(os);
}

std::istream &colvar_grid_count::read_raw(std::istream &is)
{
  return colvar_grid<size_t>::read_raw(is);
}

cvm::memory_stream &colvar_grid_count::read_raw(cvm::memory_stream &is)
{
  return colvar_grid<size_t>::read_raw(is);
}

std::ostream &colvar_grid_count::write_raw(std::ostream &os, size_t const buf_size) const
{
  return colvar_grid<size_t>::write_raw(os, buf_size);
}

cvm::memory_stream &colvar_grid_count::write_raw(cvm::memory_stream &os,
                                                 size_t const buf_size) const
{
  return colvar_grid<size_t>::write_raw(os, buf_size);
}

std::istream & colvar_grid_count::read_multicol(std::istream &is, bool add)
{
  return colvar_grid<size_t>::read_multicol(is, add);
}

int colvar_grid_count::read_multicol(std::string const &filename,
                                     std::string description,
                                     bool add)
{
  return colvar_grid<size_t>::read_multicol(filename, description, add);
}

std::ostream & colvar_grid_count::write_multicol(std::ostream &os) const
{
  return colvar_grid<size_t>::write_multicol(os);
}

int colvar_grid_count::write_multicol(std::string const &filename,
                                      std::string description) const
{
  return colvar_grid<size_t>::write_multicol(filename, description);
}

std::ostream & colvar_grid_count::write_opendx(std::ostream &os) const
{
  return colvar_grid<size_t>::write_opendx(os);
}

int colvar_grid_count::write_opendx(std::string const &filename,
                                    std::string description) const
{
  return colvar_grid<size_t>::write_opendx(filename, description);
}



colvar_grid_scalar::colvar_grid_scalar()
  : colvar_grid<cvm::real>(), samples(NULL)
{}

colvar_grid_scalar::colvar_grid_scalar(colvar_grid_scalar const &g)
  : colvar_grid<cvm::real>(g), samples(NULL)
{
}

colvar_grid_scalar::colvar_grid_scalar(std::vector<colvar *> &colvars,
                                       std::shared_ptr<const colvar_grid_params> params,
                                       bool add_extra_bin,
                                       std::string config)
  : colvar_grid<cvm::real>(colvars, 0.0, 1, add_extra_bin, params, config), samples(NULL)
{
}

colvar_grid_scalar::colvar_grid_scalar(std::string const &filename)
  : colvar_grid<cvm::real>(filename, 1),
    samples(nullptr)
{
}

colvar_grid_scalar::~colvar_grid_scalar()
{
}

std::string colvar_grid_scalar::get_state_params() const
{
  return colvar_grid<cvm::real>::get_state_params();
}

int colvar_grid_scalar::parse_params(std::string const &conf,
                                    colvarparse::Parse_Mode const parse_mode)
{
  return colvar_grid<cvm::real>::parse_params(conf, parse_mode);
}

std::istream &colvar_grid_scalar::read_restart(std::istream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

cvm::memory_stream &colvar_grid_scalar::read_restart(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

std::ostream &colvar_grid_scalar::write_restart(std::ostream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

cvm::memory_stream &colvar_grid_scalar::write_restart(cvm::memory_stream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

std::istream &colvar_grid_scalar::read_raw(std::istream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

cvm::memory_stream &colvar_grid_scalar::read_raw(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

std::ostream &colvar_grid_scalar::write_raw(std::ostream &os, size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

cvm::memory_stream &colvar_grid_scalar::write_raw(cvm::memory_stream &os,
                                                  size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

std::istream & colvar_grid_scalar::read_multicol(std::istream &is, bool add)
{
  return colvar_grid<cvm::real>::read_multicol(is, add);
}

int colvar_grid_scalar::read_multicol(std::string const &filename,
                                      std::string description,
                                      bool add)
{
  return colvar_grid<cvm::real>::read_multicol(filename, description, add);
}

std::ostream & colvar_grid_scalar::write_multicol(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_multicol(os);
}

int colvar_grid_scalar::write_multicol(std::string const &filename,
                                       std::string description) const
{
  return colvar_grid<cvm::real>::write_multicol(filename, description);
}

std::ostream & colvar_grid_scalar::write_opendx(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_opendx(os);
}

int colvar_grid_scalar::write_opendx(std::string const &filename,
                                     std::string description) const
{
  return colvar_grid<cvm::real>::write_opendx(filename, description);
}


cvm::real colvar_grid_scalar::maximum_value() const
{
  cvm::real max = data[0];
  for (size_t i = 0; i < nt; i++) {
    if (data[i] > max) max = data[i];
  }
  return max;
}


cvm::real colvar_grid_scalar::minimum_value() const
{
  cvm::real min = data[0];
  for (size_t i = 0; i < nt; i++) {
    if (data[i] < min) min = data[i];
  }
  return min;
}

cvm::real colvar_grid_scalar::minimum_pos_value() const
{
  cvm::real minpos = data[0];
  size_t i;
  for (i = 0; i < nt; i++) {
    if(data[i] > 0) {
      minpos = data[i];
      break;
    }
  }
  for (i = 0; i < nt; i++) {
    if (data[i] > 0 && data[i] < minpos) minpos = data[i];
  }
  return minpos;
}

cvm::real colvar_grid_scalar::integral() const
{
  cvm::real sum = 0.0;
  for (size_t i = 0; i < nt; i++) {
    sum += data[i];
  }
  cvm::real bin_volume = 1.0;
  for (size_t id = 0; id < widths.size(); id++) {
    bin_volume *= widths[id];
  }
  return bin_volume * sum;
}


cvm::real colvar_grid_scalar::entropy() const
{
  cvm::real sum = 0.0;
  for (size_t i = 0; i < nt; i++) {
    if (data[i] >0) {
      sum += -1.0 * data[i] * cvm::logn(data[i]);
    }
  }
  cvm::real bin_volume = 1.0;
  for (size_t id = 0; id < widths.size(); id++) {
    bin_volume *= widths[id];
  }
  return bin_volume * sum;
}

/// \brief Return the RMSD between this grid's data and another one
/// Grids must have the same dimensions.
cvm::real colvar_grid_scalar::grid_rmsd(colvar_grid_scalar const &other_grid) const
{
  if (other_grid.data.size() != this->data.size()) {
    cvm::error("Error: trying to subtract two grids with "
                "different size.\n");
    return -1.;
  }

  cvm::real sum2 = 0.0;

  if (samples && other_grid.samples) {
    for (size_t i = 0; i < data.size(); i++) {
      size_t n = samples->get_value(i);
      cvm::real us = n ? data[i] / n : 0.0;
      n = other_grid.samples->get_value(i);
      cvm::real them = n ? other_grid.data[i ] / n : 0.0;
      cvm::real d = us - them;
      sum2 += d*d;
    }
  } else {
    for (size_t i = 0; i < data.size(); i++) {
      cvm::real d = other_grid.data[i] - data[i];
      sum2 += d*d;
    }
  }

  return sqrt(sum2/this->data.size());
}


colvar_grid_gradient::colvar_grid_gradient()
  : colvar_grid<cvm::real>(), samples(NULL)
{}


// colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars, std::string config)
//   : colvar_grid<cvm::real>(colvars, 0.0, colvars.size(), false, nullptr, config), samples(NULL)
// {}

// colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars,
//                                            std::shared_ptr<colvar_grid_count> samples_in)
//   : colvar_grid<cvm::real>(colvars, 0.0, colvars.size(), false, samples_in), samples(samples_in)
// {
//   if (samples_in)
//     samples_in->has_parent_data = true;
// }

colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars,
                                           std::shared_ptr<colvar_grid_count> samples_in,
                                           std::shared_ptr<const colvar_grid_params> params,
                                           std::string config)
  : colvar_grid<cvm::real>(colvars, 0.0, colvars.size(), false, params, config), samples(samples_in)
{
  if (samples_in)
    samples_in->has_parent_data = true;
}


colvar_grid_gradient::colvar_grid_gradient(std::string const &filename)
  : colvar_grid<cvm::real>(filename, 0),
    samples(nullptr)
{
}

std::string colvar_grid_gradient::get_state_params() const
{
  return colvar_grid<cvm::real>::get_state_params();
}

int colvar_grid_gradient::parse_params(std::string const &conf,
                                       colvarparse::Parse_Mode const parse_mode)
{
  return colvar_grid<cvm::real>::parse_params(conf, parse_mode);
}

std::istream &colvar_grid_gradient::read_restart(std::istream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

cvm::memory_stream &colvar_grid_gradient::read_restart(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_restart(is);
}

std::ostream &colvar_grid_gradient::write_restart(std::ostream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

cvm::memory_stream &colvar_grid_gradient::write_restart(cvm::memory_stream &os)
{
  return colvar_grid<cvm::real>::write_restart(os);
}

std::istream &colvar_grid_gradient::read_raw(std::istream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

cvm::memory_stream &colvar_grid_gradient::read_raw(cvm::memory_stream &is)
{
  return colvar_grid<cvm::real>::read_raw(is);
}

std::ostream &colvar_grid_gradient::write_raw(std::ostream &os, size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

cvm::memory_stream &colvar_grid_gradient::write_raw(cvm::memory_stream &os,
                                                    size_t const buf_size) const
{
  return colvar_grid<cvm::real>::write_raw(os, buf_size);
}

std::istream & colvar_grid_gradient::read_multicol(std::istream &is, bool add)
{
  return colvar_grid<cvm::real>::read_multicol(is, add);
}

int colvar_grid_gradient::read_multicol(std::string const &filename,
                                        std::string description,
                                        bool add)
{
  return colvar_grid<cvm::real>::read_multicol(filename, description, add);
}

std::ostream & colvar_grid_gradient::write_multicol(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_multicol(os);
}

int colvar_grid_gradient::write_multicol(std::string const &filename,
                                         std::string description) const
{
  return colvar_grid<cvm::real>::write_multicol(filename, description);
}

std::ostream & colvar_grid_gradient::write_opendx(std::ostream &os) const
{
  return colvar_grid<cvm::real>::write_opendx(os);
}

int colvar_grid_gradient::write_opendx(std::string const &filename,
                                       std::string description) const
{
  return colvar_grid<cvm::real>::write_opendx(filename, description);
}


void colvar_grid_gradient::write_1D_integral(std::ostream &os)
{
  cvm::real bin, min, integral;
  std::vector<cvm::real> int_vals;

  os << "#       xi            A(xi)\n";

  if (cv.size() != 1) {
    cvm::error("Cannot write integral for multi-dimensional gradient grids.");
    return;
  }

  integral = 0.0;
  int_vals.push_back(0.0);
  min = 0.0;

  // correction for periodic colvars, so that the PMF is periodic
  cvm::real corr;
  if (periodic[0]) {
    corr = average();
  } else {
    corr = 0.0;
  }

  for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {

    if (samples) {
      size_t const samples_here = samples->value(ix);
      if (samples_here)
        integral += (value(ix) / cvm::real(samples_here) - corr) * cv[0]->width;
    } else {
      integral += (value(ix) - corr) * cv[0]->width;
    }

    if ( integral < min ) min = integral;
    int_vals.push_back(integral);
  }

  bin = 0.0;
  for ( int i = 0; i < nx[0]; i++, bin += 1.0 ) {
    os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
       << std::setw(cvm::cv_width)
       << std::setprecision(cvm::cv_prec)
       << int_vals[i] - min << "\n";
  }

  os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
     << std::setw(cvm::cv_width)
     << std::setprecision(cvm::cv_prec)
     << int_vals[nx[0]] - min << "\n";

  return;
}


/// \brief Return the RMSD between this grid's data and another one
/// Grids must have the same dimensions.
cvm::real colvar_grid_gradient::grid_rmsd(colvar_grid_gradient const &other_grid) const
{
  if (other_grid.multiplicity() != this->multiplicity()) {
    cvm::error("Error: trying to subtract two grids with "
                "different multiplicity.\n");
    return -1.;
  }

  if (other_grid.data.size() != this->data.size()) {
    cvm::error("Error: trying to subtract two grids with "
                "different size.\n");
    return -1.;
  }

  cvm::real sum2 = 0.0;
  std::vector<int> ix;
  size_t imult;
  for (ix = new_index(); index_ok(ix); incr(ix)) {
    for (imult = 0; imult < this->multiplicity(); imult++) {
      cvm::real d = this->value_output(ix, imult) - other_grid.value_output(ix, imult);
      sum2 += d*d;
    }
  }
  return sqrt(sum2/this->data.size());
}

