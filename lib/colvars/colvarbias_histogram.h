// -*- c++ -*-

#ifndef COLVARBIAS_HISTOGRAM_H
#define COLVARBIAS_HISTOGRAM_H

#include <vector>
#include <list>
#include <sstream>
#include <iomanip>

#include "colvarbias.h"
#include "colvargrid.h"

/// Histogram "bias" (does as the name says)
class colvarbias_histogram : public colvarbias {

public:

  colvarbias_histogram(char const *key);
  ~colvarbias_histogram();
  virtual int init(std::string const &conf);
  virtual int update();
  virtual int write_output_files();

protected:

  /// n-dim histogram
  colvar_grid_scalar *grid;
  std::vector<int> bin;
  std::string out_name, out_name_dx;
  size_t output_freq;

  /// If one or more of the variables are \link type_vector \endlink, treat them as arrays of this length
  size_t colvar_array_size;
  /// If colvar_array_size is larger than 1, weigh each one by this number before accumulating the histogram
  std::vector<cvm::real> weights;

  virtual std::istream& read_restart(std::istream&);
  virtual std::ostream& write_restart(std::ostream&);
};

#endif
