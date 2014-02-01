// -*- c++ -*-

#ifndef COLVARGRID_H
#define COLVARGRID_H

#include <iostream>
#include <iomanip>
#include <cmath>

#include "colvar.h"
#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"

/// \brief Grid of values of a function of several collective
/// variables \param T The data type
///
/// Only scalar colvars supported so far
template <class T> class colvar_grid : public colvarparse {

protected:

  /// Number of dimensions
  size_t nd;

  /// Number of points along each dimension
  std::vector<int> nx;

  /// Cumulative number of points along each dimension
  std::vector<int> nxc;

  /// \brief Multiplicity of each datum (allow the binning of
  /// non-scalar types)
  size_t mult;

  /// Total number of grid points
  size_t nt;

  /// Low-level array of values
  std::vector<T> data;

  /// Newly read data (used for count grids, when adding several grids read from disk)
  std::vector<size_t> new_data;

  /// Colvars collected in this grid
  std::vector<colvar *> cv;

  /// Do we request actual value (for extended-system colvars)?
  std::vector<bool> actual_value;

  /// Get the low-level index corresponding to an index
  inline size_t address (std::vector<int> const &ix) const
  {
    size_t addr = 0;
    for (size_t i = 0; i < nd; i++) {
      addr += ix[i]*nxc[i];
      if (cvm::debug()) {
        if (ix[i] >= nx[i])
          cvm::fatal_error ("Error: exceeding bounds in colvar_grid.\n");
      }
    }
    return addr;
  }

public:

  /// Lower boundaries of the colvars in this grid
  std::vector<colvarvalue> lower_boundaries;

  /// Upper boundaries of the colvars in this grid
  std::vector<colvarvalue> upper_boundaries;

  /// Whether some colvars are periodic
  std::vector<bool>        periodic;

  /// Whether some colvars have hard lower boundaries
  std::vector<bool>        hard_lower_boundaries;

  /// Whether some colvars have hard upper boundaries
  std::vector<bool>        hard_upper_boundaries;

  /// Widths of the colvars in this grid
  std::vector<cvm::real>   widths;

  /// True if this is a count grid related to another grid of data
  bool has_parent_data;

  /// Whether this grid has been filled with data or is still empty
  bool has_data;

  /// Return the number of colvars
  inline size_t number_of_colvars() const
  {
    return nd;
  }

  /// Return the number of points in the i-th direction, if provided, or
  /// the total number
  inline size_t number_of_points (int const icv = -1) const
  {
    if (icv < 0) {
      return nt;
    } else {
      return nx[icv];
    }
  }

  /// Get the sizes in each direction
  inline std::vector<int> const & sizes() const
  {
    return nx;
  }

  /// Set the sizes in each direction
  inline void set_sizes (std::vector<int> const &new_sizes)
  {
    nx = new_sizes;
  }

  /// Return the multiplicity of the type used
  inline size_t multiplicity() const
  {
    return mult;
  }

  /// \brief Allocate data (allow initialization also after construction)
  void create (std::vector<int> const &nx_i,
               T const &t = T(),
               size_t const &mult_i = 1)
  {
    mult = mult_i;
    nd = nx_i.size();
    nxc.resize (nd);
    nx = nx_i;

    nt = mult;
    for (int i = nd-1; i >= 0; i--) {
      if (nx_i[i] <= 0)
        cvm::fatal_error ("Error: providing an invalid number of points, "+
                          cvm::to_str (nx_i[i])+".\n");
      nxc[i] = nt;
      nt *= nx[i];
    }

    data.reserve (nt);
    data.assign (nt, t);
  }

  /// \brief Allocate data (allow initialization also after construction)
  void create()
  {
    create (this->nx, T(), this->mult);
  }

  /// \brief Reset data (in case the grid is being reused)
  void reset (T const &t = T())
  {
    data.assign (nt, t);
  }


  /// Default constructor
  colvar_grid() : has_data (false)
  {
    save_delimiters = false;
    nd = nt = 0;
  }

  /// Destructor
  virtual ~colvar_grid()
  {}

  /// \brief "Almost copy-constructor": only copies configuration
  /// parameters from another grid, but doesn't reallocate stuff;
  /// create() must be called after that;
  colvar_grid (colvar_grid<T> const &g) : has_data (false),
                                          nd (g.nd),
                                          nx (g.nx),
                                          mult (g.mult),
                                          cv (g.cv),
                                          lower_boundaries (g.lower_boundaries),
                                          upper_boundaries (g.upper_boundaries),
                                          hard_lower_boundaries (g.hard_lower_boundaries),
                                          hard_upper_boundaries (g.hard_upper_boundaries),
                                          periodic (g.periodic),
                                          widths (g.widths),
                                          actual_value (g.actual_value),
                                          data()
  {
    save_delimiters = false;
  }

  /// \brief Constructor from explicit grid sizes \param nx_i Number
  /// of grid points along each dimension \param t Initial value for
  /// the function at each point (optional) \param mult_i Multiplicity
  /// of each value
  colvar_grid (std::vector<int> const &nx_i,
               T const &t = T(),
               size_t const &mult_i = 1) : has_data (false)
  {
    save_delimiters = false;
    this->create (nx_i, t, mult_i);
  }

  /// \brief Constructor from a vector of colvars
  colvar_grid (std::vector<colvar *> const &colvars,
               T const &t = T(),
               size_t const &mult_i = 1,
               bool margin = false)
    : cv (colvars), has_data (false)
  {
    save_delimiters = false;

    std::vector<int> nx_i;

    if (cvm::debug())
      cvm::log ("Allocating a grid for "+cvm::to_str (colvars.size())+
                " collective variables.\n");

    for (size_t i =  0; i < cv.size(); i++) {

      if (cv[i]->type() != colvarvalue::type_scalar) {
        cvm::fatal_error ("Colvar grids can only be automatically "
                          "constructed for scalar variables.  "
                          "ABF and histogram can not be used; metadynamics "
                          "can be used with useGrids disabled.\n");
      }

      if (cv[i]->width <= 0.0) {
        cvm::fatal_error ("Tried to initialize a grid on a "
                          "variable with negative or zero width.\n");
      }

      if (!cv[i]->tasks[colvar::task_lower_boundary] || !cv[i]->tasks[colvar::task_upper_boundary]) {
        cvm::fatal_error ("Tried to initialize a grid on a "
                          "variable with undefined boundaries.\n");
      }

      widths.push_back (cv[i]->width);
      hard_lower_boundaries.push_back (cv[i]->hard_lower_boundary);
      hard_upper_boundaries.push_back (cv[i]->hard_upper_boundary);
      periodic.push_back (cv[i]->periodic_boundaries());

      // By default, get reported colvar value (for extended Lagrangian colvars)
      actual_value.push_back (false);

      // except if a colvar is specified twice in a row
      // then the first instance is the actual value
      if (i > 0 && cv[i-1] == cv[i]) {
        actual_value[i-1] = true;
      }

      if (margin) {
        if (periodic[i]) {
          // Shift the grid by half the bin width (values at edges instead of center of bins)
          lower_boundaries.push_back (cv[i]->lower_boundary.real_value - 0.5 * widths[i]);
          upper_boundaries.push_back (cv[i]->upper_boundary.real_value - 0.5 * widths[i]);
        } else {
          // Make this grid larger by one bin width
          lower_boundaries.push_back (cv[i]->lower_boundary.real_value - 0.5 * widths[i]);
          upper_boundaries.push_back (cv[i]->upper_boundary.real_value + 0.5 * widths[i]);
        }
      } else {
        lower_boundaries.push_back (cv[i]->lower_boundary);
        upper_boundaries.push_back (cv[i]->upper_boundary);
      }


      {
        cvm::real nbins = ( upper_boundaries[i].real_value -
                            lower_boundaries[i].real_value ) / widths[i];
        int nbins_round = (int)(nbins+0.5);

        if (std::fabs (nbins - cvm::real (nbins_round)) > 1.0E-10) {
          cvm::log ("Warning: grid interval ("+
                    cvm::to_str (lower_boundaries[i], cvm::cv_width, cvm::cv_prec)+" - "+
                    cvm::to_str (upper_boundaries[i], cvm::cv_width, cvm::cv_prec)+
                    ") is not commensurate to its bin width ("+
                    cvm::to_str (widths[i], cvm::cv_width, cvm::cv_prec)+").\n");
          upper_boundaries[i].real_value = lower_boundaries[i].real_value +
            (nbins_round * widths[i]);
        }

        if (cvm::debug())
          cvm::log ("Number of points is "+cvm::to_str ((int) nbins_round)+
                    " for the colvar no. "+cvm::to_str (i+1)+".\n");

        nx_i.push_back (nbins_round);
      }

    }

    create (nx_i, t, mult_i);
  }


  /// Wrap an index vector around periodic boundary conditions
  /// also checks validity of non-periodic indices
  inline void wrap (std::vector<int> & ix) const
  {
    for (size_t i = 0; i < nd; i++) {
      if (periodic[i]) {
        ix[i] = (ix[i] + nx[i]) % nx[i]; //to ensure non-negative result
      } else {
        if (ix[i] < 0 || ix[i] >= nx[i])
          cvm::fatal_error ("Trying to wrap illegal index vector (non-PBC): "
                            + cvm::to_str (ix));
      }
    }
  }

  /// \brief Report the bin corresponding to the current value of variable i
  inline int current_bin_scalar(int const i) const
  {
    return value_to_bin_scalar (actual_value[i] ? cv[i]->actual_value() : cv[i]->value(), i);
  }

  /// \brief Use the lower boundary and the width to report which bin
  /// the provided value is in
  inline int value_to_bin_scalar (colvarvalue const &value, const int i) const
  {
    return (int) std::floor ( (value.real_value - lower_boundaries[i].real_value) / widths[i] );
  }

  /// \brief Same as the standard version, but uses another grid definition
  inline int value_to_bin_scalar (colvarvalue const &value,
                                  colvarvalue const &new_offset,
                                  cvm::real   const &new_width) const
  {
    return (int) std::floor ( (value.real_value - new_offset.real_value) / new_width );
  }

  /// \brief Use the two boundaries and the width to report the
  /// central value corresponding to a bin index
  inline colvarvalue bin_to_value_scalar (int const &i_bin, int const i) const
  {
    return lower_boundaries[i].real_value + widths[i] * (0.5 + i_bin);
  }

  /// \brief Same as the standard version, but uses different parameters
  inline colvarvalue bin_to_value_scalar (int const &i_bin,
                                          colvarvalue const &new_offset,
                                          cvm::real const &new_width) const
  {
    return new_offset.real_value + new_width * (0.5 + i_bin);
  }

  /// Set the value at the point with index ix
  inline void set_value (std::vector<int> const &ix,
                         T const &t,
                         size_t const &imult = 0)
  {
    data[this->address (ix)+imult] = t;
    has_data = true;
  }


  /// \brief Get the binned value indexed by ix, or the first of them
  /// if the multiplicity is larger than 1
  inline T const & value (std::vector<int> const &ix,
                          size_t const &imult = 0) const
  {
    return data[this->address (ix) + imult];
  }


  /// \brief Add a constant to all elements (fast loop)
  inline void add_constant (T const &t)
  {
    for (size_t i = 0; i < nt; i++)
      data[i] += t;
    has_data = true;
  }

  /// \brief Multiply all elements by a scalar constant (fast loop)
  inline void multiply_constant (cvm::real const &a)
  {
    for (size_t i = 0; i < nt; i++)
      data[i] *= a;
  }


  /// \brief Get the bin indices corresponding to the provided values of
  /// the colvars
  inline std::vector<int> const get_colvars_index (std::vector<colvarvalue> const &values) const
  {
    std::vector<int> index = new_index();
    for (size_t i = 0; i < nd; i++) {
      index[i] = value_to_bin_scalar (values[i], i);
    }
    return index;
  }

  /// \brief Get the bin indices corresponding to the current values
  /// of the colvars
  inline std::vector<int> const get_colvars_index() const
  {
    std::vector<int> index = new_index();
    for (size_t i = 0; i < nd; i++) {
      index[i] = current_bin_scalar (i);
    }
    return index;
  }

  /// \brief Get the minimal distance (in number of bins) from the
  /// boundaries; a negative number is returned if the given point is
  /// off-grid
  inline cvm::real bin_distance_from_boundaries (std::vector<colvarvalue> const &values,
                                                 bool skip_hard_boundaries = false)
  {
    cvm::real minimum = 1.0E+16;
    for (size_t i = 0; i < nd; i++) {

      if (periodic[i]) continue;

      cvm::real dl = std::sqrt (cv[i]->dist2 (values[i], lower_boundaries[i])) / widths[i];
      cvm::real du = std::sqrt (cv[i]->dist2 (values[i], upper_boundaries[i])) / widths[i];

      if (values[i].real_value < lower_boundaries[i])
        dl *= -1.0;
      if (values[i].real_value > upper_boundaries[i])
        du *= -1.0;

      if ( ((!skip_hard_boundaries) || (!hard_lower_boundaries[i])) && (dl < minimum))
        minimum = dl;
      if ( ((!skip_hard_boundaries) || (!hard_upper_boundaries[i])) && (du < minimum))
        minimum = du;
    }

    return minimum;
  }


  /// \brief Add data from another grid of the same type
  ///
  /// Note: this function maps other_grid inside this one regardless
  /// of whether it fits or not.
  void map_grid (colvar_grid<T> const &other_grid)
  {
    if (other_grid.multiplicity() != this->multiplicity())
      cvm::fatal_error ("Error: trying to merge two grids with values of "
                        "different multiplicity.\n");

    std::vector<colvarvalue> const &gb  = this->lower_boundaries;
    std::vector<cvm::real> const &gw    = this->widths;
    std::vector<colvarvalue> const &ogb = other_grid.lower_boundaries;
    std::vector<cvm::real> const &ogw   = other_grid.widths;

    std::vector<int> ix = this->new_index();
    std::vector<int> oix = other_grid.new_index();

    if (cvm::debug())
      cvm::log ("Remapping grid...\n");
    for ( ; this->index_ok (ix); this->incr (ix)) {

      for (size_t i = 0; i < nd; i++) {
        oix[i] =
          value_to_bin_scalar (bin_to_value_scalar (ix[i], gb[i], gw[i]),
                               ogb[i],
                               ogw[i]);
      }

      if (! other_grid.index_ok (oix)) {
        continue;
      }

      for (size_t im = 0; im < mult; im++) {
        this->set_value (ix, other_grid.value (oix, im), im);
      }
    }

    has_data = true;
    if (cvm::debug())
      cvm::log ("Remapping done.\n");
  }

  /// \brief Add data from another grid of the same type, AND
  /// identical definition (boundaries, widths)
  void add_grid (colvar_grid<T> const &other_grid,
                 cvm::real scale_factor = 1.0)
  {
    if (other_grid.multiplicity() != this->multiplicity())
      cvm::fatal_error ("Error: trying to sum togetehr two grids with values of "
                        "different multiplicity.\n");
    if (scale_factor != 1.0)
      for (size_t i = 0; i < data.size(); i++) {
        data[i] += scale_factor * other_grid.data[i];
      }
    else
      // skip multiplication if possible
      for (size_t i = 0; i < data.size(); i++) {
        data[i] += other_grid.data[i];
      }
    has_data = true;
  }

  /// \brief Return the value suitable for output purposes (so that it
  /// may be rescaled or manipulated without changing it permanently)
  virtual inline T value_output (std::vector<int> const &ix,
                                 size_t const &imult = 0)
  {
    return value (ix, imult);
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (the two may be different,
  /// e.g. when using colvar_grid_count)
  virtual inline void value_input (std::vector<int> const &ix,
                                   T const &t,
                                   size_t const &imult = 0,
                                   bool add = false)
  {
    if ( add )
      data[address (ix) + imult] += t;
    else
      data[address (ix) + imult] = t;
    has_data = true;
  }

  //   /// Get the pointer to the binned value indexed by ix
  //   inline T const *value_p (std::vector<int> const &ix)
  //   {
  //     return &(data[address (ix)]);
  //   }

  /// \brief Get the index corresponding to the "first" bin, to be
  /// used as the initial value for an index in looping
  inline std::vector<int> const new_index() const
  {
    return std::vector<int> (nd, 0);
  }

  /// \brief Check that the index is within range in each of the
  /// dimensions
  inline bool index_ok (std::vector<int> const &ix) const
  {
    for (size_t i = 0; i < nd; i++) {
      if ( (ix[i] < 0) || (ix[i] >= int (nx[i])) )
        return false;
    }
    return true;
  }

  /// \brief Increment the index, in a way that will make it loop over
  /// the whole nd-dimensional array
  inline void incr (std::vector<int> &ix) const
  {
    for (int i = ix.size()-1; i >= 0; i--) {

      ix[i]++;

      if (ix[i] >= nx[i]) {

        if (i > 0) {
          ix[i] = 0;
          continue;
        } else {
          // this is the last iteration, a non-valid index is being
          // set for the outer index, which will be caught by
          // index_ok()
          ix[0] = nx[0];
          return;
        }
      } else {
        return;
      }
    }
  }

  /// \brief Write the grid parameters (number of colvars, boundaries, width and number of points)
  std::ostream & write_params (std::ostream &os)
  {
    os << "grid_parameters {\n  n_colvars " << nd << "\n";

    os << "  lower_boundaries ";
    for (size_t i = 0; i < nd; i++)
      os << " " << lower_boundaries[i];
    os << "\n";

    os << "  upper_boundaries ";
    for (size_t i = 0; i < nd; i++)
      os << " " << upper_boundaries[i];
    os << "\n";

    os << "  widths ";
    for (size_t i = 0; i < nd; i++)
      os << " " << widths[i];
    os << "\n";

    os << "  sizes ";
    for (size_t i = 0; i < nd; i++)
      os << " " << nx[i];
    os << "\n";

    os << "}\n";
    return os;
  }


  bool parse_params (std::string const &conf)
  {
    std::vector<int> old_nx = nx;
    std::vector<colvarvalue> old_lb = lower_boundaries;

    {
      size_t nd_in = 0;
      colvarparse::get_keyval (conf, "n_colvars", nd_in, nd, colvarparse::parse_silent);
      if (nd_in != nd)
        cvm::fatal_error ("Error: trying to read data for a grid "
                          "that contains a different number of colvars ("+
                          cvm::to_str (nd_in)+") than the grid defined "
                          "in the configuration file ("+cvm::to_str (nd)+
                          ").\n");
    }

    colvarparse::get_keyval (conf, "lower_boundaries",
                             lower_boundaries, lower_boundaries, colvarparse::parse_silent);

    colvarparse::get_keyval (conf, "upper_boundaries",
                             upper_boundaries, upper_boundaries, colvarparse::parse_silent);

    colvarparse::get_keyval (conf, "widths", widths, widths, colvarparse::parse_silent);

    colvarparse::get_keyval (conf, "sizes", nx, nx, colvarparse::parse_silent);

    bool new_params = false;
    for (size_t i = 0; i < nd; i++) {
      if ( (old_nx[i] != nx[i]) ||
           (std::sqrt (cv[i]->dist2 (old_lb[i],
                                     lower_boundaries[i])) > 1.0E-10) ) {
        new_params = true;
      }
    }

    // reallocate the array in case the grid params have just changed
    if (new_params) {
      data.resize (0);
      this->create (nx, T(), mult);
    }

    return true;
  }

  /// \brief Check that the grid information inside (boundaries,
  /// widths, ...) is consistent with the current setting of the
  /// colvars
  void check_consistency()
  {
    for (size_t i = 0; i < nd; i++) {
      if ( (std::sqrt (cv[i]->dist2 (cv[i]->lower_boundary,
                                     lower_boundaries[i])) > 1.0E-10) ||
           (std::sqrt (cv[i]->dist2 (cv[i]->upper_boundary,
                                     upper_boundaries[i])) > 1.0E-10) ||
           (std::sqrt (cv[i]->dist2 (cv[i]->width,
                                     widths[i])) > 1.0E-10) ) {
        cvm::fatal_error ("Error: restart information for a grid is "
                          "inconsistent with that of its colvars.\n");
      }
    }
  }


  /// \brief Check that the grid information inside (boundaries,
  /// widths, ...) is consistent with that of another grid
  void check_consistency (colvar_grid<T> const &other_grid)
  {
    for (size_t i = 0; i < nd; i++) {
      // we skip dist2(), because periodicities and the like should
      // matter: boundaries should be EXACTLY the same (otherwise,
      // map_grid() should be used)
      if ( (std::fabs (other_grid.lower_boundaries[i] -
                       lower_boundaries[i]) > 1.0E-10) ||
           (std::fabs (other_grid.upper_boundaries[i] -
                       upper_boundaries[i]) > 1.0E-10) ||
           (std::fabs (other_grid.widths[i] -
                       widths[i]) > 1.0E-10) ||
           (data.size() != other_grid.data.size()) ) {
      cvm::fatal_error ("Error: inconsistency between "
                        "two grids that are supposed to be equal, "
                        "aside from the data stored.\n");
    }
  }
}


/// \brief Write the grid data without labels, as they are
/// represented in memory
/// \param buf_size Number of values per line
  std::ostream & write_raw (std::ostream &os,
                            size_t const buf_size = 3)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    std::vector<int> ix = new_index();
    size_t count = 0;
    for ( ; index_ok (ix); incr (ix)) {
      for (size_t imult = 0; imult < mult; imult++) {
        os << " "
           << std::setw (w) << std::setprecision (p)
           << value_output (ix, imult);
        if (((++count) % buf_size) == 0)
          os << "\n";
      }
    }
    // write a final newline only if buffer is not empty
    if ((count % buf_size) != 0)
      os << "\n";

    return os;
  }

/// \brief Read data written by colvar_grid::write_raw()
std::istream & read_raw (std::istream &is)
{
  size_t const start_pos = is.tellg();

  for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix)) {
    for (size_t imult = 0; imult < mult; imult++) {
      T new_value;
      if (is >> new_value) {
        value_input (ix, new_value, imult);
      } else {
        is.clear();
        is.seekg (start_pos, std::ios::beg);
        is.setstate (std::ios::failbit);
        return is;
      }
    }
  }

  has_data = true;
  return is;
}

/// \brief To be called after colvar_grid::read_raw() returns an error
void read_raw_error()
{
  cvm::fatal_error ("Error: failed to read all of the grid points from file.  Possible explanations: grid parameters in the configuration (lowerBoundary, upperBoundary, width) are different from those in the file, or the file is corrupt/incomplete.\n");
}

/// \brief Write the grid in a format which is both human readable
/// and suitable for visualization e.g. with gnuplot
void write_multicol (std::ostream &os)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  // Data in the header: nColvars, then for each
  // xiMin, dXi, nPoints, periodic

  os << std::setw (2) << "# " << nd << "\n";
  for (size_t i = 0; i < nd; i++) {
    os << "# "
       << std::setw (10) << lower_boundaries[i]
       << std::setw (10) << widths[i]
       << std::setw (10) << nx[i] << "  "
       << periodic[i] << "\n";
  }

  for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix) ) {

    if (ix.back() == 0) {
      // if the last index is 0, add a new line to mark the new record
      os << "\n";
    }

    for (size_t i = 0; i < nd; i++) {
      os << " "
         << std::setw (w) << std::setprecision (p)
         << bin_to_value_scalar (ix[i], i);
    }
    os << " ";
    for (size_t imult = 0; imult < mult; imult++) {
      os << " "
         << std::setw (w) << std::setprecision (p)
         << value_output (ix, imult);
    }
    os << "\n";
  }
}

/// \brief Read a grid written by colvar_grid::write_multicol()
/// Adding data if add is true, replacing if false
std::istream & read_multicol (std::istream &is, bool add = false)
{
  // Data in the header: nColvars, then for each
  // xiMin, dXi, nPoints, periodic

  std::string   hash;
  cvm::real     lower, width, x;
  size_t        n, periodic;
  bool          remap;
  std::vector<T>        new_value;
  std::vector<int>      nx_read;
  std::vector<int>      bin;

  if ( cv.size() != nd ) {
    cvm::fatal_error ("Cannot read grid file: missing reference to colvars.");
  }

  if ( !(is >> hash) || (hash != "#") ) {
    cvm::fatal_error ("Error reading grid at position "+
                      cvm::to_str (is.tellg())+" in stream (read \"" + hash + "\")\n");
  }

  is >> n;
  if ( n != nd ) {
    cvm::fatal_error ("Error reading grid: wrong number of collective variables.\n");
  }

  nx_read.resize (n);
  bin.resize (n);
  new_value.resize (mult);

  if (this->has_parent_data && add) {
    new_data.resize (data.size());
  }

  remap = false;
  for (size_t i = 0; i < nd; i++ ) {
    if ( !(is >> hash) || (hash != "#") ) {
      cvm::fatal_error ("Error reading grid at position "+
                        cvm::to_str (is.tellg())+" in stream (read \"" + hash + "\")\n");
    }

    is >> lower >> width >> nx_read[i] >> periodic;


    if ( (std::fabs (lower - lower_boundaries[i].real_value) > 1.0e-10) ||
         (std::fabs (width - widths[i] ) > 1.0e-10) ||
         (nx_read[i] != nx[i]) ) {
      cvm::log ("Warning: reading from different grid definition (colvar "
                + cvm::to_str (i+1) + "); remapping data on new grid.\n");
      remap = true;
    }
  }

  if ( remap ) {
    // re-grid data
    while (is.good()) {
      bool end_of_file = false;

      for (size_t i = 0; i < nd; i++ ) {
        if ( !(is >> x) ) end_of_file = true;
        bin[i] = value_to_bin_scalar (x, i);
      }
      if (end_of_file) break;

      for (size_t imult = 0; imult < mult; imult++) {
        is >> new_value[imult];
      }

      if ( index_ok(bin) ) {
        for (size_t imult = 0; imult < mult; imult++) {
          value_input (bin, new_value[imult], imult, add);
        }
      }
    }
  } else {
    // do not re-grid the data but assume the same grid is used
    for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix) ) {
      for (size_t i = 0; i < nd; i++ ) {
        is >> x;
      }
      for (size_t imult = 0; imult < mult; imult++) {
        is >> new_value[imult];
        value_input (ix, new_value[imult], imult, add);
      }
    }
  }
  has_data = true;
  return is;
}

};



/// \brief Colvar_grid derived class to hold counters in discrete
/// n-dim colvar space
class colvar_grid_count : public colvar_grid<size_t>
{
public:

  /// Default constructor
  colvar_grid_count();

  /// Destructor
  virtual inline ~colvar_grid_count()
  {}

  /// Constructor
  colvar_grid_count (std::vector<int> const &nx_i,
                     size_t const           &def_count = 0);

  /// Constructor from a vector of colvars
  colvar_grid_count (std::vector<colvar *>  &colvars,
                     size_t const           &def_count = 0);

  /// Increment the counter at given position
  inline void incr_count (std::vector<int> const &ix)
  {
    ++(data[this->address (ix)]);
  }

  /// \brief Get the binned count indexed by ix from the newly read data
  inline size_t const & new_count (std::vector<int> const &ix,
                                   size_t const &imult = 0)
  {
    return new_data[address (ix) + imult];
  }

  /// \brief Read the grid from a restart
  std::istream & read_restart (std::istream &is);

  /// \brief Write the grid to a restart
  std::ostream & write_restart (std::ostream &os);

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual inline void value_input (std::vector<int> const &ix,
                                   size_t const &t,
                                   size_t const &imult = 0,
                                   bool add = false)
  {
    if (add) {
      data[address (ix)] += t;
      if (this->has_parent_data) {
        // save newly read data for inputting parent grid
        new_data[address (ix)] = t;
      }
    } else {
      data[address (ix)] = t;
    }
    has_data = true;
  }
};


/// Class for accumulating a scalar function on a grid
class colvar_grid_scalar : public colvar_grid<cvm::real>
{
public:

  /// \brief Provide the associated sample count by which each binned value
  /// should be divided
  colvar_grid_count *samples;

  /// Default constructor
  colvar_grid_scalar();

  /// Copy constructor (needed because of the grad pointer)
  colvar_grid_scalar (colvar_grid_scalar const &g);

  /// Destructor
  ~colvar_grid_scalar();

  /// Constructor from specific sizes arrays
  colvar_grid_scalar (std::vector<int> const &nx_i);

  /// Constructor from a vector of colvars
  colvar_grid_scalar (std::vector<colvar *> &colvars,
                      bool margin = 0);

  /// Accumulate the value
  inline void acc_value (std::vector<int> const &ix,
                         cvm::real const &new_value,
                         size_t const &imult = 0)
  {
    // only legal value of imult here is 0
    data[address (ix)] += new_value;
    if (samples)
      samples->incr_count (ix);
    has_data = true;
  }

  /// Return the gradient of the scalar field from finite differences
  inline const cvm::real * gradient_finite_diff ( const std::vector<int> &ix0 )
  {
    cvm::real A0, A1;
    std::vector<int> ix;
    if (nd != 2) cvm::fatal_error ("Finite differences available in dimension 2 only.");
    for (int n = 0; n < nd; n++) {
      ix = ix0;
      A0 = data[address (ix)];
      ix[n]++; wrap (ix);
      A1 = data[address (ix)];
      ix[1-n]++; wrap (ix);
      A1 += data[address (ix)];
      ix[n]--; wrap (ix);
      A0 += data[address (ix)];
      grad[n] = 0.5 * (A1 - A0) / widths[n];
    }
    return grad;
  }

  /// \brief Return the value of the function at ix divided by its
  /// number of samples (if the count grid is defined)
  virtual cvm::real value_output (std::vector<int> const &ix,
                                  size_t const &imult = 0)
  {
    if (imult > 0)
      cvm::fatal_error ("Error: trying to access a component "
                        "larger than 1 in a scalar data grid.\n");
    if (samples)
      return (samples->value (ix) > 0) ?
        (data[address (ix)] / cvm::real (samples->value (ix))) :
        0.0;
    else
      return data[address (ix)];
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual void value_input (std::vector<int> const &ix,
                            cvm::real const &new_value,
                            size_t const &imult = 0,
                            bool add = false)
  {
    if (imult > 0)
      cvm::fatal_error ("Error: trying to access a component "
                        "larger than 1 in a scalar data grid.\n");
    if (add) {
      if (samples)
        data[address (ix)] += new_value * samples->new_count (ix);
      else
        data[address (ix)] += new_value;
    } else {
      if (samples)
        data[address (ix)] = new_value * samples->value (ix);
      else
        data[address (ix)] = new_value;
    }
    has_data = true;
  }

  /// \brief Read the grid from a restart
  std::istream & read_restart (std::istream &is);

  /// \brief Write the grid to a restart
  std::ostream & write_restart (std::ostream &os);

  /// \brief Return the highest value
  inline cvm::real maximum_value()
  {
    cvm::real max = data[0];
    for (size_t i = 0; i < nt; i++) {
      if (data[i] > max) max = data[i];
    }
    return max;
  }

  /// \brief Return the lowest value
  inline cvm::real minimum_value()
  {
    cvm::real min = data[0];
    for (size_t i = 0; i < nt; i++) {
      if (data[i] < min) min = data[i];
    }
    return min;
  }

private:
  // gradient
  cvm::real * grad;
};



/// Class for accumulating the gradient of a scalar function on a grid
class colvar_grid_gradient : public colvar_grid<cvm::real>
{
public:

  /// \brief Provide the sample count by which each binned value
  /// should be divided
  colvar_grid_count *samples;

  /// Default constructor
  colvar_grid_gradient();

  /// Destructor
  virtual inline ~colvar_grid_gradient()
  {}

  /// Constructor from specific sizes arrays
  colvar_grid_gradient (std::vector<int> const &nx_i);

  /// Constructor from a vector of colvars
  colvar_grid_gradient (std::vector<colvar *>  &colvars);

  /// \brief Accumulate the gradient
  inline void acc_grad (std::vector<int> const &ix, cvm::real const *grads) {
    for (size_t imult = 0; imult < mult; imult++) {
      data[address (ix) + imult] += grads[imult];
    }
    if (samples)
      samples->incr_count (ix);
  }

  /// \brief Accumulate the gradient based on the force (i.e. sums the
  /// opposite of the force)
  inline void acc_force (std::vector<int> const &ix, cvm::real const *forces) {
    for (size_t imult = 0; imult < mult; imult++) {
      data[address (ix) + imult] -= forces[imult];
    }
    if (samples)
      samples->incr_count (ix);
  }

  /// \brief Return the value of the function at ix divided by its
  /// number of samples (if the count grid is defined)
  virtual inline cvm::real value_output (std::vector<int> const &ix,
                                         size_t const &imult = 0)
  {
    if (samples)
      return (samples->value (ix) > 0) ?
        (data[address (ix) + imult] / cvm::real (samples->value (ix))) :
        0.0;
    else
      return data[address (ix) + imult];
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual inline void value_input (std::vector<int> const &ix,
                                   cvm::real const &new_value,
                                   size_t const &imult = 0,
                                   bool add = false)
  {
    if (add) {
      if (samples)
        data[address (ix) + imult] += new_value * samples->new_count (ix);
      else
        data[address (ix) + imult] += new_value;
    } else {
      if (samples)
        data[address (ix) + imult] = new_value * samples->value (ix);
      else
        data[address (ix) + imult] = new_value;
    }
    has_data = true;
  }


  /// \brief Read the grid from a restart
  std::istream & read_restart (std::istream &is);

  /// \brief Write the grid to a restart
  std::ostream & write_restart (std::ostream &os);

  /// Compute and return average value for a 1D gradient grid
  inline cvm::real average()
  {
    size_t n = 0;

    if (nd != 1 || nx[0] == 0) {
      return 0.0;
    }

    cvm::real sum = 0.0;
    std::vector<int> ix = new_index();
    if (samples) {
      for ( ; index_ok (ix); incr (ix)) {
        if ( (n = samples->value (ix)) )
          sum += value (ix) / n;
      }
    } else {
      for ( ; index_ok (ix); incr (ix)) {
        sum += value (ix);
      }
    }
    return (sum / cvm::real (nx[0]));
  }

  /// \brief If the grid is 1-dimensional, integrate it and write the
  /// integral to a file
  void write_1D_integral (std::ostream &os);

};


#endif

