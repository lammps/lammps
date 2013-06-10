// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvargrid.h"


colvar_grid_count::colvar_grid_count()
  : colvar_grid<size_t>()
{}

colvar_grid_count::colvar_grid_count (std::vector<int> const &nx_i,
                                      size_t const           &def_count)
  : colvar_grid<size_t> (nx_i, def_count)
{}

colvar_grid_count::colvar_grid_count (std::vector<colvar *>  &colvars,
                                      size_t const           &def_count)
  : colvar_grid<size_t> (colvars, def_count)
{}

std::istream & colvar_grid_count::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();
  std::string key, conf;
  if ((is >> key) && (key == std::string ("grid_parameters"))) {
    is.seekg (start_pos, std::ios::beg);
    is >> colvarparse::read_block ("grid_parameters", conf);
    parse_params (conf);
  } else {
    cvm::log ("Grid parameters are missing in the restart file, using those from the configuration.\n");
    is.seekg (start_pos, std::ios::beg);
  }
  read_raw (is);
  return is;
}

std::ostream & colvar_grid_count::write_restart (std::ostream &os)
{
  write_params (os);
  write_raw (os);
  return os;
}



colvar_grid_scalar::colvar_grid_scalar()
  : colvar_grid<cvm::real>(), samples (NULL), grad (NULL)
{}

colvar_grid_scalar::colvar_grid_scalar (colvar_grid_scalar const &g)
  : colvar_grid<cvm::real> (g), samples (NULL), grad (NULL)
{
  grad = new cvm::real[nd];
}

colvar_grid_scalar::colvar_grid_scalar (std::vector<int> const &nx_i)
  : colvar_grid<cvm::real> (nx_i, 0.0, 1), samples (NULL)
{
  grad = new cvm::real[nd];
}

colvar_grid_scalar::colvar_grid_scalar (std::vector<colvar *> &colvars, bool margin)
  : colvar_grid<cvm::real> (colvars, 0.0, 1, margin), samples (NULL)
{
  grad = new cvm::real[nd];
}

colvar_grid_scalar::~colvar_grid_scalar()
{
  if (grad) {
    delete [] grad;
    grad = NULL;
  }
}

std::istream & colvar_grid_scalar::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();
  std::string key, conf;
  if ((is >> key) && (key == std::string ("grid_parameters"))) {
    is.seekg (start_pos, std::ios::beg);
    is >> colvarparse::read_block ("grid_parameters", conf);
    parse_params (conf);
  } else {
    cvm::log ("Grid parameters are missing in the restart file, using those from the configuration.\n");
    is.seekg (start_pos, std::ios::beg);
  }
  read_raw (is);
  return is;
}

std::ostream & colvar_grid_scalar::write_restart (std::ostream &os)
{
  write_params (os);
  write_raw (os);
  return os;
}



colvar_grid_gradient::colvar_grid_gradient()
  : colvar_grid<cvm::real>(), samples (NULL)
{}

colvar_grid_gradient::colvar_grid_gradient (std::vector<int> const &nx_i)
  : colvar_grid<cvm::real> (nx_i, 0.0, nx_i.size()), samples (NULL)
{}

colvar_grid_gradient::colvar_grid_gradient (std::vector<colvar *> &colvars)
  : colvar_grid<cvm::real> (colvars, 0.0, colvars.size()), samples (NULL)
{}

std::istream & colvar_grid_gradient::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();
  std::string key, conf;
  if ((is >> key) && (key == std::string ("grid_parameters"))) {
    is.seekg (start_pos, std::ios::beg);
    is >> colvarparse::read_block ("grid_parameters", conf);
    parse_params (conf);
  } else {
    cvm::log ("Grid parameters are missing in the restart file, using those from the configuration.\n");
    is.seekg (start_pos, std::ios::beg);
  }
  read_raw (is);
  return is;
}

std::ostream & colvar_grid_gradient::write_restart (std::ostream &os)
{
  write_params (os);
  write_raw (os);
  return os;
}

void colvar_grid_gradient::write_1D_integral (std::ostream &os)
{
  cvm::real bin, min, integral;
  std::vector<cvm::real> int_vals;

  os << "#       xi            A(xi)\n";

  if ( cv.size() != 1 ) {
    cvm::fatal_error ("Cannot write integral for multi-dimensional gradient grids.");
  }

  integral = 0.0;
  int_vals.push_back ( 0.0 );
  bin = 0.0;
  min = 0.0;

  // correction for periodic colvars, so that the PMF is periodic
  cvm::real corr;
  if ( periodic[0] ) {
    corr = average();
  } else {
    corr = 0.0;
  }

  for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix), bin += 1.0 ) {

    if (samples) {
      size_t const samples_here = samples->value (ix);
      if (samples_here)
        integral += (value (ix) / cvm::real (samples_here) - corr) * cv[0]->width;
    } else {
      integral += (value (ix) - corr) * cv[0]->width;
    }

    if ( integral < min ) min = integral;
    int_vals.push_back ( integral );
  }

  bin = 0.0;
  for ( int i = 0; i < nx[0]; i++, bin += 1.0 ) {
    os << std::setw (10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
       << std::setw (cvm::cv_width)
       << std::setprecision (cvm::cv_prec)
       << int_vals[i] - min << "\n";
  }

  os << std::setw (10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
     << std::setw (cvm::cv_width)
     << std::setprecision (cvm::cv_prec)
     << int_vals[nx[0]] - min << "\n";

  return;
}




// quaternion_grid::quaternion_grid (std::vector<colvar *>      const &cv_i,
//                                   std::vector<std::string>   const &grid_str)
// {
//   cv = cv_i;

//   std::istringstream is (grid_str[0]);
//   is >> grid_size;

//   min.assign (3, -1.0);
//   max.assign (3,  1.0);
//   np.assign  (3, grid_size);
//   dx.assign  (3, 2.0/(cvm::real (grid_size)));

//   // assumes a uniform grid in the three directions; change
//   // get_value() if you want to use different sizes
//   cvm::log ("Allocating quaternion grid ("+cvm::to_str (np.size())+" dimensional)...");
//   data.create (np, 0.0);
//   cvm::log ("done.\n");
//   if (cvm::debug()) cvm::log ("Grid size = "+data.size());
// }

