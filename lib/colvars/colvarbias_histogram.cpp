// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias_histogram.h"


colvarbias_histogram::colvarbias_histogram(char const *key)
  : colvarbias(key),
    grid(NULL), out_name("")
{
}


int colvarbias_histogram::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_scalar_variables);
  enable(f_cvb_history_dependent);

  size_t i;

  get_keyval(conf, "outputFile", out_name, std::string(""));
  get_keyval(conf, "outputFileDX", out_name_dx, std::string(""));
  get_keyval(conf, "outputFreq", output_freq, cvm::restart_out_freq);

  /// with VMD, this may not be an error
  // if ( output_freq == 0 ) {
  //   cvm::error("User required histogram with zero output frequency");
  // }

  colvar_array_size = 0;
  {
    bool colvar_array = false;
    get_keyval(conf, "gatherVectorColvars", colvar_array, colvar_array);

    if (colvar_array) {
      for (i = 0; i < num_variables(); i++) { // should be all vector
        if (colvars[i]->value().type() != colvarvalue::type_vector) {
          cvm::error("Error: used gatherVectorColvars with non-vector colvar.\n", INPUT_ERROR);
          return INPUT_ERROR;
        }
        if (i == 0) {
          colvar_array_size = colvars[i]->value().size();
          if (colvar_array_size < 1) {
            cvm::error("Error: vector variable has dimension less than one.\n", INPUT_ERROR);
            return INPUT_ERROR;
          }
        } else {
          if (colvar_array_size != colvars[i]->value().size()) {
            cvm::error("Error: trying to combine vector colvars of different lengths.\n", INPUT_ERROR);
            return INPUT_ERROR;
          }
        }
      }
    } else {
      for (i = 0; i < num_variables(); i++) { // should be all scalar
        if (colvars[i]->value().type() != colvarvalue::type_scalar) {
          cvm::error("Error: only scalar colvars are supported when gatherVectorColvars is off.\n", INPUT_ERROR);
          return INPUT_ERROR;
        }
      }
    }
  }

  if (colvar_array_size > 0) {
    weights.assign(colvar_array_size, 1.0);
    get_keyval(conf, "weights", weights, weights);
  }

  for (i = 0; i < num_variables(); i++) {
    colvars[i]->enable(f_cv_grid);
  }

  grid = new colvar_grid_scalar();
  grid->init_from_colvars(colvars);

  {
    std::string grid_conf;
    if (key_lookup(conf, "histogramGrid", &grid_conf)) {
      grid->parse_params(grid_conf);
      grid->check_keywords(grid_conf, "histogramGrid");
    }
  }

  return COLVARS_OK;
}


colvarbias_histogram::~colvarbias_histogram()
{
  if (grid) {
    delete grid;
    grid = NULL;
  }
}


int colvarbias_histogram::update()
{
  int error_code = COLVARS_OK;
  // update base class
  error_code |= colvarbias::update();

  if (cvm::debug()) {
    cvm::log("Updating histogram bias " + this->name);
  }

  // assign a valid bin size
  bin.assign(num_variables(), 0);

  if (out_name.size() == 0) {
    // At the first timestep, we need to assign out_name since
    // output_prefix is unset during the constructor
    if (cvm::step_relative() == 0) {
      out_name = cvm::output_prefix() + "." + this->name + ".dat";
      cvm::log("Histogram " + this->name + " will be written to file \"" + out_name + "\"");
    }
  }

  if (out_name_dx.size() == 0) {
    if (cvm::step_relative() == 0) {
      out_name_dx = cvm::output_prefix() + "." + this->name + ".dx";
      cvm::log("Histogram " + this->name + " will be written to file \"" + out_name_dx + "\"");
    }
  }

  if (colvar_array_size == 0) {
    // update indices for scalar values
    size_t i;
    for (i = 0; i < num_variables(); i++) {
      bin[i] = grid->current_bin_scalar(i);
    }

    if (grid->index_ok(bin)) {
      grid->acc_value(bin, 1.0);
    }
  } else {
    // update indices for vector/array values
    size_t iv, i;
    for (iv = 0; iv < colvar_array_size; iv++) {
      for (i = 0; i < num_variables(); i++) {
        bin[i] = grid->current_bin_scalar(i, iv);
      }

      if (grid->index_ok(bin)) {
        grid->acc_value(bin, weights[iv]);
      }
    }
  }

  if (output_freq && (cvm::step_absolute() % output_freq) == 0) {
    write_output_files();
  }

  error_code |= cvm::get_error();
  return error_code;
}


int colvarbias_histogram::write_output_files()
{
  if (!has_data) {
    // nothing to write
    return COLVARS_OK;
  }

  if (out_name.size()) {
    cvm::log("Writing the histogram file \""+out_name+"\".\n");
    cvm::backup_file(out_name.c_str());
    std::ostream *grid_os = cvm::proxy->output_stream(out_name);
    if (!grid_os) {
      return cvm::error("Error opening histogram file "+out_name+
                        " for writing.\n", FILE_ERROR);
    }
    grid->write_multicol(*grid_os);
    cvm::proxy->close_output_stream(out_name);
  }

  if (out_name_dx.size()) {
    cvm::log("Writing the histogram file \""+out_name_dx+"\".\n");
    cvm::backup_file(out_name_dx.c_str());
    std::ostream *grid_os = cvm::proxy->output_stream(out_name_dx);
    if (!grid_os) {
      return cvm::error("Error opening histogram file "+out_name_dx+
                        " for writing.\n", FILE_ERROR);
    }
    grid->write_opendx(*grid_os);
    cvm::proxy->close_output_stream(out_name_dx);
  }

  return COLVARS_OK;
}


std::istream & colvarbias_histogram::read_state_data(std::istream& is)
{
  if (! read_state_data_key(is, "grid")) {
    return is;
  }
  if (! grid->read_raw(is)) {
    return is;
  }

  return is;
}


std::ostream & colvarbias_histogram::write_state_data(std::ostream& os)
{
  std::ios::fmtflags flags(os.flags());
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  os << "grid\n";
  grid->write_raw(os, 8);
  os.flags(flags);
  return os;
}
