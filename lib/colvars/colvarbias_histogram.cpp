/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarbias_histogram.h"

/// Histogram "bias" constructor

colvarbias_histogram::colvarbias_histogram(std::string const &conf, char const *key)
  : colvarbias(conf, key),
    grid(NULL), out_name("")
{
  get_keyval(conf, "outputFile", out_name, std::string(""));
  get_keyval(conf, "outputFileDX", out_name_dx, std::string(""));
  get_keyval(conf, "outputFreq", output_freq, cvm::restart_out_freq);

  /// with VMD, this may not be an error
  // if ( output_freq == 0 ) {
  //   cvm::error("User required histogram with zero output frequency");
  // }

  colvar_array_size = 0;
  {
    size_t i;
    bool colvar_array = false;
    get_keyval(conf, "gatherVectorColvars", colvar_array, colvar_array);

    if (colvar_array) {
      for (i = 0; i < colvars.size(); i++) { // should be all vector
        if (colvars[i]->value().type() != colvarvalue::type_vector) {
          cvm::error("Error: used gatherVectorColvars with non-vector colvar.\n", INPUT_ERROR);
          return;
        }
        if (i == 0) {
          colvar_array_size = colvars[i]->value().size();
          if (colvar_array_size < 1) {
            cvm::error("Error: vector variable has dimension less than one.\n", INPUT_ERROR);
            return;
          }
        } else {
          if (colvar_array_size != colvars[i]->value().size()) {
            cvm::error("Error: trying to combine vector colvars of different lengths.\n", INPUT_ERROR);
            return;
          }
        }
      }
    } else {
      for (i = 0; i < colvars.size(); i++) { // should be all scalar
        if (colvars[i]->value().type() != colvarvalue::type_scalar) {
          cvm::error("Error: only scalar colvars are supported when gatherVectorColvars is off.\n", INPUT_ERROR);
          return;
        }
      }
    }
  }

  if (colvar_array_size > 0) {
    weights.assign(colvar_array_size, 1.0);
    get_keyval(conf, "weights", weights, weights, colvarparse::parse_silent);
  }

  grid = new colvar_grid_scalar();

  {
    std::string grid_conf;
    if (key_lookup(conf, "grid", grid_conf)) {
      grid->parse_params(grid_conf);
    } else {
      grid->init_from_colvars(colvars);
    }
  }

  cvm::log("Finished histogram setup.\n");
}

/// Destructor
colvarbias_histogram::~colvarbias_histogram()
{
  if (grid) {
    delete grid;
    grid = NULL;
  }

  if (cvm::n_histo_biases > 0)
    cvm::n_histo_biases -= 1;
}

/// Update the grid
cvm::real colvarbias_histogram::update()
{
  // update base class
  colvarbias::update();

  if (cvm::debug()) {
    cvm::log("Updating histogram bias " + this->name);
  }

  // assign a valid bin size
  bin.assign(colvars.size(), 0);

  if (out_name.size() == 0) {
    // At the first timestep, we need to assign out_name since
    // output_prefix is unset during the constructor
    if (cvm::step_relative() == 0) {
      out_name = cvm::output_prefix + "." + this->name + ".dat";
      cvm::log("Histogram " + this->name + " will be written to file \"" + out_name + "\"");
    }
  }

  if (out_name_dx.size() == 0) {
    if (cvm::step_relative() == 0) {
      out_name_dx = cvm::output_prefix + "." + this->name + ".dx";
      cvm::log("Histogram " + this->name + " will be written to file \"" + out_name_dx + "\"");
    }
  }

  if (colvar_array_size == 0) {
    // update indices for scalar values
    size_t i;
    for (i = 0; i < colvars.size(); i++) {
      bin[i] = grid->value_to_bin_scalar(colvars[i]->value(), i);
    }

    if (grid->index_ok(bin)) {
      grid->acc_value(bin, 1.0);
    }
  } else {
    // update indices for vector/array values
    size_t iv, i;
    for (iv = 0; iv < colvar_array_size; iv++) {
      for (i = 0; i < colvars.size(); i++) {
        bin[i] = grid->value_to_bin_scalar(colvars[i]->value().vector1d_value[iv], i);
      }

      if (grid->index_ok(bin)) {
        grid->acc_value(bin, weights[iv]);
      }
    }
  }

  if (output_freq && (cvm::step_absolute() % output_freq) == 0) {
    write_output_files();
  }

  return 0.0; // no bias energy for histogram
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
    cvm::ofstream grid_os(out_name.c_str());
    if (!grid_os.is_open()) {
      cvm::error("Error opening histogram file " + out_name + " for writing.\n", FILE_ERROR);
    }
    // TODO add return code here
    grid->write_multicol(grid_os);
    grid_os.close();
  }

  if (out_name_dx.size()) {
    cvm::log("Writing the histogram file \""+out_name_dx+"\".\n");
    cvm::backup_file(out_name_dx.c_str());
    cvm::ofstream grid_os(out_name_dx.c_str());
    if (!grid_os.is_open()) {
      cvm::error("Error opening histogram file " + out_name_dx + " for writing.\n", FILE_ERROR);
    }
    // TODO add return code here
    grid->write_opendx(grid_os);
    grid_os.close();
  }
  return COLVARS_OK;
}


std::istream & colvarbias_histogram::read_restart(std::istream& is)
{
  size_t const start_pos = is.tellg();

  cvm::log("Restarting collective variable histogram \""+
            this->name+"\".\n");
  std::string key, brace, conf;

  if ( !(is >> key)   || !(key == "histogram") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block("configuration", conf)) ) {
    cvm::log("Error: in reading restart configuration for histogram \""+
              this->name+"\" at position "+
              cvm::to_str(is.tellg())+" in stream.\n");
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  int id = -1;
  std::string name = "";
  if ( (colvarparse::get_keyval(conf, "name", name, std::string(""), colvarparse::parse_silent)) &&
         (name != this->name) )
    cvm::error("Error: in the restart file, the "
                      "\"histogram\" block has a wrong name: different system?\n");
  if ( (id == -1) && (name == "") ) {
    cvm::error("Error: \"histogram\" block in the restart file "
                      "has no name.\n");
  }

  if ( !(is >> key)   || !(key == "grid")) {
    cvm::error("Error: in reading restart configuration for histogram \""+
              this->name+"\" at position "+
              cvm::to_str(is.tellg())+" in stream.\n");
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }
  if (! grid->read_raw(is)) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  is >> brace;
  if (brace != "}") {
    cvm::error("Error: corrupt restart information for ABF bias \""+
                this->name+"\": no matching brace at position "+
                cvm::to_str(is.tellg())+" in the restart file.\n");
    is.setstate(std::ios::failbit);
  }
  return is;
}

std::ostream & colvarbias_histogram::write_restart(std::ostream& os)
{
  os << "histogram {\n"
     << "  configuration {\n"
     << "    name " << this->name << "\n";
  os << "  }\n";

  os << "grid\n";
  grid->write_raw(os, 8);

  os << "}\n\n";

  return os;
}
