// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// used to set the absolute path of a replica file
#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR ::_chdir
#define GETCWD ::_getcwd
#define PATHSEP "\\"
#else
#include <unistd.h>
#define CHDIR ::chdir
#define GETCWD ::getcwd
#define PATHSEP "/"
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias_meta.h"


colvarbias_meta::colvarbias_meta(char const *key)
  : colvarbias(key), colvarbias_ti(key)
{
  new_hills_begin = hills.end();
  hills_traj_os = NULL;

  ebmeta_equil_steps = 0L;
}


int colvarbias_meta::init(std::string const &conf)
{
  colvarbias::init(conf);
  colvarbias_ti::init(conf);

  enable(f_cvb_calc_pmf);

  get_keyval(conf, "hillWeight", hill_weight, 0.0);
  if (hill_weight > 0.0) {
    enable(f_cvb_apply_force);
  } else {
    cvm::error("Error: hillWeight must be provided, and a positive number.\n", INPUT_ERROR);
  }

  get_keyval(conf, "newHillFrequency", new_hill_freq, 1000);
  if (new_hill_freq > 0) {
    enable(f_cvb_history_dependent);
  }

  get_keyval(conf, "hillWidth", hill_width, cvm::sqrt(2.0 * PI) / 2.0);
  cvm::log("Half-widths of the Gaussian hills (sigma's):\n");
  for (size_t i = 0; i < num_variables(); i++) {
    cvm::log(variables(i)->name+std::string(": ")+
             cvm::to_str(0.5 * variables(i)->width * hill_width));
  }

  {
    bool b_replicas = false;
    get_keyval(conf, "multipleReplicas", b_replicas, false);
    if (b_replicas)
      comm = multiple_replicas;
    else
      comm = single_replica;
  }

  // in all cases, the first replica is this bias itself
  if (replicas.size() == 0) {
    replicas.push_back(this);
  }

  get_keyval(conf, "useGrids", use_grids, true);

  if (use_grids) {
    get_keyval(conf, "gridsUpdateFrequency", grids_freq, new_hill_freq);
    get_keyval(conf, "rebinGrids", rebin_grids, false);

    expand_grids = false;
    size_t i;
    for (i = 0; i < num_variables(); i++) {
      variables(i)->enable(f_cv_grid);
      if (variables(i)->expand_boundaries) {
        expand_grids = true;
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                 ": Will expand grids when the colvar \""+
                 variables(i)->name+"\" approaches its boundaries.\n");
      }
    }

    get_keyval(conf, "keepHills", keep_hills, false);
    if (! get_keyval(conf, "writeFreeEnergyFile", dump_fes, true))
      get_keyval(conf, "dumpFreeEnergyFile", dump_fes, true, colvarparse::parse_silent);
    if (get_keyval(conf, "saveFreeEnergyFile", dump_fes_save, false, colvarparse::parse_silent)) {
      cvm::log("Option \"saveFreeEnergyFile\" is deprecated, "
               "please use \"keepFreeEnergyFiles\" instead.");
    }
    get_keyval(conf, "keepFreeEnergyFiles", dump_fes_save, dump_fes_save);

    hills_energy           = new colvar_grid_scalar(colvars);
    hills_energy_gradients = new colvar_grid_gradient(colvars);
  } else {
    rebin_grids = false;
    keep_hills = false;
    dump_fes = false;
    dump_fes_save = false;
    dump_replica_fes = false;

    hills_energy           = NULL;
    hills_energy_gradients = NULL;
  }

  if (comm != single_replica) {

    if (expand_grids)
      cvm::fatal_error("Error: expandBoundaries is not supported when "
                       "using more than one replicas; please allocate "
                       "wide enough boundaries for each colvar"
                       "ahead of time.\n");

    if (get_keyval(conf, "dumpPartialFreeEnergyFile", dump_replica_fes, false)) {
      if (dump_replica_fes && (! dump_fes)) {
        cvm::log("Enabling \"dumpFreeEnergyFile\".\n");
      }
    }

    get_keyval(conf, "replicaID", replica_id, std::string(""));
    if (!replica_id.size())
      cvm::error("Error: replicaID must be defined "
                 "when using more than one replica.\n", INPUT_ERROR);

    get_keyval(conf, "replicasRegistry",
               replicas_registry_file,
               (this->name+".replicas.registry.txt"));

    get_keyval(conf, "replicaUpdateFrequency",
               replica_update_freq, new_hill_freq);

    if (keep_hills)
      cvm::log("Warning: in metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": keepHills with more than one replica can lead to a very "
               "large amount of input/output and slow down your calculations.  "
               "Please consider disabling it.\n");

  }

  get_keyval(conf, "writeHillsTrajectory", b_hills_traj, false);

  init_well_tempered_params(conf);
  init_ebmeta_params(conf);

  if (cvm::debug())
    cvm::log("Done initializing the metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");

  return COLVARS_OK;
}


int colvarbias_meta::init_well_tempered_params(std::string const &conf)
{
  // for well-tempered metadynamics
  get_keyval(conf, "wellTempered", well_tempered, false);
  get_keyval(conf, "biasTemperature", bias_temperature, -1.0);
  if ((bias_temperature == -1.0) && well_tempered) {
    cvm::fatal_error("Error: biasTemperature is not set.\n");
  }
  if (well_tempered) {
    cvm::log("Well-tempered metadynamics is used.\n");
    cvm::log("The bias temperature is "+cvm::to_str(bias_temperature)+".\n");
  }
  return COLVARS_OK;
}


int colvarbias_meta::init_ebmeta_params(std::string const &conf)
{
  // for ebmeta
  target_dist = NULL;
  get_keyval(conf, "ebMeta", ebmeta, false);
  if(ebmeta){
    if (use_grids && expand_grids) {
      cvm::fatal_error("Error: expandBoundaries is not supported with "
                       "ebMeta please allocate wide enough boundaries for "
                       "each colvar ahead of time and set targetdistfile "
                       "accordingly. \n");
    }
    target_dist = new colvar_grid_scalar();
    target_dist->init_from_colvars(colvars);
    std::string target_dist_file;
    get_keyval(conf, "targetDistFile", target_dist_file);
    std::ifstream targetdiststream(target_dist_file.c_str());
    target_dist->read_multicol(targetdiststream);
    cvm::real min_val = target_dist->minimum_value();
    cvm::real max_val = target_dist->maximum_value();
    if(min_val<0){
      cvm::error("Error: Target distribution of EBMetaD "
                 "has negative values!.\n", INPUT_ERROR);
    }
    cvm::real target_dist_min_val;
    get_keyval(conf, "targetDistMinVal", target_dist_min_val, 1/1000000.0);
    if(target_dist_min_val>0 && target_dist_min_val<1){
      target_dist_min_val=max_val*target_dist_min_val;
      target_dist->remove_small_values(target_dist_min_val);
    } else {
      if (target_dist_min_val==0) {
        cvm::log("NOTE: targetDistMinVal is set to zero, the minimum value of the target \n");
        cvm::log(" distribution will be set as the minimum positive value.\n");
        cvm::real min_pos_val = target_dist->minimum_pos_value();
        if(min_pos_val<=0){
          cvm::error("Error: Target distribution of EBMetaD has negative "
                     "or zero minimum positive value!.\n", INPUT_ERROR);
        }
        if(min_val==0){
          cvm::log("WARNING: Target distribution has zero values.\n");
          cvm::log("Zeros will be converted to the minimum positive value.\n");
          target_dist->remove_small_values(min_pos_val);
        }
      } else {
          cvm::error("Error: targetDistMinVal must be a value between 0 and 1!.\n", INPUT_ERROR);
      }
    }
    // normalize target distribution and multiply by effective volume = exp(differential entropy)
    target_dist->multiply_constant(1.0/target_dist->integral());
    cvm::real volume = cvm::exp(target_dist->entropy());
    target_dist->multiply_constant(volume);
    get_keyval(conf, "ebMetaEquilSteps", ebmeta_equil_steps, ebmeta_equil_steps);
  }

  return COLVARS_OK;
}


colvarbias_meta::~colvarbias_meta()
{
  colvarbias_meta::clear_state_data();
  colvarproxy *proxy = cvm::proxy;

  if (proxy->get_output_stream(replica_hills_file)) {
    proxy->close_output_stream(replica_hills_file);
  }

  if (hills_traj_os) {
    proxy->close_output_stream(hills_traj_file_name());
    hills_traj_os = NULL;
  }

  if (target_dist) {
    delete target_dist;
    target_dist = NULL;
  }
}


int colvarbias_meta::clear_state_data()
{
  if (hills_energy) {
    delete hills_energy;
    hills_energy = NULL;
  }

  if (hills_energy_gradients) {
    delete hills_energy_gradients;
    hills_energy_gradients = NULL;
  }

  hills.clear();
  hills_off_grid.clear();

  return COLVARS_OK;
}


// **********************************************************************
// Hill management member functions
// **********************************************************************

std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::create_hill(colvarbias_meta::hill const &h)
{
  hill_iter const hills_end = hills.end();
  hills.push_back(h);
  if (new_hills_begin == hills_end) {
    // if new_hills_begin is unset, set it for the first time
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {

    // also add it to the list of hills that are off-grid, which may
    // need to be computed analytically when the colvar returns
    // off-grid
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries(h.centers, true);
    if (min_dist < (3.0 * cvm::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(h);
    }
  }

  // output to trajectory (if specified)
  if (hills_traj_os) {
    *hills_traj_os << (hills.back()).output_traj();
    cvm::proxy->flush_output_stream(hills_traj_os);
  }

  has_data = true;
  return hills.end();
}


std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::delete_hill(hill_iter &h)
{
  if (cvm::debug()) {
    cvm::log("Deleting hill from the metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ", with step number "+
             cvm::to_str(h->it)+(h->replica.size() ?
                                 ", replica id \""+h->replica :
                                 "")+".\n");
  }

  if (use_grids && !hills_off_grid.empty()) {
    for (hill_iter hoff = hills_off_grid.begin();
         hoff != hills_off_grid.end(); hoff++) {
      if (*h == *hoff) {
        hills_off_grid.erase(hoff);
        break;
      }
    }
  }

  if (hills_traj_os) {
    // output to the trajectory
    *hills_traj_os << "# DELETED this hill: "
                   << (hills.back()).output_traj()
                   << "\n";
    cvm::proxy->flush_output_stream(hills_traj_os);
  }

  return hills.erase(h);
}


int colvarbias_meta::update()
{
  int error_code = COLVARS_OK;

  // update base class
  error_code |= colvarbias::update();

  // update the TI estimator (if defined)
  error_code |= colvarbias_ti::update();

  // update grid definition, if needed
  error_code |= update_grid_params();
  // add new biasing energy/forces
  error_code |= update_bias();
  // update grid content to reflect new bias
  error_code |= update_grid_data();

  if (comm != single_replica &&
      (cvm::step_absolute() % replica_update_freq) == 0) {
    // sync with the other replicas (if needed)
    error_code |= replica_share();
  }

  error_code |= calc_energy();
  error_code |= calc_forces();

  return error_code;
}


int colvarbias_meta::update_grid_params()
{
  if (use_grids) {

    std::vector<int> curr_bin = hills_energy->get_colvars_index();
    if (cvm::debug()) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": current coordinates on the grid: "+
               cvm::to_str(curr_bin)+".\n");
    }

    if (expand_grids) {
      // first of all, expand the grids, if specified
      bool changed_grids = false;
      int const min_buffer =
        (3 * (size_t) cvm::floor(hill_width)) + 1;

      std::vector<int>         new_sizes(hills_energy->sizes());
      std::vector<colvarvalue> new_lower_boundaries(hills_energy->lower_boundaries);
      std::vector<colvarvalue> new_upper_boundaries(hills_energy->upper_boundaries);

      for (size_t i = 0; i < num_variables(); i++) {

        if (! variables(i)->expand_boundaries)
          continue;

        cvm::real &new_lb   = new_lower_boundaries[i].real_value;
        cvm::real &new_ub   = new_upper_boundaries[i].real_value;
        int       &new_size = new_sizes[i];
        bool changed_lb = false, changed_ub = false;

        if (!variables(i)->hard_lower_boundary)
          if (curr_bin[i] < min_buffer) {
            int const extra_points = (min_buffer - curr_bin[i]);
            new_lb -= extra_points * variables(i)->width;
            new_size += extra_points;
            // changed offset in this direction => the pointer needs to
            // be changed, too
            curr_bin[i] += extra_points;

            changed_lb = true;
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     ": new lower boundary for colvar \""+
                     variables(i)->name+"\", at "+
                     cvm::to_str(new_lower_boundaries[i])+".\n");
          }

        if (!variables(i)->hard_upper_boundary)
          if (curr_bin[i] > new_size - min_buffer - 1) {
            int const extra_points = (curr_bin[i] - (new_size - 1) + min_buffer);
            new_ub += extra_points * variables(i)->width;
            new_size += extra_points;

            changed_ub = true;
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     ": new upper boundary for colvar \""+
                     variables(i)->name+"\", at "+
                     cvm::to_str(new_upper_boundaries[i])+".\n");
          }

        if (changed_lb || changed_ub)
          changed_grids = true;
      }

      if (changed_grids) {

        // map everything into new grids

        colvar_grid_scalar *new_hills_energy =
          new colvar_grid_scalar(*hills_energy);
        colvar_grid_gradient *new_hills_energy_gradients =
          new colvar_grid_gradient(*hills_energy_gradients);

        // supply new boundaries to the new grids

        new_hills_energy->lower_boundaries = new_lower_boundaries;
        new_hills_energy->upper_boundaries = new_upper_boundaries;
        new_hills_energy->setup(new_sizes, 0.0, 1);

        new_hills_energy_gradients->lower_boundaries = new_lower_boundaries;
        new_hills_energy_gradients->upper_boundaries = new_upper_boundaries;
        new_hills_energy_gradients->setup(new_sizes, 0.0, num_variables());

        new_hills_energy->map_grid(*hills_energy);
        new_hills_energy_gradients->map_grid(*hills_energy_gradients);

        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy = new_hills_energy;
        hills_energy_gradients = new_hills_energy_gradients;

        curr_bin = hills_energy->get_colvars_index();
        if (cvm::debug())
          cvm::log("Coordinates on the new grid: "+
                   cvm::to_str(curr_bin)+".\n");
      }
    }
  }
  return COLVARS_OK;
}


int colvarbias_meta::update_bias()
{
  // add a new hill if the required time interval has passed
  if ((cvm::step_absolute() % new_hill_freq) == 0 &&
      is_enabled(f_cvb_history_dependent)) {

    if (cvm::debug()) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
               ": adding a new hill at step "+cvm::to_str(cvm::step_absolute())+".\n");
    }

    cvm::real hills_scale=1.0;

    if (ebmeta) {
      hills_scale *= 1.0/target_dist->value(target_dist->get_colvars_index());
      if(cvm::step_absolute() <= ebmeta_equil_steps) {
        cvm::real const hills_lambda =
          (cvm::real(ebmeta_equil_steps - cvm::step_absolute())) /
          (cvm::real(ebmeta_equil_steps));
        hills_scale = hills_lambda + (1-hills_lambda)*hills_scale;
      }
    }

    if (well_tempered) {
      cvm::real hills_energy_sum_here = 0.0;
      if (use_grids) {
        std::vector<int> curr_bin = hills_energy->get_colvars_index();
        hills_energy_sum_here = hills_energy->value(curr_bin);
      } else {
        calc_hills(new_hills_begin, hills.end(), hills_energy_sum_here);
      }
      hills_scale *= cvm::exp(-1.0*hills_energy_sum_here/(bias_temperature*cvm::boltzmann()));
    }

    switch (comm) {

    case single_replica:

      create_hill(hill(hill_weight*hills_scale, colvars, hill_width));

      break;

    case multiple_replicas:
      create_hill(hill(hill_weight*hills_scale, colvars, hill_width, replica_id));
      std::ostream *replica_hills_os =
        cvm::proxy->get_output_stream(replica_hills_file);
      if (replica_hills_os) {
        *replica_hills_os << hills.back();
      } else {
        return cvm::error("Error: in metadynamics bias \""+this->name+"\""+
                          ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                          " while writing hills for the other replicas.\n", FILE_ERROR);
      }
      break;
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::update_grid_data()
{
  if ((cvm::step_absolute() % grids_freq) == 0) {
    // map the most recent gaussians to the grids
    project_hills(new_hills_begin, hills.end(),
                  hills_energy,    hills_energy_gradients);
    new_hills_begin = hills.end();

    // TODO: we may want to condense all into one replicas array,
    // including "this" as the first element
    if (comm == multiple_replicas) {
      for (size_t ir = 0; ir < replicas.size(); ir++) {
        replicas[ir]->project_hills(replicas[ir]->new_hills_begin,
                                    replicas[ir]->hills.end(),
                                    replicas[ir]->hills_energy,
                                    replicas[ir]->hills_energy_gradients);
        replicas[ir]->new_hills_begin = replicas[ir]->hills.end();
      }
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::calc_energy(std::vector<colvarvalue> const &values)
{
  size_t ir = 0;

  for (ir = 0; ir < replicas.size(); ir++) {
    replicas[ir]->bias_energy = 0.0;
  }

  std::vector<int> const curr_bin = values.size() ?
    hills_energy->get_colvars_index(values) :
    hills_energy->get_colvars_index();

  if (hills_energy->index_ok(curr_bin)) {
    // index is within the grid: get the energy from there
    for (ir = 0; ir < replicas.size(); ir++) {

      bias_energy += replicas[ir]->hills_energy->value(curr_bin);
      if (cvm::debug()) {
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                 ": current coordinates on the grid: "+
                 cvm::to_str(curr_bin)+".\n");
        cvm::log("Grid energy = "+cvm::to_str(bias_energy)+".\n");
      }
    }
  } else {
    // off the grid: compute analytically only the hills at the grid's edges
    for (ir = 0; ir < replicas.size(); ir++) {
      calc_hills(replicas[ir]->hills_off_grid.begin(),
                 replicas[ir]->hills_off_grid.end(),
                 bias_energy,
                 values);
    }
  }

  // now include the hills that have not been binned yet (starting
  // from new_hills_begin)

  for (ir = 0; ir < replicas.size(); ir++) {
    calc_hills(replicas[ir]->new_hills_begin,
               replicas[ir]->hills.end(),
               bias_energy);
    if (cvm::debug()) {
      cvm::log("Hills energy = "+cvm::to_str(bias_energy)+".\n");
    }
  }

  return COLVARS_OK;
}


int colvarbias_meta::calc_forces(std::vector<colvarvalue> const &values)
{
  size_t ir = 0, ic = 0;
  for (ir = 0; ir < replicas.size(); ir++) {
    for (ic = 0; ic < num_variables(); ic++) {
      replicas[ir]->colvar_forces[ic].reset();
    }
  }

  std::vector<int> const curr_bin = values.size() ?
    hills_energy->get_colvars_index(values) :
    hills_energy->get_colvars_index();

  if (hills_energy->index_ok(curr_bin)) {
    for (ir = 0; ir < replicas.size(); ir++) {
      cvm::real const *f = &(replicas[ir]->hills_energy_gradients->value(curr_bin));
      for (ic = 0; ic < num_variables(); ic++) {
        // the gradients are stored, not the forces
        colvar_forces[ic].real_value += -1.0 * f[ic];
      }
    }
  } else {
    // off the grid: compute analytically only the hills at the grid's edges
    for (ir = 0; ir < replicas.size(); ir++) {
      for (ic = 0; ic < num_variables(); ic++) {
        calc_hills_force(ic,
                         replicas[ir]->hills_off_grid.begin(),
                         replicas[ir]->hills_off_grid.end(),
                         colvar_forces,
                         values);
      }
    }
  }

  // now include the hills that have not been binned yet (starting
  // from new_hills_begin)

  if (cvm::debug()) {
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ": adding the forces from the other replicas.\n");
  }

  for (ir = 0; ir < replicas.size(); ir++) {
    for (ic = 0; ic < num_variables(); ic++) {
      calc_hills_force(ic,
                       replicas[ir]->new_hills_begin,
                       replicas[ir]->hills.end(),
                       colvar_forces,
                       values);
      if (cvm::debug()) {
        cvm::log("Hills forces = "+cvm::to_str(colvar_forces)+".\n");
      }
    }
  }

  return COLVARS_OK;
}



void colvarbias_meta::calc_hills(colvarbias_meta::hill_iter      h_first,
                                 colvarbias_meta::hill_iter      h_last,
                                 cvm::real                      &energy,
                                 std::vector<colvarvalue> const &colvar_values)
{
  size_t i = 0;
  std::vector<colvarvalue> curr_values(num_variables());
  for (i = 0; i < num_variables(); i++) {
    curr_values[i].type(variables(i)->value());
  }

  if (colvar_values.size()) {
    for (i = 0; i < num_variables(); i++) {
      curr_values[i] = colvar_values[i];
    }
  } else {
    for (i = 0; i < num_variables(); i++) {
      curr_values[i] = variables(i)->value();
    }
  }

  for (hill_iter h = h_first; h != h_last; h++) {

    // compute the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    for (i = 0; i < num_variables(); i++) {
      colvarvalue const &x  = curr_values[i];
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      cv_sqdev += (variables(i)->dist2(x, center)) / (half_width*half_width);
    }

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-05)
      h->value(0.0);
    } else {
      h->value(cvm::exp(-0.5*cv_sqdev));
    }
    energy += h->energy();
  }
}


void colvarbias_meta::calc_hills_force(size_t const &i,
                                       colvarbias_meta::hill_iter      h_first,
                                       colvarbias_meta::hill_iter      h_last,
                                       std::vector<colvarvalue>       &forces,
                                       std::vector<colvarvalue> const &values)
{
  // Retrieve the value of the colvar
  colvarvalue const x(values.size() ? values[i] : variables(i)->value());

  // do the type check only once (all colvarvalues in the hills series
  // were already saved with their types matching those in the
  // colvars)

  hill_iter h;
  switch (variables(i)->value().type()) {

  case colvarvalue::type_scalar:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].real_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (variables(i)->dist2_lgrad(x, center)).real_value );
    }
    break;

  case colvarvalue::type_3vector:
  case colvarvalue::type_unit3vector:
  case colvarvalue::type_unit3vectorderiv:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].rvector_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (variables(i)->dist2_lgrad(x, center)).rvector_value );
    }
    break;

  case colvarvalue::type_quaternion:
  case colvarvalue::type_quaternionderiv:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].quaternion_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (variables(i)->dist2_lgrad(x, center)).quaternion_value );
    }
    break;

  case colvarvalue::type_vector:
    for (h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].vector1d_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (variables(i)->dist2_lgrad(x, center)).vector1d_value );
    }
    break;

  case colvarvalue::type_notset:
  case colvarvalue::type_all:
  default:
    break;
  }
}


// **********************************************************************
// grid management functions
// **********************************************************************

void colvarbias_meta::project_hills(colvarbias_meta::hill_iter  h_first,
                                    colvarbias_meta::hill_iter  h_last,
                                    colvar_grid_scalar         *he,
                                    colvar_grid_gradient       *hg,
                                    bool print_progress)
{
  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
             ": projecting hills.\n");

  // TODO: improve it by looping over a small subgrid instead of the whole grid

  std::vector<colvarvalue> colvar_values(num_variables());
  std::vector<cvm::real> colvar_forces_scalar(num_variables());

  std::vector<int> he_ix = he->new_index();
  std::vector<int> hg_ix = (hg != NULL) ? hg->new_index() : std::vector<int> (0);
  cvm::real hills_energy_here = 0.0;
  std::vector<colvarvalue> hills_forces_here(num_variables(), 0.0);

  size_t count = 0;
  size_t const print_frequency = ((hills.size() >= 1000000) ? 1 : (1000000/(hills.size()+1)));

  if (hg != NULL) {

    // loop over the points of the grid
    for ( ;
          (he->index_ok(he_ix)) && (hg->index_ok(hg_ix));
          count++) {
      size_t i;
      for (i = 0; i < num_variables(); i++) {
        colvar_values[i] = hills_energy->bin_to_value_scalar(he_ix[i], i);
      }

      // loop over the hills and increment the energy grid locally
      hills_energy_here = 0.0;
      calc_hills(h_first, h_last, hills_energy_here, colvar_values);
      he->acc_value(he_ix, hills_energy_here);

      for (i = 0; i < num_variables(); i++) {
        hills_forces_here[i].reset();
        calc_hills_force(i, h_first, h_last, hills_forces_here, colvar_values);
        colvar_forces_scalar[i] = hills_forces_here[i].real_value;
      }
      hg->acc_force(hg_ix, &(colvar_forces_scalar.front()));

      he->incr(he_ix);
      hg->incr(hg_ix);

      if ((count % print_frequency) == 0) {
        if (print_progress) {
          cvm::real const progress = cvm::real(count) / cvm::real(hg->number_of_points());
          std::ostringstream os;
          os.setf(std::ios::fixed, std::ios::floatfield);
          os << std::setw(6) << std::setprecision(2)
             << 100.0 * progress
             << "% done.";
          cvm::log(os.str());
        }
      }
    }

  } else {

    // simpler version, with just the energy

    for ( ; (he->index_ok(he_ix)); ) {

      for (size_t i = 0; i < num_variables(); i++) {
        colvar_values[i] = hills_energy->bin_to_value_scalar(he_ix[i], i);
      }

      hills_energy_here = 0.0;
      calc_hills(h_first, h_last, hills_energy_here, colvar_values);
      he->acc_value(he_ix, hills_energy_here);

      he->incr(he_ix);

      count++;
      if ((count % print_frequency) == 0) {
        if (print_progress) {
          cvm::real const progress = cvm::real(count) / cvm::real(he->number_of_points());
          std::ostringstream os;
          os.setf(std::ios::fixed, std::ios::floatfield);
          os << std::setw(6) << std::setprecision(2)
             << 100.0 * progress
             << "% done.";
          cvm::log(os.str());
        }
      }
    }
  }

  if (print_progress) {
    cvm::log("100.00% done.");
  }

  if (! keep_hills) {
    hills.erase(hills.begin(), hills.end());
  }
}


void colvarbias_meta::recount_hills_off_grid(colvarbias_meta::hill_iter  h_first,
                                             colvarbias_meta::hill_iter  h_last,
                                             colvar_grid_scalar         *he)
{
  hills_off_grid.clear();

  for (hill_iter h = h_first; h != h_last; h++) {
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries(h->centers, true);
    if (min_dist < (3.0 * cvm::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(*h);
    }
  }
}



// **********************************************************************
// multiple replicas functions
// **********************************************************************


int colvarbias_meta::replica_share()
{
  colvarproxy *proxy = cvm::proxy;
  // sync with the other replicas (if needed)
  if (comm == multiple_replicas) {
    // reread the replicas registry
    update_replicas_registry();
    // empty the output buffer
    std::ostream *replica_hills_os =
      proxy->get_output_stream(replica_hills_file);
    if (replica_hills_os) {
      proxy->flush_output_stream(replica_hills_os);
    }
    read_replica_files();
  }
  return COLVARS_OK;
}


void colvarbias_meta::update_replicas_registry()
{
  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\""+
             ": updating the list of replicas, currently containing "+
             cvm::to_str(replicas.size())+" elements.\n");

  {
    // copy the whole file into a string for convenience
    std::string line("");
    std::ifstream reg_file(replicas_registry_file.c_str());
    if (reg_file.is_open()) {
      replicas_registry.clear();
      while (colvarparse::getline_nocomments(reg_file, line))
        replicas_registry.append(line+"\n");
    } else {
      cvm::error("Error: failed to open file \""+replicas_registry_file+
                 "\" for reading.\n", FILE_ERROR);
    }
  }

  // now parse it
  std::istringstream reg_is(replicas_registry);
  if (reg_is.good()) {

    std::string new_replica("");
    std::string new_replica_file("");
    while ((reg_is >> new_replica) && new_replica.size() &&
           (reg_is >> new_replica_file) && new_replica_file.size()) {

      if (new_replica == this->replica_id) {
        // this is the record for this same replica, skip it
        if (cvm::debug())
          cvm::log("Metadynamics bias \""+this->name+"\""+
                   ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                   ": skipping this replica's own record: \""+
                   new_replica+"\", \""+new_replica_file+"\"\n");
        new_replica_file.clear();
        new_replica.clear();
        continue;
      }

      bool already_loaded = false;
      for (size_t ir = 0; ir < replicas.size(); ir++) {
        if (new_replica == (replicas[ir])->replica_id) {
          // this replica was already added
          if (cvm::debug())
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                     ": skipping a replica already loaded, \""+
                     (replicas[ir])->replica_id+"\".\n");
          already_loaded = true;
          break;
        }
      }

      if (!already_loaded) {
        // add this replica to the registry
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": accessing replica \""+new_replica+"\".\n");
        replicas.push_back(new colvarbias_meta("metadynamics"));
        (replicas.back())->replica_id = new_replica;
        (replicas.back())->replica_list_file = new_replica_file;
        (replicas.back())->replica_state_file = "";
        (replicas.back())->replica_state_file_in_sync = false;

        // Note: the following could become a copy constructor?
        (replicas.back())->name = this->name;
        (replicas.back())->colvars = colvars;
        (replicas.back())->use_grids = use_grids;
        (replicas.back())->dump_fes = false;
        (replicas.back())->expand_grids = false;
        (replicas.back())->rebin_grids = false;
        (replicas.back())->keep_hills = false;
        (replicas.back())->colvar_forces = colvar_forces;

        (replicas.back())->comm = multiple_replicas;

        if (use_grids) {
          (replicas.back())->hills_energy           = new colvar_grid_scalar(colvars);
          (replicas.back())->hills_energy_gradients = new colvar_grid_gradient(colvars);
        }
        if (is_enabled(f_cvb_calc_ti_samples)) {
          (replicas.back())->enable(f_cvb_calc_ti_samples);
          (replicas.back())->colvarbias_ti::init_grids();
        }
      }
    }
  } else {
    cvm::fatal_error("Error: cannot read the replicas registry file \""+
                     replicas_registry+"\".\n");
  }

  // now (re)read the list file of each replica
  for (size_t ir = 0; ir < replicas.size(); ir++) {
    if (cvm::debug())
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ": reading the list file for replica \""+
               (replicas[ir])->replica_id+"\".\n");

    std::ifstream list_is((replicas[ir])->replica_list_file.c_str());
    std::string key;
    std::string new_state_file, new_hills_file;
    if (!(list_is >> key) ||
        !(key == std::string("stateFile")) ||
        !(list_is >> new_state_file) ||
        !(list_is >> key) ||
        !(key == std::string("hillsFile")) ||
        !(list_is >> new_hills_file)) {
      cvm::log("Metadynamics bias \""+this->name+"\""+
               ": failed to read the file \""+
               (replicas[ir])->replica_list_file+"\": will try again after "+
               cvm::to_str(replica_update_freq)+" steps.\n");
      (replicas[ir])->update_status++;
    } else {
      (replicas[ir])->update_status = 0;
      if (new_state_file != (replicas[ir])->replica_state_file) {
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": replica \""+(replicas[ir])->replica_id+
                 "\" has supplied a new state file, \""+new_state_file+
                 "\".\n");
        (replicas[ir])->replica_state_file_in_sync = false;
        (replicas[ir])->replica_state_file = new_state_file;
        (replicas[ir])->replica_hills_file = new_hills_file;
      }
    }
  }

  if (cvm::debug())
    cvm::log("Metadynamics bias \""+this->name+"\": the list of replicas contains "+
             cvm::to_str(replicas.size())+" elements.\n");
}


void colvarbias_meta::read_replica_files()
{
  // Note: we start from the 2nd replica.
  for (size_t ir = 1; ir < replicas.size(); ir++) {

    if (! (replicas[ir])->replica_state_file_in_sync) {
      // if a new state file is being read, the hills file is also new
      (replicas[ir])->replica_hills_file_pos = 0;
    }

    // (re)read the state file if necessary
    if ( (! (replicas[ir])->has_data) ||
         (! (replicas[ir])->replica_state_file_in_sync) ) {

      cvm::log("Metadynamics bias \""+this->name+"\""+
               ": reading the state of replica \""+
               (replicas[ir])->replica_id+"\" from file \""+
               (replicas[ir])->replica_state_file+"\".\n");

      std::ifstream is((replicas[ir])->replica_state_file.c_str());
      if ((replicas[ir])->read_state(is)) {
        // state file has been read successfully
        (replicas[ir])->replica_state_file_in_sync = true;
        (replicas[ir])->update_status = 0;
      }
      is.close();
    }

    // now read the hills added after writing the state file
    if ((replicas[ir])->replica_hills_file.size()) {

      if (cvm::debug())
        cvm::log("Metadynamics bias \""+this->name+"\""+
                 ": checking for new hills from replica \""+
                 (replicas[ir])->replica_id+"\" in the file \""+
                 (replicas[ir])->replica_hills_file+"\".\n");

      // read hills from the other replicas' files; for each file, resume
      // the position recorded previously

      std::ifstream is((replicas[ir])->replica_hills_file.c_str());
      if (is.is_open()) {

        // try to resume the previous position
        is.seekg((replicas[ir])->replica_hills_file_pos, std::ios::beg);
        if (!is.is_open()){
          // if fail (the file may have been overwritten), reset this
          // position
          is.clear();
          is.seekg(0, std::ios::beg);
          // reset the counter
          (replicas[ir])->replica_hills_file_pos = 0;
          // schedule to reread the state file
          (replicas[ir])->replica_state_file_in_sync = false;
          // and record the failure
          (replicas[ir])->update_status++;
          cvm::log("Failed to read the file \""+(replicas[ir])->replica_hills_file+
                   "\" at the previous position: will try again in "+
                   cvm::to_str(replica_update_freq)+" steps.\n");
        } else {

          while ((replicas[ir])->read_hill(is)) {
            //           if (cvm::debug())
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ": received a hill from replica \""+
                     (replicas[ir])->replica_id+
                     "\" at step "+
                     cvm::to_str(((replicas[ir])->hills.back()).it)+".\n");
          }
          is.clear();
          // store the position for the next read
          (replicas[ir])->replica_hills_file_pos = is.tellg();
          if (cvm::debug())
            cvm::log("Metadynamics bias \""+this->name+"\""+
                     ": stopped reading file \""+(replicas[ir])->replica_hills_file+
                     "\" at position "+
                     cvm::to_str((replicas[ir])->replica_hills_file_pos)+".\n");

          // test whether this is the end of the file
          is.seekg(0, std::ios::end);
          if (is.tellg() > (replicas[ir])->replica_hills_file_pos+1) {
            (replicas[ir])->update_status++;
          } else {
            (replicas[ir])->update_status = 0;
          }
        }

      } else {
        cvm::log("Failed to read the file \""+(replicas[ir])->replica_hills_file+
                 "\": will try again in "+
                 cvm::to_str(replica_update_freq)+" steps.\n");
        (replicas[ir])->update_status++;
        // cvm::fatal_error ("Error: cannot read from file \""+
        //                   (replicas[ir])->replica_hills_file+"\".\n");
      }
      is.close();
    }

    size_t const n_flush = (replica_update_freq/new_hill_freq + 1);
    if ((replicas[ir])->update_status > 3*n_flush) {
      // TODO: suspend the calculation?
      cvm::log("WARNING: in metadynamics bias \""+this->name+"\""+
               " failed to read completely the output of replica \""+
               (replicas[ir])->replica_id+
               "\" after more than "+
               cvm::to_str((replicas[ir])->update_status * replica_update_freq)+
               " steps.  Ensure that it is still running.\n");
    }
  }
}


int colvarbias_meta::set_state_params(std::string const &state_conf)
{
  std::string new_replica = "";
  if (colvarparse::get_keyval(state_conf, "replicaID", new_replica,
                              std::string(""), colvarparse::parse_silent) &&
      (new_replica != this->replica_id)) {
    cvm::error("Error: in the state file, the "
               "\"metadynamics\" block has a different replicaID: different system?\n",
               INPUT_ERROR);
    return INPUT_ERROR;
  }

  return COLVARS_OK;
}


std::istream & colvarbias_meta::read_state_data(std::istream& is)
{
  bool grids_from_restart_file = use_grids;

  if (use_grids) {

    if (expand_grids) {
      // the boundaries of the colvars may have been changed; TODO:
      // this reallocation is only for backward-compatibility, and may
      // be deleted when grid_parameters (i.e. colvargrid's own
      // internal reallocation) has kicked in
      delete hills_energy;
      delete hills_energy_gradients;
      hills_energy = new colvar_grid_scalar(colvars);
      hills_energy_gradients = new colvar_grid_gradient(colvars);
    }

    colvar_grid_scalar   *hills_energy_backup = NULL;
    colvar_grid_gradient *hills_energy_gradients_backup = NULL;

    if (has_data) {
      if (cvm::debug())
        cvm::log("Backupping grids for metadynamics bias \""+
                 this->name+"\""+
                 ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");
      hills_energy_backup           = hills_energy;
      hills_energy_gradients_backup = hills_energy_gradients;
      hills_energy                  = new colvar_grid_scalar(colvars);
      hills_energy_gradients        = new colvar_grid_gradient(colvars);
    }

    size_t const hills_energy_pos = is.tellg();
    std::string key;
    if (!(is >> key)) {
      if (hills_energy_backup != NULL) {
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy           = hills_energy_backup;
        hills_energy_gradients = hills_energy_gradients_backup;
      }
      is.clear();
      is.seekg(hills_energy_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    } else if (!(key == std::string("hills_energy")) ||
               !(hills_energy->read_restart(is))) {
      is.clear();
      is.seekg(hills_energy_pos, std::ios::beg);
      grids_from_restart_file = false;
      if (!rebin_grids) {
        if (hills_energy_backup == NULL)
          cvm::fatal_error("Error: couldn't read the free energy grid for metadynamics bias \""+
                           this->name+"\""+
                           ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                           "; if useGrids was off when the state file was written, "
                           "enable rebinGrids now to regenerate the grids.\n");
        else {
          if (comm == single_replica)
            cvm::log("Error: couldn't read the free energy grid for metadynamics bias \""+
                     this->name+"\".\n");
          delete hills_energy;
          delete hills_energy_gradients;
          hills_energy           = hills_energy_backup;
          hills_energy_gradients = hills_energy_gradients_backup;
          is.setstate(std::ios::failbit);
          return is;
        }
      }
    }

    size_t const hills_energy_gradients_pos = is.tellg();
    if (!(is >> key)) {
      if (hills_energy_backup != NULL)  {
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy           = hills_energy_backup;
        hills_energy_gradients = hills_energy_gradients_backup;
      }
      is.clear();
      is.seekg(hills_energy_gradients_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    } else if (!(key == std::string("hills_energy_gradients")) ||
               !(hills_energy_gradients->read_restart(is))) {
      is.clear();
      is.seekg(hills_energy_gradients_pos, std::ios::beg);
      grids_from_restart_file = false;
      if (!rebin_grids) {
        if (hills_energy_backup == NULL)
          cvm::fatal_error("Error: couldn't read the free energy gradients grid for metadynamics bias \""+
                           this->name+"\""+
                           ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                           "; if useGrids was off when the state file was written, "
                           "enable rebinGrids now to regenerate the grids.\n");
        else {
          if (comm == single_replica)
            cvm::log("Error: couldn't read the free energy gradients grid for metadynamics bias \""+
                     this->name+"\".\n");
          delete hills_energy;
          delete hills_energy_gradients;
          hills_energy           = hills_energy_backup;
          hills_energy_gradients = hills_energy_gradients_backup;
          is.setstate(std::ios::failbit);
          return is;
        }
      }
    }

    if (cvm::debug())
      cvm::log("Successfully read new grids for bias \""+
               this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+"\n");

    if (hills_energy_backup != NULL) {
      // now that we have successfully updated the grids, delete the
      // backup copies
      if (cvm::debug())
        cvm::log("Deallocating the older grids.\n");

      delete hills_energy_backup;
      delete hills_energy_gradients_backup;
    }
  }

  bool const existing_hills = !hills.empty();
  size_t const old_hills_size = hills.size();
  hill_iter old_hills_end = hills.end();
  hill_iter old_hills_off_grid_end = hills_off_grid.end();

  // read the hills explicitly written (if there are any)
  while (read_hill(is)) {
    if (cvm::debug())
      cvm::log("Read a previously saved hill under the "
               "metadynamics bias \""+
               this->name+"\", created at step "+
               cvm::to_str((hills.back()).it)+".\n");
  }
  is.clear();
  new_hills_begin = hills.end();
  if (grids_from_restart_file) {
    if (hills.size() > old_hills_size)
      cvm::log("Read "+cvm::to_str(hills.size())+
               " hills in addition to the grids.\n");
  } else {
    if (!hills.empty())
      cvm::log("Read "+cvm::to_str(hills.size())+" hills.\n");
  }

  if (rebin_grids) {

    // allocate new grids (based on the new boundaries and widths just
    // read from the configuration file), and project onto them the
    // grids just read from the restart file

    colvar_grid_scalar   *new_hills_energy =
      new colvar_grid_scalar(colvars);
    colvar_grid_gradient *new_hills_energy_gradients =
      new colvar_grid_gradient(colvars);

    if (!grids_from_restart_file || (keep_hills && !hills.empty())) {
      // if there are hills, recompute the new grids from them
      cvm::log("Rebinning the energy and forces grids from "+
               cvm::to_str(hills.size())+" hills (this may take a while)...\n");
      project_hills(hills.begin(), hills.end(),
                    new_hills_energy, new_hills_energy_gradients, true);
      cvm::log("rebinning done.\n");

    } else {
      // otherwise, use the grids in the restart file
      cvm::log("Rebinning the energy and forces grids "
               "from the grids in the restart file.\n");
      new_hills_energy->map_grid(*hills_energy);
      new_hills_energy_gradients->map_grid(*hills_energy_gradients);
    }

    delete hills_energy;
    delete hills_energy_gradients;
    hills_energy = new_hills_energy;
    hills_energy_gradients = new_hills_energy_gradients;

    // assuming that some boundaries have expanded, eliminate those
    // off-grid hills that aren't necessary any more
    if (!hills.empty())
      recount_hills_off_grid(hills.begin(), hills.end(), hills_energy);
  }

  if (use_grids) {
    if (!hills_off_grid.empty()) {
      cvm::log(cvm::to_str(hills_off_grid.size())+" hills are near the "
               "grid boundaries: they will be computed analytically "
               "and saved to the state files.\n");
    }
  }

  colvarbias_ti::read_state_data(is);

  if (cvm::debug())
    cvm::log("colvarbias_meta::read_restart() done\n");

  if (existing_hills) {
    hills.erase(hills.begin(), old_hills_end);
    hills_off_grid.erase(hills_off_grid.begin(), old_hills_off_grid_end);
  }

  has_data = true;

  if (comm != single_replica) {
    read_replica_files();
  }

  return is;
}


std::istream & colvarbias_meta::read_hill(std::istream &is)
{
  if (!is) return is; // do nothing if failbit is set

  size_t const start_pos = is.tellg();

  std::string data;
  if ( !(is >> read_block("hill", data)) ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  cvm::step_number h_it;
  get_keyval(data, "step", h_it, 0L, parse_silent);
  if (h_it <= state_file_step) {
    if (cvm::debug())
      cvm::log("Skipping a hill older than the state file for metadynamics bias \""+
               this->name+"\""+
               ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+"\n");
    return is;
  }

  cvm::real h_weight;
  get_keyval(data, "weight", h_weight, hill_weight, parse_silent);

  std::vector<colvarvalue> h_centers(num_variables());
  for (size_t i = 0; i < num_variables(); i++) {
    h_centers[i].type(variables(i)->value());
  }
  {
    // it is safer to read colvarvalue objects one at a time;
    // TODO: change this it later
    std::string centers_input;
    key_lookup(data, "centers", &centers_input);
    std::istringstream centers_is(centers_input);
    for (size_t i = 0; i < num_variables(); i++) {
      centers_is >> h_centers[i];
    }
  }

  std::vector<cvm::real> h_widths(num_variables());
  get_keyval(data, "widths", h_widths,
             std::vector<cvm::real>(num_variables(), (cvm::sqrt(2.0 * PI) / 2.0)),
             parse_silent);

  std::string h_replica = "";
  if (comm != single_replica) {
    get_keyval(data, "replicaID", h_replica, replica_id, parse_silent);
    if (h_replica != replica_id)
      cvm::fatal_error("Error: trying to read a hill created by replica \""+h_replica+
                       "\" for replica \""+replica_id+
                       "\"; did you swap output files?\n");
  }

  hill_iter const hills_end = hills.end();
  hills.push_back(hill(h_it, h_weight, h_centers, h_widths, h_replica));
  if (new_hills_begin == hills_end) {
    // if new_hills_begin is unset, set it for the first time
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {
    // add this also to the list of hills that are off-grid, which will
    // be computed analytically
    cvm::real const min_dist =
      hills_energy->bin_distance_from_boundaries((hills.back()).centers, true);
    if (min_dist < (3.0 * cvm::floor(hill_width)) + 1.0) {
      hills_off_grid.push_back(hills.back());
    }
  }

  has_data = true;
  return is;
}


int colvarbias_meta::setup_output()
{
  output_prefix = cvm::output_prefix();
  if (cvm::main()->num_biases_feature(colvardeps::f_cvb_calc_pmf) > 1) {
    // if this is not the only free energy integrator, append
    // this bias's name, to distinguish it from the output of the other
    // biases producing a .pmf file
    output_prefix += ("."+this->name);
  }

  if (comm == multiple_replicas) {

    // TODO: one may want to specify the path manually for intricated filesystems?
    char *pwd = new char[3001];
    if (GETCWD(pwd, 3000) == NULL)
      cvm::fatal_error("Error: cannot get the path of the current working directory.\n");
    replica_list_file =
      (std::string(pwd)+std::string(PATHSEP)+
       this->name+"."+replica_id+".files.txt");
    // replica_hills_file and replica_state_file are those written
    // by the current replica; within the mirror biases, they are
    // those by another replica
    replica_hills_file =
      (std::string(pwd)+std::string(PATHSEP)+
       cvm::output_prefix()+".colvars."+this->name+"."+replica_id+".hills");
    replica_state_file =
      (std::string(pwd)+std::string(PATHSEP)+
       cvm::output_prefix()+".colvars."+this->name+"."+replica_id+".state");
    delete[] pwd;

    // now register this replica

    // first check that it isn't already there
    bool registered_replica = false;
    std::ifstream reg_is(replicas_registry_file.c_str());
    if (reg_is.is_open()) {  // the file may not be there yet
      std::string existing_replica("");
      std::string existing_replica_file("");
      while ((reg_is >> existing_replica) && existing_replica.size() &&
             (reg_is >> existing_replica_file) && existing_replica_file.size()) {
        if (existing_replica == replica_id) {
          // this replica was already registered
          replica_list_file = existing_replica_file;
          reg_is.close();
          registered_replica = true;
          break;
        }
      }
      reg_is.close();
    }

    // if this replica was not included yet, we should generate a
    // new record for it: but first, we write this replica's files,
    // for the others to read

    // open the "hills" buffer file
    reopen_replica_buffer_file();

    // write the state file (so that there is always one available)
    write_replica_state_file();

    // schedule to read the state files of the other replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }

    // if we're running without grids, use a growing list of "hills" files
    // otherwise, just one state file and one "hills" file as buffer
    std::ostream *list_os =
      cvm::proxy->output_stream(replica_list_file,
                                (use_grids ? std::ios_base::trunc :
                                 std::ios_base::app));
    if (!list_os) {
      return cvm::get_error();
    }
    *list_os << "stateFile " << replica_state_file << "\n";
    *list_os << "hillsFile " << replica_hills_file << "\n";
    cvm::proxy->close_output_stream(replica_list_file);

    // finally, add a new record for this replica to the registry
    if (! registered_replica) {
      std::ostream *reg_os =
        cvm::proxy->output_stream(replicas_registry_file,
                                  std::ios::app);
      if (!reg_os) {
        return cvm::get_error();
      }
      *reg_os << replica_id << " " << replica_list_file << "\n";
      cvm::proxy->close_output_stream(replicas_registry_file);
    }
  }

  if (b_hills_traj) {
    if (!hills_traj_os) {
      hills_traj_os = cvm::proxy->output_stream(hills_traj_file_name());
      if (!hills_traj_os) return cvm::get_error();
    }
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::string const colvarbias_meta::hills_traj_file_name() const
{
  return std::string(cvm::output_prefix()+
                     ".colvars."+this->name+
                     ( (comm != single_replica) ?
                       ("."+replica_id) :
                       ("") )+
                     ".hills.traj");
}


std::string const colvarbias_meta::get_state_params() const
{
  std::ostringstream os;
  if (this->comm != single_replica)
    os << "replicaID " << this->replica_id << "\n";
  return (colvarbias::get_state_params() + os.str());
}


std::ostream & colvarbias_meta::write_state_data(std::ostream& os)
{
  if (use_grids) {

    // this is a very good time to project hills, if you haven't done
    // it already!
    project_hills(new_hills_begin, hills.end(),
                  hills_energy,    hills_energy_gradients);
    new_hills_begin = hills.end();

    // write down the grids to the restart file
    os << "  hills_energy\n";
    hills_energy->write_restart(os);
    os << "  hills_energy_gradients\n";
    hills_energy_gradients->write_restart(os);
  }

  if ( (!use_grids) || keep_hills ) {
    // write all hills currently in memory
    for (std::list<hill>::const_iterator h = this->hills.begin();
         h != this->hills.end();
         h++) {
      os << *h;
    }
  } else {
    // write just those that are near the grid boundaries
    for (std::list<hill>::const_iterator h = this->hills_off_grid.begin();
         h != this->hills_off_grid.end();
         h++) {
      os << *h;
    }
  }

  colvarbias_ti::write_state_data(os);
  return os;
}


int colvarbias_meta::write_state_to_replicas()
{
  int error_code = COLVARS_OK;
  if (comm != single_replica) {
    error_code |= write_replica_state_file();
    error_code |= reopen_replica_buffer_file();
    // schedule to reread the state files of the other replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }
  }
  return error_code;
}


int colvarbias_meta::write_output_files()
{
  colvarbias_ti::write_output_files();
  if (dump_fes) {
    write_pmf();
  }
  return COLVARS_OK;
}


void colvarbias_meta::write_pmf()
{
  // allocate a new grid to store the pmf
  colvar_grid_scalar *pmf = new colvar_grid_scalar(*hills_energy);
  pmf->setup();

  if ((comm == single_replica) || (dump_replica_fes)) {
    // output the PMF from this instance or replica
    pmf->reset();
    pmf->add_grid(*hills_energy);

    if (ebmeta) {
      int nt_points=pmf->number_of_points();
      for (int i = 0; i < nt_points; i++) {
         cvm:: real pmf_val=0.0;
         cvm:: real target_val=target_dist->value(i);
         if (target_val>0) {
           pmf_val=pmf->value(i);
           pmf_val=pmf_val+cvm::temperature() * cvm::boltzmann() * std::log(target_val);
         }
         pmf->set_value(i,pmf_val);
      }
    }

    cvm::real const max = pmf->maximum_value();
    pmf->add_constant(-1.0 * max);
    pmf->multiply_constant(-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + cvm::temperature()) / bias_temperature;
      pmf->multiply_constant(well_temper_scale);
    }
    {
      std::string const fes_file_name(this->output_prefix +
                                      ((comm != single_replica) ? ".partial" : "") +
                                      (dump_fes_save ?
                                       "."+cvm::to_str(cvm::step_absolute()) : "") +
                                      ".pmf");
      cvm::proxy->backup_file(fes_file_name);
      std::ostream *fes_os = cvm::proxy->output_stream(fes_file_name);
      pmf->write_multicol(*fes_os);
      cvm::proxy->close_output_stream(fes_file_name);
    }
  }

  if (comm != single_replica) {
    // output the combined PMF from all replicas
    pmf->reset();
    // current replica already included in the pools of replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      pmf->add_grid(*(replicas[ir]->hills_energy));
    }

    if (ebmeta) {
      int nt_points=pmf->number_of_points();
      for (int i = 0; i < nt_points; i++) {
         cvm:: real pmf_val=0.0;
         cvm:: real target_val=target_dist->value(i);
         if (target_val>0) {
           pmf_val=pmf->value(i);
           pmf_val=pmf_val+cvm::temperature() * cvm::boltzmann() * std::log(target_val);
         }
         pmf->set_value(i,pmf_val);
      }
    }

    cvm::real const max = pmf->maximum_value();
    pmf->add_constant(-1.0 * max);
    pmf->multiply_constant(-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + cvm::temperature()) / bias_temperature;
      pmf->multiply_constant(well_temper_scale);
    }
    std::string const fes_file_name(this->output_prefix +
                                    (dump_fes_save ?
                                     "."+cvm::to_str(cvm::step_absolute()) : "") +
                                    ".pmf");
    cvm::proxy->backup_file(fes_file_name);
    std::ostream *fes_os = cvm::proxy->output_stream(fes_file_name);
    pmf->write_multicol(*fes_os);
    cvm::proxy->close_output_stream(fes_file_name);
  }

  delete pmf;
}



int colvarbias_meta::write_replica_state_file()
{
  colvarproxy *proxy = cvm::proxy;

  if (cvm::debug()) {
    cvm::log("Writing replica state file for bias \""+name+"\"\n");
  }

  int error_code = COLVARS_OK;

  // Write to temporary state file
  std::string const tmp_state_file(replica_state_file+".tmp");
  error_code |= proxy->remove_file(tmp_state_file);
  std::ostream *rep_state_os = cvm::proxy->output_stream(tmp_state_file);
  if (rep_state_os) {
    if (!write_state(*rep_state_os)) {
      error_code |= cvm::error("Error: in writing to temporary file \""+
                               tmp_state_file+"\".\n", FILE_ERROR);
    }
  }
  error_code |= proxy->close_output_stream(tmp_state_file);

  error_code |= proxy->rename_file(tmp_state_file, replica_state_file);

  return error_code;
}


int colvarbias_meta::reopen_replica_buffer_file()
{
  int error_code = COLVARS_OK;
  colvarproxy *proxy = cvm::proxy;
  if (proxy->get_output_stream(replica_hills_file) != NULL) {
    error_code |= proxy->close_output_stream(replica_hills_file);
  }
  error_code |= proxy->remove_file(replica_hills_file);
  std::ostream *replica_hills_os = proxy->output_stream(replica_hills_file);
  if (replica_hills_os) {
    replica_hills_os->setf(std::ios::scientific, std::ios::floatfield);
  } else {
    error_code |= FILE_ERROR;
  }
  return error_code;
}


std::string colvarbias_meta::hill::output_traj()
{
  std::ostringstream os;
  os.setf(std::ios::fixed, std::ios::floatfield);
  os << std::setw(cvm::it_width) << it << " ";

  os.setf(std::ios::scientific, std::ios::floatfield);

  size_t i;
  os << "  ";
  for (i = 0; i < centers.size(); i++) {
    os << " ";
    os << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)  << centers[i];
  }

  os << "  ";
  for (i = 0; i < widths.size(); i++) {
    os << " ";
    os << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width) << widths[i];
  }

  os << "  ";
  os << std::setprecision(cvm::en_prec)
     << std::setw(cvm::en_width) << W << "\n";

  return os.str();
}


std::ostream & operator << (std::ostream &os, colvarbias_meta::hill const &h)
{
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << "hill {\n";
  os << "  step " << std::setw(cvm::it_width) << h.it << "\n";
  os << "  weight   "
     << std::setprecision(cvm::en_prec)
     << std::setw(cvm::en_width)
     << h.W << "\n";

  if (h.replica.size())
    os << "  replicaID  " << h.replica << "\n";

  size_t i;
  os << "  centers ";
  for (i = 0; i < (h.centers).size(); i++) {
    os << " "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << h.centers[i];
  }
  os << "\n";

  os << "  widths  ";
  for (i = 0; i < (h.widths).size(); i++) {
    os << " "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << h.widths[i];
  }
  os << "\n";

  os << "}\n";

  return os;
}
