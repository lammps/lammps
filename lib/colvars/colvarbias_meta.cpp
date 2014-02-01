#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// used to set the absolute path of a replica file
#if defined (WIN32) && !defined(__CYGWIN__)
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


#include "colvar.h"
#include "colvarbias_meta.h"


colvarbias_meta::colvarbias_meta()
  : colvarbias(),
    state_file_step (0),
    new_hills_begin (hills.end())
{
}


colvarbias_meta::colvarbias_meta (std::string const &conf, char const *key)
  : colvarbias (conf, key),
    state_file_step (0),
    new_hills_begin (hills.end())
{
  if (cvm::n_abf_biases > 0)
    cvm::log ("Warning: running ABF and metadynamics together is not recommended unless applyBias is off for ABF.\n");

  get_keyval (conf, "hillWeight", hill_weight, 0.01);
  if (hill_weight == 0.0)
    cvm::log ("Warning: hillWeight has been set to zero, "
              "this bias will have no effect.\n");

  get_keyval (conf, "newHillFrequency", new_hill_freq, 1000);

  get_keyval (conf, "hillWidth", hill_width, std::sqrt (2.0 * PI) / 2.0);

  {
    bool b_replicas = false;
    get_keyval (conf, "multipleReplicas", b_replicas, false);
    if (b_replicas)
      comm = multiple_replicas;
    else
      comm = single_replica;
  }

  get_keyval (conf, "useGrids", use_grids, true);

  if (use_grids) {
    get_keyval (conf, "gridsUpdateFrequency", grids_freq, new_hill_freq);
    get_keyval (conf, "rebinGrids", rebin_grids, false);

    expand_grids = false;
    for (size_t i = 0; i < colvars.size(); i++) {
      if (colvars[i]->expand_boundaries) {
        expand_grids = true;
        cvm::log ("Metadynamics bias \""+this->name+"\""+
                  ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                  ": Will expand grids when the colvar \""+
                  colvars[i]->name+"\" approaches its boundaries.\n");
      }
    }

    get_keyval (conf, "keepHills", keep_hills, false);
    if (! get_keyval (conf, "writeFreeEnergyFile", dump_fes, true))
      get_keyval (conf, "dumpFreeEnergyFile", dump_fes, true, colvarparse::parse_silent);
    get_keyval (conf, "saveFreeEnergyFile", dump_fes_save, false);

    for (size_t i = 0; i < colvars.size(); i++) {
      colvars[i]->enable (colvar::task_grid);
    }

    hills_energy           = new colvar_grid_scalar   (colvars);
    hills_energy_gradients = new colvar_grid_gradient (colvars);
  } else {
    rebin_grids = false;
    keep_hills = false;
    dump_fes = false;
    dump_fes_save = false;
    dump_replica_fes = false;
  }

  if (comm != single_replica) {

    if (expand_grids)
      cvm::fatal_error ("Error: expandBoundaries is not supported when "
                        "using more than one replicas; please allocate "
                        "wide enough boundaries for each colvar"
                        "ahead of time.\n");

    if (get_keyval (conf, "dumpPartialFreeEnergyFile", dump_replica_fes, false)) {
      if (dump_replica_fes && (! dump_fes)) {
        cvm::log ("Enabling \"dumpFreeEnergyFile\".\n");
      }
    }

    get_keyval (conf, "replicaID", replica_id, std::string (""));
    if (!replica_id.size())
      cvm::fatal_error ("Error: replicaID must be defined "
                        "when using more than one replica.\n");

    get_keyval (conf, "replicasRegistry",
                replicas_registry_file,
                (this->name+".replicas.registry.txt"));

    get_keyval (conf, "replicaUpdateFrequency",
                replica_update_freq, new_hill_freq);

    if (keep_hills)
      cvm::log ("Warning: in metadynamics bias \""+this->name+"\""+
                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                ": keepHills with more than one replica can lead to a very "
                "large amount input/output and slow down your calculations.  "
                "Please consider disabling it.\n");


    {
      // TODO: one may want to specify the path manually for intricated filesystems?
      char *pwd = new char[3001];
      if (GETCWD (pwd, 3000) == NULL)
        cvm::fatal_error ("Error: cannot get the path of the current working directory.\n");
      replica_list_file =
        (std::string (pwd)+std::string (PATHSEP)+
         this->name+"."+replica_id+".files.txt");
      // replica_hills_file and replica_state_file are those written
      // by the current replica; within the mirror biases, they are
      // those by another replica
      replica_hills_file =
        (std::string (pwd)+std::string (PATHSEP)+
         cvm::output_prefix+".colvars."+this->name+"."+replica_id+".hills");
      replica_state_file =
        (std::string (pwd)+std::string (PATHSEP)+
         cvm::output_prefix+".colvars."+this->name+"."+replica_id+".state");
      delete pwd;
    }

    // now register this replica

    // first check that it isn't already there
    bool registered_replica = false;
    std::ifstream reg_is (replicas_registry_file.c_str());
    if (reg_is.good()) {  // the file may not be there yet
      std::string existing_replica ("");
      std::string existing_replica_file ("");
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
    replica_hills_os.open (replica_hills_file.c_str());
    if (!replica_hills_os.good())
      cvm::fatal_error ("Error: in opening file \""+
                        replica_hills_file+"\" for writing.\n");
    replica_hills_os.setf (std::ios::scientific, std::ios::floatfield);

    // write the state file (so that there is always one available)
    write_replica_state_file();
    // schedule to read the state files of the other replicas
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }

    // if we're running without grids, use a growing list of "hills" files
    // otherwise, just one state file and one "hills" file as buffer
    std::ofstream list_os (replica_list_file.c_str(),
                           (use_grids ? std::ios::trunc : std::ios::app));
    if (! list_os.good())
      cvm::fatal_error ("Error: in opening file \""+
                        replica_list_file+"\" for writing.\n");
    list_os << "stateFile " << replica_state_file << "\n";
    list_os << "hillsFile " << replica_hills_file << "\n";
    list_os.close();

    // finally, if add a new record for this replica to the registry
    if (! registered_replica) {
      std::ofstream reg_os (replicas_registry_file.c_str(), std::ios::app);
      if (! reg_os.good())
        cvm::fatal_error ("Error: in opening file \""+
                          replicas_registry_file+"\" for writing.\n");
      reg_os << replica_id << " " << replica_list_file << "\n";
      reg_os.close();
    }
  }

  get_keyval (conf, "writeHillsTrajectory", b_hills_traj, false);
  if (b_hills_traj) {
    std::string const traj_file_name (cvm::output_prefix+
                                      ".colvars."+this->name+
                                      ( (comm != single_replica) ?
                                        ("."+replica_id) :
                                        ("") )+
                                      ".hills.traj");
    hills_traj_os.open (traj_file_name.c_str());
    if (!hills_traj_os.good())
      cvm::fatal_error ("Error: in opening hills output file \"" +
                        traj_file_name + "\".\n");
  }

  // for well-tempered metadynamics
  get_keyval (conf, "wellTempered", well_tempered, false);
  get_keyval (conf, "biasTemperature", bias_temperature, -1.0);
  if ((bias_temperature == -1.0) && well_tempered) {
    cvm::fatal_error ("Error: biasTemperature is not set.\n");
  }
  if (well_tempered) {
    cvm::log("Well-tempered metadynamics is used.\n");
    cvm::log("The bias temperature is "+cvm::to_str(bias_temperature)+".\n");
  }

  if (cvm::debug())
    cvm::log ("Done initializing the metadynamics bias \""+this->name+"\""+
              ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");

  save_delimiters = false;
}


colvarbias_meta::~colvarbias_meta()
{
  if (hills_energy) {
    delete hills_energy;
    hills_energy = NULL;
  }

  if (hills_energy_gradients) {
    delete hills_energy_gradients;
    hills_energy_gradients = NULL;
  }

  if (replica_hills_os.good())
    replica_hills_os.close();

  if (hills_traj_os.good())
    hills_traj_os.close();

  if (cvm::n_meta_biases > 0)
    cvm::n_meta_biases -= 1;
}



// **********************************************************************
// Hill management member functions
// **********************************************************************

std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::create_hill (colvarbias_meta::hill const &h)
{
  hill_iter const hills_end = hills.end();
  hills.push_back (h);
  if (new_hills_begin == hills_end) {
    // if new_hills_begin is unset, set it for the first time
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {

    // also add it to the list of hills that are off-grid, which may
    // need to be computed analytically when the colvar returns
    // off-grid
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries (h.centers, true);
    if (min_dist < (3.0 * std::floor (hill_width)) + 1.0) {
      hills_off_grid.push_back (h);
    }
  }

  // output to trajectory (if specified)
  if (hills_traj_os.good()) {
    hills_traj_os << (hills.back()).output_traj();
    if (cvm::debug()) {
      hills_traj_os.flush();
    }
  }

  has_data = true;
  return hills.end();
}


std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::delete_hill (hill_iter &h)
{
  if (cvm::debug()) {
    cvm::log ("Deleting hill from the metadynamics bias \""+this->name+"\""+
              ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
              ", with step number "+
              cvm::to_str (h->it)+(h->replica.size() ?
                                   ", replica id \""+h->replica :
                                   "")+".\n");
  }

  if (use_grids && hills_off_grid.size()) {
    for (hill_iter hoff = hills_off_grid.begin();
         hoff != hills_off_grid.end(); hoff++) {
      if (*h == *hoff) {
        hills_off_grid.erase (hoff);
        break;
      }
    }
  }

  if (hills_traj_os.good()) {
    // output to the trajectory
    hills_traj_os << "# DELETED this hill: "
                  << (hills.back()).output_traj()
                  << "\n";
    if (cvm::debug())
      hills_traj_os.flush();
  }

  return hills.erase (h);
}


cvm::real colvarbias_meta::update()
{
  if (cvm::debug())
    cvm::log ("Updating the metadynamics bias \""+this->name+"\""+
              ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");

  if (use_grids) {

    std::vector<int> curr_bin = hills_energy->get_colvars_index();

    if (expand_grids) {

      // first of all, expand the grids, if specified
      if (cvm::debug())
        cvm::log ("Metadynamics bias \""+this->name+"\""+
                  ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                  ": current coordinates on the grid: "+
                  cvm::to_str (curr_bin)+".\n");

      bool changed_grids = false;
      int const min_buffer =
        (3 * (size_t) std::floor (hill_width)) + 1;

      std::vector<int>         new_sizes (hills_energy->sizes());
      std::vector<colvarvalue> new_lower_boundaries (hills_energy->lower_boundaries);
      std::vector<colvarvalue> new_upper_boundaries (hills_energy->upper_boundaries);

      for (size_t i = 0; i < colvars.size(); i++) {

        if (! colvars[i]->expand_boundaries)
          continue;

        cvm::real &new_lb   = new_lower_boundaries[i].real_value;
        cvm::real &new_ub   = new_upper_boundaries[i].real_value;
        int       &new_size = new_sizes[i];
        bool changed_lb = false, changed_ub = false;

        if (!colvars[i]->hard_lower_boundary)
          if (curr_bin[i] < min_buffer) {
            int const extra_points = (min_buffer - curr_bin[i]);
            new_lb -= extra_points * colvars[i]->width;
            new_size += extra_points;
            // changed offset in this direction => the pointer needs to
            // be changed, too
            curr_bin[i] += extra_points;

            changed_lb = true;
            cvm::log ("Metadynamics bias \""+this->name+"\""+
                      ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                      ": new lower boundary for colvar \""+
                      colvars[i]->name+"\", at "+
                      cvm::to_str (new_lower_boundaries[i])+".\n");
          }

        if (!colvars[i]->hard_upper_boundary)
          if (curr_bin[i] > new_size - min_buffer - 1) {
            int const extra_points = (curr_bin[i] - (new_size - 1) + min_buffer);
            new_ub += extra_points * colvars[i]->width;
            new_size += extra_points;

            changed_ub = true;
            cvm::log ("Metadynamics bias \""+this->name+"\""+
                      ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                      ": new upper boundary for colvar \""+
                      colvars[i]->name+"\", at "+
                      cvm::to_str (new_upper_boundaries[i])+".\n");
          }

        if (changed_lb || changed_ub)
          changed_grids = true;
      }

      if (changed_grids) {

        // map everything into new grids

        colvar_grid_scalar *new_hills_energy =
          new colvar_grid_scalar (*hills_energy);
        colvar_grid_gradient *new_hills_energy_gradients =
          new colvar_grid_gradient (*hills_energy_gradients);

        // supply new boundaries to the new grids

        new_hills_energy->lower_boundaries = new_lower_boundaries;
        new_hills_energy->upper_boundaries = new_upper_boundaries;
        new_hills_energy->create (new_sizes, 0.0, 1);

        new_hills_energy_gradients->lower_boundaries = new_lower_boundaries;
        new_hills_energy_gradients->upper_boundaries = new_upper_boundaries;
        new_hills_energy_gradients->create (new_sizes, 0.0, colvars.size());

        new_hills_energy->map_grid (*hills_energy);
        new_hills_energy_gradients->map_grid (*hills_energy_gradients);

        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy = new_hills_energy;
        hills_energy_gradients = new_hills_energy_gradients;

        curr_bin = hills_energy->get_colvars_index();
        if (cvm::debug())
          cvm::log ("Coordinates on the new grid: "+
                    cvm::to_str (curr_bin)+".\n");
      }
    }
  }

  // add a new hill if the required time interval has passed
  if ((cvm::step_absolute() % new_hill_freq) == 0) {

    if (cvm::debug())
      cvm::log ("Metadynamics bias \""+this->name+"\""+
                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                ": adding a new hill at step "+cvm::to_str (cvm::step_absolute())+".\n");

    switch (comm) {

    case single_replica:
      if (well_tempered) {
        std::vector<int> curr_bin = hills_energy->get_colvars_index();
        cvm::real const hills_energy_sum_here = hills_energy->value(curr_bin);
        cvm::real const exp_weight = std::exp(-hills_energy_sum_here/(bias_temperature*cvm::boltzmann()));
        create_hill (hill ((hill_weight*exp_weight), colvars, hill_width));
      } else {
        create_hill (hill (hill_weight, colvars, hill_width));
      }
      break;

    case multiple_replicas:
      if (well_tempered) {
        std::vector<int> curr_bin = hills_energy->get_colvars_index();
        cvm::real const hills_energy_sum_here = hills_energy->value(curr_bin);
        cvm::real const exp_weight = std::exp(-hills_energy_sum_here/(bias_temperature*cvm::boltzmann()));
        create_hill (hill ((hill_weight*exp_weight), colvars, hill_width, replica_id));
      } else {
        create_hill (hill (hill_weight, colvars, hill_width, replica_id));
      }
      if (replica_hills_os.good()) {
        replica_hills_os << hills.back();
      } else {
        cvm::fatal_error ("Error: in metadynamics bias \""+this->name+"\""+
                          ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                          " while writing hills for the other replicas.\n");
      }
      break;
    }
  }

  // sync with the other replicas (if needed)
  if (comm != single_replica) {

    // reread the replicas registry
    if ((cvm::step_absolute() % replica_update_freq) == 0) {
      update_replicas_registry();
      // empty the output buffer
      replica_hills_os.flush();

      read_replica_files();
    }
  }

  // calculate the biasing energy and forces
  bias_energy = 0.0;
  for (size_t ir = 0; ir < colvars.size(); ir++) {
    colvar_forces[ir].reset();
  }
  if (comm == multiple_replicas)
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      replicas[ir]->bias_energy = 0.0;
      for (size_t ic = 0; ic < colvars.size(); ic++) {
        replicas[ir]->colvar_forces[ic].reset();
      }
    }

  if (use_grids) {

    // get the forces from the grid

    if ((cvm::step_absolute() % grids_freq) == 0) {
      // map the most recent gaussians to the grids
      project_hills (new_hills_begin, hills.end(),
                     hills_energy,    hills_energy_gradients);
      new_hills_begin = hills.end();

      // TODO: we may want to condense all into one replicas array,
      // including "this" as the first element
      if (comm == multiple_replicas) {
        for (size_t ir = 0; ir < replicas.size(); ir++) {
          replicas[ir]->project_hills (replicas[ir]->new_hills_begin,
                                       replicas[ir]->hills.end(),
                                       replicas[ir]->hills_energy,
                                       replicas[ir]->hills_energy_gradients);
          replicas[ir]->new_hills_begin = replicas[ir]->hills.end();
        }
      }
    }

    std::vector<int> curr_bin = hills_energy->get_colvars_index();
    if (cvm::debug())
      cvm::log ("Metadynamics bias \""+this->name+"\""+
                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                ": current coordinates on the grid: "+
                cvm::to_str (curr_bin)+".\n");

    if (hills_energy->index_ok (curr_bin)) {

      // within the grid: add the energy and the forces from there

      bias_energy += hills_energy->value (curr_bin);
      for (size_t ic = 0; ic < colvars.size(); ic++) {
        cvm::real const *f = &(hills_energy_gradients->value (curr_bin));
        colvar_forces[ic].real_value += -1.0 * f[ic];
        // the gradients are stored, not the forces
      }
      if (comm == multiple_replicas)
        for (size_t ir = 0; ir < replicas.size(); ir++) {
          bias_energy += replicas[ir]->hills_energy->value (curr_bin);
          cvm::real const *f = &(replicas[ir]->hills_energy_gradients->value (curr_bin));
          for (size_t ic = 0; ic < colvars.size(); ic++) {
            colvar_forces[ic].real_value += -1.0 * f[ic];
          }
        }

    } else {

      // off the grid: compute analytically only the hills at the grid's edges

      calc_hills (hills_off_grid.begin(), hills_off_grid.end(), bias_energy);
      for (size_t ic = 0; ic < colvars.size(); ic++) {
        calc_hills_force (ic, hills_off_grid.begin(), hills_off_grid.end(), colvar_forces);
      }

      if (comm == multiple_replicas)
        for (size_t ir = 0; ir < replicas.size(); ir++) {
          calc_hills (replicas[ir]->hills_off_grid.begin(),
                      replicas[ir]->hills_off_grid.end(),
                      bias_energy);
          for (size_t ic = 0; ic < colvars.size(); ic++) {
            calc_hills_force (ic,
                              replicas[ir]->hills_off_grid.begin(),
                              replicas[ir]->hills_off_grid.end(),
                              colvar_forces);
          }
        }
    }
  }

  // now include the hills that have not been binned yet (starting
  // from new_hills_begin)

  calc_hills (new_hills_begin, hills.end(), bias_energy);
  for (size_t ic = 0; ic < colvars.size(); ic++) {
    calc_hills_force (ic, new_hills_begin, hills.end(), colvar_forces);
  }

  if (cvm::debug())
    cvm::log ("Hills energy = "+cvm::to_str (bias_energy)+
              ", hills forces = "+cvm::to_str (colvar_forces)+".\n");

  if (cvm::debug())
    cvm::log ("Metadynamics bias \""+this->name+"\""+
              ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
              ": adding the forces from the other replicas.\n");

  if (comm == multiple_replicas)
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      calc_hills (replicas[ir]->new_hills_begin,
                  replicas[ir]->hills.end(),
                  bias_energy);
      for (size_t ic = 0; ic < colvars.size(); ic++) {
        calc_hills_force (ic,
                          replicas[ir]->new_hills_begin,
                          replicas[ir]->hills.end(),
                          colvar_forces);
      }
      if (cvm::debug())
        cvm::log ("Hills energy = "+cvm::to_str (bias_energy)+
                  ", hills forces = "+cvm::to_str (colvar_forces)+".\n");
    }

  return bias_energy;
}


void colvarbias_meta::calc_hills (colvarbias_meta::hill_iter      h_first,
                                  colvarbias_meta::hill_iter      h_last,
                                  cvm::real                      &energy,
                                  std::vector<colvarvalue> const &colvar_values)
{
  std::vector<colvarvalue> curr_values (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    curr_values[i].type (colvars[i]->type());
  }

  if (colvar_values.size()) {
    for (size_t i = 0; i < colvars.size(); i++) {
      curr_values[i] = colvar_values[i];
    }
  } else {
    for (size_t i = 0; i < colvars.size(); i++) {
      curr_values[i] = colvars[i]->value();
    }
  }

  for (hill_iter h = h_first; h != h_last; h++) {

    // compute the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    for (size_t i = 0; i < colvars.size(); i++) {
      colvarvalue const &x  = curr_values[i];
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      cv_sqdev += (colvars[i]->dist2 (x, center)) / (half_width*half_width);
    }

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-05)
      h->value (0.0);
    } else {
      h->value (std::exp (-0.5*cv_sqdev));
    }
    energy += h->energy();
  }
}


void colvarbias_meta::calc_hills_force (size_t const &i,
                                        colvarbias_meta::hill_iter      h_first,
                                        colvarbias_meta::hill_iter      h_last,
                                        std::vector<colvarvalue>       &forces,
                                        std::vector<colvarvalue> const &values)
{
  // Retrieve the value of the colvar
  colvarvalue x (values.size() ? values[i].type() : colvars[i]->type());
  x = (values.size() ? values[i] : colvars[i]->value());

  // do the type check only once (all colvarvalues in the hills series
  // were already saved with their types matching those in the
  // colvars)

  switch (colvars[i]->type()) {

  case colvarvalue::type_scalar:
    for (hill_iter h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].real_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (colvars[i]->dist2_lgrad (x, center)).real_value );
    }
    break;

  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    for (hill_iter h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].rvector_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (colvars[i]->dist2_lgrad (x, center)).rvector_value );
    }
    break;

  case colvarvalue::type_quaternion:
    for (hill_iter h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].quaternion_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (colvars[i]->dist2_lgrad (x, center)).quaternion_value );
    }
    break;

  case colvarvalue::type_notset:
    break;
  }
}


// **********************************************************************
// grid management functions
// **********************************************************************

void colvarbias_meta::project_hills (colvarbias_meta::hill_iter  h_first,
                                     colvarbias_meta::hill_iter  h_last,
                                     colvar_grid_scalar         *he,
                                     colvar_grid_gradient       *hg,
                                     bool print_progress)
{
  if (cvm::debug())
    cvm::log ("Metadynamics bias \""+this->name+"\""+
              ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
              ": projecting hills.\n");

  // TODO: improve it by looping over a small subgrid instead of the whole grid

  std::vector<colvarvalue> colvar_values (colvars.size());
  std::vector<cvm::real> colvar_forces_scalar (colvars.size());

  std::vector<int> he_ix = he->new_index();
  std::vector<int> hg_ix = (hg != NULL) ? hg->new_index() : std::vector<int> (0);
  cvm::real hills_energy_here = 0.0;
  std::vector<colvarvalue> hills_forces_here (colvars.size(), 0.0);

  size_t count = 0;
  size_t const print_frequency = ((hills.size() >= 1000000) ? 1 : (1000000/(hills.size()+1)));

  if (hg != NULL) {

    // loop over the points of the grid
    for ( ;
          (he->index_ok (he_ix)) && (hg->index_ok (hg_ix));
          count++) {

      for (size_t i = 0; i < colvars.size(); i++) {
        colvar_values[i] = hills_energy->bin_to_value_scalar (he_ix[i], i);
      }

      // loop over the hills and increment the energy grid locally
      hills_energy_here = 0.0;
      calc_hills (h_first, h_last, hills_energy_here, colvar_values);
      he->acc_value (he_ix, hills_energy_here);

      for (size_t i = 0; i < colvars.size(); i++) {
        hills_forces_here[i].reset();
        calc_hills_force (i, h_first, h_last, hills_forces_here, colvar_values);
        colvar_forces_scalar[i] = hills_forces_here[i].real_value;
      }
      hg->acc_force (hg_ix, &(colvar_forces_scalar.front()));

      he->incr (he_ix);
      hg->incr (hg_ix);

      if ((count % print_frequency) == 0) {
        if (print_progress) {
          cvm::real const progress = cvm::real (count) / cvm::real (hg->number_of_points());
          std::ostringstream os;
          os.setf (std::ios::fixed, std::ios::floatfield);
          os << std::setw (6) << std::setprecision (2)
             << 100.0 * progress
             << "% done.";
          cvm::log (os.str());
        }
      }
    }

  } else {

    // simpler version, with just the energy

    for ( ; (he->index_ok (he_ix)); ) {

      for (size_t i = 0; i < colvars.size(); i++) {
        colvar_values[i] = hills_energy->bin_to_value_scalar (he_ix[i], i);
      }

      hills_energy_here = 0.0;
      calc_hills (h_first, h_last, hills_energy_here, colvar_values);
      he->acc_value (he_ix, hills_energy_here);

      he->incr (he_ix);

      count++;
      if ((count % print_frequency) == 0) {
        if (print_progress) {
          cvm::real const progress = cvm::real (count) / cvm::real (he->number_of_points());
          std::ostringstream os;
          os.setf (std::ios::fixed, std::ios::floatfield);
          os << std::setw (6) << std::setprecision (2)
             << 100.0 * progress
             << "% done.";
          cvm::log (os.str());
        }
      }
    }
  }

  if (print_progress) {
    cvm::log ("100.00% done.");
  }

  if (! keep_hills) {
    hills.erase (hills.begin(), hills.end());
  }
}


void colvarbias_meta::recount_hills_off_grid (colvarbias_meta::hill_iter  h_first,
                                              colvarbias_meta::hill_iter  h_last,
                                              colvar_grid_scalar         *he)
{
  hills_off_grid.clear();

  for (hill_iter h = h_first; h != h_last; h++) {
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries (h->centers, true);
    if (min_dist < (3.0 * std::floor (hill_width)) + 1.0) {
      hills_off_grid.push_back (*h);
    }
  }
}



// **********************************************************************
// multiple replicas functions
// **********************************************************************


void colvarbias_meta::update_replicas_registry()
{
  if (cvm::debug())
    cvm::log ("Metadynamics bias \""+this->name+"\""+
              ": updating the list of replicas, currently containing "+
              cvm::to_str (replicas.size())+" elements.\n");

  {
    // copy the whole file into a string for convenience
    std::string line ("");
    std::ifstream reg_file (replicas_registry_file.c_str());
    if (reg_file.good()) {
      replicas_registry.clear();
      while (colvarparse::getline_nocomments (reg_file, line))
        replicas_registry.append (line+"\n");
    } else {
      cvm::fatal_error ("Error: failed to open file \""+replicas_registry_file+
                        "\" for reading.\n");
    }
  }

  // now parse it
  std::istringstream reg_is (replicas_registry);
  if (reg_is.good()) {

    std::string new_replica ("");
    std::string new_replica_file ("");
    while ((reg_is >> new_replica) && new_replica.size() &&
           (reg_is >> new_replica_file) && new_replica_file.size()) {

      if (new_replica == this->replica_id) {
        // this is the record for this same replica, skip it
        if (cvm::debug())
          cvm::log ("Metadynamics bias \""+this->name+"\""+
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
            cvm::log ("Metadynamics bias \""+this->name+"\""+
                      ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                      ": skipping a replica already loaded, \""+
                      (replicas[ir])->replica_id+"\".\n");
          already_loaded = true;
          break;
        }
      }

      if (!already_loaded) {
        // add this replica to the registry
        cvm::log ("Metadynamics bias \""+this->name+"\""+
                  ": accessing replica \""+new_replica+"\".\n");
        replicas.push_back (new colvarbias_meta());
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
          (replicas.back())->hills_energy           = new colvar_grid_scalar   (colvars);
          (replicas.back())->hills_energy_gradients = new colvar_grid_gradient (colvars);
        }
      }
    }

    // continue the cycle
    new_replica_file = "";
    new_replica = "";
  } else {
    cvm::fatal_error ("Error: cannot read the replicas registry file \""+
                      replicas_registry+"\".\n");
  }

  // now (re)read the list file of each replica
  for (size_t ir = 0; ir < replicas.size(); ir++) {
    if (cvm::debug())
      cvm::log ("Metadynamics bias \""+this->name+"\""+
                ": reading the list file for replica \""+
                (replicas[ir])->replica_id+"\".\n");

    std::ifstream list_is ((replicas[ir])->replica_list_file.c_str());
    std::string key;
    std::string new_state_file, new_hills_file;
    if (!(list_is >> key) ||
        !(key == std::string ("stateFile")) ||
        !(list_is >> new_state_file) ||
        !(list_is >> key) ||
        !(key == std::string ("hillsFile")) ||
        !(list_is >> new_hills_file)) {
      cvm::log ("Metadynamics bias \""+this->name+"\""+
                ": failed to read the file \""+
                (replicas[ir])->replica_list_file+"\": will try again after "+
                cvm::to_str (replica_update_freq)+" steps.\n");
      (replicas[ir])->update_status++;
    } else {
      (replicas[ir])->update_status = 0;
      if (new_state_file != (replicas[ir])->replica_state_file) {
        cvm::log ("Metadynamics bias \""+this->name+"\""+
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
    cvm::log ("Metadynamics bias \""+this->name+"\": the list of replicas contains "+
              cvm::to_str (replicas.size())+" elements.\n");
}


void colvarbias_meta::read_replica_files()
{
  for (size_t ir = 0; ir < replicas.size(); ir++) {

    if (! (replicas[ir])->replica_state_file_in_sync) {
      // if a new state file is being read, the hills file is also new
      (replicas[ir])->replica_hills_file_pos = 0;
    }

    // (re)read the state file if necessary
    if ( (! (replicas[ir])->has_data) ||
         (! (replicas[ir])->replica_state_file_in_sync) ) {

      cvm::log ("Metadynamics bias \""+this->name+"\""+
                ": reading the state of replica \""+
                (replicas[ir])->replica_id+"\" from file \""+
                (replicas[ir])->replica_state_file+"\".\n");

      std::ifstream is ((replicas[ir])->replica_state_file.c_str());
      if (! (replicas[ir])->read_restart (is)) {
        cvm::log ("Reading from file \""+(replicas[ir])->replica_state_file+
                  "\" failed or incomplete: will try again in "+
                  cvm::to_str (replica_update_freq)+" steps.\n");
      } else {
        // state file has been read successfully
        (replicas[ir])->replica_state_file_in_sync = true;
        (replicas[ir])->update_status = 0;
      }
      is.close();
    }

    // now read the hills added after writing the state file
    if ((replicas[ir])->replica_hills_file.size()) {

      if (cvm::debug())
        cvm::log ("Metadynamics bias \""+this->name+"\""+
                  ": checking for new hills from replica \""+
                  (replicas[ir])->replica_id+"\" in the file \""+
                  (replicas[ir])->replica_hills_file+"\".\n");

      // read hills from the other replicas' files; for each file, resume
      // the position recorded previously

      std::ifstream is ((replicas[ir])->replica_hills_file.c_str());
      if (is.good()) {

        // try to resume the previous position
        is.seekg ((replicas[ir])->replica_hills_file_pos, std::ios::beg);
        if (!is.good()){
          // if fail (the file may have been overwritten), reset this
          // position
          is.clear();
          is.seekg (0, std::ios::beg);
          // reset the counter
          (replicas[ir])->replica_hills_file_pos = 0;
          // schedule to reread the state file
          (replicas[ir])->replica_state_file_in_sync = false;
          // and record the failure
          (replicas[ir])->update_status++;
          cvm::log ("Failed to read the file \""+(replicas[ir])->replica_hills_file+
                    "\" at the previous position: will try again in "+
                    cvm::to_str (replica_update_freq)+" steps.\n");
        } else {

          while ((replicas[ir])->read_hill (is)) {
            //           if (cvm::debug())
              cvm::log ("Metadynamics bias \""+this->name+"\""+
                        ": received a hill from replica \""+
                        (replicas[ir])->replica_id+
                        "\" at step "+
                        cvm::to_str (((replicas[ir])->hills.back()).it)+".\n");
          }
          is.clear();
          // store the position for the next read
          (replicas[ir])->replica_hills_file_pos = is.tellg();
          if (cvm::debug())
            cvm::log ("Metadynamics bias \""+this->name+"\""+
                      ": stopped reading file \""+(replicas[ir])->replica_hills_file+
                      "\" at position "+
                      cvm::to_str ((replicas[ir])->replica_hills_file_pos)+".\n");

          // test whether this is the end of the file
          is.seekg (0, std::ios::end);
          if (is.tellg() > (replicas[ir])->replica_hills_file_pos+1) {
            (replicas[ir])->update_status++;
          } else {
            (replicas[ir])->update_status = 0;
          }
        }

      } else {
        cvm::log ("Failed to read the file \""+(replicas[ir])->replica_hills_file+
                  "\": will try again in "+
                  cvm::to_str (replica_update_freq)+" steps.\n");
        (replicas[ir])->update_status++;
        // cvm::fatal_error ("Error: cannot read from file \""+
        //                   (replicas[ir])->replica_hills_file+"\".\n");
      }
      is.close();
    }

    size_t const n_flush = (replica_update_freq/new_hill_freq + 1);
    if ((replicas[ir])->update_status > 3*n_flush) {
      // TODO: suspend the calculation?
      cvm::log ("WARNING: in metadynamics bias \""+this->name+"\""+
                " failed to read completely the output of replica \""+
                (replicas[ir])->replica_id+
                "\" after more than "+
                cvm::to_str ((replicas[ir])->update_status * replica_update_freq)+
                " steps.  Ensure that it is still running.\n");
    }
  }
}


// **********************************************************************
// input functions
// **********************************************************************


std::istream & colvarbias_meta::read_restart (std::istream& is)
{
  size_t const start_pos = is.tellg();

  if (comm == single_replica) {
    // if using a multiple replicas scheme, output messages
    // are printed before and after calling this function
    cvm::log ("Restarting metadynamics bias \""+this->name+"\""+
              ".\n");
  }
  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == "metadynamics") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block ("configuration", conf)) ) {

    if (comm == single_replica)
      cvm::log ("Error: in reading restart configuration for metadynamics bias \""+
                this->name+"\""+
                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                (replica_state_file_in_sync ? ("at position "+
                                               cvm::to_str (start_pos)+
                                               " in the state file") : "")+".\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  std::string name = "";
  if ( colvarparse::get_keyval (conf, "name", name,
                                std::string (""), colvarparse::parse_silent) &&
       (name != this->name) )
    cvm::fatal_error ("Error: in the restart file, the "
                      "\"metadynamics\" block has a different name: different system?\n");

  if (name.size() == 0) {
    cvm::fatal_error ("Error: \"metadynamics\" block within the restart file "
                      "has no identifiers.\n");
  }

  if (comm != single_replica) {
    std::string replica = "";
    if (colvarparse::get_keyval (conf, "replicaID", replica,
                                 std::string (""), colvarparse::parse_silent) &&
        (replica != this->replica_id))
      cvm::fatal_error ("Error: in the restart file, the "
                        "\"metadynamics\" block has a different replicaID: different system?\n");

    colvarparse::get_keyval (conf, "step", state_file_step,
                             cvm::step_absolute(), colvarparse::parse_silent);
  }

  bool grids_from_restart_file = use_grids;

  if (use_grids) {

    if (expand_grids) {
      // the boundaries of the colvars may have been changed; TODO:
      // this reallocation is only for backward-compatibility, and may
      // be deleted when grid_parameters (i.e. colvargrid's own
      // internal reallocation) has kicked in
      delete hills_energy;
      delete hills_energy_gradients;
      hills_energy = new colvar_grid_scalar (colvars);
      hills_energy_gradients = new colvar_grid_gradient (colvars);
    }

    colvar_grid_scalar   *hills_energy_backup = NULL;
    colvar_grid_gradient *hills_energy_gradients_backup = NULL;

    if (has_data) {
      if (cvm::debug())
        cvm::log ("Backupping grids for metadynamics bias \""+
                  this->name+"\""+
                  ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+".\n");
      hills_energy_backup           = hills_energy;
      hills_energy_gradients_backup = hills_energy_gradients;
      hills_energy                  = new colvar_grid_scalar (colvars);
      hills_energy_gradients        = new colvar_grid_gradient (colvars);
    }

    size_t const hills_energy_pos = is.tellg();
    if (!(is >> key)) {
      if (hills_energy_backup != NULL) {
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy           = hills_energy_backup;
        hills_energy_gradients = hills_energy_gradients_backup;
      }
      is.clear();
      is.seekg (hills_energy_pos, std::ios::beg);
      is.setstate (std::ios::failbit);
      return is;
    } else if (!(key == std::string ("hills_energy")) ||
               !(hills_energy->read_restart (is))) {
      is.clear();
      is.seekg (hills_energy_pos, std::ios::beg);
      grids_from_restart_file = false;
      if (!rebin_grids) {
        if (hills_energy_backup == NULL)
          cvm::fatal_error ("Error: couldn't read the free energy grid for metadynamics bias \""+
                            this->name+"\""+
                            ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                            "; if useGrids was off when the state file was written, "
                            "enable rebinGrids now to regenerate the grids.\n");
        else {
          if (comm == single_replica)
            cvm::log ("Error: couldn't read the free energy grid for metadynamics bias \""+
                      this->name+"\".\n");
          delete hills_energy;
          delete hills_energy_gradients;
          hills_energy           = hills_energy_backup;
          hills_energy_gradients = hills_energy_gradients_backup;
          is.setstate (std::ios::failbit);
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
      is.seekg (hills_energy_gradients_pos, std::ios::beg);
      is.setstate (std::ios::failbit);
      return is;
    } else if (!(key == std::string ("hills_energy_gradients")) ||
               !(hills_energy_gradients->read_restart (is))) {
      is.clear();
      is.seekg (hills_energy_gradients_pos, std::ios::beg);
      grids_from_restart_file = false;
      if (!rebin_grids) {
        if (hills_energy_backup == NULL)
          cvm::fatal_error ("Error: couldn't read the free energy gradients grid for metadynamics bias \""+
                            this->name+"\""+
                            ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
                            "; if useGrids was off when the state file was written, "
                            "enable rebinGrids now to regenerate the grids.\n");
        else {
          if (comm == single_replica)
            cvm::log ("Error: couldn't read the free energy gradients grid for metadynamics bias \""+
                      this->name+"\".\n");
          delete hills_energy;
          delete hills_energy_gradients;
          hills_energy           = hills_energy_backup;
          hills_energy_gradients = hills_energy_gradients_backup;
          is.setstate (std::ios::failbit);
          return is;
        }
      }
    }

    if (cvm::debug())
      cvm::log ("Successfully read new grids for bias \""+
                this->name+"\""+
                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+"\n");

    if (hills_energy_backup != NULL) {
      // now that we have successfully updated the grids, delete the
      // backup copies
      if (cvm::debug())
        cvm::log ("Deallocating the older grids.\n");

      delete hills_energy_backup;
      delete hills_energy_gradients_backup;
    }
  }

  bool const existing_hills = (hills.size() > 0);
  size_t const old_hills_size = hills.size();
  hill_iter old_hills_end = hills.end();
  hill_iter old_hills_off_grid_end = hills_off_grid.end();

  // read the hills explicitly written (if there are any)
  while (read_hill (is)) {
    if (cvm::debug())
      cvm::log ("Read a previously saved hill under the "
                "metadynamics bias \""+
                this->name+"\", created at step "+
                cvm::to_str ((hills.back()).it)+".\n");
  }
  is.clear();
  new_hills_begin = hills.end();
  if (grids_from_restart_file) {
    if (hills.size() > old_hills_size)
      cvm::log ("Read "+cvm::to_str (hills.size())+
                " hills in addition to the grids.\n");
  } else {
    if (hills.size())
      cvm::log ("Read "+cvm::to_str (hills.size())+" hills.\n");
  }

  if (rebin_grids) {

    // allocate new grids (based on the new boundaries and widths just
    // read from the configuration file), and project onto them the
    // grids just read from the restart file

    colvar_grid_scalar   *new_hills_energy =
      new colvar_grid_scalar (colvars);
    colvar_grid_gradient *new_hills_energy_gradients =
      new colvar_grid_gradient (colvars);

    if (!grids_from_restart_file || (keep_hills && hills.size())) {
      // if there are hills, recompute the new grids from them
      cvm::log ("Rebinning the energy and forces grids from "+
                cvm::to_str (hills.size())+" hills (this may take a while)...\n");
      project_hills (hills.begin(), hills.end(),
                     new_hills_energy, new_hills_energy_gradients, true);
      cvm::log ("rebinning done.\n");

    } else {
      // otherwise, use the grids in the restart file
      cvm::log ("Rebinning the energy and forces grids "
                "from the grids in the restart file.\n");
      new_hills_energy->map_grid (*hills_energy);
      new_hills_energy_gradients->map_grid (*hills_energy_gradients);
    }

    delete hills_energy;
    delete hills_energy_gradients;
    hills_energy = new_hills_energy;
    hills_energy_gradients = new_hills_energy_gradients;

    // assuming that some boundaries have expanded, eliminate those
    // off-grid hills that aren't necessary any more
    if (hills.size())
      recount_hills_off_grid (hills.begin(), hills.end(), hills_energy);
  }

  if (use_grids) {
    if (hills_off_grid.size()) {
      cvm::log (cvm::to_str (hills_off_grid.size())+" hills are near the "
                "grid boundaries: they will be computed analytically "
                "and saved to the state files.\n");
    }
  }

  is >> brace;
  if (brace != "}") {
    cvm::log ("Incomplete restart information for metadynamics bias \""+
              this->name+"\""+
              ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+
              ": no closing brace at position "+
              cvm::to_str (is.tellg())+" in the file.\n");
    is.setstate (std::ios::failbit);
    return is;
  }

  if (cvm::debug())
    cvm::log ("colvarbias_meta::read_restart() done\n");

  if (existing_hills) {
    hills.erase (hills.begin(), old_hills_end);
    hills_off_grid.erase (hills_off_grid.begin(), old_hills_off_grid_end);
  }

  has_data = true;

  if (comm != single_replica) {
    read_replica_files();
  }

  return is;
}


std::istream & colvarbias_meta::read_hill (std::istream &is)
{
  if (!is) return is; // do nothing if failbit is set

  size_t const start_pos = is.tellg();

  std::string data;
  if ( !(is >> read_block ("hill", data)) ) {
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  size_t h_it;
  get_keyval (data, "step", h_it, 0, parse_silent);
  if (h_it <= state_file_step) {
    if (cvm::debug())
      cvm::log ("Skipping a hill older than the state file for metadynamics bias \""+
                this->name+"\""+
                ((comm != single_replica) ? ", replica \""+replica_id+"\"" : "")+"\n");
    return is;
  }

  cvm::real h_weight;
  get_keyval (data, "weight", h_weight, hill_weight, parse_silent);

  std::vector<colvarvalue> h_centers (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    h_centers[i].type ((colvars[i]->value()).type());
  }
  {
    // it is safer to read colvarvalue objects one at a time;
    // TODO: change this it later
    std::string centers_input;
    key_lookup (data, "centers", centers_input);
    std::istringstream centers_is (centers_input);
    for (size_t i = 0; i < colvars.size(); i++) {
      centers_is >> h_centers[i];
    }
  }

  std::vector<cvm::real> h_widths (colvars.size());
  get_keyval (data, "widths", h_widths,
              std::vector<cvm::real> (colvars.size(), (std::sqrt (2.0 * PI) / 2.0)),
              parse_silent);

  std::string h_replica = "";
  if (comm != single_replica) {
    get_keyval (data, "replicaID", h_replica, replica_id, parse_silent);
    if (h_replica != replica_id)
      cvm::fatal_error ("Error: trying to read a hill created by replica \""+h_replica+
                        "\" for replica \""+replica_id+
                        "\"; did you swap output files?\n");
  }

  hill_iter const hills_end = hills.end();
  hills.push_back (hill (h_it, h_weight, h_centers, h_widths, h_replica));
  if (new_hills_begin == hills_end) {
    // if new_hills_begin is unset, set it for the first time
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {
    // add this also to the list of hills that are off-grid, which will
    // be computed analytically
    cvm::real const min_dist =
      hills_energy->bin_distance_from_boundaries ((hills.back()).centers, true);
    if (min_dist < (3.0 * std::floor (hill_width)) + 1.0) {
      hills_off_grid.push_back (hills.back());
    }
  }

  has_data = true;
  return is;
}




// **********************************************************************
// output functions
// **********************************************************************

std::ostream & colvarbias_meta::write_restart (std::ostream& os)
{
  os << "metadynamics {\n"
     << "  configuration {\n"
     << "    step " << cvm::step_absolute() << "\n"
     << "    name " << this->name << "\n";
  if (this->comm != single_replica)
    os << "    replicaID " << this->replica_id << "\n";
  os << "  }\n\n";

  if (use_grids) {

    // this is a very good time to project hills, if you haven't done
    // it already!
    project_hills (new_hills_begin, hills.end(),
                   hills_energy,    hills_energy_gradients);
    new_hills_begin = hills.end();

    // write down the grids to the restart file
    os << "  hills_energy\n";
    hills_energy->write_restart (os);
    os << "  hills_energy_gradients\n";
    hills_energy_gradients->write_restart (os);
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

  os << "}\n\n";

  if (comm != single_replica) {
    write_replica_state_file();
    // schedule to reread the state files of the other replicas (they
    // have also rewritten them)
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      (replicas[ir])->replica_state_file_in_sync = false;
    }
  }

  if (dump_fes) {
    write_pmf();
  }

  return os;
}


void colvarbias_meta::write_pmf()
{
  // allocate a new grid to store the pmf
  colvar_grid_scalar *pmf = new colvar_grid_scalar (*hills_energy);
  pmf->create();

  std::string fes_file_name_prefix (cvm::output_prefix);

  if ((cvm::n_meta_biases > 1) || (cvm::n_abf_biases > 0)) {
    // if this is not the only free energy integrator, append
    // this bias's name, to distinguish it from the output of the other
    // biases producing a .pmf file
    // TODO: fix for ABF with updateBias == no
    fes_file_name_prefix += ("."+this->name);
  }

  if ((comm == single_replica) || (dump_replica_fes)) {
    // output the PMF from this instance or replica
    pmf->reset();
    pmf->add_grid (*hills_energy);
    cvm::real const max = pmf->maximum_value();
    pmf->add_constant (-1.0 * max);
    pmf->multiply_constant (-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + cvm::temperature()) / bias_temperature;
      pmf->multiply_constant (well_temper_scale);
    }
    {
      std::string const fes_file_name (fes_file_name_prefix +
                                       ((comm != single_replica) ? ".partial" : "") +
                                       (dump_fes_save ?
                                        "."+cvm::to_str (cvm::step_absolute()) : "") +
                                       ".pmf");
      cvm::backup_file (fes_file_name.c_str());
      std::ofstream fes_os (fes_file_name.c_str());
      pmf->write_multicol (fes_os);
      fes_os.close();
    }
  }
  if (comm != single_replica) {
    // output the combined PMF from all replicas
    pmf->reset();
    pmf->add_grid (*hills_energy);
    for (size_t ir = 0; ir < replicas.size(); ir++) {
      pmf->add_grid (*(replicas[ir]->hills_energy));
    }
    cvm::real const max = pmf->maximum_value();
    pmf->add_constant (-1.0 * max);
    pmf->multiply_constant (-1.0);
    if (well_tempered) {
      cvm::real const well_temper_scale = (bias_temperature + cvm::temperature()) / bias_temperature;
      pmf->multiply_constant (well_temper_scale);
    }
    std::string const fes_file_name (fes_file_name_prefix +
                                     (dump_fes_save ?
                                      "."+cvm::to_str (cvm::step_absolute()) : "") +
                                     ".pmf");
    cvm::backup_file (fes_file_name.c_str());
    std::ofstream fes_os (fes_file_name.c_str());
    pmf->write_multicol (fes_os);
    fes_os.close();
  }

  delete pmf;
}



void colvarbias_meta::write_replica_state_file()
{
  // write down also the restart for the other replicas: TODO: this
  // is duplicated code, that could be cleaned up later
  cvm::backup_file (replica_state_file.c_str());
  std::ofstream rep_state_os (replica_state_file.c_str());
  if (!rep_state_os.good())
    cvm::fatal_error ("Error: in opening file \""+
                      replica_state_file+"\" for writing.\n");

  rep_state_os.setf (std::ios::scientific, std::ios::floatfield);
  rep_state_os << "\n"
               << "metadynamics {\n"
               << "  configuration {\n"
               << "    name " << this->name << "\n"
               << "    step " << cvm::step_absolute() << "\n";
  if (this->comm != single_replica) {
    rep_state_os << "    replicaID " << this->replica_id << "\n";
  }
  rep_state_os << "  }\n\n";
  rep_state_os << "  hills_energy\n";
  rep_state_os << std::setprecision (cvm::cv_prec)
               << std::setw (cvm::cv_width);
  hills_energy->write_restart (rep_state_os);
  rep_state_os << "  hills_energy_gradients\n";
  rep_state_os << std::setprecision (cvm::cv_prec)
               << std::setw (cvm::cv_width);
  hills_energy_gradients->write_restart (rep_state_os);

  if ( (!use_grids) || keep_hills ) {
    // write all hills currently in memory
    for (std::list<hill>::const_iterator h = this->hills.begin();
         h != this->hills.end();
         h++) {
      rep_state_os << *h;
    }
  } else {
    // write just those that are near the grid boundaries
    for (std::list<hill>::const_iterator h = this->hills_off_grid.begin();
         h != this->hills_off_grid.end();
         h++) {
      rep_state_os << *h;
    }
  }

  rep_state_os << "}\n\n";
  rep_state_os.close();

  // reopen the hills file
  replica_hills_os.close();
  replica_hills_os.open (replica_hills_file.c_str());
  if (!replica_hills_os.good())
    cvm::fatal_error ("Error: in opening file \""+
                      replica_hills_file+"\" for writing.\n");
  replica_hills_os.setf (std::ios::scientific, std::ios::floatfield);
}

std::string colvarbias_meta::hill::output_traj()
{
  std::ostringstream os;
  os.setf (std::ios::fixed, std::ios::floatfield);
  os << std::setw (cvm::it_width) << it << " ";

  os.setf (std::ios::scientific, std::ios::floatfield);

  os << "  ";
  for (size_t i = 0; i < centers.size(); i++) {
    os << " ";
    os << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)  << centers[i];
  }

  os << "  ";
  for (size_t i = 0; i < widths.size(); i++) {
    os << " ";
    os << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width) << widths[i];
  }

  os << "  ";
  os << std::setprecision (cvm::en_prec)
     << std::setw (cvm::en_width) << W << "\n";

  return os.str();
}


std::ostream & operator << (std::ostream &os, colvarbias_meta::hill const &h)
{
  os.setf (std::ios::scientific, std::ios::floatfield);

  os << "hill {\n";
  os << "  step " << std::setw (cvm::it_width) << h.it << "\n";
  os << "  weight   "
     << std::setprecision (cvm::en_prec)
     << std::setw (cvm::en_width)
     << h.W << "\n";

  if (h.replica.size())
    os << "  replicaID  " << h.replica << "\n";

  os << "  centers ";
  for (size_t i = 0; i < (h.centers).size(); i++) {
    os << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << h.centers[i];
  }
  os << "\n";

  os << "  widths  ";
  for (size_t i = 0; i < (h.widths).size(); i++) {
    os << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << h.widths[i];
  }
  os << "\n";

  os << "}\n";

  return os;
}
