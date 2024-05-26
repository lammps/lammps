// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_META_H
#define COLVARBIAS_META_H

#include <vector>
#include <list>
#include <iosfwd>

#include "colvarbias.h"
#include "colvargrid.h"

/// Metadynamics bias (implementation of \link colvarbias \endlink)
class colvarbias_meta
  : public virtual colvarbias,
    public virtual colvarbias_ti
{

public:

  /// Communication between different replicas
  enum Communication {
    /// One replica (default)
    single_replica,
    /// Hills added concurrently by several replicas
    multiple_replicas
  };

  /// Communication between different replicas
  Communication comm;

  colvarbias_meta(char const *key);
  virtual ~colvarbias_meta();

  virtual int init(std::string const &conf);
  virtual int init_replicas_params(std::string const &conf);
  virtual int init_well_tempered_params(std::string const &conf);
  virtual int init_ebmeta_params(std::string const &conf);

  virtual int clear_state_data();

  virtual int update();
  virtual int update_grid_params();
  virtual int update_bias();
  virtual int update_grid_data();
  virtual int replica_share();

  virtual int calc_energy(std::vector<colvarvalue> const *values);
  virtual int calc_forces(std::vector<colvarvalue> const *values);

  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &state_conf);
  virtual std::ostream & write_state_data(std::ostream &os);
  virtual std::istream & read_state_data(std::istream &os);

  virtual int setup_output();
  virtual int write_output_files();
  virtual void write_pmf();
  virtual int write_state_to_replicas();

  class hill;
  typedef std::list<hill>::iterator hill_iter;

protected:

  /// Width of a hill in number of grid points
  ///
  /// The local width of each collective variable, multiplied by this
  /// number, provides the hill width along that direction
  cvm::real hill_width;

  /// The sigma parameters of the Gaussian hills
  std::vector<cvm::real> colvar_sigmas;

  /// \brief Number of simulation steps between two hills
  size_t     new_hill_freq;

  /// Write the hill logfile
  bool b_hills_traj;

  /// Name of the hill logfile
  std::string const hills_traj_file_name() const;

  /// \brief List of hills used on this bias (total); if a grid is
  /// employed, these don't need to be updated at every time step
  std::list<hill> hills;

  /// \brief Iterator to the first of the "newest" hills (when using
  /// grids, those who haven't been mapped yet)
  hill_iter new_hills_begin;

  /// \brief List of hills used on this bias that are on the boundary
  /// edges; these are updated regardless of whether hills are used
  std::list<hill> hills_off_grid;

  /// \brief Same as new_hills_begin, but for the off-grid ones
  hill_iter new_hills_off_grid_begin;

  /// Regenerate the hills_off_grid list
  void recount_hills_off_grid(hill_iter h_first, hill_iter h_last,
                               colvar_grid_scalar *ge);

  /// Read a hill from a file
  std::istream & read_hill(std::istream &is);

  /// \brief Add a new hill; if a .hills trajectory is written,
  /// write it there; if there is more than one replica, communicate
  /// it to the others
  std::list<hill>::const_iterator add_hill(hill const &h);

  /// \brief Remove a previously saved hill (returns an iterator for
  /// the next hill in the list)
  std::list<hill>::const_iterator delete_hill(hill_iter &h);

  /// \brief Calculate the values of the hills, incrementing
  /// bias_energy
  virtual void calc_hills(hill_iter  h_first,
                          hill_iter  h_last,
                          cvm::real &energy,
                          std::vector<colvarvalue> const *values);

  /// \brief Calculate the forces acting on the i-th colvar,
  /// incrementing colvar_forces[i]; must be called after calc_hills
  /// each time the values of the colvars are changed
  virtual void calc_hills_force(size_t const &i,
                                hill_iter h_first,
                                hill_iter h_last,
                                std::vector<colvarvalue> &forces,
                                std::vector<colvarvalue> const *values);


  /// Height of new hills
  cvm::real  hill_weight;

  /// \brief Bin the hills on grids of energy and forces, and use them
  /// to force the colvars (as opposed to deriving the hills analytically)
  bool       use_grids;

  /// \brief Rebin the hills upon restarting
  bool       rebin_grids;

  /// \brief Should the grids be expanded if necessary?
  bool       expand_grids;

  /// \brief How often the hills should be projected onto the grids
  size_t     grids_freq;

  /// Keep hills in the restart file (e.g. to accurately rebin later)
  bool       keep_hills;

  /// value of keepHills saved in the most recent restart file
  bool restart_keep_hills;

  /// \brief Dump the free energy surface (.pmf file) every restartFrequency
  bool       dump_fes;

  /// \brief Dump the free energy surface (.pmf file) every restartFrequency
  /// using only the hills from this replica (only applicable to more than one replica)
  bool       dump_replica_fes;

  /// \brief Dump the free energy surface files at different
  /// time steps, appending the step number to each file
  bool       dump_fes_save;

  /// \brief Whether to use well-tempered metadynamics
  bool       well_tempered;

  /// \brief Biasing temperature in well-tempered metadynamics
  cvm::real  bias_temperature;

  /// Ensemble-biased metadynamics (EBmeta) flag
  bool       ebmeta;

  /// Target distribution for EBmeta
  colvar_grid_scalar* target_dist;

  /// Number of equilibration steps for EBmeta
  cvm::step_number ebmeta_equil_steps;


  /// \brief Try to read the restart information by allocating new
  /// grids before replacing the current ones (used e.g. in
  /// multiple_replicas)
  bool       safely_read_restart;

  /// Hill energy, cached on a grid
  colvar_grid_scalar    *hills_energy;

  /// Hill forces, cached on a grid
  colvar_grid_gradient  *hills_energy_gradients;

  /// \brief Project the selected hills onto grids
  void project_hills(hill_iter h_first, hill_iter h_last,
                      colvar_grid_scalar *ge, colvar_grid_gradient *gf,
                      bool print_progress = false);


  // Multiple Replicas variables and functions

  /// \brief Identifier for this replica
  std::string            replica_id;

  /// \brief File containing the paths to the output files from this replica
  std::string            replica_file_name;

  /// \brief Read the existing replicas on registry
  virtual int update_replicas_registry();

  /// \brief Read new data from replicas' files
  virtual int read_replica_files();

  /// Write full state information to be read by other replicas
  virtual int write_replica_state_file();

  /// Call this after write_replica_state_file()
  virtual int reopen_replica_buffer_file();

  /// \brief Additional, "mirror" metadynamics biases, to collect info
  /// from the other replicas
  ///
  /// These are supposed to be synchronized by reading data from the
  /// other replicas, and not be modified by the "local" replica
  std::vector<colvarbias_meta *> replicas;

  /// \brief Frequency at which data the "mirror" biases are updated
  size_t                 replica_update_freq;

  /// List of replicas (and their output list files): contents are
  /// copied into replicas_registry for convenience
  std::string            replicas_registry_file;
  /// List of replicas (and their output list files)
  std::string            replicas_registry;
  /// List of files written by this replica
  std::string            replica_list_file;

  /// Hills energy and gradients written specifically for other
  /// replica (in addition to its own restart file)
  std::string            replica_state_file;
  /// Whether a mirror bias has read the latest version of its state file
  bool                   replica_state_file_in_sync;

  /// If there was a failure reading one of the files (because they
  /// are not complete), this counter is incremented
  size_t                 update_status;

  /// Explicit hills communicated between replicas
  ///
  /// This file becomes empty after replica_state_file is rewritten
  std::string            replica_hills_file;

  /// Position within replica_hills_file (when reading it)
  std::streampos         replica_hills_file_pos;

};




/// \brief A hill for the metadynamics bias
class colvarbias_meta::hill {

protected:

  /// Time step at which this hill was added
  cvm::step_number it;

  /// Value of the hill function (ranges between 0 and 1)
  cvm::real hill_value;

  /// Scale factor, which could be modified at runtime (default: 1)
  cvm::real sW;

  /// Maximum height in energy of the hill
  cvm::real W;

  /// Centers of the hill in the collective variable space
  std::vector<colvarvalue> centers;

  /// Half-widths of the hill in the collective variable space
  std::vector<cvm::real> sigmas;

  /// Identity of the replica who added this hill
  std::string replica;

public:

  friend class colvarbias_meta;

  /// Constructor of a hill object
  /// \param it Step number at which the hill was added
  /// \param W Weight of the hill (energy units)
  /// \param cv_values Array of collective variable values
  /// \param cv_sigmas Array of collective variable values
  /// \param replica ID of the replica that creates the hill (optional)
  hill(cvm::step_number it, cvm::real W,
       std::vector<colvarvalue> const &cv_values,
       std::vector<cvm::real> const &cv_sigmas,
       std::string const &replica = "");

  /// Copy constructor
  hill(colvarbias_meta::hill const &h);

  /// Destructor
  ~hill();

  /// Assignment operator
  hill & operator = (colvarbias_meta::hill const &h);

  /// Get the energy
  inline cvm::real energy()
  {
    return W * sW * hill_value;
  }

  /// Get the energy using another hill weight
  inline cvm::real energy(cvm::real const &new_weight)
  {
    return new_weight * sW * hill_value;
  }

  /// Get the current hill value
  inline cvm::real const &value()
  {
    return hill_value;
  }

  /// Set the hill value as specified
  inline void value(cvm::real const &new_value)
  {
    hill_value = new_value;
  }

  /// Get the weight
  inline cvm::real weight()
  {
    return W * sW;
  }

  /// Scale the weight with this factor (by default 1.0 is used)
  inline void scale(cvm::real const &new_scale_fac)
  {
    sW = new_scale_fac;
  }

  /// Get the center of the hill
  inline std::vector<colvarvalue> & center()
  {
    return centers;
  }

  /// Get the i-th component of the center
  inline colvarvalue & center(size_t const &i)
  {
    return centers[i];
  }

  /// Comparison operator
  inline friend bool operator < (hill const &h1, hill const &h2)
  {
    if (h1.it < h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator <= (hill const &h1, hill const &h2)
  {
    if (h1.it <= h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator > (hill const &h1, hill const &h2)
  {
    if (h1.it > h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator >= (hill const &h1, hill const &h2)
  {
    if (h1.it >= h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator == (hill const &h1, hill const &h2)
  {
    if ( (h1.it >= h2.it) && (h1.replica == h2.replica) ) return true;
    else return false;
  }

  /// Represent the hill ina string suitable for a trajectory file
  std::string output_traj();

  /// Write the hill to an output stream
  friend std::ostream & operator << (std::ostream &os, hill const &h);

};


#endif
