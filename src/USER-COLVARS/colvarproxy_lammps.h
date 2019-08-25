// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_LAMMPS_H
#define COLVARPROXY_LAMMPS_H

#include "colvarproxy_lammps_version.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"

#include "lammps.h"
#include "domain.h"
#include "force.h"
#include "update.h"

#include <string>
#include <vector>
#include <iostream>

/* struct for packed data communication of coordinates and forces. */
struct commdata {
  int tag,type;
  double x,y,z,m,q;
};

inline std::ostream & operator<< (std::ostream &out, const commdata &cd)
{
  out << " (" << cd.tag << "/" << cd.type << ": "
      << cd.x << ", " << cd.y << ", " << cd.z << ") ";
  return out;
}

/// \brief Communication between colvars and LAMMPS
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_lammps : public colvarproxy {

  // LAMMPS specific data objects and flags
 protected:

  // pointers to LAMMPS class instances
  class LAMMPS_NS::LAMMPS *_lmp;
  class LAMMPS_NS::RanPark *_random;

  // state of LAMMPS properties
  double t_target, my_timestep, my_boltzmann, my_angstrom;
  double bias_energy;
  int  restart_every;
  int  previous_step;

  bool first_timestep;
  bool total_force_requested;
  bool do_exit;

  // std::vector<int>          colvars_atoms;
  // std::vector<size_t>       colvars_atoms_ncopies;
  // std::vector<struct commdata> positions;
  // std::vector<struct commdata> total_forces;
  // std::vector<struct commdata> applied_forces;
  // std::vector<struct commdata> previous_applied_forces;

  std::vector<int>          atoms_types;

  MPI_Comm inter_comm;     // MPI comm with 1 root proc from each world
  int inter_me, inter_num; // rank for the inter replica comm

 public:
  friend class cvm::atom;
  colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp, const char *,
                     const char *, const int, const double, MPI_Comm);
  virtual ~colvarproxy_lammps();
  void init(const char*);
  int setup();

 // disable default and copy constructor
 private:
  colvarproxy_lammps() {};
  colvarproxy_lammps(const colvarproxy_lammps &) {};

  // methods for lammps to move data or trigger actions in the proxy
 public:
  void set_temperature(double t) { t_target = t; };
  bool total_forces_enabled() const { return total_force_requested; };
  bool total_forces_same_step() const { return true; };
  bool want_exit() const { return do_exit; };

  // perform colvars computation. returns biasing energy
  double compute();

  // dump status to string
  void serialize_status(std::string &);

  // set status from string
  bool deserialize_status(std::string &);

  // Write files expected from Colvars (called by post_run())
  void write_output_files();

  // read additional config from file
  void add_config_file(char const *config_filename);

  // read additional config from string
  void add_config_string(const std::string &config);

  // implementation of pure methods from base class
 public:

  inline cvm::real unit_angstrom() { return my_angstrom; };
  inline cvm::real boltzmann() { return my_boltzmann; };
  inline cvm::real temperature() { return t_target; };
  inline cvm::real dt() { return my_timestep; }; // return _lmp->update->dt * _lmp->force->femtosecond; };

  inline size_t restart_frequency() { return restart_every; };

  void add_energy(cvm::real energy) { bias_energy += energy; };
  void request_total_force(bool yesno) { total_force_requested = yesno; };

  void log(std::string const &message);
  void error(std::string const &message);
  void fatal_error(std::string const &message);

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                 cvm::atom_pos const &pos2) const;

  int backup_file(char const *filename);

  cvm::real rand_gaussian(void) { return _random->gaussian(); };

  int init_atom(int atom_number);
  int check_atom_id(int atom_number);

  inline std::vector<int> *modify_atom_types() { return &atoms_types; }

  // implementation of optional methods from base class
 public:

  // Multi-replica support
  // Indicate if multi-replica support is available and active
  virtual bool replica_enabled() { return (inter_comm != MPI_COMM_NULL); }

  // Index of this replica
  virtual int replica_index() { return inter_me; }

  // Total number of replica
  virtual int replica_num() { return inter_num; }

  // Synchronize replica
  virtual void replica_comm_barrier();

  // Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  // Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);
};

#endif

