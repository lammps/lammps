// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_LAMMPS_H
#define COLVARPROXY_LAMMPS_H

#include "colvarproxy_lammps_version.h"    // IWYU pragma: export

#include "colvarmodule.h"
#include "colvarproxy.h"

// forward declarations

namespace LAMMPS_NS {
class LAMMPS;
class RanPark;
}    // namespace LAMMPS_NS

/// \brief Communication between colvars and LAMMPS
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_lammps : public colvarproxy {

  // LAMMPS specific data objects and flags
 protected:
  // pointers to LAMMPS class instances
  LAMMPS_NS::LAMMPS *_lmp;
  LAMMPS_NS::RanPark *_random;

  // state of LAMMPS properties
  double bias_energy;
  cvm::step_number previous_step;

  bool first_timestep;
  bool do_exit;

  std::vector<int> atoms_types;

 public:
  colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp);
  ~colvarproxy_lammps() override;

  // disable default and copy constructor
  colvarproxy_lammps() = delete;
  colvarproxy_lammps(const colvarproxy_lammps &) = delete;

  void init();

  /// Set the internal seed used by \link rand_gaussian() \endlink
  void set_random_seed(int seed);

  int setup() override;

  // methods for lammps to move data or trigger actions in the proxy
 public:
  bool total_forces_enabled() const override { return total_force_requested; };
  // Total forces are saved at end of step, only processed at the next step
  bool total_forces_same_step() const override { return false; };
  bool want_exit() const { return do_exit; };

  // perform colvars computation. returns biasing energy
  double compute();

  // Request to set the units used internally by Colvars
  int set_unit_system(std::string const &units_in, bool check_only) override;

  /// Convert a command-line argument to string
  char const *script_obj_to_str(unsigned char *obj);

  /// Convert a command-line argument to a vector of strings
  std::vector<std::string> script_obj_to_str_vector(unsigned char *obj);

  void add_energy(cvm::real energy) override { bias_energy += energy; };
  void request_total_force(bool yesno) override { total_force_requested = yesno; };

  void log(std::string const &message) override;
  void error(std::string const &message) override;

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                 cvm::atom_pos const &pos2) const override;

  cvm::real rand_gaussian() override;

  int init_atom(int atom_number) override;
  int init_atom(cvm::residue_id const &residue, std::string const &atom_name,
                std::string const &segment_id) override
  {
    return colvarproxy::init_atom(residue, atom_name, segment_id);
  }

  int check_atom_id(int atom_number) override;
  int check_atom_id(cvm::residue_id const &residue, std::string const &atom_name,
                    std::string const &segment_id) override
  {
    return colvarproxy::check_atom_id(residue, atom_name, segment_id);
  }

  inline std::vector<int> *modify_atom_types() { return &atoms_types; }
};

#endif
