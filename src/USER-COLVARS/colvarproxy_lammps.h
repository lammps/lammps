#ifndef COLVARPROXY_LAMMPS_H
#define COLVARPROXY_LAMMPS_H

#include "colvarmodule.h"
#include "colvarproxy.h"

#include "lammps.h"
#include "domain.h"
#include "force.h"
#include "random_park.h"
#include "update.h"

#include <string>
#include <vector>


/// \brief Communication between colvars and LAMMPS 
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_lammps : public colvarproxy {

  // LAMMPS specific data objects and flags
 protected:
  /// pointer to a LAMMPS class instance
  class LAMMPS_NS::LAMMPS *_lmp;
  class LAMMPS_NS::RanPark *_random;

  std::string input_prefix_str, output_prefix_str, restart_prefix_str;
  double t_target;
  int  restart_every;
  int  previous_step;

  bool first_timestep;
  bool system_force_requested;

  std::vector<int>          colvars_atoms;
  std::vector<size_t>       colvars_atoms_ncopies;
  std::vector<cvm::rvector> positions;
  std::vector<cvm::rvector> total_forces;
  std::vector<cvm::rvector> applied_forces;

 public:
  friend class cvm::atom;
  colvarproxy_lammps (LAMMPS_NS::LAMMPS *lmp, const char *,
		      const char *, const char *, const int);
  virtual ~colvarproxy_lammps();

 // disable default and copy constructor
 private:
  colvarproxy_lammps() {};
  colvarproxy_lammps(const colvarproxy_lammps &) {};

  // methods for lammps to push data or trigger action in the proxy
 public:

  // update temperature info
  void set_temperature(const double tt) { t_target = tt; };

  // perform colvars computation
  void compute();

  // implementation of pure methods from base class
 public:

  inline cvm::real unit_angstrom() { return _lmp->force->angstrom; };
  cvm::real boltzmann() { return _lmp->force->boltz; };
  cvm::real temperature() { return t_target; };
  cvm::real dt() { return _lmp->update->dt * _lmp->force->femtosecond; };

  inline std::string input_prefix() { return input_prefix_str; };
  inline std::string output_prefix() { return output_prefix_str; };
  inline std::string restart_output_prefix() { return restart_prefix_str; };
  inline size_t restart_frequency() { return restart_every; };

  void add_energy (cvm::real energy) { /* XXX */ ; };
  void request_system_force (bool yesno) { system_force_requested = yesno; };

  void log (std::string const &message);
  void fatal_error (std::string const &message);
  void exit (std::string const &message);

  cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2);
  cvm::real position_dist2 (cvm::atom_pos const &pos1,
                            cvm::atom_pos const &pos2);
  void select_closest_image (cvm::atom_pos &pos,
                             cvm::atom_pos const &ref_pos);


  void load_atoms(char const *filename,
                  std::vector<cvm::atom> &atoms,
                  std::string const pdb_field,
                  double const pdb_field_value = 0.0);

  void load_coords(char const *filename,
                   std::vector<cvm::atom_pos> &pos,
                   const std::vector<int> &indices,
                   std::string const pdb_field,
                   double const pdb_field_value = 0.0);

  void backup_file(char const *filename);

  cvm::real rand_gaussian(void) { return _random->gaussian(); };

};

#endif


// Emacs
// Local Variables:
// mode: C++
// End:
