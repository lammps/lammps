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
#include <iostream>

/* struct for packed data communication of coordinates and forces. */
struct commdata { 
  int tag,type; 
  double x,y,z; 
};

inline std::ostream & operator<< (std::ostream &out, const commdata &cd)
{
  out << " (" << cd.tag << "/" << cd.type << ": " 
      << cd.x << ", " << cd.y << ", " << cd.z << ") ";
  return out;
};

/// \brief Communication between colvars and LAMMPS 
/// (implementation of \link colvarproxy \endlink)
class colvarproxy_lammps : public colvarproxy {

  // LAMMPS specific data objects and flags
 protected:

  // pointers to LAMMPS class instances
  class LAMMPS_NS::LAMMPS *_lmp;
  class LAMMPS_NS::RanPark *_random;

  // pointers to LAMMPS provided storage
  const int *_typemap;

  // state of LAMMPS properties
  double t_target;
  double bias_energy;
  int  restart_every;
  int  previous_step;

  std::string input_prefix_str;
  std::string output_prefix_str;
  std::string restart_prefix_str;

  bool first_timestep;
  bool system_force_requested;
  bool do_exit;

  std::vector<int>          colvars_atoms;
  std::vector<size_t>       colvars_atoms_ncopies;
  std::vector<struct commdata> positions;
  std::vector<struct commdata> total_forces;
  std::vector<struct commdata> applied_forces;

 public:
  friend class cvm::atom;
  colvarproxy_lammps (LAMMPS_NS::LAMMPS *lmp, const char *, const char *,
		      const char *, const int, const double, const int *);
  virtual ~colvarproxy_lammps();

 // disable default and copy constructor
 private:
  colvarproxy_lammps() {};
  colvarproxy_lammps(const colvarproxy_lammps &) {};

  // methods for lammps to move data or trigger actions in the proxy
 public:
  void set_temperature(double t) { t_target = t; };
  bool need_system_forces() const { return  system_force_requested; };
  bool want_exit() const { return do_exit; };
  std::vector<int> *             get_tags()   { return &colvars_atoms; };
  std::vector<struct commdata> * get_coords() { return &positions; };
  std::vector<struct commdata> * get_forces() { return &applied_forces; };
  std::vector<struct commdata> * get_oforce() { return &total_forces; };

  // initialize atom structure
  int init_lammps_atom(const int &, cvm::atom *);

  // perform colvars computation. returns biasing energy
  double compute();

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

  void add_energy (cvm::real energy) { bias_energy = energy; };
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
