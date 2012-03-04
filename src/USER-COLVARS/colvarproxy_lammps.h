#ifndef COLVARPROXY_LAMMPS_H
#define COLVARPROXY_LAMMPS_H

#include "colvarmodule.h"
#include "colvarproxy.h"

#include "lammps.h"
#include "force.h"
#include "update.h"

#include <string>
#include <vector>


/// \brief Communication between colvars and LAMMPS (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_lammps : public colvarproxy {

 private:
  class LAMMPS_NS::LAMMPS *_lmp;

  double thermostat_temperature;
  
#if 0
  /// Pointer to the NAMD simulation input object
  SimParameters const *simparams;
  /// Pointer to the NAMD boundary conditions object
  Lattice const       *lattice;


  /// NAMD-style PRNG object
  Random random;
#endif

  std::string input_prefix_str, output_prefix_str, restart_output_prefix_str;
  size_t      restart_frequency_s;
  size_t      previous_NAMD_step;
  bool        first_timestep;
  bool        system_force_requested;

  std::vector<int>          colvars_atoms;
  std::vector<size_t>       colvars_atoms_ncopies;
  std::vector<cvm::rvector> positions;
  std::vector<cvm::rvector> total_forces;
  std::vector<cvm::rvector> applied_forces;

//  size_t init_namd_atom (AtomID const &aid);

//  SubmitReduction *reduction;
public:

  friend class cvm::atom;

  colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp);
  virtual ~colvarproxy_lammps();

#if 0
  /// \brief Reimplements GlobalMaster member function, to be called
  /// at every step
  void calculate();

  void add_energy (cvm::real energy);
  void request_system_force (bool yesno);
  void log (std::string const &message);
  void fatal_error (std::string const &message);
  void exit (std::string const &message);
#endif

  inline cvm::real unit_angstrom()
  {
    return _lmp->force->angstrom;
  }

  cvm::real boltzmann()
  {
    return _lmp->force->boltz;
  }

  cvm::real temperature()
  {
    return thermostat_temperature;
  }

  cvm::real dt()
  {
    return _lmp->update->dt;
  }

  inline std::string input_prefix()
  {
    return input_prefix_str;
  }
  inline std::string restart_output_prefix()
  {
    return restart_output_prefix_str;
  }
  inline std::string output_prefix()
  {
    return output_prefix_str;
  }
  inline size_t restart_frequency()
  {
    return restart_frequency_s;
  }

#if 0
  cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2);
  cvm::real position_dist2 (cvm::atom_pos const &pos1,
                            cvm::atom_pos const &pos2);

  void select_closest_image (cvm::atom_pos &pos,
                             cvm::atom_pos const &ref_pos);


  void load_atoms (char const *filename,
                   std::vector<cvm::atom> &atoms,
                   std::string const pdb_field,
                   double const pdb_field_value = 0.0);

  void load_coords (char const *filename,
                    std::vector<cvm::atom_pos> &pos,
                    const std::vector<int> &indices,
                    std::string const pdb_field,
                    double const pdb_field_value = 0.0);

  void backup_file (char const *filename);

  cvm::real rand_gaussian (void)
  {
    return random.gaussian();
  }

#endif
};

#if 0
inline cvm::rvector colvarproxy_namd::position_distance (cvm::atom_pos const &pos1,
                                                         cvm::atom_pos const &pos2)
{
  Position const p1 (pos1.x, pos1.y, pos1.z);
  Position const p2 (pos2.x, pos2.y, pos2.z);
  // return p2 - p1
  Vector const d = this->lattice->delta (p2, p1);
  return cvm::rvector (d.x, d.y, d.z);
}


inline void colvarproxy_namd::select_closest_image (cvm::atom_pos &pos,
                                                    cvm::atom_pos const &ref_pos)
{
  Position const p (pos.x, pos.y, pos.z);
  Position const rp (ref_pos.x, ref_pos.y, ref_pos.z);
  ScaledPosition const srp = this->lattice->scale (rp);
  Position const np = this->lattice->nearest (p, srp);
  pos.x = np.x;
  pos.y = np.y;
  pos.z = np.z;
}


inline cvm::real colvarproxy_namd::position_dist2 (cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2)
{
  Lattice const *l = this->lattice;
  Vector const p1 (pos1.x, pos1.y, pos1.z);
  Vector const p2 (pos2.x, pos2.y, pos2.z);
  Vector const d = l->delta (p1, p2);
  return cvm::real (d.x*d.x + d.y*d.y + d.z*d.z);
}

#endif

#endif


// Emacs
// Local Variables:
// mode: C++
// End:
