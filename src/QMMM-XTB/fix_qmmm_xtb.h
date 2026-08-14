/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qmmm/xtb,FixQMMMXTB);
// clang-format on
#else

#ifndef LMP_FIX_QMMM_XTB_H
#define LMP_FIX_QMMM_XTB_H

#include "fix.h"

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

class PPPMTIP4PXTB;
class PPPMXTB;
class QMMMXTBEwald;

class FixQMMMXTB : public Fix {
 public:
  FixQMMMXTB(class LAMMPS *, int, char **);
  ~FixQMMMXTB() override;

  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void setup(int) override;
  void pre_force(int) override;
  void post_force(int) override;
  void min_setup(int) override;
  void min_pre_force(int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double memory_usage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  enum class PairCoulombMapping { FULL, MM_ONLY, INVALID };

  struct MMPoint {
    tagint tag;
    int atomic_number;
    double x[3];
    double charge;
    int nforce;
    tagint force_tags[3];
    double force_weights[3];
  };

  // Long-range Coulomb pair style captured for the MM-only/full reference
  // subtraction.  It may also contain LJ terms, which cancel between the two
  // reference evaluations and therefore remain unchanged in the final force.
  class Pair *pair_long;
  // True when hybrid/overlay pair_coeff mappings already restrict the
  // Coulomb-only sub-style to exactly the MM-MM type pairs.
  bool pair_coulomb_mm_only;
  // TIP4P parent atoms can only be checked after communication has acquired
  // the bonded H ghosts using the enlarged TIP4P communication cutoff.
  bool tip4p_qm_group_validated;
  std::vector<int> elements;
  int xtb_method, qm_charge, qm_uhf, maxiter;
  double cutoff, accuracy, electronic_temperature, mm_hardness;
  double image_alpha;
  std::array<int, 3> image_kmax;
  int image_ksqmax;

  int nqm;
  std::vector<tagint> qm_tags;
  std::unordered_map<tagint, int> qm_index;
  std::vector<int> qm_atomic_numbers;
  std::vector<double> qm_x, qm_x_wrapped, qm_charge_scf;
  std::vector<double> qm_gradient, qm_force_correction;
  std::vector<MMPoint> mm_points;
  std::vector<double> mm_gradient, mm_force_correction;

  std::vector<double> entry_force;
  std::vector<double> pair_mm_force, pair_full_force, qmqm_kspace_force;
  std::array<double, 6> pair_mm_virial, pair_full_virial, qmqm_kspace_virial;
  double pair_mm_energy, pair_full_energy, mm_kspace_energy;
  double qm_energy, energy_correction;

  PPPMXTB *pppm_xtb;
  PPPMTIP4PXTB *pppm_tip4p_xtb;
  std::unique_ptr<QMMMXTBEwald> image_ewald;
  bool adapter_active;

  void gather_qm_atoms(bool);
  void gather_mm_points();
  int get_charge_site(int, double *, int *, double *);
  void compute_group_potential(double *, int, int, bool);
  void validate_tip4p_qm_group();
  PairCoulombMapping classify_pair_coulomb_mapping(char *);
  void set_qm_charges(double);
  void restore_qm_charges();
  void save_entry_forces();
  void restore_entry_forces();
  void clear_forces();
  void capture_pair(std::vector<double> &, double &, std::array<double, 6> &);
  void capture_kspace(std::vector<double> &, double &, std::array<double, 6> &);
  void compute_mm_shift(std::vector<double> &);
  void run_xtb(const std::vector<double> &, const std::vector<double> &);
  void build_periodic_forces();
};

}    // namespace LAMMPS_NS

#endif
#endif
