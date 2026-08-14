/* -*- c++ -*- ----------------------------------------------------------
   pair_style template/offset

   Adds a constant energy offset to the system for every distinct match of a
   molecule template, decoupled from fix msevb.  Because the offset is a per-
   configuration constant it contributes energy but no forces.  Matching reuses
   the msevb superimpose graph-isomorphism kernel.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(template/offset,PairTemplateOffset);
// clang-format on
#else

#ifndef LMP_PAIR_TEMPLATE_OFFSET_H
#define LMP_PAIR_TEMPLATE_OFFSET_H

#include "fix_msevb_superimpose.h"    // RefTopo
#include "pair.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

class PairTemplateOffset : public Pair {
 public:
  PairTemplateOffset(class LAMMPS *);
  ~PairTemplateOffset() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  double cut_global;

  // One registered (molecule template, offset) pair.
  struct Entry {
    std::string mol_id;    // molecule-command id
    int mol_index;         // index into atom->molecules[]
    double offset;         // energy added per distinct match (energy units)
  };
  std::vector<Entry> entries;

  // Global 1-2 topology snapshot (by global tag) so template matches that span
  // MPI ranks are found; rebuilt each compute() call.
  RefTopo gref;
  std::vector<int> g_type;
  std::vector<int> g_nspecial;
  std::vector<tagint> g_special;

  void allocate();
  void build_global_topology();
  int count_matches(class Molecule *mol);
};

}    // namespace LAMMPS_NS

#endif
#endif
