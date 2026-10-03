/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_VECTOR_TIP4P_H
#define LMP_ELECTRODE_VECTOR_TIP4P_H

#include "electrode_vector.h"

namespace LAMMPS_NS {

class ElectrodeVectorTIP4P : public ElectrodeVector {
 public:
  ElectrodeVectorTIP4P(class LAMMPS *, int, int, double, bool);
  ~ElectrodeVectorTIP4P() override = default;

  void setup(class Pair *, class NeighList *, bool) override;
  void charge_force(int, const double *) override;
  int charge_force_virial(int, const double *, double *, int *) override;
  double charge_force_alpha() const override;

 protected:
  void charge_position(int, double *) override;
  double pair_cutsq(int, int) const override;

 private:
  double qdist;
  double cut_coul;
  double alpha;
  int typeO;
  int typeH;

  void find_M(int, int &, int &, double *);
};

}    // namespace LAMMPS_NS

#endif
