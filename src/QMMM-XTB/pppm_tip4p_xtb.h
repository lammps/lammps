/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/tip4p/xtb,PPPMTIP4PXTB);
// clang-format on
#else

#ifndef LMP_PPPM_TIP4P_XTB_H
#define LMP_PPPM_TIP4P_XTB_H

#include "pppm_tip4p.h"

namespace LAMMPS_NS {

template <class> class QMMMXTBPPPMHelper;

class PPPMTIP4PXTB : public PPPMTIP4P {
  friend class QMMMXTBPPPMHelper<PPPMTIP4PXTB>;

 public:
  explicit PPPMTIP4PXTB(class LAMMPS *);
  void compute_group_potential(double *, int, int, bool);
  int get_charge_site(int, double *, int *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
