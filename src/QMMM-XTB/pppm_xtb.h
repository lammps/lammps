/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/xtb,PPPMXTB);
// clang-format on
#else

#ifndef LMP_PPPM_XTB_H
#define LMP_PPPM_XTB_H

#include "pppm.h"

namespace LAMMPS_NS {

template <class> class QMMMXTBPPPMHelper;

class PPPMXTB : public PPPM {
  friend class QMMMXTBPPPMHelper<PPPMXTB>;

 public:
  explicit PPPMXTB(class LAMMPS *);
  void compute_group_potential(double *, int, int, bool);
  int get_charge_site(int, double *, int *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
