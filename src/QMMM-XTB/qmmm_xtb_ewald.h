/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_QMMM_XTB_EWALD_H
#define LMP_QMMM_XTB_EWALD_H

#include <array>
#include <vector>

namespace LAMMPS_NS {

// Small direct-Ewald operator used inside the xTB SCC loop.  PME remains the
// production long-range solver for the MM system; this dense operator is only
// for the usually small QM-image response that changes every SCC iteration.
class QMMMXTBEwald {
 public:
  void setup(const std::array<double, 3> &, double, const std::array<int, 3> &, int);
  void response(const std::vector<double> &, std::vector<double> &) const;
  double energy(const std::vector<double> &, const std::vector<double> &) const;
  void add_force(const std::vector<double> &, const std::vector<double> &,
                 std::vector<double> &) const;

  static void add_erf_pair_force(const double *, const double *, double, double, double, double *,
                                 double *);

 private:
  struct KTerm {
    double k[3];
    double coefficient;
  };

  std::array<double, 3> box_{};
  double alpha_ = 0.0;
  double volume_ = 0.0;
  std::vector<KTerm> kterms_;
};

}    // namespace LAMMPS_NS

#endif
