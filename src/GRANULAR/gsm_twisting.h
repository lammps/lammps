/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef GSM_CLASS
// clang-format off
GSMStyle(none,
         GSMTwistingNone,
         TWISTING);

GSMStyle(marshall,
         GSMTwistingMarshall,
         TWISTING);

GSMStyle(sds,
         GSMTwistingSDS,
         TWISTING);
// clang-format on
#else

#ifndef GSM_TWISTING_H_
#define GSM_TWISTING_H_

#include "gsm.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GSMTwisting : public GSM {
 public:
  GSMTwisting(class GranularModel *, class LAMMPS *);
  virtual ~GSMTwisting() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual void calculate_forces() = 0;
};

/* ---------------------------------------------------------------------- */

class GSMTwistingNone : public GSMTwisting {
 public:
  GSMTwistingNone(class GranularModel *, class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class GSMTwistingMarshall : public GSMTwisting {
 public:
  GSMTwistingMarshall(class GranularModel *, class LAMMPS *);
  void calculate_forces();
  void init();
 protected:
  double k_tang, mu_tang;
};

/* ---------------------------------------------------------------------- */

class GSMTwistingSDS : public GSMTwisting {
 public:
  GSMTwistingSDS(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double k, mu, damp;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GSM_TWISTING_H_ */
#endif /*GSM_CLASS_H_ */
