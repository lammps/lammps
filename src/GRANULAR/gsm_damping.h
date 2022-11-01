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
         GSMDampingNone,
         DAMPING);

GSMStyle(velocity,
         GSMDampingVelocity,
         DAMPING);

GSMStyle(mass_velocity,
         GSMDampingMassVelocity,
         DAMPING);

GSMStyle(viscoelastic,
         GSMDampingViscoelastic,
         DAMPING);

GSMStyle(tsuji,
         GSMDampingTsuji,
         DAMPING);
// clang-format on
#else

#ifndef GSM_DAMPING_H_
#define GSM_DAMPING_H_

#include "gsm.h"
#include "pointers.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GSMDamping : public GSM {
 public:
  GSMDamping(class GranularModel *, class LAMMPS *);
  ~GSMDamping() {};
  virtual void coeffs_to_local() {};
  virtual void mix_coeffs(double*, double*) {};
  virtual void init();
  virtual double calculate_forces() = 0;
  double damp_prefactor; // Used by tangential models
 protected:
  double damp;
};

/* ---------------------------------------------------------------------- */

class GSMDampingNone : public GSMDamping {
 public:
  GSMDampingNone(class GranularModel *, class LAMMPS *);
  void init() override {};
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GSMDampingVelocity : public GSMDamping {
 public:
  GSMDampingVelocity(class GranularModel *, class LAMMPS *);
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GSMDampingMassVelocity : public GSMDamping {
 public:
  GSMDampingMassVelocity(class GranularModel *, class LAMMPS *);
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GSMDampingViscoelastic : public GSMDamping {
 public:
  GSMDampingViscoelastic(class GranularModel *, class LAMMPS *);
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GSMDampingTsuji : public GSMDamping {
 public:
  GSMDampingTsuji(class GranularModel *, class LAMMPS *);
  void init() override;
  double calculate_forces();
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GSM_DAMPING_H_ */
#endif /*GSM_CLASS_H_ */
