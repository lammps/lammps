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
         GSMTangentialNone,
         TANGENTIAL);

GSMStyle(linear_nohistory,
         GSMTangentialLinearNoHistory,
         TANGENTIAL);

GSMStyle(linear_history,
         GSMTangentialLinearHistory,
         TANGENTIAL);

GSMStyle(linear_history_classic,
         GSMTangentialLinearHistoryClassic,
         TANGENTIAL);

GSMStyle(mindlin_classic,
         GSMTangentialMindlinClassic,
         TANGENTIAL);

GSMStyle(mindlin,
         GSMTangentialMindlin,
         TANGENTIAL);

GSMStyle(mindlin/force,
         GSMTangentialMindlinForce,
         TANGENTIAL);

GSMStyle(mindlin_rescale,
         GSMTangentialMindlinRescale,
         TANGENTIAL);

GSMStyle(mindlin_rescale/force,
         GSMTangentialMindlinRescaleForce,
         TANGENTIAL);
// clang-format on
#else

#ifndef GSM_TANGENTIAL_H_
#define GSM_TANGENTIAL_H_

#include "gsm.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GSMTangential : public GSM {
 public:
  GSMTangential(class GranularModel *, class LAMMPS *);
  virtual ~GSMTangential() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual void calculate_forces() = 0;
  double k, damp, mu; // Used by Marshall twisting model
};

/* ---------------------------------------------------------------------- */

class GSMTangentialNone : public GSMTangential {
 public:
  GSMTangentialNone(class GranularModel *, class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class GSMTangentialLinearNoHistory : public GSMTangential {
 public:
  GSMTangentialLinearNoHistory(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double xt;
};

/* ---------------------------------------------------------------------- */

class GSMTangentialLinearHistory : public GSMTangential {
 public:
  GSMTangentialLinearHistory(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double xt;
};

/* ---------------------------------------------------------------------- */

class GSMTangentialLinearHistoryClassic : public GSMTangentialLinearHistory {
 public:
  GSMTangentialLinearHistoryClassic(class GranularModel *, class LAMMPS *);
  void calculate_forces();
 protected:
  double xt;
};

/* ---------------------------------------------------------------------- */

class GSMTangentialMindlinClassic : public GSMTangentialLinearHistoryClassic {
 public:
  GSMTangentialMindlinClassic(class GranularModel *, class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class GSMTangentialMindlin : public GSMTangential {
 public:
  GSMTangentialMindlin(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
  void calculate_forces();
 protected:
  int mindlin_rescale, mindlin_force;
  double xt;
};

/* ---------------------------------------------------------------------- */

class GSMTangentialMindlinForce : public GSMTangentialMindlin {
 public:
  GSMTangentialMindlinForce(class GranularModel *, class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class GSMTangentialMindlinRescale : public GSMTangentialMindlin {
 public:
  GSMTangentialMindlinRescale(class GranularModel *, class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class GSMTangentialMindlinRescaleForce : public GSMTangentialMindlin {
 public:
  GSMTangentialMindlinRescaleForce(class GranularModel *, class LAMMPS *);
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GSM_TANGENTIAL_H_ */
#endif /*GSM_CLASS_H_ */
