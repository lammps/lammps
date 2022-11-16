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

#ifdef GRAN_SUB_MOD_CLASS
// clang-format off
GranSubModStyle(none,
         GranSubModTangentialNone,
         TANGENTIAL);

GranSubModStyle(linear_nohistory,
         GranSubModTangentialLinearNoHistory,
         TANGENTIAL);

GranSubModStyle(linear_history,
         GranSubModTangentialLinearHistory,
         TANGENTIAL);

GranSubModStyle(linear_history_classic,
         GranSubModTangentialLinearHistoryClassic,
         TANGENTIAL);

GranSubModStyle(mindlin_classic,
         GranSubModTangentialMindlinClassic,
         TANGENTIAL);

GranSubModStyle(mindlin,
         GranSubModTangentialMindlin,
         TANGENTIAL);

GranSubModStyle(mindlin/force,
         GranSubModTangentialMindlinForce,
         TANGENTIAL);

GranSubModStyle(mindlin_rescale,
         GranSubModTangentialMindlinRescale,
         TANGENTIAL);

GranSubModStyle(mindlin_rescale/force,
         GranSubModTangentialMindlinRescaleForce,
         TANGENTIAL);
// clang-format on
#else

#ifndef GRAN_SUB_MOD_TANGENTIAL_H
#define GRAN_SUB_MOD_TANGENTIAL_H

#include "gran_sub_mod.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GranSubModTangential : public GranSubMod {
 public:
  GranSubModTangential(class GranularModel *, class LAMMPS *);
  virtual ~GranSubModTangential() {};
  virtual void calculate_forces() = 0;
  double k, damp, mu; // Used by Marshall twisting model
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialNone : public GranSubModTangential {
 public:
  GranSubModTangentialNone(class GranularModel *, class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialLinearNoHistory : public GranSubModTangential {
 public:
  GranSubModTangentialLinearNoHistory(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double xt;
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialLinearHistory : public GranSubModTangential {
 public:
  GranSubModTangentialLinearHistory(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double xt;
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialLinearHistoryClassic : public GranSubModTangentialLinearHistory {
 public:
  GranSubModTangentialLinearHistoryClassic(class GranularModel *, class LAMMPS *);
  void calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialMindlinClassic : public GranSubModTangentialLinearHistoryClassic {
 public:
  GranSubModTangentialMindlinClassic(class GranularModel *, class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialMindlin : public GranSubModTangential {
 public:
  GranSubModTangentialMindlin(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
  void calculate_forces();
 protected:
  int mindlin_rescale, mindlin_force;
  double xt;
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialMindlinForce : public GranSubModTangentialMindlin {
 public:
  GranSubModTangentialMindlinForce(class GranularModel *, class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialMindlinRescale : public GranSubModTangentialMindlin {
 public:
  GranSubModTangentialMindlinRescale(class GranularModel *, class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class GranSubModTangentialMindlinRescaleForce : public GranSubModTangentialMindlin {
 public:
  GranSubModTangentialMindlinRescaleForce(class GranularModel *, class LAMMPS *);
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GRAN_SUB_MOD_TANGENTIAL_H */
#endif /*GRAN_SUB_MOD_CLASS_H */
