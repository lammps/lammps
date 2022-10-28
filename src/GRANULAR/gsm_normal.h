/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef GSM_CLASS
// clang-format off
GSMStyle(none,
         GSMNormalNone,
         NORMAL);

GSMStyle(hooke,
         GSMNormalHooke,
         NORMAL);

GSMStyle(hertz,
         GSMNormalHertz,
         NORMAL);

GSMStyle(hertz/material,
         GSMNormalHertzMaterial,
         NORMAL);

GSMStyle(dmt,
         GSMNormalDMT,
         NORMAL);

GSMStyle(jkr,
         GSMNormalJKR,
         NORMAL);
// clang-format on
#else

#ifndef GSM_NORMAL_H_
#define GSM_NORMAL_H_

#include "gsm.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GSMNormal : public GSM {
 public:
  GSMNormal(class GranularModel *, class LAMMPS *);
  ~GSMNormal() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual bool touch();
  virtual double pulloff_distance(double, double);
  virtual double calculate_area();
  virtual void set_knfac() = 0;
  virtual double calculate_forces() = 0;
  virtual void set_fncrit();
  double damp;  // Vestigial argument needed by damping
  double Emod, poiss;
  double Fncrit, Fne, knfac;
  int material_properties;
};

/* ---------------------------------------------------------------------- */

class GSMNormalNone : public GSMNormal {
 public:
  GSMNormalNone(class GranularModel *, class LAMMPS *);
  void set_knfac() {};
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GSMNormalHooke : public GSMNormal {
 public:
  GSMNormalHooke(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void set_knfac();
  double calculate_forces();
 protected:
  double k;
};

/* ---------------------------------------------------------------------- */

class GSMNormalHertz : public GSMNormal {
 public:
  GSMNormalHertz(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void set_knfac();
  double calculate_forces();
 protected:
  double k;
};

/* ---------------------------------------------------------------------- */

class GSMNormalHertzMaterial : public GSMNormalHertz {
 public:
  GSMNormalHertzMaterial(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
};

/* ---------------------------------------------------------------------- */

class GSMNormalDMT : public GSMNormal {
 public:
  GSMNormalDMT(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
  void set_knfac();
  double calculate_forces();
  void set_fncrit() override;
 protected:
  double k, cohesion;
  double F_pulloff;
};

/* ---------------------------------------------------------------------- */

class GSMNormalJKR : public GSMNormal {
 public:
  GSMNormalJKR(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
  bool touch() override;
  double pulloff_distance(double, double) override;
  double calculate_area() override;
  void set_knfac();
  double calculate_forces();
  void set_fncrit() override;
 protected:
  double k, cohesion;
  double Emix, F_pulloff;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GSM_NORMAL_H_ */
#endif /*GSM_CLASS_H_ */
