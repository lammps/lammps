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

#ifdef BOND_CLASS
// clang-format off
BondStyle(bpm/peri,BondBPMPeri);
// clang-format on
#else

#ifndef LMP_BOND_BPM_PERI_H
#define LMP_BOND_BPM_PERI_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMPeri : public BondBPM {
 public:
  BondBPMPeri(class LAMMPS *);
  ~BondBPMPeri() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  // peridynamic constitutive models (internal ids; lower-case keywords in input)
  enum { PMB };

  int *model;             // per bond type: constitutive model id
  double *c;              // PMB micromodulus c = 18 K / (pi delta^4)
  double *cut;            // horizon delta (per bond type)
  double *s00, *alpha;    // critical-stretch bond-break parameters (#984 rule)

  // internal per-atom storage (fix property/atom), created on demand:
  // s0   = diagnostic critical stretch, smin = min (most compressive) stretch
  char *id_fix_property_peri;
  int index_s0, index_smin;

  void allocate();
  void store_data() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
