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

#ifdef FIX_CLASS
// clang-format off
FixStyle(pimd/nvt/bosonic,FixPIMDBNVT);
// clang-format on
#else

#ifndef LMP_FIX_PIMD_NVT_BOSONIC_H
#define LMP_FIX_PIMD_NVT_BOSONIC_H

#include "fix_pimd_nvt.h"

namespace LAMMPS_NS {

class FixPIMDBNVT : public FixPIMDNVT {
 public:
  FixPIMDBNVT(class LAMMPS *, int, char **);
  ~FixPIMDBNVT() override;

 protected:
  bool parse_bosonic_keyword(int, char **, int &);
  void spring_force() override;
  void compute_spring_energy() override;
  void compute_t_prim() override;

 private:
  class BosonicExchange *bosonic_exchange;
  double **f_tag_order;
  int nbosons;
};

}    // namespace LAMMPS_NS

#endif
#endif
