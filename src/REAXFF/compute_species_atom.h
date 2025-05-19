// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Labo0ratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(SPECIES/ATOM,ComputeSpeciesAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPECIES_ATOM_H
#define LMP_COMPUTE_SPECIES_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSpeciesAtom : public Compute {
 public:
  ComputeSpeciesAtom(class LAMMPS *, int, char **);
  ~ComputeSpeciesAtom() override;
  void init() override {}
  void compute_peratom() override;
  double memory_usage() override;

 protected:
  int nvalues;
  int nmax;
  double *buf;
  double *vbuf;

  typedef void (ComputeSpeciesAtom::*FnPtrPack)(int);
  FnPtrPack *pack_choice;

  void pack_q(int);
  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);

  void pack_abo01(int);
  void pack_abo02(int);
  void pack_abo03(int);
  void pack_abo04(int);
  void pack_abo05(int);
  void pack_abo06(int);
  void pack_abo07(int);
  void pack_abo08(int);
  void pack_abo09(int);
  void pack_abo10(int);
  void pack_abo11(int);
  void pack_abo12(int);
  void pack_abo13(int);
  void pack_abo14(int);
  void pack_abo15(int);
  void pack_abo16(int);
  void pack_abo17(int);
  void pack_abo18(int);
  void pack_abo19(int);
  void pack_abo20(int);
  void pack_abo21(int);
  void pack_abo22(int);
  void pack_abo23(int);
  void pack_abo24(int);

  unsigned int get_property_type(FnPtrPack) const;

  class PairReaxFF *reaxff;
};

}    // namespace LAMMPS_NS

#endif
#endif
