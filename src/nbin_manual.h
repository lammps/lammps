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

#ifndef LMP_NBIN_MANUAL_H
#define LMP_NBIN_MANUAL_H

#include "nbin.h"

namespace LAMMPS_NS {

class NBinManual : public NBin {
 public:
  NBinManual(class LAMMPS *);

  void bin_atoms_setup(int) override;
  void setup_bins(int) override;
  void bin_atoms() override;
  double memory_usage() override;

  void bin_custom_setup(double **, int);
  void bin_custom(double **, int);
  void assign_neighbor_info(double);

  int *binhead_custom;     // index of first custom element in each bin
  int *bins_custom;        // index of next custom element in same bin
  int *custom2bin;         // bin assignment for each custom element
 protected:

  int maxcustom;    // size of bins array
  int maxbin_custom;     // size of binhead array
};

}    // namespace LAMMPS_NS

#endif
