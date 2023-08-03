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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(hma,ComputeHMA);
// clang-format on
#else

#ifndef LMP_COMPUTE_HMA_H
#define LMP_COMPUTE_HMA_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHMA : public Compute {
 public:
  ComputeHMA(class LAMMPS *, int, char **);
  ~ComputeHMA() override;
  void setup() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_vector() override;
  void set_arrays(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;

 private:
  int nmax;
  int atomsingroup;
  char *id_fix;
  char *id_temp;
  double finaltemp;
  class FixStoreAtom *fix;
  double boltz, nktv2p, inv_volume;
  double deltaPcap;
  double virial_compute(int);
  static double sumVirial(int n, double *v)
  {
    double x = 0;
    for (int i = 0; i < n; i++) x += v[i];
    return x;
  }
  int computeU, computeP, computeCv;
  class NeighList *list;    // half neighbor list
  double **deltaR;
  int returnAnharmonic;
  double uLat, pLat;
};

}    // namespace LAMMPS_NS

#endif
#endif
