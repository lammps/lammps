/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(HMA,ComputeHMA)

#else

#ifndef LMP_COMPUTE_HMA_H
#define LMP_COMPUTE_HMA_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHMA : public Compute {
 public:
  ComputeHMA(class LAMMPS *, int, char **);
  ~ComputeHMA();
  void setup();
  void init();
  void init_list(int, class NeighList *);
  void compute_vector();
  void set_arrays(int);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 private:
  int nmax;
  int atomsingroup;
  char *id_fix;
  char *id_temp;
  double finaltemp;
  class FixStore *fix;
  double boltz, nktv2p, inv_volume;
  double deltaPcap;
  double virial_compute(int);
  static double sumVirial(int n, double* v) {
    double x = 0;
    for (int i=0; i<n; i++) x += v[i];
    return x;
  }
  int computeU, computeP, computeCv;
  class NeighList *list; // half neighbor list
  double **deltaR;
  int returnAnharmonic;
  double uLat, pLat;
};

}

#endif
#endif

 
