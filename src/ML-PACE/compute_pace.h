/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pace,ComputePACE);
// clang-format on
#else

#ifndef LMP_COMPUTE_PACE_H
#define LMP_COMPUTE_PACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePACE : public Compute {
 public:
  ComputePACE(class LAMMPS *, int, char **);
  ~ComputePACE();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int natoms, nmax, size_peratom, lastcol;
  int nvalues, yoffset, zoffset;
  int ndims_peratom, ndims_force, ndims_virial;
  double **cutsq;
  class NeighList *list;
  double **pace, **paceall;
  double **pace_peratom;
  int *map;    // map types to [0,nelements)
  int bikflag, bik_rows, dgradflag, dgrad_rows;
  double cutmax;

  Compute *c_pe;
  Compute *c_virial;
  std::string id_virial;

  void dbdotr_compute();
  struct ACECimpl *acecimpl;
};

}    // namespace LAMMPS_NS

#endif
#endif
