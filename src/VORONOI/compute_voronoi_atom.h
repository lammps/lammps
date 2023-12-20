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
ComputeStyle(voronoi/atom,ComputeVoronoi);
// clang-format on
#else

#ifndef LMP_COMPUTE_VORONOI_H
#define LMP_COMPUTE_VORONOI_H

#include "compute.h"

namespace voro {
class container;
class container_poly;
class voronoicell_neighbor;
}    // namespace voro

namespace LAMMPS_NS {

class ComputeVoronoi : public Compute {
 public:
  ComputeVoronoi(class LAMMPS *, int, char **);
  ~ComputeVoronoi() override;
  void init() override;
  void compute_peratom() override;
  void compute_vector() override;
  void compute_local() override;
  double memory_usage() override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  voro::container *con_mono;
  voro::container_poly *con_poly;

  void buildCells();
  void checkOccupation();
  void loopCells();
  void processCell(voro::voronoicell_neighbor &, int);

  int nmax, rmax, maxedge, sgroupbit;
  char *radstr;
  double fthresh, ethresh;
  double **voro;
  double *edge, *sendvector, *rfield;
  enum { VOROSURF_NONE, VOROSURF_ALL, VOROSURF_GROUP } surface;
  bool onlyGroup, occupation;

  tagint *tags, oldmaxtag;
  int *occvec, *sendocc, *lroot, *lnext, lmax, oldnatoms, oldnall;
  int faces_flag, nfaces, nfacesmax;
  double **faces;
};

}    // namespace LAMMPS_NS

#endif
#endif
