/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// dummy interface to USER-CUDA
// used when USER-CUDA is not installed

#ifndef LMP_ACCELERATOR_H
#define LMP_ACCELERATOR_H

#include "comm.h"
#include "modify.h"
#include "verlet.h"

namespace LAMMPS_NS {

class Cuda {
 public:
  int cuda_exists;
  int oncpu;
  int neighbor_decide_by_integrator;
  
  Cuda(class LAMMPS *) {cuda_exists = 0;}
  ~Cuda() {}
  void setDevice(class LAMMPS *) {}
  void accelerator(int, char **) {}
  void evsetup_eatom_vatom(int, int) {}
  void downloadAll() {}
  void uploadAll() {}
};

class CommCuda : public Comm {
 public:
 CommCuda(class LAMMPS *lmp) : Comm(lmp) {}
  ~CommCuda() {}
};

class DomainCuda : public Domain {
 public:
 DomainCuda(class LAMMPS *lmp) : Domain(lmp) {}
  ~DomainCuda() {}
};

class NeighborCuda : public Neighbor {
 public:
 NeighborCuda(class LAMMPS *lmp) : Neighbor(lmp) {}
  ~NeighborCuda() {}
};

class ModifyCuda : public Modify {
 public:
 ModifyCuda(class LAMMPS *lmp) : Modify(lmp) {}
  ~ModifyCuda() {}
};
 
class VerletCuda : public Verlet {
 public:
 VerletCuda(class LAMMPS *lmp, int narg, char **arg) : Verlet(lmp,narg,arg) {}
  ~VerletCuda() {}
};

}

#endif
