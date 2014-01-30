/* -------------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/split,VerletSplit)

#else

#ifndef LMP_VERLET_SPLIT_H
#define LMP_VERLET_SPLIT_H

#include "verlet.h"

namespace LAMMPS_NS {

class VerletSplit : public Verlet {
 public:
  VerletSplit(class LAMMPS *, int, char **);
  ~VerletSplit();
  void init();
  void setup();
  void setup_minimal(int);
  void run(int);
  bigint memory_usage();

 private:
  int master;                        // 1 if an Rspace proc, 0 if Kspace
  int me_block;                      // proc ID within Rspace/Kspace block
  int ratio;                         // ratio of Rspace procs to Kspace procs
  int *qsize,*qdisp,*xsize,*xdisp;   // MPI gather/scatter params for block comm
  MPI_Comm block;                    // communicator within one block
  int tip4p_flag;                    // 1 if PPPM/tip4p so do extra comm

  double **f_kspace;                 // copy of Kspace forces on Rspace procs
  int maxatom;

  void rk_setup();
  void r2k_comm();
  void k2r_comm();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Verlet/split requires 2 partitions

See the -partition command-line switch.

E: Verlet/split requires Rspace partition size be multiple of Kspace partition size

This is so there is an equal number of Rspace processors for every
Kspace processor.

E: Verlet/split requires Rspace partition layout be multiple of Kspace partition layout in each dim

This is controlled by the processors command.

W: No Kspace calculation with verlet/split

The 2nd partition performs a kspace calculation so the kspace_style
command must be used.

E: Verlet/split does not yet support TIP4P

This is a current limitation.

*/
