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

#ifndef KOKKOS_LMP_H
#define KOKKOS_LMP_H

#include "pointers.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

class KokkosLMP : protected Pointers {
 public:
  int kokkos_exists;
  int neighflag;
  int neighflag_qeq;
  int neighflag_qeq_set;
  int exchange_comm_classic;
  int forward_comm_classic;
  int reverse_comm_classic;
  int exchange_comm_on_host;
  int forward_comm_on_host;
  int reverse_comm_on_host;
  int num_threads,ngpu;
  int numa;
  int auto_sync;
  int gpu_direct_flag;
  int team_flag;

  KokkosLMP(class LAMMPS *, int, char **);
  ~KokkosLMP();
  void accelerator(int, char **);
  int neigh_count(int);

  template<class DeviceType>
  int need_dup()
  {
    int value = 0;

    if (neighflag == HALFTHREAD)
      value = NeedDup<HALFTHREAD,DeviceType>::value;

    return value;
  }

 private:
  static void my_signal_handler(int);
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid Kokkos command-line args

Self-explanatory.  See Section 2.7 of the manual for details.

E: GPUs are requested but Kokkos has not been compiled for CUDA

Recompile Kokkos with CUDA support to use GPUs.

E: Kokkos has been compiled for CUDA but no GPUs are requested

One or more GPUs must be used when Kokkos is compiled for CUDA.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Must use Kokkos half/thread or full neighbor list with threads or GPUs

Using Kokkos half-neighbor lists with threading is not allowed.

E: Must use KOKKOS package option 'neigh full' with 'team on'

The 'team on' option requires a full neighbor list

*/
