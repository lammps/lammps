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

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/gpu,PPPMGPU)

#else

#ifndef LMP_PPPM_GPU_H
#define LMP_PPPM_GPU_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMGPU : public PPPM {
 public:
  PPPMGPU(class LAMMPS *, int, char **);
  virtual ~PPPMGPU();
  void init();
  void setup();
  void compute(int, int);
  int timing_1d(int, double &);
  int timing_3d(int, double &);
  double memory_usage();

 protected:
  FFT_SCALAR ***density_brick_gpu, ***vd_brick;
  bool kspace_split, im_real_space;
  int old_nlocal;
  double poisson_time;

  void brick2fft();
  virtual void poisson_ik();

  void pack_forward(int, FFT_SCALAR *, int, int *);
  void unpack_forward(int, FFT_SCALAR *, int, int *);
  void pack_reverse(int, FFT_SCALAR *, int, int *);
  void unpack_reverse(int, FFT_SCALAR *, int, int *);

  FFT_SCALAR ***create_3d_offset(int, int, int, int, int, int, const char *,
                                 FFT_SCALAR *, int);
  void destroy_3d_offset(FFT_SCALAR ***, int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot currently use pppm/gpu with fix balance.

Self-explanatory.

E: Cannot (yet) do analytic differentiation with pppm/gpu

This is a current restriction of this command.

E: Cannot use order greater than 8 with pppm/gpu.

Self-explanatory.

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor.  This is likely for one of two
reasons, both of them bad.  First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors.  This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command.  The safest settings are "delay 0
every 1 check yes".  Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

*/
