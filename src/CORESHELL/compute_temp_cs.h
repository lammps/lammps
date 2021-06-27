/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/cs,ComputeTempCS);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_CS_H
#define LMP_COMPUTE_TEMP_CS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempCS : public Compute {
 public:
  ComputeTempCS(class LAMMPS *, int, char **);
  ~ComputeTempCS();
  void init();
  void setup();
  double compute_scalar();
  void compute_vector();
  double memory_usage();

  void remove_bias(int, double *);
  void remove_bias_all();
  void reapply_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 private:
  int groupbit_c, groupbit_s;
  int nshells;
  int firstflag;
  int maxatom;
  int cgroup, sgroup;

  double tfactor;
  double **vint;

  char *id_fix;
  class FixStore *fix;

  void dof_compute();
  void vcm_pairs();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute temp/cs used when bonds are not allowed

This compute only works on pairs of bonded particles.

E: Cannot find specified group ID for core particles

Self-explanatory.

E: Cannot find specified group ID for shell particles

Self-explanatory.

E: Compute temp/cs requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Number of core atoms != number of shell atoms

There must be a one-to-one pairing of core and shell atoms.

E: Core/shell partner atom not found

Could not find one of the atoms in the bond pair.

E: Core/shell partners were not all found

Could not find or more atoms in the bond pairs.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
