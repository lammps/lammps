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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/cs,ComputeTempCS)

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
  int groupbit_c,groupbit_s;
  int nshells;
  int firstflag;
  int maxatom;
  int cgroup,sgroup;

  int fix_dof;
  double tfactor;
  double **vint;

  char *id_fix;
  class FixStore *fix;

  void dof_compute();
  void vcm_pairs();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Option mol of compute temp/vcm requires molecular atom style

Self-explanatory.

E: Option prop of compute temp/vcm requires one set of parameters 
added by the property/atom fix

Self-explanatory.

E: Fix property/atom vector must contain only intergers to assign 
sub-ID property

Self-explanatory.

E: Specified sub-ID property does not exist or has not been created 
by the property/atom fix

Self-explanatory. Usually this means that the specified fix 
property/atom ID does not match the ID stated in the compute temp/vcm.

E: Molecule count changed in compute com/temp/molecule

Number of molecules must remain constant over time.

E: Sub-ID count changed in compute vcm/temp

Number of Sub-ID groups must remain constant over time.

W: Atom with sub-ID = 0 included in compute group

Self-explanatory. A sub-ID with value 0 will be counted as a normal sub-ID 
and not left out of by the compute treatment. Therefore a sub-ID of 0 is to 
be avoided.

E: Too many sub-ID groups for compute

Self-explanatory.

W: More than 2 atoms specified with the same sub-ID, in the case 
of a core-shell model simulation only core and shell should share the same ID

Self-explanatory.

*/
