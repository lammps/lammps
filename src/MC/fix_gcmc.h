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

#ifdef FIX_CLASS

FixStyle(gcmc,FixGCMC)

#else

#ifndef LMP_FIX_GCMC_H
#define LMP_FIX_GCMC_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixGCMC : public Fix {
 public:
  FixGCMC(class LAMMPS *, int, char **);
  ~FixGCMC();
  int setmask();
  void init();
  void pre_exchange();
  void attempt_move();
  void attempt_deletion();
  void attempt_insertion();
  double energy(int, double *);
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int ntype,nevery,seed;
  int ncycles,nexchanges,nmcmoves;
  int ngas;           // # of gas molecules (or atoms) on all procs 
  int ngas_local;     // # of gas molecules (or atoms) on this proc 
  int ngas_before;    // # of gas molecules (or atoms) on procs < this proc
  int molflag;        // 0 = atomic, 1 = molecular system

  double nmove_attempts;   
  double nmove_successes;  
  double ndel_attempts;    
  double ndel_successes;   
  double ninsert_attempts; 
  double ninsert_successes;
  
  int nmax;
  double reservoir_temperature;
  double chemical_potential;
  double displace;
  double beta,zz,sigma,volume;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *sublo,*subhi;
  int *local_gas_list;                           
  double **cutsq;
  class Pair *pair;
 
  class RanPark *random_equal;
  class RanPark *random_unequal;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid atom type in fix GCMC command

The atom type specified in the GCMC command does not exist.

E: Cannot do GCMC on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.

W: Fix GCMC may delete atom with non-zero molecule ID

This is probably an error, since you should not delete only one atom
of a molecule. The GCMC molecule exchange feature does not yet work.

E: Fix GCMC molecule command requires atom attribute molecule

Should not choose the GCMC molecule feature if no molecules are being
simulated. The general molecule flag is off, but GCMC's molecule flag
is on.

E: Fix GCMC molecule feature does not yet work

Fix GCMC cannot (yet) be used to exchange molecules, only atoms.

E: Fix GCMC incompatible with given pair_style

Some pair_styles do not provide single-atom energies, which are needed
by fix GCMC.

*/
