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

#ifdef FIX_CLASS

FixStyle(integration/spin,FixIntegrationSpin)

#else

#ifndef LMP_FIX_INTEGRATION_SPIN_H
#define LMP_FIX_INTEGRATION_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixIntegrationSpin : public Fix {
	
 public:
  FixIntegrationSpin(class LAMMPS *, int, char **);
  virtual ~FixIntegrationSpin();
  int setmask();
  void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

  void ComputeInteractionsSpin(int);	// compute and advance single spin functions 
  void AdvanceSingleSpin(int, double);

  void sectoring();			// sectoring operation functions 
  int coords2sector(double *);

  void setup_pre_neighbor();
  void pre_neighbor();

 protected:
  int extra;
  int mpi_flag;			//mpi_flag =  if parallel algorithm
  int mech_flag; 		// mech_flag = 0 if spins only
  				// mech_flag = 1 if spin-lattice calc. 

  double dtv,dtf,dts;		// velocity, force, and spin timesteps
  
  int magpair_flag;		// magnetic pair flags
  int exch_flag;
  int soc_neel_flag, soc_dmi_flag;
  int me_flag;
  int magforce_flag;		// magnetic force flags
  int zeeman_flag, aniso_flag;
  int maglangevin_flag;		// magnetic langevin flags
  int tdamp_flag, temp_flag;

  // pointers to magnetic interaction classes

  class PairHybrid *lockhybrid;    
  class PairSpinExchange *lockpairspinexchange;
  class PairSpinSocNeel *lockpairspinsocneel;
  class PairSpinSocDmi *lockpairspinsocdmi;
  class PairSpinMe *lockpairspinme;
  class FixForceSpin *lockforcespin;
  class FixLangevinSpin *locklangevinspin; 

  int nsectors;			// sectoring variables
  double *rsec;

  // stacking variables for sectoring algorithm
  
  int *stack_head;	// index of first atom in backward_stacks  
  int *stack_foot;	// index of first atom in forward_stacks
  int *backward_stacks;	// index of next atom in backward stack
  int *forward_stacks;	// index of next atom in forward stack

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix integration/spin command

Self-explanatory.  Check the input script syntax and compare to the 
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix integration/spin requires spin attribute mumag

An atom/spin style with this attribute is needed.

E: Illegal sectoring operation

The number of processes does not match the size of the system. 
See the documentation of the sectoring method.

*/
