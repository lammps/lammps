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
  int coords2sector(double x[3]);

 protected:
  int extra, mpi_flag;

  double dtv,dtf,dts;		// velocity, force, and spin timesteps
  
  int magpair_flag;		// magnetic pair flags
  int soc_flag, exch_flag;
  int magforce_flag;		// magnetic force flags
  int zeeman_flag, aniso_flag;
  int maglangevin_flag;		// magnetic langevin flags
  int tdamp_flag, temp_flag;

  // pointers to magnetic interaction classes

  class PairHybrid *lockhybrid;    
  class PairSpinExchange *lockpairspinexchange;
  class PairSpinSocNeel *lockpairspinsocneel;
  class FixForceSpin *lockforcespin;
  class FixLangevinSpin *locklangevinspin; 

  int nsectors;			// sectoring variables
  double *rsec;
  int *k, **adv_list;

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
