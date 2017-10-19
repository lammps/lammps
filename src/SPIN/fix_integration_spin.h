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

  // compute and advance single spin
  void ComputeInteractionsSpin(int);   
  void AdvanceSingleSpin(int, double, double **, double **);

  // sectoring operations
  void sectoring(); 
  int coords2sector(double *);

 protected:
  int extra, mpi_flag;

  // vel., force, and spin timesteps
  double dtv,dtf,dts;
  
  // mag. interaction flags
  int magpair_flag;
  int soc_flag;
  int exch_flag;
  int magforce_flag;
  int zeeman_flag, aniso_flag;
  int maglangevin_flag;
  int tdamp_flag, temp_flag;

  // pointers to interaction classes
  class PairHybrid *lockhybrid; 
  class PairSpin *lockpairspin;
  class PairSpinExchange *lockpairspinexchange;
  class PairSpinSocNeel *lockpairspinsocneel;
  class FixForceSpin *lockforcespin;
  class FixLangevinSpin *locklangevinspin; 

  // temporary variables
  double *xi, *rij;
  double *spi, *spj;
  double *fmi, *fmj; 
 
  // sectoring variables
  int nsectors;
  int *sec, *seci;
  double *rsec;

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
