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

/* ----------------------------------------------------------------------
   Contributing authors: 
   James Larentzos and Timothy I. Mattox (Engility Corporation)

   Martin Lisal (Institute of Chemical Process Fundamentals 
   of the Czech Academy of Sciences and J. E. Purkinje University)

   John Brennan, Joshua Moore and William Mattson (Army Research Lab)

   Please cite the related publications:
   J. P. Larentzos, J. K. Brennan, J. D. Moore, M. Lisal, W. D. Mattson,
   "Parallel implementation of isothermal and isoenergetic Dissipative
   Particle Dynamics using Shardlow-like splitting algorithms", 
   Computer Physics Communications, 2014, 185, pp 1987--1998.

   M. Lisal, J. K. Brennan, J. Bonet Avalos, "Dissipative particle dynamics
   at isothermal, isobaric, isoenergetic, and isoenthalpic conditions using
   Shardlow-like splitting algorithms", Journal of Chemical Physics, 2011,
   135, 204105.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(shardlow,FixShardlow)

#else

#ifndef LMP_FIX_SHARDLOW_H
#define LMP_FIX_SHARDLOW_H

#include "fix.h"

namespace LAMMPS_NS {

class FixShardlow : public Fix {
 public:
  FixShardlow(class LAMMPS *, int, char **);
  virtual ~FixShardlow() {}
  int setmask();
  virtual void init_list(int,class NeighList *);
  virtual void setup(int);
  virtual void setup_pre_force(int);
  virtual void initial_integrate(int);

  void setup_pre_neighbor();
  void pre_neighbor();

 protected:
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

  class PairDPDfdt *pairDPD;
  class PairDPDfdtEnergy *pairDPDE;
  double **dvSSA;

  private:
  class NeighList *list;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Must use dpd/fdt pair_style with fix shardlow

Self-explanatory.

E: Must use pair_style dpd/fdt or dpd/fdt/energy with fix shardlow

E: A deterministic integrator must be specified after fix shardlow in input 
file (e.g. fix nve or fix nph).

Self-explanatory.

E: Cannot use constant temperature integration routines with DPD

Self-explanatory.  Must use deterministic integrators such as nve or nph

E: Fix shardlow does not yet support triclinic geometries

Self-explanatory.

E:  Shardlow algorithm requires sub-domain length > 2*(rcut+skin). Either 
reduce the number of processors requested, or change the cutoff/skin

The Shardlow splitting algorithm requires the size of the sub-domain lengths 
to be are larger than twice the cutoff+skin.  Generally, the domain decomposition 
is dependant on the number of processors requested.

*/
