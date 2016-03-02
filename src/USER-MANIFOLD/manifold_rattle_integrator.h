/* ----------------------------------------------------------------------
   Lammps - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the USER-MANIFOLD package.

   Copyright (2013-2015) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License v.2

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterised by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */

#ifndef LMP_MANIFOLD_RATTLE_INTEGRATOR_H
#define LMP_MANIFOLD_RATTLE_INTEGRATOR_H

#include "manifold.h"
#include "pointers.h"

namespace LAMMPS_NS {

  // A general util function for handling args.
  void shift_args_back( int argi, char ***arg_ptr, int delta, int &narg );

  class ManifoldRattleIntegratorBase : public Pointers
  {
   public:
	  
    struct statistics {

      statistics() : x_iters(0), v_iters(0), x_iters_per_atom(0),
                     v_iters_per_atom(0), natoms(0), dofs_removed(0),
                     last_out(0) {}
      double x_iters, v_iters;
      double x_iters_per_atom;
      double v_iters_per_atom;
      int natoms;
      int dofs_removed;
      bigint last_out;
    };

    ManifoldRattleIntegratorBase( LAMMPS*, int &, char *** );
    ~ManifoldRattleIntegratorBase();
    
    virtual void init() = 0;
    virtual void update_params() = 0;
    virtual double memory_usage() = 0;

    
    virtual void nve_x_rattle(int, int);
    virtual void nve_v_rattle(int, int);
    virtual int nparams() = 0;

    void init_manifold( int, char **, manifold * );

    void print_stats( const char * );
    void reset_dt();

    int parse_arg( int , int &, char *** );
    int is_legal_keyword( const char * );

    int dof(int,int);

    int nevery;

   protected:
    double dtv, dtf;
    double tolerance;
    int max_iter;

    char **tstrs;
    int *tvars;
    int *tstyle;
    int *is_var;

    statistics stats;
    int update_style;
    int nvars;

   private:
    virtual void rattle_manifold_x( double*, double*, double*, double, double, int = 0 ) = 0;
    virtual void rattle_manifold_v( double*, double*, double*, double, int = 0 ) = 0;
  };


  class ManifoldRattleIntegrator : public ManifoldRattleIntegratorBase
  {
   public:
    ManifoldRattleIntegrator(LAMMPS *, int&, char *** );
    virtual ~ManifoldRattleIntegrator();

    virtual void init();
    virtual void update_params();
    virtual double memory_usage();

    virtual int nparams(){ return ptr_m->nparams(); }

	  
   private:
    void update_var_params();
    virtual void rattle_manifold_x( double*, double*, double*, double, double, int = 0 );
    virtual void rattle_manifold_v( double*, double*, double*, double, int = 0);

    manifold *ptr_m;

  };

}


#endif // LMP_MANIFOLD_RATTLE_INTEGRATOR_H
