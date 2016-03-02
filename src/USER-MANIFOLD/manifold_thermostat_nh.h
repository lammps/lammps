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

#ifndef LMP_MANIFOLD_THERMOSTAT_NH_H
#define LMP_MANIFOLD_THERMOSTAT_NH_H

#include "pointers.h"

namespace LAMMPS_NS {

  class ManifoldThermostatNH : protected Pointers {
   public:

    ManifoldThermostatNH(LAMMPS *, int&, char ***, const char *,
                         int, int, int );
    ~ManifoldThermostatNH();

    int parse_arg( int , int &, char *** );
    int is_legal_keyword( const char * );
    double memory_usage();
    
    void print_stats( const char * );
    void init();
    void setup(int);
    void compute_temp_target();
    void nhc_temp_integrate();
    void nh_v_temp();
    void reset_dt();

    int nevery;

    void print_stuff();
    void update_t_current();

   private:
    // Fix NVT needs more params for Nose-Hoover chain
    double dthalf, dt4, dt8;

    char *id_temp;
    class Compute* temperature;
    double t_start,t_stop, t_period;
    double t_current,t_target,ke_target;
    double t_freq, drag, tdrag_factor;
    double boltz,nktv2p,tdof;
    double *eta,*eta_dot;
    double *eta_dotdot;
    double *eta_mass;
    int mtchain;
    double factor_eta;
    int which, got_temp;

    const char *fix_id;

    int igroup, groupbit;
	  
  };

}


#endif // LMP_MANIFOLD_THERMOSTAT_NH_H
