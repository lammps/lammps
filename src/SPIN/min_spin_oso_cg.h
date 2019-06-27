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

#ifdef MINIMIZE_CLASS

MinimizeStyle(spin_oso_cg, MinSpinOSO_CG)

#else

#ifndef LMP_MIN_SPIN_OSO_CG_H
#define LMP_MIN_SPIN_OSO_CG_H

#include "min.h"

namespace LAMMPS_NS {

class MinSpinOSO_CG : public Min {

public:
    MinSpinOSO_CG(class LAMMPS *);  //?
    ~MinSpinOSO_CG() {}  //?
    void init();
    void setup_style();
    int modify_param(int, char **);
    void reset_vectors();
    int iterate(int);
    double evaluate_dt();
    void advance_spins();
    double fmnorm_sqr();
    void calc_gradient(double);
    void calc_search_direction(int);

private:
    // global and spin timesteps

    double dt;
    double dts;

    double alpha_damp;            // damping for spin minimization
    double discrete_factor;       // factor for spin timestep evaluation

    double *spvec;               // variables for atomic dof, as 1d vector
    double *fmvec;               // variables for atomic dof, as 1d vector

    double *g_old;  // gradient vector at previous iteration
    double *g;  // gradient vector
    double *p;  // search direction vector

    bigint last_negative;
};

}

#endif
#endif
