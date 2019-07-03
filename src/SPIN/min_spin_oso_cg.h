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

MinimizeStyle(spin/oso_cg, MinSpinOSO_CG)

#else

#ifndef LMP_MIN_SPIN_OSO_CG_H
#define LMP_MIN_SPIN_OSO_CG_H

#include "min.h"

namespace LAMMPS_NS {

class MinSpinOSO_CG : public Min {

public:
    MinSpinOSO_CG(class LAMMPS *);
    virtual ~MinSpinOSO_CG();
    void init();
    void setup_style();
    int modify_param(int, char **);
    void reset_vectors();
    int iterate(int);
    double evaluate_dt();
    void advance_spins();
    double fmnorm_sqr();
    void calc_gradient(double);
    void calc_search_direction();

private:
    // global and spin timesteps

    int nlocal_max;		// max value of nlocal (for size of lists)
    double dt;
    double dts;

    double alpha_damp;		// damping for spin minimization
    double discrete_factor;	// factor for spin timestep evaluation

    double *spvec;		// variables for atomic dof, as 1d vector
    double *fmvec;		// variables for atomic dof, as 1d vector

    double *g_old;  		// gradient vector at previous iteration
    double *g_cur;  		// current gradient vector
    double *p_s;  		// search direction vector
    int local_iter;  // number of times we call search_direction

    void vm3(const double *, const double *, double *);
    void rodrigues_rotation(const double *, double *);

    bigint last_negative;
};

}

#endif
#endif
