/* -*- c++ -*- ----------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov
 
 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.
 
 See the README file in the top-level LAMMPS directory.
 
Contributing authors: Steven E. Strong and Joel D. Eaves
 Joel.Eaves@Colorado.edu
 ------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(gaussFlow,FixGaussFlow)

#else

#ifndef LMP_FIX_GAUSSFLOW_H
#define LMP_FIX_GAUSSFLOW_H

#include "fix.h"

namespace LAMMPS_NS {
    
    class FixGaussFlow : public Fix {
    public:
        FixGaussFlow(class LAMMPS *, int, char **);
        int setmask();
        double compute_scalar();
        double compute_vector(int n);
        void post_force(int);
        void setup(int);
        
    protected:
        int dimension;
        bool flow[3];       //flag if each direction is conserved
        double a_app[3];    //applied acceleration
        double mTot;        //total mass of constrained group
        double f_tot[3];    //total applied force
        double peAdded;     //total added energy per proc
        double pe_tot;      //total added energy
        bool force_flag;    //if force has been computed this timestep already
        double dt;          //timestep
        bool workflag;      //if calculate work done by fix
        
    };
    
}

#endif
#endif
