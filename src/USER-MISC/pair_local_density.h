/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
   pair_LocalDensity written by:
   Tanmoy Sanyal and M. Scott Shell from UC Santa Barbara
   David Rosenberger: TU Darmstadt
-------------------------------------------------------------------------*/   
    

#ifdef PAIR_CLASS

PairStyle(local/density,PairLocalDensity)

#else

#ifndef LMP_PAIR_LOCAL_DENSITY_H
#define LMP_PAIR_LOCAL_DENSITY_H

#include "pair.h"


namespace LAMMPS_NS {

class PairLocalDensity : public Pair {
  public:
    PairLocalDensity(class LAMMPS *);
    virtual ~PairLocalDensity();
    virtual void compute(int, int);
    void settings(int, char **);
    virtual void coeff(int, char **);
    void init_style();
    double init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);
    
    virtual int pack_comm(int, int *, double *, int, int *);
    virtual void unpack_comm(int, int, double *);
    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);
    double memory_usage();


  protected:
    //------------------------------------------------------------------------
    //This information is read from the tabulated input file
    
    int nLD, nrho;                          // number of LD types
    int **a, **b;                           // central and neigh atom filters
    double *uppercut, *lowercut;            // upper and lower cutoffs
    double *uppercutsq, *lowercutsq;        // square of above cutoffs
    double *c0, *c2, *c4, *c6;              // coeffs for indicator function
    double *rho_min, *rho_max, *delta_rho;  // min, max & grid-size for LDs
    double **rho, **frho;                   // LD and LD function tables
    
    //------------------------------------------------------------------------
    
    double ***frho_spline; // splined LD potentials
    double cutmax;          // max cutoff for all elements
    double cutforcesq;      // square of global upper cutoff
    
    int nmax;               // max size of per-atom arrays
    double **localrho;     // per-atom LD
    double **fp;           // per-atom LD potential function derivative
    
    void allocate();
    
	// read tabulated input file
	void parse_file(char *);
    
	// convert array to spline
	void array2spline();
	
	// cubic spline interpolation
    void interpolate_cbspl(int, double, double *, double **);
};

}

#endif
#endif
