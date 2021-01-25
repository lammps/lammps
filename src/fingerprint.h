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

/*  ----------------------------------------------------------------------
   Contributing author: Christopher Barrett (MSU) barrett@me.msstate.edu
    ----------------------------------------------------------------------*/

#ifndef FINGERPRINT_H_
#define FINGERPRINT_H_

#include "pair_rann.h"

namespace LAMMPS_NS {

	class Fingerprint {
	public:
		Fingerprint(PairRANN *);
		virtual ~Fingerprint();
		virtual bool parse_values(char*,char*);
		virtual void write_values(FILE *);
		virtual void init(int*,int);
		virtual void allocate();
		void init_screen(int);
		virtual void compute_fingerprint(double*,double*,double*,double*,int,int,double*,double*,double*,int*,int,int*);//no screen,no spin
		virtual void compute_fingerprint(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,bool*,int,int,double*,double*,double*,int*,int,int*);//screen
		virtual void compute_fingerprint(double*,double*,double*,double*,double*,double*,double*,int,int,double*,double*,double*,int*,int,int*);//spin
		virtual void compute_fingerprint(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,bool*,int,int,double*,double*,double*,int*,int,int*);//spin,screen
		virtual int get_length();
		virtual double cutofffunction(double,double, double);
		virtual void generate_rinvssqrttable();
		bool spin;
		bool screen;
		int n_body_type;//i-j vs. i-j-k vs. i-j-k-l, etc.
		bool empty;
		bool fullydefined;
		int startingneuron;
	    int id;//based on ordering of fingerprints listed for i-j in potential file
	    const char *style;
	    int* atomtypes;
	    double *rinvsqrttable;
	    double rc;
	    PairRANN *pair;
	};

}


#endif /* FINGERPRINT_H_ */
