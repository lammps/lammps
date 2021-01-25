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

#ifdef FINGERPRINT_CLASS

FingerprintStyle(bondscreened,Fingerprint_bondscreened)

#else

#ifndef FINGERPRINT_BONDSCREENED_H_
#define FINGERPRINT_BONDSCREENED_H_

#include "fingerprint.h"

namespace LAMMPS_NS {

class Fingerprint_bondscreened : public Fingerprint {
 public:
	Fingerprint_bondscreened(PairRANN *);
	~Fingerprint_bondscreened();
	bool parse_values(char*,char*);
	void write_values(FILE *);
	void init(int*,int);
	void allocate();
	void compute_fingerprint(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,bool*,int,int,double*,double*,double*,int*,int,int*);
	void do3bodyfeatureset_doubleneighborloop(double *,double *,double *,double *,double*,double*,double*,double*,double*,double*,double*,bool*,int,int,double*,double*,double*,int*,int,int*);
	void do3bodyfeatureset_singleneighborloop(double *,double *,double *,double *,double*,double*,double*,double*,double*,double*,double*,bool*,int,int,double*,double*,double*,int*,int,int*);
	int factorial(int);
	void generate_exp_cut_table();
	void generate_coefficients();
	int get_length();

	double *expcuttable;
    double *dfctable;
    double dr;
    double *alpha_k;
    double re;
    int **coeff;
    int **coeffx;
    int **coeffy;
    int **coeffz;
    int k;
    int m;
	int **Mf;

};


}


#endif
#endif /* FINGERPRINT_BOND_H_ */
