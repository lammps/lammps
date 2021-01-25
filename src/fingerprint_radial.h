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

FingerprintStyle(radial,Fingerprint_radial)

#else

#ifndef FINGERPRINT_RADIAL_H_
#define FINGERPRINT_RADIAL_H_

#include "fingerprint.h"

namespace LAMMPS_NS {

class Fingerprint_radial : public Fingerprint {
 public:
	Fingerprint_radial(PairRANN *);
	~Fingerprint_radial();
	bool parse_values(char*,char*);
	void write_values(FILE *);
	void init(int*,int);
	void allocate();
	void compute_fingerprint(double*,double*,double*,double*,int,int,double*,double*,double*,int*,int,int*);
	int get_length();

	double *radialtable;
    double *dfctable;
    double dr;
    double *alpha;
    double re;
    int n;//highest term
    int o;//lowest term

};


}

#endif
#endif /* FINGERPRINT_RADIAL_H_ */
