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
   Contributing authors: Christopher Barrett (MSU) barrett@me.msstate.edu
    ----------------------------------------------------------------------*/

#include "fingerprint.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>


using namespace LAMMPS_NS;

Fingerprint::Fingerprint(PairRANN *pair)
{
	spin = false;
	screen = false;
	empty = true;
	fullydefined = false;
	n_body_type = 0;
	style = "empty";
	this->pair = pair;
}

Fingerprint::~Fingerprint(){

}

bool Fingerprint::parse_values(char *word,char *line1){
	return false;
}

void Fingerprint::init(int *i,int id){

}

void Fingerprint::allocate(){

}

void Fingerprint::compute_fingerprint(double *features,double *dfeaturesx,double *dfeaturesy,double * dfeaturesz, int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int *jl){

}

void Fingerprint::compute_fingerprint(double *features,double *dfeaturesx,double *dfeaturesy,double * dfeaturesz,double *Sik, double *dSikx, double*dSiky, double *dSikz, double *dSijkx, double *dSijky, double *dSijkz, bool *Bij, int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int *jl){

}

void Fingerprint::compute_fingerprint(double *features,double *dfeaturesx,double *dfeaturesy,double * dfeaturesz,double *sx, double *sy, double *sz, int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int *jl){

}

void Fingerprint::compute_fingerprint(double *features,double *dfeaturesx,double *dfeaturesy,double * dfeaturesz,double *sx, double *sy, double *sz, double *Sik, double *dSikx, double*dSiky, double *dSikz, double *dSijkx, double *dSijky, double *dSijkz, bool *Bij, int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int *jl){

}

void Fingerprint::write_values(FILE *fid){

}

int Fingerprint::get_length(){
	return 0;
}

//Smooth cutoff, goes from 1 to zero over the interval rc-dr to rc. Same as MEAM uses. Used by generateradialtable and generatexpcuttable.
double Fingerprint::cutofffunction(double r,double rc, double dr){
	double out;
	if (r < (rc -dr))out=1;
	else if (r>rc)out=0;
	else {
		out = pow(1-pow(1-(rc-r)/dr,4.0),2.0);
	}
	return out;
}

void Fingerprint::generate_rinvssqrttable(){
	int buf = 5;
	int m;
	double r1;
	double cutmax = pair->cutmax;
	int res = pair->res;
	rinvsqrttable = new double[res+buf];
	for (m=0;m<(res+buf);m++){
		r1 = cutmax*cutmax*(double)(m)/(double)(res);
		rinvsqrttable[m] = 1/sqrt(r1);
	}
}




