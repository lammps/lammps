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

#include <math.h>
#include "activation_sigI.h"

using namespace LAMMPS_NS;

Activation_sigI::Activation_sigI(PairRANN *pair) : Activation(pair){
	empty = false;
	style = "sigI";
}

double Activation_sigI::activation_function(double in){
	if (in>34)return in;
	return 0.1*in + 0.9*log(exp(in) + 1);
}

double Activation_sigI::dactivation_function(double in){
	if (in>34)return 1;
	return 0.1 + 0.9/(exp(in)+1)*exp(in);
}

double Activation_sigI::ddactivation_function(double in){
	if (in>34)return 0;
	return 0.9*exp(in)/(exp(in)+1)/(exp(in)+1);;
}
