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
#include "activation_linear.h"
#include "activation.h"


using namespace LAMMPS_NS;

Activation_linear::Activation_linear(PairRANN *pair) : Activation(pair){
	empty = false;
	style = "linear";
}

double Activation_linear::activation_function(double A)
{
	return A;
}

double Activation_linear::dactivation_function(double A)
{
	return 1.0;
}

double Activation_linear::ddactivation_function(double){
	return 0.0;
}
