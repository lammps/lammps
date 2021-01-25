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

#ifndef ACTIVATION_H_
#define ACTIVATION_H_

#include "pair_rann.h"

namespace LAMMPS_NS {

	class Activation {
	public:
		Activation(class PairRANN *);
		virtual ~Activation();
		virtual double activation_function(double);
		virtual double dactivation_function(double);
		virtual double ddactivation_function(double);
		bool empty;
		const char *style;
	};
}



#endif /* ACTIVATION_H_ */
