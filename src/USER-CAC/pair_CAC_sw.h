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

#ifdef PAIR_CLASS

PairStyle(CAC/sw,PairCACSW)

#else

#ifndef LMP_PAIR_SW_CAC_H
#define LMP_PAIR_SW_CAC_H

//#include "asa_user.h"
#include "pair.h"
#include "pair_CAC.h"
#include <stdint.h>


namespace LAMMPS_NS {

class PairCACSW : public PairCAC {
 public:
  PairCACSW(class LAMMPS *);
  virtual ~PairCACSW();

  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
 // virtual double single(int, int, int, int, double, double, double, double &);
  //double LJEOS(int);
  /*
  struct asa_objective_struct;
  typedef asa_objective_struct asa_objective;

  struct asacg_parm_struct;
  typedef asacg_parm_struct asacg_parm;

  struct asa_parm_struct;
  typedef asa_parm_struct asa_parm;
  */
	struct Param {
		double epsilon, sigma;
		double littlea, lambda, gamma, costheta;
		double biga, bigb;
		double powerp, powerq;
		double tol;
		double cut, cutsq;
		double sigma_gamma, lambda_epsilon, lambda_epsilon2;
		double c1, c2, c3, c4, c5, c6;
		int ielement, jelement, kelement;
	};

 protected:
	
  


  
    int neigh_nodes_per_element;

   




	//potential params
	
	int nelements;                // # of unique elements
	char **elements;              // names of unique elements
	int ***elem2param;            // mapping from element triplets to parameters
	int *map;                     // mapping from atom types to elements
	int nparams;                  // # of stored parameter sets
	int maxparam;                 // max # of parameter sets
	Param *params;                // parameter set for an I-J-K interaction

    
  double **cut;
  double **inner_neighbor_coords;
  double **outer_neighbor_coords;
  int *inner_neighbor_types;
  int *outer_neighbor_types;
 
	

	

  void allocate();
  void read_file(char *);
  virtual void setup_params();
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, Param *, Param *, double, double, double *, double *,
	  double *, double *, int, double &);
  //double density_map(double);
  
  
  

  
 
  void force_densities(int, double, double, double, double, double
	  &fx, double &fy, double &fz);
  
};

}

#endif
#endif
