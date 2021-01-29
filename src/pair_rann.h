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
   	   	   	   	   	     Doyl Dickel (MSU) doyl@me.msstate.edu
    ----------------------------------------------------------------------*/
/*
“The research described and the resulting data presented herein, unless
otherwise noted, was funded under PE 0602784A, Project T53 "Military
Engineering Applied Research", Task 002 under Contract No. W56HZV-17-C-0095,
managed by the U.S. Army Combat Capabilities Development Command (CCDC) and
the Engineer Research and Development Center (ERDC).  The work described in
this document was conducted at CAVS, MSU.  Permission was granted by ERDC
to publish this information. Any opinions, findings and conclusions or
recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the United States Army.​”

DISTRIBUTION A. Approved for public release; distribution unlimited. OPSEC#4918
 */

#ifdef PAIR_CLASS

PairStyle(rann,PairRANN)

#else

#ifndef LMP_PAIR_RANN
#define LMP_PAIR_RANN

#define MAXLINE 1024

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "pair.h"
#include <map>
#include <string>

namespace LAMMPS_NS {

class PairRANN : public Pair {
 public:

  //inherited functions
  PairRANN(class LAMMPS *);
  ~PairRANN();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void init_list(int , NeighList *);
  void errorf(const char*);
  int count_words(char *);
  //black magic for modular fingerprints and activations
  class Activation ***activation;
  class Fingerprint ***fingerprints;
  typedef Fingerprint *(*FingerprintCreator)(PairRANN *);
  typedef Activation *(*ActivationCreator)(PairRANN *);
  typedef std::map<std::string,FingerprintCreator> FingerprintCreatorMap;
  typedef std::map<std::string,ActivationCreator> ActivationCreatorMap;
  FingerprintCreatorMap *fingerprint_map;
  ActivationCreatorMap *activation_map;
  Fingerprint * create_fingerprint(const char *);
  Activation * create_activation(const char *);

  //global variables
  int nelements;                // # of elements (distinct from LAMMPS atom types since multiple atom types can be mapped to one element)
  int nelementsp;				// nelements+1
  char **elements;              // names of elements
  char **elementsp;				// names of elements with "all" appended as the last "element"
  double *mass;                 // mass of each element
  double cutmax;				// max radial distance for neighbor lists
  int *map;                     // mapping from atom types to elements
  int *fingerprintcount;		// static variable used in initialization
  int *fingerprintlength;       // # of input neurons defined by fingerprints of each element.
  int *fingerprintperelement;   // # of fingerprints for each element
  bool doscreen;//screening is calculated if any defined fingerprint uses it
  bool allscreen;//all fingerprints use screening so screened neighbors can be completely ignored
  bool dospin;
  int res;//Resolution of function tables for cubic interpolation.
  int memguess;
	double *screening_min;
	double *screening_max;
	bool **weightdefined;
	bool **biasdefined;

	struct Simulation{
	  	int *id;
		bool forces;
		bool spins;
		double **x;
		double **f;
		double **s;
		double box[3][3];
		double origin[3];
		double **features;
		double **dfx;
		double **dfy;
		double **dfz;
		double **dsx;
		double **dsy;
		double **dsz;
		int *ilist,*numneigh,**firstneigh,*type,inum,gnum;
	};
	Simulation *sims;

  struct NNarchitecture{
	  int layers;
	  int *dimensions;//vector of length layers with entries for neurons per layer
	  double **Weights;
	  double **Biases;
	  int *activations;//unused
	  int maxlayer;//longest layer (for memory allocation)
  };
  NNarchitecture *net;//array of networks, 1 for each element.

 private:
  template <typename T> static Fingerprint *fingerprint_creator(PairRANN *);
  template <typename T> static Activation *activation_creator(PairRANN *);
  //new functions
  void allocate(char **);//called after reading element list, but before reading the rest of the potential
  void read_file(char *);//read potential file
  void read_atom_types(char **,char *);
  void read_mass(char **,char *);
  void read_fpe(char**,char *);//fingerprints per element. Count total fingerprints defined for each 1st element in element combinations
  void read_fingerprints(char **,int,char *);
  void read_fingerprint_constants(char **,int,char *);
  void read_network_layers(char**,char*);//include input and output layer (hidden layers + 2)
  void read_layer_size(char**,char*);
  void read_weight(char**,char*,FILE*);//weights should be formatted as properly shaped matrices
  void read_bias(char**,char*,FILE*);//biases should be formatted as properly shaped vectors
  void read_activation_functions(char**,char*);
  void read_screening(char**,int, char*);
  bool check_potential();//after finishing reading potential file
  void propagateforward(double *,double *,double *,double *,double *,double **,double **,int,int,int*);//called by compute to get force and energy
  void propagateforwardspin(double *,double *,double *,double *,double *,double *,double *,double *,double **,double **,double**,int,int,int*);//called by compute to get force and energy
  void screen(double*,double*,double*,double*,double*,double*,double*,bool*,int,int,double*,double*,double*,int *,int);
  void cull_neighbor_list(double *,double *,double *,int *,int *,int *,int,int);
  void screen_neighbor_list(double *,double *,double *,int *,int *,int *,int,int,bool*,double*,double*,double*,double*,double*,double*,double*);
};

}

#endif
#endif



