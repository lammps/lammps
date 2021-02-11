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


#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include <map>


namespace LAMMPS_NS {

  namespace RANN {
    //forward declarations
    class Activation;
    class Fingerprint;
  }

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
    void errorf(const char*,int,const char*);
    int factorial(int);
    //black magic for modular fingerprints and activations
    class RANN::Activation ***activation;
    class RANN::Fingerprint ***fingerprints;
    typedef RANN::Fingerprint *(*FingerprintCreator)(PairRANN *);
    typedef RANN::Activation *(*ActivationCreator)(PairRANN *);
    typedef std::map<std::string,FingerprintCreator> FingerprintCreatorMap;
    typedef std::map<std::string,ActivationCreator> ActivationCreatorMap;
    FingerprintCreatorMap *fingerprint_map;
    ActivationCreatorMap *activation_map;
    RANN::Fingerprint * create_fingerprint(const char *);
    RANN::Activation * create_activation(const char *);

    //global variables
    int nelements;                // # of elements (distinct from LAMMPS atom types since multiple atom types can be mapped to one element)
    int nelementsp;        // nelements+1
    char **elements;              // names of elements
    char **elementsp;        // names of elements with "all" appended as the last "element"
    double *mass;                 // mass of each element
    double cutmax;        // max radial distance for neighbor lists
    int *map;                     // mapping from atom types to elements
    int *fingerprintcount;    // static variable used in initialization
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
    int nmax1;
    int nmax2;
    int fmax;
    int fnmax;
    //memory actively written to during each compute:
    double *xn,*yn,*zn,*Sik,*dSikx,*dSiky,*dSikz,*dSijkx,*dSijky,*dSijkz,*sx,*sy,*sz,**dSijkxc,**dSijkyc,**dSijkzc,*dfeaturesx,*dfeaturesy,*dfeaturesz,*features;
    double *layer,*sum,*sum1,**dlayerx,**dlayery,**dlayerz,**dlayersumx,**dlayersumy,**dlayersumz;
    double **dsx,**dsy,**dsz,**dssumx,**dssumy,**dssumz;
    int *tn,*jl;
    bool *Bij;

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

    struct NNarchitecture{
    int layers;
    int *dimensions;//vector of length layers with entries for neurons per layer
    double **Weights;
    double **Biases;
    int *activations;//unused
    int maxlayer;//longest layer (for memory allocation)
    };

    Simulation *sims;
    NNarchitecture *net;//array of networks, 1 for each element.

   private:
    template <typename T> static RANN::Fingerprint *fingerprint_creator(PairRANN *);
    template <typename T> static RANN::Activation *activation_creator(PairRANN *);
    //new functions
    void allocate(std::vector<std::string>);//called after reading element list, but before reading the rest of the potential
    void read_file(char *);//read potential file
    void strip_comments(char *,FILE *fp);
    void read_atom_types(std::vector<std::string>,std::vector<std::string>);
    void read_mass(std::vector<std::string>,std::vector<std::string>);
    void read_fpe(std::vector<std::string>,std::vector<std::string>);//fingerprints per element. Count total fingerprints defined for each 1st element in element combinations
    void read_fingerprints(std::vector<std::string>,std::vector<std::string>);
    void read_fingerprint_constants(std::vector<std::string>,std::vector<std::string>);
    void read_network_layers(std::vector<std::string>,std::vector<std::string>);//include input and output layer (hidden layers + 2)
    void read_layer_size(std::vector<std::string>,std::vector<std::string>);
    void read_weight(std::vector<std::string>,std::vector<std::string>,FILE*);//weights should be formatted as properly shaped matrices
    void read_bias(std::vector<std::string>,std::vector<std::string>,FILE*);//biases should be formatted as properly shaped vectors
    void read_activation_functions(std::vector<std::string>,std::vector<std::string>);
    void read_screening(std::vector<std::string>,std::vector<std::string>);
    bool check_potential();//after finishing reading potential file
    void propagateforward(double *,double **,double **,int,int);//called by compute to get force and energy
    void propagateforwardspin(double *,double **,double **,double**,int,int);//called by compute to get force and energy
    void screen(int,int,int);
    void cull_neighbor_list(int *,int,int);
    void screen_neighbor_list(int *,int,int);
  };

}

#endif
#endif



