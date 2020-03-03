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

PairStyle(cac/eam/interp,PairCACEAMInterp)

#else

#ifndef LMP_PAIR_EAM_CAC_INTERP_H
#define LMP_PAIR_EAM_CAC_INTERP_H

#include "pair_cac.h"

namespace LAMMPS_NS {

class PairCACEAMInterp : public PairCAC {
 public:
	 friend class FixSemiGrandCanonicalMC;   // Alex Stukowski option
  // potentials as array data

  int nrho, nr;
  int nfrho, nrhor, nz2r;
  double **frho, **rhor, **z2r;
  int *type2frho, **type2rhor, **type2z2r;

	// potentials in spline form used for force computation

  double dr, rdr, drho, rdrho, rhomax;
  double ***rhor_spline, ***frho_spline, ***z2r_spline;

  PairCACEAMInterp(class LAMMPS *);
  virtual ~PairCACEAMInterp();
  
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  
  virtual void *extract(const char *, int &);
  virtual void swap_eam(double *, double **);

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
  double cutforcesq;
  double **scale;

	// per-atom arrays
  double *rho, *fp;
  double density;
  double **quad_electron_densities;
  int max_density;
 
  virtual void allocate();
  virtual void read_file(char *);
  virtual void array2spline();
  virtual void file2array();

  //EAM variables and structures

  // potentials as file data

  int *map;                   // which element each atom type maps to

  struct Funcfl {
	  char *file;
	  int nrho, nr;
	  double drho, dr, cut, mass;
	  double *frho, *rhor, *zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
	  char **elements;
	  int nelements, nrho, nr;
	  double drho, dr, cut;
	  double *mass;
	  double **frho, **rhor, ***z2r;
  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, ***rhor, ***z2r;
  };
  Fs *fs;

  virtual void interpolate(int, double, double *, double **);
  virtual void grab(FILE *, int, double *);

  //further CAC functions 
  virtual void force_densities(int, double, double, double, double, double
    &fx, double &fy, double &fz);
  virtual void pre_force_densities();
  virtual void compute_electron_densities(int);
  virtual void quad_electron_density(int, double, double, double);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

E: Invalid EAM potential file

UNDOCUMENTED

*/