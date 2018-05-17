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

PairStyle(gw,PairGW)

#else

#ifndef LMP_PAIR_GW_H
#define LMP_PAIR_GW_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGW : public Pair {
 public:
  PairGW(class LAMMPS *);
  virtual ~PairGW();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    int ielement,jelement,kelement;
    int powermint;
    double Z_i,Z_j;
    double ZBLcut,ZBLexpscale;
  };

  Param *params;                // parameter set for an I-J-K interaction
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to paramegw
  int *map;                     // mapping from atom types to elements
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets

  int **pages;                     // neighbor list pages
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int maxpage;                     // # of pages currently allocated
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom


  int *GW_numneigh;             // # of pair neighbors for each atom
  int **GW_firstneigh;          // ptr to 1st neighbor of each atom

  void GW_neigh();
  void add_pages(int howmany = 1);

  void allocate();
  virtual void read_file(char *);
  void setup_params();
  virtual void repulsive(Param *, double, double &, int, double &);
  double zeta(Param *, double, double, double *, double *);
  virtual void force_zeta(Param *, double, double, double &,
                          double &, int, double &);
  void attractive(Param *, double, double, double, double *, double *,
                  double *, double *, double *);

  double gw_fc(double, Param *);
  double gw_fc_d(double, Param *);
  virtual double gw_fa(double, Param *);
  virtual double gw_fa_d(double, Param *);
  double gw_bij(double, Param *);
  double gw_bij_d(double, Param *);

  void gw_zetaterm_d(double, double *, double, double *, double,
                               double *, double *, double *, Param *);
  void costheta_d(double *, double, double *, double,
                  double *, double *, double *);

  // inlined functions for efficiency

  inline double gw_gijk(const double costheta,
                          const Param * const param) const {
    const double gw_c = param->c * param->c;
    const double gw_d = param->d * param->d;
    const double hcth = param->h - costheta;

          //printf("gw_gijk: gw_c=%f gw_d=%f hcth=%f=%f-%f\n", gw_c, gw_d, hcth, param->h, costheta);

    return param->gamma*(1.0 + gw_c/gw_d - gw_c / (gw_d + hcth*hcth));
  }

  inline double gw_gijk_d(const double costheta,
                            const Param * const param) const {
    const double gw_c = param->c * param->c;
    const double gw_d = param->d * param->d;
    const double hcth = param->h - costheta;
    const double numerator = -2.0 * gw_c * hcth;
    const double denominator = 1.0/(gw_d + hcth*hcth);
    return param->gamma*numerator*denominator*denominator;
  }

  inline double vec3_dot(const double x[3], const double y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  inline void vec3_add(const double x[3], const double y[3],
                       double * const z) const {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }

  inline void vec3_scale(const double k, const double x[3],
                         double y[3]) const {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }

  inline void vec3_scaleadd(const double k, const double x[3],
                            const double y[3], double * const z) const {
    z[0] = k*x[0]+y[0];
    z[1] = k*x[1]+y[1];
    z[2] = k*x[2]+y[2];
  }
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

E: Pair style GW requires atom IDs

This is a requirement to use the GW potential.

E: Pair style GW requires newton pair on

See the newton command.  This is a restriction to use the GW
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open GW potential file %s

The specified GW potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in GW potential file

Incorrect number of words per line in the potential file.

E: Illegal GW parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/
