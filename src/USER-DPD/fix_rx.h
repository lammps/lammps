/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(rx,FixRX)

#else

#ifndef LMP_FIX_RX_H
#define LMP_FIX_RX_H

#include "fix.h"

typedef int (*fnptr)(double, const double *, double *, void *);
namespace LAMMPS_NS {

enum { ODE_LAMMPS_RK4 };

class FixRX : public Fix {
 public:
  FixRX(class LAMMPS *, int, char **);
  ~FixRX();
  int setmask();
  void post_constructor();
  virtual void init();
  virtual void setup_pre_force(int);
  virtual void pre_force(int);

 protected:
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

  double tmpArg;

  int *mol2param;               // mapping from molecule to parameters
  int nreactions;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  struct Param {
    double cp;
    int ispecies;
    char *name;      // names of unique molecules and interaction type
  };
  Param *params;                // parameter set for an I-J-K interaction

  int nspecies;
  void read_file(char *);
  void setupParams();
  double *Arr, *nArr, *Ea, *tempExp;
  double **stoich, **stoichReactants, **stoichProducts;
  double *kR;
  void rk4(int);

  class PairDPDfdtEnergy *pairDPDE;
  double *dpdThetaLocal;
  double *sumWeights;
  void computeLocalTemperature();
  int localTempFlag,wtFlag,odeIntegrationFlag;
  double sigFactor;

  int rhs(double, const double *, double *, void *);

 private:
  char *kineticsFile;
  char *id_fix_species,*id_fix_species_old;
  class FixPropertyAtom *fix_species,*fix_species_old;

  // ODE Parameters
  int minSteps; //!< Minimum # of steps for the ODE solver(s).

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E:  fix rx cannot be combined with fix property/atom

Self-explanatory

E:  Cannot open rx file %s

Self-explanatory

E:  Exceeded the maximum number of species permitted in fix rx

Reduce the number of species in the fix rx reaction kinetics file

E:  There are no rx species specified.

Self-explanatory

E:  Must use pair_style dpd/fdt/energy with fix rx.

Self-explanatory

E:  fix rx requires fix eos/table/rx to be specified.

Self-explanatory

E:  Missing parameters in reaction kinetic equation.

Self-explanatory

E:  Potential file has duplicate entry.

Self-explanatory

E:  Computed concentration in RK4 solver is < -1.0e-10.

Self-explanatory:  Adjust settings for the RK4 solver. 

*/
