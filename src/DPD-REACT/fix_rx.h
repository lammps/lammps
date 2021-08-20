/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(rx,FixRX);
// clang-format on
#else

#ifndef LMP_FIX_RX_H
#define LMP_FIX_RX_H

#include "fix.h"

typedef int (*fnptr)(double, const double *, double *, void *);
namespace LAMMPS_NS {

enum { ODE_LAMMPS_RK4, ODE_LAMMPS_RKF45 };

class FixRX : public Fix {
 public:
  FixRX(class LAMMPS *, int, char **);
  ~FixRX();
  int setmask();
  void post_constructor();
  virtual void init();
  void init_list(int, class NeighList *);
  virtual void setup_pre_force(int);
  virtual void pre_force(int);

 protected:
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

  class NeighList *list;

  double tmpArg;

  int *mol2param;    // mapping from molecule to parameters
  int nreactions;    // # of stored parameter sets
  int maxparam;      // max # of parameter sets
  struct Param {
    double cp;
    int ispecies;
    char *name;    // names of unique molecules and interaction type
  };
  Param *params;    // parameter set for an I-J-K interaction

  int nspecies;
  void read_file(char *);
  void setupParams();
  double *Arr, *nArr, *Ea, *tempExp;
  double **stoich, **stoichReactants, **stoichProducts;
  double *kR;

  //!< Classic Runge-Kutta 4th-order stepper.
  void rk4(int, double *, void *);

  //!< Runge-Kutta-Fehlberg ODE Solver.
  void rkf45(int, double *, void *, int ode_counter[]);

  //!< Runge-Kutta-Fehlberg ODE stepper function.
  void rkf45_step(const int neq, const double h, double y[], double y_out[], double rwk[], void *);

  //!< Initial step size estimation for the Runge-Kutta-Fehlberg ODE solver.
  int rkf45_h0(const int neq, const double t, const double t_stop, const double hmin,
               const double hmax, double &h0, double y[], double rwk[], void *v_params);

  class PairDPDfdtEnergy *pairDPDE;
  double *dpdThetaLocal;
  double *sumWeights;
  void computeLocalTemperature();
  int localTempFlag, wtFlag, odeIntegrationFlag;
  double sigFactor;

  int rhs(double, const double *, double *, void *);
  int rhs_dense(double, const double *, double *, void *);

  // User-defined data container needed in rhs.
  struct UserRHSData {
    double *kFor;
    double *rxnRateLaw;
  };

  // Sparse stoichiometric matrix storage format and methods.
  bool useSparseKinetics;
  //SparseKinetics sparseKinetics;
  void initSparse(void);
  int rhs_sparse(double, const double *, double *, void *) const;

  int sparseKinetics_maxReactants;    //<! Max # of reactants species in any reaction
  int sparseKinetics_maxProducts;     //<! Max # of product species in any reaction
  int sparseKinetics_maxSpecies;    //<! Max # of species (maxReactants + maxProducts) in any reaction

  //! Objects to hold the stoichiometric coefficients using a sparse matrix
  //! format. Enables a sparse formulation for the reaction rates:
  //! \f${\omega}_i = \Pi_{j=1}^{NS_i} K^{f}_i [x_j]^{\nu^{'}_{ij}} -
  //!                                  K^{r}_i x_j^{\nu^{''}_{ij}}\f$.
  double **sparseKinetics_nu;    //<! Stoichiometric matrix with FLT values.
  int **sparseKinetics_nuk;    //<! Index (base-0) of species ... this is the column sparse matrix.
  int **sparseKinetics_inu;    //<! Stoichiometric matrix with integral values.

  bool *
      sparseKinetics_isIntegralReaction;    //<! Flag indicating if a reaction has integer stoichiometric values.

  // ODE Parameters
  int minSteps;             //!< Minimum # of steps for the ODE solver(s).
  int maxIters;             //!< Maximum # of iterations for the ODE solver(s).
  double relTol, absTol;    //!< Relative and absolute tolerances for the ODE solver(s).

  // ODE Diagnostics
  //int nSteps; //!< # of accepted steps taken over all atoms.
  //int nIters; //!< # of attempted steps for all atoms.
  //int nFuncs; //!< # of RHS evaluations for all atoms.
  //int nFails; //!< # of ODE systems that failed (for some reason).

  int diagnosticFrequency;    //!< Frequency (LMP steps) that run-time diagnostics will be printed to the log.
  enum { numDiagnosticCounters = 5 };
  enum { StepSum = 0, FuncSum, TimeSum, AtomSum, CountSum };
  double diagnosticCounter[numDiagnosticCounters];

  int *diagnosticCounterPerODE[numDiagnosticCounters];

  //!< ODE Solver diagnostics.
  void odeDiagnostics(void);

 protected:
  char *kineticsFile;
  char *id_fix_species, *id_fix_species_old;
  class FixPropertyAtom *fix_species, *fix_species_old;
  int restartFlag;
};

}    // namespace LAMMPS_NS

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

W:  in FixRX::pre_force, ODE solver failed for %d atoms.

Self-explanatory

E:  Missing parameters in reaction kinetic equation.

Self-explanatory

E:  Potential file has duplicate entry.

Self-explanatory

E:  Computed concentration in RK4 (RKF45) solver is < -1.0e-10.

Self-explanatory:  Adjust settings for the RK4 solver.

*/
