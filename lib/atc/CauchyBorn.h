#ifndef CAUCHYBORN_H
#define CAUCHYBORN_H
#include "CbPotential.h"
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

// This file provides all routines necessary for computing the Cauchy-Born
// stresses.

namespace ATC {
  // forward declaration of CbPotential
  class CbPotential;
  class AtomCluster;

  /**
   *  @class  StressAtIP 
   *  @brief  Class for passing the vector of stresses at quadrature points 
   *          Done by storing the quadrature point and providing indexing 
   */

  class StressAtIP {
    INDEX q;    //*> Quadrature point.
    DENS_MAT_VEC &S;  //*> Stress components at all quad points.
  public:
    //* Constructor - sets stress reference and current quadrature point.
    StressAtIP(DENS_MAT_VEC &s, INDEX _q=0) : q(_q), S(s) {}
    //* Indexing operator, gets stress components for the quadrature point.
    double& operator()(INDEX i, INDEX j) const { return S[i](q, j); }
    //* Sets quadrature point to a new index.
    void set_quadrature_point(INDEX qp) { q = qp; }
    //* Operator that outputs the stress
    friend std::ostream& operator<<(std::ostream& o, const StressAtIP& s)
    { o << "stress\n"; 
      o << s(0,0) << " " << s(0,1) << " " << s(0,2) << "\n";
      o << s(1,0) << " " << s(1,1) << " " << s(1,2) << "\n";
      o << s(2,0) << " " << s(2,1) << " " << s(2,2) << "\n";
      return o;
    }
  };

  /**
   *  @class  StressArgs
   *  @brief  Class for storing parameters needed for computing the Cauchy-Born stress
   */

  class StressArgs {
  public:
    StressArgs(AtomCluster &v, CbPotential *p, double kB, double hbar, double T)
      : vac(v), potential(p), boltzmann_constant(kB), planck_constant(hbar),
      temperature(T) {}
    AtomCluster &vac;
    CbPotential *potential;
    double boltzmann_constant;
    double planck_constant;
    double temperature;
  };

  /**
   *  @class  PairParam 
   *  @brief  Class for storing parameters used in pairwise stress computations
   */
  struct PairParam {
    PairParam(const DENS_VEC &_R, double _d) : R(_R), d(_d), di(1.0/_d) {}
    const DENS_VEC &R;     //*> Reference bond vector.
    DENS_VEC r;            //*> Current bond vector.
    double d, di;          //*> Current bond length and its inverse.
    double phi_r;          //*> First derivative of pairwise term.
    double phi_rr;         //*> Second derivative of pairwise term.
    double phi_rrr;        //*> Third derivative of pairwise term.
    double rho_r;
    double rho_rr;
    double rho_rrr;
    double F_p;
    double F_pp;
    double F_ppp;
  };

  //* for EAM, calculate electron density
  double cb_electron_density(const StressArgs &args);
  //* Compute stress, from virtual atom cluster, potential
  //* temperature (can be 0K), and StressAtIP object to write to.
  void cb_stress(const StressArgs &args, StressAtIP &s, double *F=0);
  //* Computes the elastic energy (free or potential if T=0).
  double cb_energy(const StressArgs &args);
  //* Computes the entropic energy 
  double cb_entropic_energy(const StressArgs &args);
  //* Auxiliary functions for cb_stress
  //@{
  //* Computes the stress contribution given the pairwise parameters.
  void pairwise_stress(const PairParam &p, StressAtIP &s);
  //* Computes the stress contribution given the embedding parameters.
  void embedding_stress(const PairParam &p, StressAtIP &s);
  //* Computes the pairwise thermal components for the stress and free energy.
  void pairwise_thermal(const PairParam &p, DENS_MAT &D, DENS_MAT_VEC *dDdF=0);
  //* Computes the embedding thermal components for the stress and free energy.
  void embedding_thermal(const PairParam &p, DENS_MAT &D, DENS_MAT &L0, DENS_MAT_VEC *dDdF=0);
  //* Last stage of the pairwise finite-T Cauchy-Born stress computation.
  void thermal_end(const DENS_MAT_VEC &DF, const DENS_MAT &D, const DENS_MAT &F,
                   const double &T, const double &kb, StressAtIP &s, double *F_w=0); 
  //* Returns the stretch tensor and its derivative with respect to C (R C-G). 
  void stretch_tensor_derivative(const DENS_VEC &C, DENS_VEC &U, DENS_MAT &dU);
  //@}

  //* Testing functions (to be removed when all CB code is completed)
  //@{
  //* Computes the dynamical matrix (TESTING FUNCTION)
  DENS_MAT compute_dynamical_matrix(const StressArgs &args);
  //* Computes the determinant of the dynamical matrix (TESTING FUNCTION)
  double compute_detD(const StressArgs &args);
  //* Computes the derivative of the dynamical matrix (TESTING FUNCTION)
  DENS_MAT_VEC compute_dynamical_derivative(StressArgs &args);
  //@}
}
#endif
