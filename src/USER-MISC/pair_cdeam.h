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

PairStyle(eam/cd,PairCDEAM_OneSite)
PairStyle(eam/cd/old,PairCDEAM_TwoSite)

#else

#ifndef LMP_PAIR_CDEAM_H
#define LMP_PAIR_CDEAM_H

#include "pair_eam_alloy.h"

namespace LAMMPS_NS {

class PairCDEAM : public PairEAMAlloy
{
public:
  /// Constructor.
  PairCDEAM(class LAMMPS*, int cdeamVersion);

  /// Destructor.
  virtual ~PairCDEAM();

  /// Calculates the energies and forces for all atoms in the system.
  virtual void compute(int, int);

  /// Parses the pair_coeff command parameters for this pair style.
  void coeff(int, char **);

  /// This is for MPI communication with neighbor nodes.
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

  /// Reports the memory usage of this pair style to LAMMPS.
  double memory_usage();

  /// Parses the coefficients of the h polynomial from the end of the EAM file.
  void read_h_coeff(char* filename);

 public:
  // The public interface exposed by this potential class.

  // Evaluates the h(x) polynomial for a given local concentration x.
  inline double evalH(double x) const {
    double v = 0.0;
    for(int i = nhcoeff-1; i >= 1; i--) {
      v = (v + hcoeff[i]) * x;
    }
    return v + hcoeff[0];
  };

  // Calculates the derivative of the h(x) polynomial.
  inline double evalHprime(double x) const {
    double v = 0.0;
    for(int i = nhcoeff-1; i >= 2; i--) {
      v = (v + (double)i * hcoeff[i]) * x;
    }
    return v + hcoeff[1];
  };

  // We have two versions of the CD-EAM potential. The user can choose which one he wants to use:
  //
  // Version 1 - One-site concentration: The h(x_i) function depends only on the concentration at the atomic site i.
  // This is a new version with a slight modification to the formula. It happens to be computationally more efficient.
  // It has been published in the MSMSE 2009 paper.
  //
  // Version 2 - Two-site concentration: The h(x_ij) function depends on the concentrations at two atomic sites i and j.
  // This is the original version from the 2005 PRL paper.
  int cdeamVersion;

  // Per-atom arrays

  // The partial density of B atoms at each atom site.
  double *rhoB;

  // The intermediate value D_i for each atom.
  // The meaning of these values depend on the version of the CD-EAM potential used:
  //
  //   For the one-site concentration version these are the v_i values defined in equation (21)
  //   of the MSMSE paper.
  //
  //   For the old two-site concentration version these are the M_i values defined in equation (12)
  //   of the MSMSE paper.
  double *D_values;

  // The atom type index that is considered to be the A species in the alloy.
  int speciesA;

  // The atom type index that is considered to be the B species in the alloy.
  int speciesB;

 protected:

  // Evaluation functions:

  // This structure specifies an entry in one of the EAM spline tables
  // and the corresponding floating point part.
  typedef struct {
    int m;
    double p;
  } EAMTableIndex;

  // Converts a radius value to an index value to be used in a spline table lookup.
  inline EAMTableIndex radiusToTableIndex(double r) const {
    EAMTableIndex index;
    index.p = r*rdr + 1.0;
    index.m = static_cast<int>(index.p);
    index.m = index.m <= (nr-1) ? index.m : (nr-1);
    index.p -= index.m;
    index.p = index.p <= 1.0 ? index.p : 1.0;
    return index;
  };

  // Converts a density value to an index value to be used in a spline table lookup.
  inline EAMTableIndex rhoToTableIndex(double rho) const {
    EAMTableIndex index;
    index.p = rho*rdrho + 1.0;
    index.m = static_cast<int>(index.p);
    index.m = index.m <= (nrho-1) ? index.m : (nrho-1);
    index.p -= index.m;
    index.p = index.p <= 1.0 ? index.p : 1.0;
    return index;
  };

  // Computes the derivative of rho(r)
  inline double RhoPrimeOfR(const EAMTableIndex& index, int itype, int jtype) const {
    const double* coeff = rhor_spline[type2rhor[itype][jtype]][index.m];
    return (coeff[0]*index.p + coeff[1])*index.p + coeff[2];
  };

  // Computes rho(r)
  inline double RhoOfR(const EAMTableIndex& index, int itype, int jtype) const {
    const double* coeff = rhor_spline[type2rhor[itype][jtype]][index.m];
    return ((coeff[3]*index.p + coeff[4])*index.p + coeff[5])*index.p + coeff[6];
  };

  // Computes the derivative of F(rho)
  inline double FPrimeOfRho(const EAMTableIndex& index, int itype) const {
    const double* coeff = frho_spline[type2frho[itype]][index.m];
    return (coeff[0]*index.p + coeff[1])*index.p + coeff[2];
  };

  // Computes F(rho)
  inline double FofRho(const EAMTableIndex& index, int itype) const {
    const double* coeff = frho_spline[type2frho[itype]][index.m];
    return ((coeff[3]*index.p + coeff[4])*index.p + coeff[5])*index.p + coeff[6];
  };

  // Computes the derivative of z2(r)
  inline double Z2PrimeOfR(const EAMTableIndex& index, int itype, int jtype) const {
    const double* coeff = z2r_spline[type2z2r[itype][jtype]][index.m];
    return (coeff[0]*index.p + coeff[1])*index.p + coeff[2];
  };

  // Computes z2(r)
  inline double Z2OfR(const EAMTableIndex& index, int itype, int jtype) const {
    const double* coeff = z2r_spline[type2z2r[itype][jtype]][index.m];
    return ((coeff[3]*index.p + coeff[4])*index.p + coeff[5])*index.p + coeff[6];
  };

  // Computes pair potential V_ij(r).
  inline double PhiOfR(const EAMTableIndex& index, int itype, int jtype, const double oneOverR) const {
    // phi = pair potential energy
    // z2 = phi * r
    const double* coeff = z2r_spline[type2z2r[itype][jtype]][index.m];
    const double z2 = ((coeff[3]*index.p + coeff[4])*index.p + coeff[5])*index.p + coeff[6];
    return z2 * oneOverR;
  };

  // Computes pair potential V_ij(r) and its derivative.
  inline double PhiOfR(const EAMTableIndex& index, int itype, int jtype, const double oneOverR, double& phid) const {
    // phi = pair potential energy
    // phip = phi'
    // z2 = phi * r
    // z2p = (phi * r)' = (phi' r) + phi
    const double* coeff = z2r_spline[type2z2r[itype][jtype]][index.m];
    const double z2p = (coeff[0]*index.p + coeff[1])*index.p + coeff[2];
    const double z2 = ((coeff[3]*index.p + coeff[4])*index.p + coeff[5])*index.p + coeff[6];
    const double phi = z2 * oneOverR;
    phid = z2p * oneOverR - phi * oneOverR;
    return phi;
  };

  // Parameters

  // h() polynomial function coefficients
  double* hcoeff;
  // The number of coefficients in the polynomial.
  int nhcoeff;

  // This specifies the calculation stage to let the pack/unpack communication routines know
  // which data to exchange.
  int communicationStage;
};

/// The one-site concentration formulation of CD-EAM.
 class PairCDEAM_OneSite : public PairCDEAM
   {
   public:
     /// Constructor.
     PairCDEAM_OneSite(class LAMMPS* lmp) : PairEAM(lmp), PairCDEAM(lmp, 1) {}
   };

 /// The two-site concentration formulation of CD-EAM.
 class PairCDEAM_TwoSite : public PairCDEAM
   {
   public:
     /// Constructor.
     PairCDEAM_TwoSite(class LAMMPS* lmp) : PairEAM(lmp), PairCDEAM(lmp, 2) {}
   };

}

#endif
#endif
