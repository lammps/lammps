#include "CauchyBorn.h"
#include "VoigtOperations.h"
#include "CBLattice.h"
#include "CbPotential.h"

using voigt3::to_voigt;

namespace ATC {
  //============================================================================
  // Computes the electron density for EAM potentials
  //============================================================================
  double cb_electron_density(const StressArgs &args )
  {
    double e_density = 0.0;
    for (INDEX a=0; a<args.vac.size(); a++) {
      PairParam pair(args.vac.R(a), args.vac.bond_length(a));
      e_density += args.potential->rho(pair.d);
    }
    return e_density;
  }
  //============================================================================
  // Computes the stress at a quadrature point
  //============================================================================
  void cb_stress(const StressArgs &args, StressAtIP &s, double *F)
  {
    const double &T = args.temperature;
    const bool finite_temp = T > 0.0;
    DENS_MAT D;         // dynamical matrix (finite temp)
    DENS_MAT_VEC dDdF;    // derivative of dynamical matrix (finite temp)
    double e_density(0.),embed(0.),embed_p(0.),embed_pp(0.),embed_ppp(0.);
    DENS_VEC l0;
    DENS_MAT L0;
    DENS_MAT_VEC M0;

    // If temperature is nonzero then allocate space for 
    // dynamical matrix and its derivative with respect to F.
    if (finite_temp)  {
      D.reset(3,3);
      dDdF.assign(6, DENS_MAT(3,3));
      M0.assign(3, DENS_MAT(3,3));
      L0.reset(3,3);
      l0.reset(3);
    }

    if (F) *F = 0.0; 

    // if using EAM potential, calculate embedding function and derivatives
    if (args.potential->terms.embedding) {
      
      for (INDEX a=0; a<args.vac.size(); a++) {
        PairParam pair(args.vac.R(a), args.vac.bond_length(a));
        e_density += args.potential->rho(pair.d);
        pair.r = args.vac.r(a);
        pair.rho_r = args.potential->rho_r(pair.d);
        pair.rho_rr = args.potential->rho_rr(pair.d);
        if (finite_temp) {
          l0 += pair.r*pair.di*pair.rho_r;
          DENS_MAT rR = tensor_product(pair.r, pair.R);
          L0.add_scaled(rR, pair.di*pair.rho_r);
          DENS_MAT rr = tensor_product(pair.r, pair.r);
          rr *= pair.di*pair.di*(pair.rho_rr - pair.di*pair.rho_r); 
          diagonal(rr) += pair.di*pair.rho_r;
          for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
               for (int L = 0; L < 3; L++) {
                 M0[i](k,L) += rr(i,k)*args.vac.R(a)(L);
               }
            }
          }
        }
      }
      embed = args.potential->F(e_density);  // "F" in usual EAM symbology
      embed_p = args.potential->F_p(e_density);
      embed_pp = args.potential->F_pp(e_density);
      embed_ppp = args.potential->F_ppp(e_density);
      if (F) *F += embed;
      if (finite_temp) {
        const DENS_MAT ll = tensor_product(l0, l0);
        D.add_scaled(ll, embed_pp);
        const DENS_VEC llvec = to_voigt(ll);
        for (int v = 0; v < 6; v++) {
          dDdF[v].add_scaled(L0, embed_ppp*llvec(v));
        }
        dDdF[0].add_scaled(M0[0], 2*embed_pp*l0(0));
        dDdF[1].add_scaled(M0[1], 2*embed_pp*l0(1));
        dDdF[2].add_scaled(M0[2], 2*embed_pp*l0(2));
        dDdF[3].add_scaled(M0[1], embed_pp*l0(2));
        dDdF[3].add_scaled(M0[2], embed_pp*l0(1));
        dDdF[4].add_scaled(M0[0], embed_pp*l0(2));
        dDdF[4].add_scaled(M0[2], embed_pp*l0(0));
        dDdF[5].add_scaled(M0[0], embed_pp*l0(1));
        dDdF[5].add_scaled(M0[1], embed_pp*l0(0));
      }
    }

    // Loop on all cluster atoms (origin atom not included).
    for (INDEX a=0; a<args.vac.size(); a++) {
      PairParam pair(args.vac.R(a), args.vac.bond_length(a));
      if (args.potential->terms.pairwise) { 
        if (F) *F += 0.5*args.potential->phi(pair.d);
        pair.phi_r = args.potential->phi_r(pair.d);
        pairwise_stress(pair, s);
      }
      if (args.potential->terms.embedding) {
        pair.F_p = embed_p;
        pair.rho_r = args.potential->rho_r(pair.d);
        embedding_stress(pair, s); 
      }

      if (finite_temp) {  // Compute finite T terms.
        pair.r = args.vac.r(a);
        if (args.potential->terms.pairwise) {
          pair.phi_rr  = args.potential->phi_rr(pair.d);
          pair.phi_rrr = args.potential->phi_rrr(pair.d);
          pairwise_thermal(pair, D, &dDdF);
        }
        if (args.potential->terms.embedding) {
          pair.rho_rr = args.potential->rho_rr(pair.d);
          pair.rho_rrr = args.potential->rho_rrr(pair.d);
          pair.F_pp = embed_pp;
          pair.F_ppp = embed_ppp;
          embedding_thermal(pair,D,L0,&dDdF);
        }
      }
      // if has three-body terms ...       TODO compute three-body terms
    }

    // Finish finite temperature Cauchy-Born.
    if (finite_temp) {
      const DENS_MAT &F = args.vac.deformation_gradient();
      thermal_end(dDdF, D, F, T, args.boltzmann_constant, s);
    }
  }
  //===========================================================================
  // Computes the elastic energy (free or potential if T=0).
  //===========================================================================
  double cb_energy(const StressArgs &args)
  {
    const double &T = args.temperature;
    bool finite_temp = (T > 0.0);
    //const bool finite_temp = T > 0.0;
    DENS_MAT D;         // dynamical matrix (finite temp)
    double e_density,embed,embed_p(0.),embed_pp(0.),embed_ppp(0.);
    DENS_VEC l0;
    DENS_MAT L0;
    DENS_MAT_VEC M0;

    // If temperature is nonzero then allocate space for dynamical matrix.
    if (finite_temp)  {
      D.reset(3,3);
      l0.reset(3);
    }

    double F = 0.0;
    // Do pairwise terms, loop on all cluster atoms (origin atom not included).
    // if using EAM potential, calculate embedding function and derivatives
    if (args.potential->terms.embedding) {
      e_density = 0.0;
      for (INDEX a=0; a<args.vac.size(); a++) {
        PairParam pair(args.vac.R(a), args.vac.bond_length(a));
        e_density += args.potential->rho(pair.d);
        pair.r = args.vac.r(a);
        if (finite_temp) {
          l0 += pair.r*pair.di*pair.rho_r;
        }
      }
      embed = args.potential->F(e_density);
      embed_p = args.potential->F_p(e_density);
      embed_pp = args.potential->F_pp(e_density);
      embed_ppp = args.potential->F_ppp(e_density);
      F += embed;
      if (finite_temp) {
        const DENS_MAT ll = tensor_product(l0, l0);
        D.add_scaled(ll, embed_pp);
      }
    }

    for (INDEX a=0; a<args.vac.size(); a++) {
      PairParam pair(args.vac.R(a), args.vac.bond_length(a));
      if (args.potential->terms.pairwise) { 
        F += 0.5*args.potential->phi(pair.d);
      }

      if (finite_temp) {  // Compute finite T terms.
        pair.r = args.vac.r(a);
        if (args.potential->terms.pairwise) {
          pair.phi_r = args.potential->phi_r(pair.d);
          pair.phi_rr  = args.potential->phi_rr(pair.d);    
          pair.phi_rrr = args.potential->phi_rrr(pair.d);    
          pairwise_thermal(pair, D);
        }
        if (args.potential->terms.embedding) { 
          pair.rho_r = args.potential->rho_r(pair.d);
          pair.rho_rr = args.potential->rho_rr(pair.d);
          pair.rho_rrr = args.potential->rho_rrr(pair.d);
          pair.F_p = embed_p;
          pair.F_pp = embed_pp;
          pair.F_ppp = embed_ppp;
          embedding_thermal(pair,D,L0);
        }
      }
      // if has three-body terms ...       TODO compute three-body terms
    }  
    // Finish finite temperature Cauchy-Born.
    const double kB = args.boltzmann_constant;
    const double hbar = args.planck_constant;
    if (finite_temp) {
      F += kB*T*log(pow(hbar/(kB*T),3.0)*sqrt(det(D)));
    }
    //if (finite_temp) F += 0.5*args.boltzmann_constant*T*log(det(D));
    return F;
  } 
  //===========================================================================
  // Computes the entropic energy TS (minus c_v T)
  //===========================================================================
  double cb_entropic_energy(const StressArgs &args)
  {
    const double &T = args.temperature;
    DENS_MAT D(3,3);         // dynamical matrix (finite temp)
    double e_density,embed_p(0.),embed_pp(0.),embed_ppp(0.);
    DENS_VEC l0(3);
    DENS_MAT L0;
    DENS_MAT_VEC M0;

    // if using EAM potential, calculate embedding function and derivatives
    if (args.potential->terms.embedding) {
      e_density = 0.0;
      for (INDEX a=0; a<args.vac.size(); a++) {
        PairParam pair(args.vac.R(a), args.vac.bond_length(a));
        e_density += args.potential->rho(pair.d);
        pair.r = args.vac.r(a);
        l0 += pair.r*pair.di*pair.rho_r;
        //DENS_MAT rR = tensor_product(pair.r, pair.R);
        //L0.add_scaled(rR, pair.di*args.potential->rho_r(pair.d));
      }
      //embed = args.potential->F(e_density);
      embed_p = args.potential->F_p(e_density);
      embed_pp = args.potential->F_pp(e_density);
      embed_ppp = args.potential->F_ppp(e_density);
      const DENS_MAT ll = tensor_product(l0, l0);
      D.add_scaled(ll, embed_pp);
    }

    // Compute the dynamical matrix
    // Loop on all cluster atoms (origin atom not included).
    for (INDEX a=0; a<args.vac.size(); a++) {
      // Compute pairwise terms needed for pairwise_stress.
      PairParam pair(args.vac.R(a), args.vac.bond_length(a));
      pair.r = args.vac.r(a);
      if (args.potential->terms.pairwise) {
        pair.phi_r = args.potential->phi_r(pair.d);
        pair.phi_rr  = args.potential->phi_rr(pair.d);
        pair.phi_rrr = args.potential->phi_rrr(pair.d);
        pairwise_thermal(pair, D);
      }
      if (args.potential->terms.embedding) {
        pair.rho_r = args.potential->rho_r(pair.d);
        pair.rho_rr = args.potential->rho_rr(pair.d);
        pair.rho_rrr = args.potential->rho_rrr(pair.d);
        pair.F_p = embed_p;
        pair.F_pp = embed_pp;
        pair.F_ppp = embed_ppp;
        embedding_thermal(pair,D,L0);
      }
    }
    // Finish finite temperature Cauchy-Born.
    const double kB = args.boltzmann_constant;
    const double hbar = args.planck_constant;;
    double F = kB*T*log(pow(hbar/(kB*T),3.0)*sqrt(det(D)));
    return F;
  } 
  //===========================================================================
  // Computes the stress contribution given the pairwise parameters.
  //===========================================================================
  inline void pairwise_stress(const PairParam &p, StressAtIP &s)
  {
    for (INDEX i=0; i<p.R.size(); i++)
      for (INDEX j=i; j<p.R.size(); j++)
        s(i,j) += 0.5*p.di * p.phi_r * p.R(i) * p.R(j);
  }

  //===========================================================================
  // Computes the stress contribution given the embedding parameters.
  //===========================================================================
  inline void embedding_stress(const PairParam &p, StressAtIP &s)
  {
    for (INDEX i=0; i<p.R.size(); i++)
      for (INDEX j=i; j<p.R.size(); j++)
        s(i,j) += p.di * p.F_p * p.rho_r * p.R(i) * p.R(j);
  }

  //===========================================================================
  // Computes the pairwise thermal components for the stress
  //===========================================================================
  void pairwise_thermal(const PairParam &p, DENS_MAT &D, DENS_MAT_VEC *dDdF)
  {
    const double di2 = p.di*p.di;
    const double g   = p.di*p.phi_r;  
    const double g_d = p.di*p.phi_rr - p.di*g;  // units (energy / length^3)
    const double f   = di2 * (p.phi_rr - g);    // units (energy / length^4)
    const double f_d = di2*(p.phi_rrr-g_d) - 2.0*p.di*f;

    // compute needed tensor products of r and R 
    const DENS_MAT rr = tensor_product(p.r, p.r);

    // compute the dynamical matrix 
    D.add_scaled(rr, f);
    diagonal(D) += g;
    
    if (!dDdF) return;  // skip derivative
    const double gp_r = g_d*p.di;
    const double fp_r = f_d*p.di;
    const double fr[] = {f*p.r(0), f*p.r(1), f*p.r(2)};    
    const DENS_MAT rR = tensor_product(p.r, p.R);

    DENS_MAT_VEC &dD = *dDdF;
  
    // compute first term in A.13
    dD[0].add_scaled(rR, fp_r*rr(0,0) + gp_r);
    dD[1].add_scaled(rR, fp_r*rr(1,1) + gp_r);
    dD[2].add_scaled(rR, fp_r*rr(2,2) + gp_r);
    dD[3].add_scaled(rR, fp_r*rr(1,2));
    dD[4].add_scaled(rR, fp_r*rr(0,2));
    dD[5].add_scaled(rR, fp_r*rr(0,1));

    // compute second term in A.13
    for (INDEX L=0; L<p.R.size(); L++) {
      dD[0](0,L) += p.R[L] * 2.0*fr[0];
      dD[1](1,L) += p.R[L] * 2.0*fr[1];
      dD[2](2,L) += p.R[L] * 2.0*fr[2];
      dD[3](1,L) += p.R[L] * fr[2];
      dD[3](2,L) += p.R[L] * fr[1];
      dD[4](0,L) += p.R[L] * fr[2];
      dD[4](2,L) += p.R[L] * fr[0];
      dD[5](0,L) += p.R[L] * fr[1];
      dD[5](1,L) += p.R[L] * fr[0];
    }
  }

  //===========================================================================
  // Computes the embedding thermal components for the stress
  //===========================================================================
  void embedding_thermal(const PairParam &p, DENS_MAT &D, DENS_MAT &L0, DENS_MAT_VEC *dDdF)
  {
    const double di = p.di;
    const double di2 = p.di*p.di;
    const double di3 = p.di*p.di*p.di;
    const double x = p.F_pp*p.rho_r*p.rho_r + 2*p.F_p*p.rho_rr;
    const double z = di*(2*p.F_p*p.rho_r);
    const double y = di2*(x-z);

    // compute needed tensor products of r and R 
    const DENS_MAT rr = tensor_product(p.r, p.r);

    // compute the dynamical matrix 
    D.add_scaled(rr, y);
    diagonal(D) += z;
    
    if (!dDdF) return;  // skip derivative
    DENS_MAT_VEC &dD = *dDdF;
    const DENS_MAT rR = tensor_product(p.r, p.R);
    double rho_term1 = p.rho_rr - di*p.rho_r;
    double rho_term2 = p.rho_r*rho_term1;
    double rho_term3 = p.rho_rrr - 3*di*p.rho_rr + 3*di2*p.rho_r;
    const double a = di2*2*p.F_p*rho_term1;
    const double b = di2*(p.F_ppp*p.rho_r*p.rho_r + 2*p.F_pp*rho_term1);
    const double c = di3*(2*p.F_pp*rho_term2 + 2*p.F_p*rho_term3);
    const double w = di2*p.F_pp*p.rho_r*p.rho_r;

    //first add terms that multiply rR
    dD[0].add_scaled(rR, a + c*rr(0,0));
    dD[1].add_scaled(rR, a + c*rr(1,1));
    dD[2].add_scaled(rR, a + c*rr(2,2));
    dD[3].add_scaled(rR, c*rr(1,2));
    dD[4].add_scaled(rR, c*rr(0,2));
    dD[5].add_scaled(rR, c*rr(0,1));

    //add terms that multiply L0
    dD[0].add_scaled(L0, di*2*p.F_pp*p.rho_r + b*rr(0,0));
    dD[1].add_scaled(L0, di*2*p.F_pp*p.rho_r + b*rr(1,1));
    dD[2].add_scaled(L0, di*2*p.F_pp*p.rho_r + b*rr(2,2));
    dD[3].add_scaled(L0, b*rr(1,2));
    dD[4].add_scaled(L0, b*rr(0,2));
    dD[5].add_scaled(L0, b*rr(0,1));
  
    //add remaining term 
    const double aw = a + w;
    const double awr[] = {aw*p.r(0), aw*p.r(1), aw*p.r(2)};    
    for (INDEX L=0; L<p.R.size(); L++) {
      dD[0](0,L) += 2*awr[0]*p.R[L];
      dD[1](1,L) += 2*awr[1]*p.R[L];
      dD[2](2,L) += 2*awr[2]*p.R[L];
      dD[3](2,L) += awr[1]*p.R[L];
      dD[3](1,L) += awr[2]*p.R[L];
      dD[4](2,L) += awr[0]*p.R[L];
      dD[4](0,L) += awr[2]*p.R[L];
      dD[5](1,L) += awr[0]*p.R[L];
      dD[5](0,L) += awr[1]*p.R[L];
    }
  }

  //===========================================================================
  // Last stage of the pairwise finite-T Cauchy-Born stress computation.
  //===========================================================================
  inline void thermal_end(const DENS_MAT_VEC &DF, // dynamical matrix derivative
                          const DENS_MAT &D,    // dynamical matrix
                          const DENS_MAT &F,    // deformation gradient
                          const double &T,      // temperature
                          const double &kb,     // boltzmann constant
                          StressAtIP &s,        // output stress (-)
                          double* F_w)          // output free energy (optional)
  {
    DENS_MAT c = adjugate(D), dd(3,3);
    dd.add_scaled(DF[0], c(0,0));
    dd.add_scaled(DF[1], c(1,1));
    dd.add_scaled(DF[2], c(2,2));
    dd.add_scaled(DF[3], c(1,2) + c(2,1));
    dd.add_scaled(DF[4], c(0,2) + c(2,0));
    dd.add_scaled(DF[5], c(0,1) + c(1,0));

    const double detD = det(D);
    const double factor = 0.5*kb*T/detD;
    // converts from PK1 to PK2 
    dd = inv(F)*dd; 
    for (INDEX i=0; i<3; i++)
      for (INDEX j=i; j<3; j++)
        s(i,j) += factor * dd(i,j);

    // If f_W is not NULL then append thermal contribution.
    if (F_w) *F_w += 0.5*kb*T*log(detD);
  }
  //============================================================================
  // Returns the stretch tensor and its derivative with respect to C (R C-G). 
  //============================================================================
  void stretch_tensor_derivative(const DENS_VEC &C, DENS_VEC &U, DENS_MAT &dU)
  {
    // Compute the invariants of C
    const DENS_VEC C2(voigt3::dsymm(C,C));
    const double Ic   = voigt3::tr(C); 
    const double IIc  = 0.5*(Ic*Ic - voigt3::tr(C2));
    const double IIIc = voigt3::det(C);
    const DENS_VEC  I = voigt3::eye(3);

    // Compute the derivatives of the invarants of C
    DENS_VEC dIc   ( I );
    DENS_VEC dIIc  ( Ic*dIc - C );
    DENS_VEC dIIIc ( voigt3::inv(C) * IIIc );
    for (INDEX i=3; i<6; i++) {
      dIIc(i)  *= 2.0;
      dIIIc(i) *= 2.0;
    }

    // Check if C is an isotropic tensor (simple case)
    const double k = Ic*Ic - 3.0*IIc;
    const DENS_VEC dk (2.0*Ic*dIc - 3.0*dIIc);
    if (k < 1e-8) {
      const double lambda = sqrt((1.0/3.0)*Ic);
      const double dlambda = 0.5/(3.0*lambda);
      U  = I*lambda;
      dU = tensor_product(dIc*dlambda, dIc); // may not be correct
      return;
    }
   
    // Find the largest eigenvalue of C
    const double L = Ic*(Ic*Ic - 4.5*IIc) + 13.5*IIIc;
    DENS_VEC dL( (3.0*Ic*Ic-4.5*IIc)*dIc ); 
    dL.add_scaled(dIIc,  -4.5*Ic);
    dL.add_scaled(dIIIc, 13.5); 
    const double kpow  = pow(k,-1.5);
    const double dkpow = -1.5*kpow/k;
    const double phi   = acos(L*kpow); // phi - good
    
    // temporary factors for dphi
    const double d1 = -1.0/sqrt(1.0-L*L*kpow*kpow);
    const double d2 = d1*kpow;
    const double d3 = d1*L*dkpow;
    const DENS_VEC dphi (d2*dL + d3*dk);

    const double sqrt_k=sqrt(k), cos_p3i=cos((1.0/3.0)*phi);
    const double lam2  = (1.0/3.0)*(Ic + 2.0*sqrt_k*cos_p3i);

    DENS_VEC dlam2 = (1.0/3.0)*dIc;    
    dlam2.add_scaled(dk, (1.0/3.0)*cos_p3i/sqrt_k);
    dlam2.add_scaled(dphi, (-2.0/9.0)*sqrt_k*sin((1.0/3.0)*phi));
    const double lambda = sqrt(lam2);
    const DENS_VEC dlambda = (0.5/lambda)*dlam2;

    // Compute the invariants of U   
    const double IIIu  = sqrt(IIIc);
    const DENS_VEC dIIIu (0.5/IIIu*dIIIc);

    const double Iu  = lambda + sqrt(-lam2 + Ic + 2.0*IIIu/lambda);
    const double invrt = 1.0/(Iu-lambda);
    DENS_VEC dIu(dlambda);   dIu *= 1.0 + invrt*(-lambda - IIIu/lam2);
    dIu.add_scaled(dIc, 0.5*invrt);
    dIu.add_scaled(dIIIu, invrt/lambda);

    const double IIu  = 0.5*(Iu*Iu - Ic);
    const DENS_VEC dIIu ( Iu*dIu - 0.5*dIc );

    // Compute U and its derivatives
    const double fct = 1.0/(Iu*IIu-IIIu);
    DENS_VEC dfct = dIu;  dfct *= IIu;
    dfct.add_scaled(dIIu, Iu);
    dfct -= dIIIu;
    dfct *= -fct*fct;
    
    U = voigt3::eye(3, Iu*IIIu);
    U.add_scaled(C, Iu*Iu-IIu);
    U -= C2;
  
    DENS_MAT da = tensor_product(I, dIu);  da *= IIIu;
    da.add_scaled(tensor_product(I, dIIIu), Iu);
    da += tensor_product(C, 2.0*Iu*dIu-dIIu);
    da.add_scaled(eye<double>(6,6),Iu*Iu-IIu);
    da -= voigt3::derivative_of_square(C);
 
    dU = tensor_product(U, dfct);
    dU.add_scaled(da, fct);
    U *= fct; 
  }
  //============================================================================
  // Computes the dynamical matrix (TESTING FUNCTION)
  //============================================================================
  DENS_MAT compute_dynamical_matrix(const StressArgs &args)
  {
    DENS_MAT D(3,3);
    for (INDEX a=0; a<args.vac.size(); a++) {
      PairParam pair(args.vac.R(a), args.vac.r(a).norm()); 
      pair.phi_r = args.potential->phi_r(pair.d);
      pair.r = args.vac.r(a);
      pair.phi_rr  = args.potential->phi_rr(pair.d);    
      pair.phi_rrr = args.potential->phi_rrr(pair.d);    
      pairwise_thermal(pair, D);
    }
    return D;
  }
  //============================================================================
  // Computes the determinant of the dynamical matrix (TESTING FUNCTION)
  //============================================================================
  double compute_detD(const StressArgs &args)
  {
    return det(compute_dynamical_matrix(args));
  }
  //============================================================================
  // Computes the derivative of the dynamical matrix (TESTING FUNCTION)
  //============================================================================
  DENS_MAT_VEC compute_dynamical_derivative(StressArgs &args)
  {
    const double EPS = 1.0e-6;
    DENS_MAT_VEC dDdF(6, DENS_MAT(3,3));
    for (INDEX i=0; i<3; i++)  {
      for (INDEX j=0; j<3; j++)  {
        // store original F
        const double Fij = args.vac.F_(i,j);
        args.vac.F_(i,j) = Fij + EPS;
        DENS_MAT Da = compute_dynamical_matrix(args);
        args.vac.F_(i,j) = Fij - EPS;
        DENS_MAT Db = compute_dynamical_matrix(args);
        args.vac.F_(i,j) = Fij;
  
        dDdF[0](i,j) = (Da(0,0)-Db(0,0))*(0.5/EPS);
        dDdF[1](i,j) = (Da(1,1)-Db(1,1))*(0.5/EPS);
        dDdF[2](i,j) = (Da(2,2)-Db(2,2))*(0.5/EPS);
        dDdF[3](i,j) = (Da(1,2)-Db(1,2))*(0.5/EPS);
        dDdF[4](i,j) = (Da(0,2)-Db(0,2))*(0.5/EPS);
        dDdF[5](i,j) = (Da(0,1)-Db(0,1))*(0.5/EPS);
      }
    }
    return dDdF;
  }
}
