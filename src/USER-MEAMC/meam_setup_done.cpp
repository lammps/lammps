extern "C" {
#include "meam.h"
#include <math.h>

void alloyparams();
void compute_pair_meam();
double phi_meam(double r,int a, int b);
void compute_reference_density();
void get_shpfcn(double *s /* s(3) */, lattice_t latt);
void get_tavref(double *t11av,double *t21av,double *t31av,double *t12av,double *t22av,double *t32av, double t11,double t21,double t31,double t12,double t22,double t32, double r,int a,int b,lattice_t latt);
void get_Zij(int *Zij, lattice_t latt);
void get_Zij2(int *Zij2, double *a, double*S, lattice_t latt,double cmin,double cmax);
void get_sijk(double C,int i,int j,int k, double *sijk);
void get_densref(double r,int a,int b,double *rho01,double *rho11,double *rho21,double *rho31, double *rho02,double *rho12,double *rho22,double *rho32);
double zbl(double r, int z1, int z2);
double erose(double r,double re,double alpha,double Ec,double repuls,double attrac,int form);
void interpolate_meam(int ind);
double compute_phi(double rij, int elti, int eltj);

// in meam_dens_init
void fcut(double xi, double *fc);
// in meam_dens_final
void G_gam(double Gamma,int ibar,double gsmooth_factor, double *G, int *errorflag);

// Declaration in pair_meam.h:
//
// void meam_setup_done(double *)
//
// Call from pair_meam.cpp:
//
// meam_setup_done(&cutmax)
//

      void meam_setup_done_(double *cutmax)
      {
      int nv2, nv3, m, n, p;

//     Force cutoff
      meam_data.cutforce = meam_data.rc_meam;
      meam_data.cutforcesq = meam_data.cutforce*meam_data.cutforce;

//     Pass cutoff back to calling program
      *cutmax = meam_data.cutforce;

//     Augment t1 term
      for (int i=1; i<=maxelt; i++)
        meam_data.t1_meam[i] = meam_data.t1_meam[i] + meam_data.augt1 * 3.0/5.0 * meam_data.t3_meam[i];

//     Compute off-diagonal alloy parameters
      alloyparams();

// indices and factors for Voight notation
      nv2 = 1;
      nv3 = 1;
      for(m = 1; m<=3; m++) {
        for(n = m; n<=3; n++) {
          meam_data.vind2D[m][n] = nv2;
          meam_data.vind2D[n][m] = nv2;
          nv2 = nv2+1;
          for (p = n; p<=3; p++) {
            meam_data.vind3D[m][n][p] = nv3;
            meam_data.vind3D[m][p][n] = nv3;
            meam_data.vind3D[n][m][p] = nv3;
            meam_data.vind3D[n][p][m] = nv3;
            meam_data.vind3D[p][m][n] = nv3;
            meam_data.vind3D[p][n][m] = nv3;
            nv3 = nv3+1;
          }
        }
      }

      meam_data.v2D[1] = 1;
      meam_data.v2D[2] = 2;
      meam_data.v2D[3] = 2;
      meam_data.v2D[4] = 1;
      meam_data.v2D[5] = 2;
      meam_data.v2D[6] = 1;

      meam_data.v3D[1] = 1;
      meam_data.v3D[2] = 3;
      meam_data.v3D[3] = 3;
      meam_data.v3D[4] = 3;
      meam_data.v3D[5] = 6;
      meam_data.v3D[6] = 3;
      meam_data.v3D[7] = 1;
      meam_data.v3D[8] = 3;
      meam_data.v3D[9] = 3;
      meam_data.v3D[10] = 1;

      nv2 = 1;
      for(m = 1; m<=meam_data.neltypes; m++) {
        for(n = m; n<=meam_data.neltypes; n++) {
          meam_data.eltind[m][n] = nv2;
          meam_data.eltind[n][m] = nv2;
          nv2 = nv2+1;
        }
      }

//     Compute background densities for reference structure
      compute_reference_density();

//     Compute pair potentials and setup arrays for interpolation
      meam_data.nr = 1000;
      meam_data.dr = 1.1*meam_data.rc_meam/meam_data.nr;
      compute_pair_meam();
    }

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// Fill off-diagonal alloy parameters
      void alloyparams(void)
      {
      
      int i,j,k;
      double eb;

// Loop over pairs
      for(i = 1; i<=meam_data.neltypes; i++) {
        for(j = 1;i<=meam_data.neltypes; i++) {
// Treat off-diagonal pairs
// If i>j, set all equal to i<j case (which has aready been set,
// here or in the input file)
          if (i > j) {
            meam_data.re_meam[i][j] = meam_data.re_meam[j][i];
            meam_data.Ec_meam[i][j] = meam_data.Ec_meam[j][i];
            meam_data.alpha_meam[i][j] = meam_data.alpha_meam[j][i];
            meam_data.lattce_meam[i][j] = meam_data.lattce_meam[j][i];
            meam_data.nn2_meam[i][j] = meam_data.nn2_meam[j][i];
// If i<j and term is unset, use default values (e.g. mean of i-i and j-j)
          } else if (j > i) {
            if (iszero(meam_data.Ec_meam[i][j])) {
              if (meam_data.lattce_meam[i][j]==L12)
                meam_data.Ec_meam[i][j] = (3*meam_data.Ec_meam[i][i]+meam_data.Ec_meam[j][j])/4.0 - meam_data.delta_meam[i][j];
              else if (meam_data.lattce_meam[i][j]==C11) {
                if (meam_data.lattce_meam[i][i]==DIA)
                  meam_data.Ec_meam[i][j] = (2*meam_data.Ec_meam[i][i]+meam_data.Ec_meam[j][j])/3.0 - meam_data.delta_meam[i][j];
                else
                  meam_data.Ec_meam[i][j] = (meam_data.Ec_meam[i][i]+2*meam_data.Ec_meam[j][j])/3.0 - meam_data.delta_meam[i][j];
              } else
                meam_data.Ec_meam[i][j] = (meam_data.Ec_meam[i][i]+meam_data.Ec_meam[j][j])/2.0 - meam_data.delta_meam[i][j];
            }
            if (iszero(meam_data.alpha_meam[i][j]))
              meam_data.alpha_meam[i][j] = (meam_data.alpha_meam[i][i]+meam_data.alpha_meam[j][j])/2.0;
            if (iszero(meam_data.re_meam[i][j]))
              meam_data.re_meam[i][j] = (meam_data.re_meam[i][i]+meam_data.re_meam[j][j])/2.0;
          }
        }
      }

// Cmin[i][k][j] is symmetric in i-j, but not k.  For all triplets
// where i>j, set equal to the i<j element.  Likewise for Cmax.
      for(i = 2;i<=meam_data.neltypes;i++){
        for(j = 1;j<=i-1;j++){
          for(k = 1;k<=meam_data.neltypes;k++){
          meam_data.Cmin_meam[i][j][k] = meam_data.Cmin_meam[j][i][k];
          meam_data.Cmax_meam[i][j][k] = meam_data.Cmax_meam[j][i][k];
          }
        }
      }
        

// ebound gives the squared distance such that, for rik2 or rjk2>ebound,
// atom k definitely lies outside the screening function ellipse (so
// there is no need to calculate its effects).  Here, compute it for all
// triplets [i][j][k] so that ebound[i][j] is the maximized over k
      for(i = 2;i<=meam_data.neltypes;i++){
        for(j = 1;j<=meam_data.neltypes;j++){
          for(k = 1;k<=meam_data.neltypes;k++){
            eb = (meam_data.Cmax_meam[i][j][k]*meam_data.Cmax_meam[i][j][k]) / (4.0*(meam_data.Cmax_meam[i][j][k]-1.0));
            meam_data.ebound_meam[i][j] = max(meam_data.ebound_meam[i][j],eb);
          }
        }
      }
    }

//-----------------------------------------------------------------------
// compute MEAM pair potential for each pair of element types
//

      void compute_pair_meam(void)
      {

      double r/*ununsed:, temp*/;
      int j,a,b,nv2;
      double astar,frac,phizbl;
      int n,nmax,Z1,Z2;
      double arat,rarat,scrn,scrn2;
      double phiaa,phibb/*unused:,phitmp*/;
      double C,s111,s112,s221,S11,S22;

// check for previously allocated arrays and free them
      if(allocated(meam_data.phir)) deallocate(meam_data.phir);
      if(allocated(meam_data.phirar)) deallocate(meam_data.phirar);
      if(allocated(meam_data.phirar1)) deallocate(meam_data.phirar1);
      if(allocated(meam_data.phirar2)) deallocate(meam_data.phirar2);
      if(allocated(meam_data.phirar3)) deallocate(meam_data.phirar3);
      if(allocated(meam_data.phirar4)) deallocate(meam_data.phirar4);
      if(allocated(meam_data.phirar5)) deallocate(meam_data.phirar5);
      if(allocated(meam_data.phirar6)) deallocate(meam_data.phirar6);

// allocate memory for array that defines the potential
      allocate_2d(meam_data.phir,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);

// allocate coeff memory

      allocate_2d(meam_data.phirar,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);
      allocate_2d(meam_data.phirar1,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);
      allocate_2d(meam_data.phirar2,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);
      allocate_2d(meam_data.phirar3,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);
      allocate_2d(meam_data.phirar4,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);
      allocate_2d(meam_data.phirar5,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);
      allocate_2d(meam_data.phirar6,meam_data.nr,(meam_data.neltypes*(meam_data.neltypes+1))/2);

// loop over pairs of element types
      nv2 = 0;
      for(a=1; a<=meam_data.neltypes; a++) {
        for(b=a; b<=meam_data.neltypes; b++) {
          nv2 = nv2 + 1;

// loop over r values and compute
          for(j=1; j<=meam_data.nr; j++) {
            r = (j-1)*meam_data.dr;

            arr2(meam_data.phir,j,nv2) = phi_meam(r,a,b);

// if using second-nearest neighbor, solve recursive problem
// (see Lee and Baskes, PRB 62(13):8564 eqn.(21))
            if (meam_data.nn2_meam[a][b]==1) {
              get_Zij(&Z1,meam_data.lattce_meam[a][b]);
              get_Zij2(&Z2,&arat,&scrn,meam_data.lattce_meam[a][b],meam_data.Cmin_meam[a][a][b],meam_data.Cmax_meam[a][a][b]);

//     The B1, B2,  and L12 cases with NN2 have a trick to them; we need to
//     compute the contributions from second nearest neighbors, like a-a
//     pairs, but need to include NN2 contributions to those pairs as
//     well.
              if (meam_data.lattce_meam[a][b]==B1 || meam_data.lattce_meam[a][b]==B2 || meam_data.lattce_meam[a][b]==L12) {
                rarat = r*arat;

//               phi_aa
                phiaa = phi_meam(rarat,a,a);
                get_Zij(&Z1,meam_data.lattce_meam[a][a]);
                get_Zij2(&Z2,&arat,&scrn,meam_data.lattce_meam[a][a], meam_data.Cmin_meam[a][a][a],meam_data.Cmax_meam[a][a][a]);
                nmax = 10;
                if (scrn > 0.0) {
                  for(n=1; n<=nmax; n++) {
                    phiaa = phiaa + pow((-Z2*scrn/Z1),n) * phi_meam(rarat*pow(arat,n),a,a);
                  }
                }

//               phi_bb
                phibb = phi_meam(rarat,b,b);
                get_Zij(&Z1,meam_data.lattce_meam[b][b]);
                get_Zij2(&Z2,&arat,&scrn,meam_data.lattce_meam[b][b], meam_data.Cmin_meam[b][b][b],meam_data.Cmax_meam[b][b][b]);
                nmax = 10;
                if (scrn > 0.0) {
                  for(n=1; n<=nmax; n++) {
                    phibb = phibb + pow((-Z2*scrn/Z1),n) * phi_meam(rarat*pow(arat,n),b,b);
                  }
                }

                if (meam_data.lattce_meam[a][b]==B1 || meam_data.lattce_meam[a][b]==B2) {
//     Add contributions to the B1 or B2 potential
                  get_Zij(&Z1,meam_data.lattce_meam[a][b]);
                  get_Zij2(&Z2,&arat,&scrn,meam_data.lattce_meam[a][b], meam_data.Cmin_meam[a][a][b],meam_data.Cmax_meam[a][a][b]);
                  arr2(meam_data.phir,j,nv2) = arr2(meam_data.phir,j,nv2) - Z2*scrn/(2*Z1) * phiaa;
                  get_Zij2(&Z2,&arat,&scrn2,meam_data.lattce_meam[a][b], meam_data.Cmin_meam[b][b][a],meam_data.Cmax_meam[b][b][a]);
                  arr2(meam_data.phir,j,nv2) = arr2(meam_data.phir,j,nv2) - Z2*scrn2/(2*Z1) * phibb;

                } else if (meam_data.lattce_meam[a][b]==L12) {
//     The L12 case has one last trick; we have to be careful to compute
//     the correct screening between 2nd-neighbor pairs.  1-1
//     second-neighbor pairs are screened by 2 type 1 atoms and two type
//     2 atoms.  2-2 second-neighbor pairs are screened by 4 type 1
//     atoms.
                  C = 1.0;
                  get_sijk(C,a,a,a,&s111);
                  get_sijk(C,a,a,b,&s112);
                  get_sijk(C,b,b,a,&s221);
                  S11 = s111 * s111 * s112 * s112;
                  S22 = pow(s221,4);
                  arr2(meam_data.phir,j,nv2) = arr2(meam_data.phir,j,nv2) - 0.75*S11*phiaa - 0.25*S22*phibb;

                }

              } else {
                nmax = 10;
                for(n=1; n<=nmax; n++) {
                  arr2(meam_data.phir,j,nv2) = arr2(meam_data.phir,j,nv2) + pow((-Z2*scrn/Z1),n) * phi_meam(r*pow(arat,n),a,b);
                }
              }

            }

// For Zbl potential:
// if astar <= -3
//   potential is zbl potential
// else if -3 < astar < -1
//   potential is linear combination with zbl potential
// endif
            if (meam_data.zbl_meam[a][b]==1) {
              astar = meam_data.alpha_meam[a][b] * (r/meam_data.re_meam[a][b] - 1.0);
              if (astar <= -3.0)
                arr2(meam_data.phir,j,nv2) = zbl(r,meam_data.ielt_meam[a],meam_data.ielt_meam[b]);
              else if (astar > -3.0 && astar < -1.0) {
                fcut(1-(astar+1.0)/(-3.0+1.0),&frac);
                phizbl = zbl(r,meam_data.ielt_meam[a],meam_data.ielt_meam[b]);
                arr2(meam_data.phir,j,nv2) = frac*arr2(meam_data.phir,j,nv2) + (1-frac)*phizbl;
              }
            }

          }

// call interpolation
          interpolate_meam(nv2);

        }
      }

      }


//----------------------------------------------------------------------c
// Compute MEAM pair potential for distance r, element types a and b
//
      double phi_meam(double r,int a, int b)
      {
      /*unused:double a1,a2,a12;*/
      double t11av,t21av,t31av,t12av,t22av,t32av;
      double G1,G2,s1[3+1],s2[3+1]/*,s12[3+1]*/,rho0_1,rho0_2;
      double Gam1,Gam2,Z1,Z2;
      double rhobar1,rhobar2,F1,F2;
      double rho01,rho11,rho21,rho31;
      double rho02,rho12,rho22,rho32;
      double scalfac,phiaa,phibb;
      double Eu;
      double arat,scrn/*unused:,scrn2*/;
      int Z12, errorflag;
      int n,nmax,Z1nn,Z2nn;
      lattice_t latta/*unused:,lattb*/;
      double rho_bkgd1, rho_bkgd2;
      
      double phi_m = 0.0;

// Equation numbers below refer to:
//   I. Huang et.al., Modelling simul. Mater. Sci. Eng. 3:615

// get number of neighbors in the reference structure
//   Nref[i][j] = # of i's neighbors of type j
      get_Zij(&Z12,meam_data.lattce_meam[a][b]);

      get_densref(r,a,b,&rho01,&rho11,&rho21,&rho31,&rho02,&rho12,&rho22,&rho32);

// if densities are too small, numerical problems may result; just return zero
      if (rho01<=1e-14 && rho02<=1e-14)
        return 0.0;

// calculate average weighting factors for the reference structure
      if (meam_data.lattce_meam[a][b]==C11) {
        if (meam_data.ialloy==2) {
          t11av = meam_data.t1_meam[a];
          t12av = meam_data.t1_meam[b];
          t21av = meam_data.t2_meam[a];
          t22av = meam_data.t2_meam[b];
          t31av = meam_data.t3_meam[a];
          t32av = meam_data.t3_meam[b];
        } else {
          scalfac = 1.0/(rho01+rho02);
          t11av = scalfac*(meam_data.t1_meam[a]*rho01 + meam_data.t1_meam[b]*rho02);
          t12av = t11av;
          t21av = scalfac*(meam_data.t2_meam[a]*rho01 + meam_data.t2_meam[b]*rho02);
          t22av = t21av;
          t31av = scalfac*(meam_data.t3_meam[a]*rho01 + meam_data.t3_meam[b]*rho02);
          t32av = t31av;
        }
      } else {
// average weighting factors for the reference structure, eqn. I.8
         get_tavref(&t11av,&t21av,&t31av,&t12av,&t22av,&t32av, meam_data.t1_meam[a],meam_data.t2_meam[a],meam_data.t3_meam[a], meam_data.t1_meam[b],meam_data.t2_meam[b],meam_data.t3_meam[b], r,a,b,meam_data.lattce_meam[a][b]);
      }

// for c11b structure, calculate background electron densities
      if (meam_data.lattce_meam[a][b]==C11) {
         latta = meam_data.lattce_meam[a][a];
         if (latta==DIA) {
            rhobar1 = pow(((Z12/2)*(rho02+rho01)),2) + t11av*pow((rho12-rho11),2) + t21av/6.0*pow(rho22+rho21,2)+ 121.0/40.0*t31av*pow((rho32-rho31),2);
            rhobar1 = sqrt(rhobar1);
            rhobar2 = pow(Z12*rho01,2) + 2.0/3.0*t21av*pow(rho21,2);
            rhobar2 = sqrt(rhobar2);
         } else {
            rhobar2 = pow(((Z12/2)*(rho01+rho02)),2) + t12av*pow((rho11-rho12),2) + t22av/6.0*pow(rho21+rho22,2) + 121.0/40.0*t32av*pow((rho31-rho32),2);
            rhobar2 = sqrt(rhobar2);
            rhobar1 = pow(Z12*rho02,2) + 2.0/3.0*t22av*pow(rho22,2);
            rhobar1 = sqrt(rhobar1);
         }
      } else {
// for other structures, use formalism developed in Huang's paper
//
//     composition-dependent scaling, equation I.7
//     If using mixing rule for t, apply to reference structure; else
//     use precomputed values
        if (meam_data.mix_ref_t==1) {
          Z1 = meam_data.Z_meam[a];
          Z2 = meam_data.Z_meam[b];
          if (meam_data.ibar_meam[a]<=0)
            G1 = 1.0;
          else {
            get_shpfcn(s1,meam_data.lattce_meam[a][a]);
            Gam1 = (s1[1]*t11av+s1[2]*t21av+s1[3]*t31av)/(Z1*Z1);
            G_gam(Gam1,meam_data.ibar_meam[a],meam_data.gsmooth_factor,&G1,&errorflag);
          }
          if (meam_data.ibar_meam[b]<=0)
            G2 = 1.0;
          else {
            get_shpfcn(s2,meam_data.lattce_meam[b][b]);
            Gam2 = (s2[1]*t12av+s2[2]*t22av+s2[3]*t32av)/(Z2*Z2);
            G_gam(Gam2,meam_data.ibar_meam[b],meam_data.gsmooth_factor,&G2,&errorflag);
          }
          rho0_1 = meam_data.rho0_meam[a]*Z1*G1;
          rho0_2 = meam_data.rho0_meam[b]*Z2*G2;
        }
        Gam1 = (t11av*rho11+t21av*rho21+t31av*rho31);
        if (rho01 < 1.0e-14)
          Gam1 = 0.0;
        else
          Gam1 = Gam1/(rho01*rho01);

        Gam2 = (t12av*rho12+t22av*rho22+t32av*rho32);
        if (rho02 < 1.0e-14)
          Gam2 = 0.0;
        else
          Gam2 = Gam2/(rho02*rho02);

        G_gam(Gam1,meam_data.ibar_meam[a],meam_data.gsmooth_factor,&G1,&errorflag);
        G_gam(Gam2,meam_data.ibar_meam[b],meam_data.gsmooth_factor,&G2,&errorflag);
        if (meam_data.mix_ref_t==1) {
          rho_bkgd1 = rho0_1;
          rho_bkgd2 = rho0_2;
        } else {
          if (meam_data.bkgd_dyn==1) {
            rho_bkgd1 = meam_data.rho0_meam[a]*meam_data.Z_meam[a];
            rho_bkgd2 = meam_data.rho0_meam[b]*meam_data.Z_meam[b];
          } else {
            rho_bkgd1 = meam_data.rho_ref_meam[a];
            rho_bkgd2 = meam_data.rho_ref_meam[b];
          }
        }
        rhobar1 = rho01/rho_bkgd1*G1;
        rhobar2 = rho02/rho_bkgd2*G2;

      }

// compute embedding functions, eqn I.5
      if (iszero(rhobar1))
        F1 = 0.0;
      else {
        if (meam_data.emb_lin_neg==1 && rhobar1<=0)
          F1 = -meam_data.A_meam[a]*meam_data.Ec_meam[a][a]*rhobar1;
        else
          F1 = meam_data.A_meam[a]*meam_data.Ec_meam[a][a]*rhobar1*log(rhobar1);
      }
      if (iszero(rhobar2))
        F2 = 0.0;
      else {
        if (meam_data.emb_lin_neg==1  &&  rhobar2<=0)
          F2 = -meam_data.A_meam[b]*meam_data.Ec_meam[b][b]*rhobar2;
        else
          F2 = meam_data.A_meam[b]*meam_data.Ec_meam[b][b]*rhobar2*log(rhobar2);
      }

// compute Rose function, I.16
      Eu = erose(r,meam_data.re_meam[a][b],meam_data.alpha_meam[a][b], meam_data.Ec_meam[a][b],meam_data.repuls_meam[a][b],meam_data.attrac_meam[a][b],meam_data.erose_form);

// calculate the pair energy
      if (meam_data.lattce_meam[a][b]==C11) {
        latta = meam_data.lattce_meam[a][a];
        if (latta==DIA) {
          phiaa = phi_meam(r,a,a);
          phi_m = (3*Eu - F2 - 2*F1 - 5*phiaa)/Z12;
        } else {
          phibb = phi_meam(r,b,b);
          phi_m = (3*Eu - F1 - 2*F2 - 5*phibb)/Z12;
        }
      } else if (meam_data.lattce_meam[a][b]==L12) {
        phiaa = phi_meam(r,a,a);
//       account for second neighbor a-a potential here...
        get_Zij(&Z1nn,meam_data.lattce_meam[a][a]);
        get_Zij2(&Z2nn,&arat,&scrn,meam_data.lattce_meam[a][a],meam_data.Cmin_meam[a][a][a],meam_data.Cmax_meam[a][a][a]);
        nmax = 10;
        if (scrn > 0.0) {
          for(n=1; n<=nmax; n++) {
            phiaa = phiaa + pow((-Z2nn*scrn/Z1nn),n) * phi_meam(r*pow(arat,n),a,a);
          }
        }
        phi_m = Eu/3.0 - F1/4.0 - F2/12.0 - phiaa;
         
      } else {
//
// potential is computed from Rose function and embedding energy
         phi_m = (2*Eu - F1 - F2)/Z12;
//
      }

// if r = 0, just return 0
      if (iszero(r)) {
        phi_m = 0.0;
      }
      
      return phi_m;
    }

//----------------------------------------------------------------------c
// Compute background density for reference structure of each element
      void compute_reference_density(void)
      {
      int a,Z,Z2,errorflag;
      double gam,Gbar,shp[3+1];
      double rho0,rho0_2nn,arat,scrn;

// loop over element types
      for(a=1; a<=meam_data.neltypes; a++) {
        Z = (int)meam_data.Z_meam[a];
        if (meam_data.ibar_meam[a]<=0)
          Gbar = 1.0;
        else {
          get_shpfcn(shp,meam_data.lattce_meam[a][a]);
          gam = (meam_data.t1_meam[a]*shp[1]+meam_data.t2_meam[a]*shp[2]+meam_data.t3_meam[a]*shp[3])/(Z*Z);
          G_gam(gam,meam_data.ibar_meam[a],meam_data.gsmooth_factor,&Gbar,&errorflag);
        }

//     The zeroth order density in the reference structure, with
//     equilibrium spacing, is just the number of first neighbors times
//     the rho0_meam coefficient...
        rho0 = meam_data.rho0_meam[a]*Z;

//     ...unless we have unscreened second neighbors, in which case we
//     add on the contribution from those (accounting for partial
//     screening)
        if (meam_data.nn2_meam[a][a]==1) {
          get_Zij2(&Z2,&arat,&scrn,meam_data.lattce_meam[a][a],meam_data.Cmin_meam[a][a][a],meam_data.Cmax_meam[a][a][a]);
          rho0_2nn = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta0_meam[a]*(arat-1));
          rho0 = rho0 + Z2*rho0_2nn*scrn;
        }

        meam_data.rho_ref_meam[a] = rho0*Gbar;
      }
      }

//----------------------------------------------------------------------c
// Shape factors for various configurations
      void get_shpfcn(double *s /* s(3) */, lattice_t latt)
      {
      if (latt==FCC || latt==BCC || latt==B1 || latt==B2) {
        s[1] = 0.0;
        s[2] = 0.0;
        s[3] = 0.0;
      } else if (latt==HCP) {
        s[1] = 0.0;
        s[2] = 0.0;
        s[3] = 1.0/3.0;
      } else if (latt==DIA) {
        s[1] = 0.0;
        s[2] = 0.0;
        s[3] = 32.0/9.0;
      } else if (latt==DIM) {
        s[1] = 1.0;
        s[2] = 2.0/3.0;
//        s(3) = 1.d0
        s[3] = 0.40;
      } else {
        s[1] = 0.0;
//        call error('Lattice not defined in get_shpfcn.')
      }
      }
      
//------------------------------------------------------------------------------c
// Average weighting factors for the reference structure
      void get_tavref(double *t11av,double *t21av,double *t31av,double *t12av,double *t22av,double *t32av, double t11,double t21,double t31,double t12,double t22,double t32, double r,int a,int b,lattice_t latt)
      {
      double rhoa01,rhoa02,a1,a2,rho01/*,rho02*/;

//     For ialloy = 2, no averaging is done
      if (meam_data.ialloy==2) {
          *t11av = t11;
          *t21av = t21;
          *t31av = t31;
          *t12av = t12;
          *t22av = t22;
          *t32av = t32;
      } else {
        if (latt==FCC || latt==BCC || latt==DIA || latt==HCP || latt==B1 || latt==DIM || latt==B2) {
//     all neighbors are of the opposite type
          *t11av = t12;
          *t21av = t22;
          *t31av = t32;
          *t12av = t11;
          *t22av = t21;
          *t32av = t31;
        } else {
          a1  = r/meam_data.re_meam[a][a] - 1.0;
          a2  = r/meam_data.re_meam[b][b] - 1.0;
          rhoa01 = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta0_meam[a]*a1);
          rhoa02 = meam_data.rho0_meam[b]*fm_exp(-meam_data.beta0_meam[b]*a2);
          if (latt==L12) {
            rho01 = 8*rhoa01 + 4*rhoa02;
            *t11av = (8*t11*rhoa01 + 4*t12*rhoa02)/rho01;
            *t12av = t11;
            *t21av = (8*t21*rhoa01 + 4*t22*rhoa02)/rho01;
            *t22av = t21;
            *t31av = (8*t31*rhoa01 + 4*t32*rhoa02)/rho01;
            *t32av = t31;
          } else {
//      call error('Lattice not defined in get_tavref.')
          }
        }
      }
      }
      
//------------------------------------------------------------------------------c
// Number of neighbors for the reference structure
      void get_Zij(int *Zij, lattice_t latt)
      {
      if (latt==FCC)
        *Zij = 12;
      else if (latt==BCC)
        *Zij = 8;
      else if (latt==HCP)
        *Zij = 12;
      else if (latt==B1)
        *Zij = 6;
      else if (latt==DIA)
        *Zij = 4;
      else if (latt==DIM)
        *Zij = 1;
      else if (latt==C11)
        *Zij = 10;
      else if (latt==L12)
        *Zij = 12;
      else if (latt==B2)
        *Zij = 8;
      else {
//        call error('Lattice not defined in get_Zij.')
      }
      }
     
//------------------------------------------------------------------------------c
// Zij2 = number of second neighbors, a = distance ratio R1/R2, and S = second
// neighbor screening function for lattice type "latt"

      void get_Zij2(int *Zij2, double *a, double*S, lattice_t latt,double cmin,double cmax)
      {
        
      double /*rratio,*/C,x,sijk;
      int numscr = 0;

      if (latt==BCC) {
        *Zij2 = 6;
        *a = 2.0/sqrt(3.0);
        numscr = 4;
      } else if (latt==FCC) {
        *Zij2 = 6;
        *a = sqrt(2.0);
        numscr = 4;
      } else if (latt==DIA) {
        *Zij2 = 0;
        *a = sqrt(8.0/3.0);
        numscr = 4;
        if (cmin < 0.500001) {
//          call error('can not do 2NN MEAM for dia')
        }
      } else if (latt==HCP) {
        *Zij2 = 6;
        *a = sqrt(2.0);
        numscr = 4;
      } else if (latt==B1) {
        *Zij2 = 12;
        *a = sqrt(2.0);
        numscr = 2;
      } else if (latt==L12) {
        *Zij2 = 6;
        *a = sqrt(2.0);
        numscr = 4;
      } else if (latt==B2) {
        *Zij2 = 6;
        *a = 2.0/sqrt(3.0);
        numscr = 4;
      } else if (latt==DIM) {
//        this really shouldn't be allowed; make sure screening is zero
         *Zij2 = 0;
         *a = 1;
         *S = 0;
         return;
      } else {
//        call error('Lattice not defined in get_Zij2.')
      }

// Compute screening for each first neighbor
      C = 4.0/(*a * *a) - 1.0;
      x = (C-cmin)/(cmax-cmin);
      fcut(x,&sijk);
// There are numscr first neighbors screening the second neighbors
      *S = pow(sijk,numscr);
      
      }
      


//------------------------------------------------------------------------------c
      void get_sijk(double C,int i,int j,int k, double *sijk)
      {
      double x;
      x = (C-meam_data.Cmin_meam[i][j][k])/(meam_data.Cmax_meam[i][j][k]-meam_data.Cmin_meam[i][j][k]);
      fcut(x,sijk);
      }

//------------------------------------------------------------------------------c
// Calculate density functions, assuming reference configuration
      void get_densref(double r,int a,int b,double *rho01,double *rho11,double *rho21,double *rho31, double *rho02,double *rho12,double *rho22,double *rho32)
      {
      double a1,a2;
      double s[3+1];
      lattice_t lat;
      int  Zij1nn,Zij2nn;
      double rhoa01nn,rhoa02nn;
      double rhoa01,rhoa11,rhoa21,rhoa31;
      double rhoa02,rhoa12,rhoa22,rhoa32;
      double arat,scrn,denom;
      double C,s111,s112,s221,S11,S22;

      a1  = r/meam_data.re_meam[a][a] - 1.0;
      a2  = r/meam_data.re_meam[b][b] - 1.0;

      rhoa01 = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta0_meam[a]*a1);
      rhoa11 = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta1_meam[a]*a1);
      rhoa21 = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta2_meam[a]*a1);
      rhoa31 = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta3_meam[a]*a1);
      rhoa02 = meam_data.rho0_meam[b]*fm_exp(-meam_data.beta0_meam[b]*a2);
      rhoa12 = meam_data.rho0_meam[b]*fm_exp(-meam_data.beta1_meam[b]*a2);
      rhoa22 = meam_data.rho0_meam[b]*fm_exp(-meam_data.beta2_meam[b]*a2);
      rhoa32 = meam_data.rho0_meam[b]*fm_exp(-meam_data.beta3_meam[b]*a2);

      lat = meam_data.lattce_meam[a][b];

      *rho11 = 0.0;
      *rho21 = 0.0;
      *rho31 = 0.0;
      *rho12 = 0.0;
      *rho22 = 0.0;
      *rho32 = 0.0;

      get_Zij(&Zij1nn,lat);

      if (lat==FCC) {
        *rho01 = 12.0*rhoa02;
        *rho02 = 12.0*rhoa01;
      } else if (lat==BCC) {
        *rho01 = 8.0*rhoa02;
        *rho02 = 8.0*rhoa01;
      } else if (lat==B1) {
        *rho01 = 6.0*rhoa02;
        *rho02 = 6.0*rhoa01;
      } else if (lat==DIA) {
        *rho01 = 4.0*rhoa02;
        *rho02 = 4.0*rhoa01;
        *rho31 = 32.0/9.0*rhoa32*rhoa32;
        *rho32 = 32.0/9.0*rhoa31*rhoa31;
      } else if (lat==HCP) {
        *rho01 = 12*rhoa02;
        *rho02 = 12*rhoa01;
        *rho31 = 1.0/3.0*rhoa32*rhoa32;
        *rho32 = 1.0/3.0*rhoa31*rhoa31;
      } else if (lat==DIM) {
        get_shpfcn(s,DIM);
        *rho01 = rhoa02;
        *rho02 = rhoa01;
        *rho11 = s[1]*rhoa12*rhoa12;
        *rho12 = s[1]*rhoa11*rhoa11;
        *rho21 = s[2]*rhoa22*rhoa22;
        *rho22 = s[2]*rhoa21*rhoa21;
        *rho31 = s[3]*rhoa32*rhoa32;
        *rho32 = s[3]*rhoa31*rhoa31;
      } else if (lat==C11) {
        *rho01 = rhoa01;
        *rho02 = rhoa02;
        *rho11 = rhoa11;
        *rho12 = rhoa12;
        *rho21 = rhoa21;
        *rho22 = rhoa22;
        *rho31 = rhoa31;
        *rho32 = rhoa32;
      } else if (lat==L12) {
        *rho01 = 8*rhoa01 + 4*rhoa02;
        *rho02 = 12*rhoa01;
        if (meam_data.ialloy==1) {
          *rho21 = 8./3.*pow(rhoa21*meam_data.t2_meam[a]-rhoa22*meam_data.t2_meam[b],2);
          denom = 8*rhoa01*pow(meam_data.t2_meam[a],2) + 4*rhoa02*pow(meam_data.t2_meam[b],2);
          if (denom > 0.)
            *rho21 = *rho21/denom * *rho01;
        } else
          *rho21 = 8./3.*(rhoa21-rhoa22)*(rhoa21-rhoa22);
      } else if (lat==B2) {
        *rho01 = 8.0*rhoa02;
        *rho02 = 8.0*rhoa01;
      } else {
//        call error('Lattice not defined in get_densref.')
      }

      if (meam_data.nn2_meam[a][b]==1) {

        get_Zij2(&Zij2nn,&arat,&scrn,lat,meam_data.Cmin_meam[a][a][b],meam_data.Cmax_meam[a][a][b]);

        a1 = arat*r/meam_data.re_meam[a][a] - 1.0;
        a2 = arat*r/meam_data.re_meam[b][b] - 1.0;

        rhoa01nn = meam_data.rho0_meam[a]*fm_exp(-meam_data.beta0_meam[a]*a1);
        rhoa02nn = meam_data.rho0_meam[b]*fm_exp(-meam_data.beta0_meam[b]*a2);

        if (lat==L12) {
//     As usual, L12 thinks it's special; we need to be careful computing
//     the screening functions
          C = 1.0;
          get_sijk(C,a,a,a,&s111);
          get_sijk(C,a,a,b,&s112);
          get_sijk(C,b,b,a,&s221);
          S11 = s111 * s111 * s112 * s112;
          S22 = pow(s221,4);
          *rho01 = *rho01 + 6*S11*rhoa01nn;
          *rho02 = *rho02 + 6*S22*rhoa02nn;

        } else {
//     For other cases, assume that second neighbor is of same type,
//     first neighbor may be of different type

          *rho01 = *rho01 + Zij2nn*scrn*rhoa01nn;

//     Assume Zij2nn and arat don't depend on order, but scrn might
          get_Zij2(&Zij2nn,&arat,&scrn,lat,meam_data.Cmin_meam[b][b][a],meam_data.Cmax_meam[b][b][a]);
          *rho02 = *rho02 + Zij2nn*scrn*rhoa02nn;
        }
      }
      }

//---------------------------------------------------------------------
// Compute ZBL potential
//
      double zbl(double r, int z1, int z2)
      {
      int i;
      const double c[]={0.028171,0.28022,0.50986,0.18175};
      const double d[]={0.20162,0.40290,0.94229,3.1998};
      const double azero=0.4685;
      const double cc=14.3997;
      double a,x;
// azero = (9pi^2/128)^1/3 (0.529) Angstroms
      a = azero/(pow(z1,0.23)+pow(z2,0.23));
      double result = 0.0;
      x = r/a;
      for(i=0; i<=3; i++) {
        result = result + c[i]*fm_exp(-d[i]*x);
      }
      if (r > 0.0) result = result*z1*z2/r*cc;
      return result;
      }

//---------------------------------------------------------------------
// Compute Rose energy function, I.16
//
      double erose(double r,double re,double alpha,double Ec,double repuls,double attrac,int form)
      {
      double astar,a3;
      double result = 0.0;

      if (r > 0.0) {
        astar = alpha * (r/re - 1.0);
        a3 = 0.0;
        if (astar>=0)
          a3 = attrac;
        else if (astar < 0)
          a3 = repuls;

        if (form==1)
          result = -Ec*(1+astar+(-attrac+repuls/r)* pow(astar,3))*fm_exp(-astar);
        else if (form==2)
          result = -Ec * (1 +astar + a3*pow(astar,3))*fm_exp(-astar);
        else
          result = -Ec * (1+ astar + a3*pow(astar,3)/(r/re))*fm_exp(-astar);
       
      }
      return result;
      }

// -----------------------------------------------------------------------

      void interpolate_meam(int ind)
      {
      int j;
      double drar;

// map to coefficient space

      meam_data.nrar = meam_data.nr;
      drar = meam_data.dr;
      meam_data.rdrar = 1.0/drar;

// phir interp
      for(j=1; j<=meam_data.nrar; j++) {
        arr2(meam_data.phirar,j,ind) = arr2(meam_data.phir,j,ind);
      }
      arr2(meam_data.phirar1,1,ind) = arr2(meam_data.phirar,2,ind)-arr2(meam_data.phirar,1,ind);
      arr2(meam_data.phirar1,2,ind) = 0.5*(arr2(meam_data.phirar,3,ind)-arr2(meam_data.phirar,1,ind));
      arr2(meam_data.phirar1,meam_data.nrar-1,ind) = 0.5*(arr2(meam_data.phirar,meam_data.nrar,ind)      -arr2(meam_data.phirar,meam_data.nrar-2,ind));
      arr2(meam_data.phirar1,meam_data.nrar,ind) = 0.0;
      for(j=3; j<=meam_data.nrar-2; j++) {
        arr2(meam_data.phirar1,j,ind) = ((arr2(meam_data.phirar,j-2,ind)-arr2(meam_data.phirar,j+2,ind)) +        8.0*(arr2(meam_data.phirar,j+1,ind)-arr2(meam_data.phirar,j-1,ind)))/12.;
      }

      for(j=1; j<=meam_data.nrar-1; j++) {
        arr2(meam_data.phirar2,j,ind) = 3.0*(arr2(meam_data.phirar,j+1,ind)-arr2(meam_data.phirar,j,ind)) -        2.0*arr2(meam_data.phirar1,j,ind) - arr2(meam_data.phirar1,j+1,ind);
        arr2(meam_data.phirar3,j,ind) = arr2(meam_data.phirar1,j,ind) + arr2(meam_data.phirar1,j+1,ind) -        2.0*(arr2(meam_data.phirar,j+1,ind)-arr2(meam_data.phirar,j,ind));
      }
      arr2(meam_data.phirar2,meam_data.nrar,ind) = 0.0;
      arr2(meam_data.phirar3,meam_data.nrar,ind) = 0.0;

      for(j=1; j<=meam_data.nrar; j++) {
        arr2(meam_data.phirar4,j,ind) = arr2(meam_data.phirar1,j,ind)/drar;
        arr2(meam_data.phirar5,j,ind) = 2.0*arr2(meam_data.phirar2,j,ind)/drar;
        arr2(meam_data.phirar6,j,ind) = 3.0*arr2(meam_data.phirar3,j,ind)/drar;
      }  
           
      }

//---------------------------------------------------------------------
// Compute Rose energy function, I.16
//
      double compute_phi(double rij, int elti, int eltj)
      {
      double pp;
      int ind, kk;

      ind = meam_data.eltind[elti][eltj];
      pp = rij*meam_data.rdrar + 1.0;
      kk = (int)pp;
      kk = min(kk,meam_data.nrar-1);
      pp = pp - kk;
      pp = min(pp,1.0);
      double result = ((arr2(meam_data.phirar3,kk,ind)*pp + arr2(meam_data.phirar2,kk,ind))*pp + arr2(meam_data.phirar1,kk,ind))*pp + arr2(meam_data.phirar,kk,ind);

      return result;
      }


}
