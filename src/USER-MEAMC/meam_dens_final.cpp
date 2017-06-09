extern "C"{
#include "meam.h"
#include <math.h>

      void G_gam(double Gamma, int ibar, double gsmooth_factor, double *G, int *errorflag);
      void dG_gam(double Gamma, int ibar, double gsmooth_factor,double *G, double *dG);

// in meam_setup_done
      void get_shpfcn(double *s /* s(3) */, lattice_t latt);

// Extern "C" declaration has the form:
//
//  void meam_dens_final_(int *, int *, int *, int *, int *, double *, double *,
//                int *, int *, int *,
//		 double *, double *, double *, double *, double *, double *,
//		 double *, double *, double *, double *, double *, double *,
//		 double *, double *, double *, double *, double *, int *);
//
// Call from pair_meam.cpp has the form:
//
//  meam_dens_final_(&nlocal,&nmax,&eflag_either,&eflag_global,&eflag_atom,
//            &eng_vdwl,eatom,ntype,type,fmap,
//	     &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],
//	     &arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],gamma,dgamma1,
//	     dgamma2,dgamma3,rho,rho0,rho1,rho2,rho3,frhop,&errorflag);
//

      void meam_dens_final_(int *nlocal, int *nmax,
           int *eflag_either, int *eflag_global, int *eflag_atom, double *eng_vdwl, double *eatom,
           int *ntype, int *type, int *fmap,
           double *Arho1, double *Arho2, double *Arho2b, double *Arho3, double *Arho3b, double *t_ave, double *tsq_ave,
           double *Gamma, double *dGamma1, double *dGamma2, double *dGamma3,
           double *rho, double *rho0, double *rho1, double *rho2, double *rho3, double *fp, int *errorflag)
      {
        
      arrdim2v(Arho1,3,*nmax)
      arrdim2v(Arho2,6,*nmax)
      arrdim2v(Arho3,10,*nmax)
      arrdim2v(Arho3b,3,*nmax)
      arrdim2v(t_ave,3,*nmax)
      arrdim2v(tsq_ave,3,*nmax)

      int i, elti;
      int m;
      double  rhob, G, dG, Gbar, dGbar, gam, shp[3+1], Z;
      double  B, denom, rho_bkgd;

//     Complete the calculation of density

      for(i=1; i<=*nlocal; i++) {
        elti = arr1v(fmap,arr1v(type,i));
        if (elti > 0) {
          arr1v(rho1,i) = 0.0;
          arr1v(rho2,i) = -1.0/3.0*arr1v(Arho2b,i)*arr1v(Arho2b,i);
          arr1v(rho3,i) = 0.0;
          for(m=1; m<=3; m++) {
            arr1v(rho1,i) = arr1v(rho1,i) + arr2v(Arho1,m,i)*arr2v(Arho1,m,i);
            arr1v(rho3,i) = arr1v(rho3,i) - 3.0/5.0*arr2v(Arho3b,m,i)*arr2v(Arho3b,m,i);
          }
          for(m=1; m<=6; m++) {
            arr1v(rho2,i) = arr1v(rho2,i) + meam_data.v2D[m]*arr2v(Arho2,m,i)*arr2v(Arho2,m,i);
          }
          for(m=1; m<=10; m++) {
            arr1v(rho3,i) = arr1v(rho3,i) + meam_data.v3D[m]*arr2v(Arho3,m,i)*arr2v(Arho3,m,i);
          }

          if( arr1v(rho0,i)  >  0.0 ) {
            if (meam_data.ialloy==1) {
              arr2v(t_ave,1,i) = arr2v(t_ave,1,i)/arr2v(tsq_ave,1,i);
              arr2v(t_ave,2,i) = arr2v(t_ave,2,i)/arr2v(tsq_ave,2,i);
              arr2v(t_ave,3,i) = arr2v(t_ave,3,i)/arr2v(tsq_ave,3,i);
            } else if (meam_data.ialloy==2) {
              arr2v(t_ave,1,i) = meam_data.t1_meam[elti];
              arr2v(t_ave,2,i) = meam_data.t2_meam[elti];
              arr2v(t_ave,3,i) = meam_data.t3_meam[elti];
            } else {
              arr2v(t_ave,1,i) = arr2v(t_ave,1,i)/arr1v(rho0,i);
              arr2v(t_ave,2,i) = arr2v(t_ave,2,i)/arr1v(rho0,i);
              arr2v(t_ave,3,i) = arr2v(t_ave,3,i)/arr1v(rho0,i);
            }
          }

          arr1v(Gamma,i) = arr2v(t_ave,1,i)*arr1v(rho1,i) + arr2v(t_ave,2,i)*arr1v(rho2,i) + arr2v(t_ave,3,i)*arr1v(rho3,i);

          if( arr1v(rho0,i)  >  0.0 ) {
            arr1v(Gamma,i) = arr1v(Gamma,i)/(arr1v(rho0,i)*arr1v(rho0,i));
          }

          Z = meam_data.Z_meam[elti];

          G_gam(arr1v(Gamma,i),meam_data.ibar_meam[elti], meam_data.gsmooth_factor,&G,errorflag);
          if (*errorflag!=0) return;
          get_shpfcn(shp,meam_data.lattce_meam[elti][elti]);
          if (meam_data.ibar_meam[elti]<=0) {
            Gbar = 1.0;
            dGbar = 0.0;
          } else {
            if (meam_data.mix_ref_t==1) {
              gam = (arr2v(t_ave,1,i)*shp[1]+arr2v(t_ave,2,i)*shp[2] +arr2v(t_ave,3,i)*shp[3])/(Z*Z);
            } else {
              gam = (meam_data.t1_meam[elti]*shp[1]+meam_data.t2_meam[elti]*shp[2] +meam_data.t3_meam[elti]*shp[3])/(Z*Z);
            }
            G_gam(gam,meam_data.ibar_meam[elti],meam_data.gsmooth_factor,&Gbar,errorflag);
          }
          arr1v(rho,i) = arr1v(rho0,i) * G;

          if (meam_data.mix_ref_t==1) {
            if (meam_data.ibar_meam[elti]<=0) {
              Gbar = 1.0;
              dGbar = 0.0;
            } else {
              gam = (arr2v(t_ave,1,i)*shp[1]+arr2v(t_ave,2,i)*shp[2] +arr2v(t_ave,3,i)*shp[3])/(Z*Z);
              dG_gam(gam,meam_data.ibar_meam[elti],meam_data.gsmooth_factor, &Gbar,&dGbar);
            }
            rho_bkgd = meam_data.rho0_meam[elti]*Z*Gbar;
          } else {
            if (meam_data.bkgd_dyn==1) {
              rho_bkgd = meam_data.rho0_meam[elti]*Z;
            } else {
              rho_bkgd = meam_data.rho_ref_meam[elti];
            }
          }
          rhob = arr1v(rho,i)/rho_bkgd;
          denom = 1.0/rho_bkgd;

          dG_gam(arr1v(Gamma,i),meam_data.ibar_meam[elti],meam_data.gsmooth_factor,&G,&dG);

          arr1v(dGamma1,i) = (G - 2*dG*arr1v(Gamma,i))*denom;

          if( !iszero(arr1v(rho0,i)) ) {
            arr1v(dGamma2,i) = (dG/arr1v(rho0,i))*denom;
          } else {
            arr1v(dGamma2,i) = 0.0;
          }

//     dGamma3 is nonzero only if we are using the "mixed" rule for
//     computing t in the reference system (which is not correct, but
//     included for backward compatibility
          if (meam_data.mix_ref_t==1) {
            arr1v(dGamma3,i) = arr1v(rho0,i)*G*dGbar/(Gbar*Z*Z)*denom;
          } else {
            arr1v(dGamma3,i) = 0.0;
          }

          B = meam_data.A_meam[elti]*meam_data.Ec_meam[elti][elti];

          if( !iszero(rhob) ) {
            if (meam_data.emb_lin_neg==1 && rhob<=0) {
              arr1v(fp,i) = -B;
            } else {
              arr1v(fp,i) = B*(log(rhob)+1.0);
            }
            if (*eflag_either!=0) {
              if (*eflag_global!=0) {
                if (meam_data.emb_lin_neg==1 && rhob<=0) {
                  *eng_vdwl = *eng_vdwl - B*rhob;
                } else {
                  *eng_vdwl = *eng_vdwl + B*rhob*log(rhob);
                }
              }
              if (*eflag_atom!=0) {
                if (meam_data.emb_lin_neg==1 && rhob<=0) {
                  arr1v(eatom,i) = arr1v(eatom,i) - B*rhob;
                } else {
                  arr1v(eatom,i) = arr1v(eatom,i) + B*rhob*log(rhob);
                }
              }
            }
          } else {
            if (meam_data.emb_lin_neg==1) {
              arr1v(fp,i) = -B;
            } else {
              arr1v(fp,i) = B;
            }
          }
        }
      }

      }

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      void G_gam(double Gamma, int ibar, double gsmooth_factor, double *G, int *errorflag)
      {
//     Compute G(Gamma) based on selection flag ibar:
//   0 => G = sqrt(1+Gamma)
//   1 => G = exp(Gamma/2)
//   2 => not implemented
//   3 => G = 2/(1+exp(-Gamma))
//   4 => G = sqrt(1+Gamma)
//  -5 => G = +-sqrt(abs(1+Gamma))
      
      double gsmooth_switchpoint;
      if (ibar==0 || ibar==4) {
        gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor+1);
        if (Gamma < gsmooth_switchpoint) {
//         e.g. gsmooth_factor is 99, {:
//         gsmooth_switchpoint = -0.99
//         G = 0.01*(-0.99/Gamma)**99
          *G = 1/(gsmooth_factor+1) * pow((gsmooth_switchpoint/Gamma),gsmooth_factor);
          *G = sqrt(*G);
        } else {
          *G = sqrt(1.0+Gamma);
        }
      } else if (ibar==1) {
        *G = fm_exp(Gamma/2.0);
      } else if (ibar==3) {
        *G = 2.0/(1.0+exp(-Gamma));
      } else if (ibar==-5) {
        if ((1.0+Gamma)>=0) {
           *G = sqrt(1.0+Gamma);
        } else {
           *G = -sqrt(-1.0-Gamma);
        }
      } else {
         *errorflag = 1;
      }
      }

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      void dG_gam(double Gamma, int ibar, double gsmooth_factor,double *G, double *dG)
      {
// Compute G(Gamma) and dG(gamma) based on selection flag ibar:
//   0 => G = sqrt(1+Gamma)
//   1 => G = fm_exp(Gamma/2)
//   2 => not implemented
//   3 => G = 2/(1+fm_exp(-Gamma))
//   4 => G = sqrt(1+Gamma)
//  -5 => G = +-sqrt(abs(1+Gamma))
      double gsmooth_switchpoint;
      
      if (ibar==0 || ibar==4) {
        gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor+1);
        if (Gamma < gsmooth_switchpoint) {
//         e.g. gsmooth_factor is 99, {:
//         gsmooth_switchpoint = -0.99
//         G = 0.01*(-0.99/Gamma)**99
          *G = 1/(gsmooth_factor+1) * pow((gsmooth_switchpoint/Gamma), gsmooth_factor);
          *G = sqrt(*G);
          *dG = -gsmooth_factor* *G/(2.0*Gamma);
        } else {
          *G = sqrt(1.0+Gamma);
          *dG = 1.0/(2.0* *G);
        }
      } else if (ibar==1) {
        *G = fm_exp(Gamma/2.0);
        *dG = *G/2.0;
      } else if (ibar==3) {
        *G = 2.0/(1.0+fm_exp(-Gamma));
        *dG = *G*(2.0-*G)/2;
      } else if (ibar==-5) {
        if ((1.0+Gamma)>=0) {
           *G = sqrt(1.0+Gamma);
           *dG = 1.0/(2.0* *G);
        } else {
           *G = -sqrt(-1.0-Gamma);
           *dG = -1.0/(2.0* *G);
        }
      }
      }

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

}