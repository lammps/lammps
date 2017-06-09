#include "meam.h"
#include <math.h>

//
// declaration in pair_meam.h:
//
//  void meam_setup_global(int *, int *, double *, int *, double *, double *,
//			 double *, double *, double *, double *, double *,
//			 double *, double *, double *, double *, double *,
//			 double *, double *, int *);
//
// call in pair_meam.cpp:
//
//  meam_setup_global(&nelements,lat,z,ielement,atwt,alpha,b0,b1,b2,b3,
//		    alat,esub,asub,t0,t1,t2,t3,rozero,ibar);
//
//

      void meam_setup_global_(int *nelt, int *lat, double *z, int *ielement, double *atwt, double *alpha,
          double *b0, double *b1, double *b2, double *b3, double *alat, double *esub, double *asub,
          double *t0, double *t1, double *t2, double *t3, double *rozero, int *ibar)
      {
      
      int i;
      double tmplat[maxelt];

      meam_data.neltypes = *nelt;

      for(i=1; i<=*nelt; i++) {
         if (arr1v(lat,i)==0)
            meam_data.lattce_meam[i][i] = FCC;
         else if (arr1v(lat,i)==1)
            meam_data.lattce_meam[i][i] = BCC;
         else if (arr1v(lat,i)==2)
            meam_data.lattce_meam[i][i] = HCP;
         else if (arr1v(lat,i)==3)
            meam_data.lattce_meam[i][i] = DIM;
         else if (arr1v(lat,i)==4)
            meam_data.lattce_meam[i][i] = DIA;
         else {
//           unknown
         }

         meam_data.Z_meam[i] = arr1v(z,i);
         meam_data.ielt_meam[i] = arr1v(ielement,i);
         meam_data.alpha_meam[i][i] = arr1v(alpha,i);
         meam_data.beta0_meam[i] = arr1v(b0,i);
         meam_data.beta1_meam[i] = arr1v(b1,i);
         meam_data.beta2_meam[i] = arr1v(b2,i);
         meam_data.beta3_meam[i] = arr1v(b3,i);
         tmplat[i] = arr1v(alat,i);
         meam_data.Ec_meam[i][i] = arr1v(esub,i);
         meam_data.A_meam[i] = arr1v(asub,i);
         meam_data.t0_meam[i] = arr1v(t0,i);
         meam_data.t1_meam[i] = arr1v(t1,i);
         meam_data.t2_meam[i] = arr1v(t2,i);
         meam_data.t3_meam[i] = arr1v(t3,i);
         meam_data.rho0_meam[i] = arr1v(rozero,i);
         meam_data.ibar_meam[i] = arr1v(ibar,i);

         if (meam_data.lattce_meam[i][i]==FCC)
            meam_data.re_meam[i][i] = tmplat[i]/sqrt(2.0);
         else if (meam_data.lattce_meam[i][i]==BCC)
            meam_data.re_meam[i][i] = tmplat[i]*sqrt(3.0)/2.0;
         else if (meam_data.lattce_meam[i][i]==HCP)
            meam_data.re_meam[i][i] = tmplat[i];
         else if (meam_data.lattce_meam[i][i]==DIM)
            meam_data.re_meam[i][i] = tmplat[i];
         else if (meam_data.lattce_meam[i][i]==DIA)
            meam_data.re_meam[i][i] = tmplat[i]*sqrt(3.0)/4.0;
         else {
//           error
         }

      }


// Set some defaults
      meam_data.rc_meam = 4.0;
      meam_data.delr_meam = 0.1;
      setall2d(meam_data.attrac_meam, 0.0);
      setall2d(meam_data.repuls_meam, 0.0);
      setall3d(meam_data.Cmax_meam, 2.8);
      setall3d(meam_data.Cmin_meam, 2.0);
      setall2d(meam_data.ebound_meam, pow(2.8,2)/(4.0*(2.8-1.0)));
      setall2d(meam_data.delta_meam, 0.0);
      setall2d(meam_data.nn2_meam, 0);
      setall2d(meam_data.zbl_meam, 1);
      meam_data.gsmooth_factor = 99.0;
      meam_data.augt1 = 1;
      meam_data.ialloy = 0;
      meam_data.mix_ref_t = 0;
      meam_data.emb_lin_neg = 0;
      meam_data.bkgd_dyn = 0;
      meam_data.erose_form = 0;

      }


