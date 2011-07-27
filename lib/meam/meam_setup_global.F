c
c declaration in pair_meam.h:
c
c  void meam_setup_global(int *, int *, double *, int *, double *, double *,
c			 double *, double *, double *, double *, double *,
c			 double *, double *, double *, double *, double *,
c			 double *, double *, int *);
c
c call in pair_meam.cpp:
c
c  meam_setup_global(&nelements,lat,z,ielement,atwt,alpha,b0,b1,b2,b3,
c		    alat,esub,asub,t0,t1,t2,t3,rozero,ibar);
c
c

      subroutine meam_setup_global(nelt, lat, z, ielement, atwt, alpha,
     $     b0, b1, b2, b3, alat, esub, asub, 
     $     t0, t1, t2, t3, rozero, ibar)

      use meam_data
      implicit none

      integer nelt, lat, ielement, ibar
      real*8  z, atwt, alpha, b0, b1, b2, b3
      real*8  alat, esub, asub, t0, t1, t2, t3
      real*8  rozero

      dimension lat(nelt), ielement(nelt), ibar(nelt)
      dimension z(nelt), atwt(nelt), alpha(nelt)
      dimension b0(nelt), b1(nelt), b2(nelt), b3(nelt)
      dimension alat(nelt), esub(nelt), asub(nelt)
      dimension t0(nelt), t1(nelt), t2(nelt), t3(nelt), rozero(nelt)

      integer i
      real*8 tmplat(maxelt)

      neltypes = nelt

      do i = 1,nelt

         if (lat(i).eq.0) then
            lattce_meam(i,i) = 'fcc'
         else if (lat(i).eq.1) then
            lattce_meam(i,i) = 'bcc'
         else if (lat(i).eq.2) then
            lattce_meam(i,i) = 'hcp'
         else if (lat(i).eq.3) then
            lattce_meam(i,i) = 'dim'
         else if (lat(i).eq.4) then
            lattce_meam(i,i) = 'dia'
         else
c           unknown
         endif

         Z_meam(i) = z(i)
         ielt_meam(i) = ielement(i)
         alpha_meam(i,i) = alpha(i)
         beta0_meam(i) = b0(i)
         beta1_meam(i) = b1(i)
         beta2_meam(i) = b2(i)
         beta3_meam(i) = b3(i)
         tmplat(i) = alat(i)
         Ec_meam(i,i) = esub(i)
         A_meam(i) = asub(i)
         t0_meam(i) = t0(i)
         t1_meam(i) = t1(i)
         t2_meam(i) = t2(i)
         t3_meam(i) = t3(i)
         rho0_meam(i) = rozero(i)
         ibar_meam(i) = ibar(i)

         if (lattce_meam(i,i).eq.'fcc') then
            re_meam(i,i) = tmplat(i)/sqrt(2.d0)
         elseif (lattce_meam(i,i).eq.'bcc') then
            re_meam(i,i) = tmplat(i)*sqrt(3.d0)/2.d0
         elseif (lattce_meam(i,i).eq.'hcp') then
            re_meam(i,i) = tmplat(i)
         elseif (lattce_meam(i,i).eq.'dim') then
            re_meam(i,i) = tmplat(i)
         elseif (lattce_meam(i,i).eq.'dia') then
            re_meam(i,i) = tmplat(i)*sqrt(3.d0)/4.d0
         else
c           error
         endif

      enddo


c Set some defaults
      rc_meam = 4.0
      delr_meam = 0.1
      attrac_meam(:,:) = 0.0
      repuls_meam(:,:) = 0.0
      Cmax_meam(:,:,:) = 2.8
      Cmin_meam(:,:,:) = 2.0
      ebound_meam(:,:) = (2.8d0**2)/(4.d0*(2.8d0-1.d0))
      delta_meam(:,:) = 0.0
      nn2_meam(:,:) = 0
      zbl_meam(:,:) = 1
      gsmooth_factor = 99.0
      augt1 = 1
      ialloy = 0
      mix_ref_t = 0
      emb_lin_neg = 0
      bkgd_dyn = 0
      erose_form = 0

      return
      end


