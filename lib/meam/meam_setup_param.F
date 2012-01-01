c     
c     Declaration in pair_meam.h:
c     
c     void meam_setup_param(int *, double *, int *, int *, int *);
c     
c     Call in pair_meam.cpp
c     
c     meam_setup_param(&which,&value,&nindex,index,&errorflag);
c     
c     
c     
c     The "which" argument corresponds to the index of the "keyword" array
c     in pair_meam.cpp:
c     
c     0 = Ec_meam
c     1 = alpha_meam
c     2 = rho0_meam
c     3 = delta_meam
c     4 = lattce_meam
c     5 = attrac_meam
c     6 = repuls_meam
c     7 = nn2_meam
c     8 = Cmin_meam
c     9 = Cmax_meam
c     10 = rc_meam
c     11 = delr_meam
c     12 = augt1
c     13 = gsmooth_factor
c     14 = re_meam
c     15 = ialloy
c     16 = mixture_ref_t
c     17 = erose_form
c     18 = zbl_meam
c     19 = emb_lin_neg
c     20 = bkgd_dyn

      subroutine meam_setup_param(which, value, nindex, 
     $     index, errorflag)

      use meam_data
      implicit none

      integer which, nindex, index(3), errorflag
      real*8  value
      integer i1, i2

      errorflag = 0

c     0 = Ec_meam
      if (which.eq.0) then
        Ec_meam(index(1),index(2)) = value

c     1 = alpha_meam
      else if (which.eq.1) then
        alpha_meam(index(1),index(2)) = value
        
c     2 = rho0_meam
      else if (which.eq.2) then
        rho0_meam(index(1)) = value

c     3 = delta_meam
      else if (which.eq.3) then
        delta_meam(index(1),index(2)) = value

c     4 = lattce_meam
      else if (which.eq.4) then
        if (value.eq.0) then
          lattce_meam(index(1),index(2)) = "fcc"
        else if (value.eq.1) then
          lattce_meam(index(1),index(2)) = "bcc"
        else if (value.eq.2) then
          lattce_meam(index(1),index(2)) = "hcp"
        else if (value.eq.3) then
          lattce_meam(index(1),index(2)) = "dim"
        else if (value.eq.4) then
          lattce_meam(index(1),index(2)) = "dia"
        else if (value.eq.5) then
          lattce_meam(index(1),index(2)) = 'b1'
        else if (value.eq.6) then
          lattce_meam(index(1),index(2)) = 'c11'
        else if (value.eq.7) then
          lattce_meam(index(1),index(2)) = 'l12'
        else if (value.eq.8) then
          lattce_meam(index(1),index(2)) = 'b2'
        endif

c     5 = attrac_meam
      else if (which.eq.5) then
        attrac_meam(index(1),index(2)) = value

c     6 = repuls_meam
      else if (which.eq.6) then
        repuls_meam(index(1),index(2)) = value

c     7 = nn2_meam
      else if (which.eq.7) then
        i1 = min(index(1),index(2))
        i2 = max(index(1),index(2))
        nn2_meam(i1,i2) = value

c     8 = Cmin_meam
      else if (which.eq.8) then
        Cmin_meam(index(1),index(2),index(3)) = value

c     9 = Cmax_meam
      else if (which.eq.9) then
        Cmax_meam(index(1),index(2),index(3)) = value

c     10 = rc_meam
      else if (which.eq.10) then
        rc_meam = value

c     11 = delr_meam
      else if (which.eq.11) then
        delr_meam = value

c     12 = augt1
      else if (which.eq.12) then
        augt1 = value

c     13 = gsmooth
      else if (which.eq.13) then
        gsmooth_factor = value

c     14 = re_meam
      else if (which.eq.14) then
        re_meam(index(1),index(2)) = value

c     15 = ialloy
      else if (which.eq.15) then
        ialloy = value

c     16 = mixture_ref_t
      else if (which.eq.16) then
        mix_ref_t = value

c     17 = erose_form
      else if (which.eq.17) then
        erose_form = value

c     18 = zbl_meam
      else if (which.eq.18) then
        i1 = min(index(1),index(2))
        i2 = max(index(1),index(2))
        zbl_meam(i1,i2) = value

c     19 = emb_lin_neg
      else if (which.eq.19) then
        emb_lin_neg = value

c     20 = bkgd_dyn
      else if (which.eq.20) then
        bkgd_dyn = value

      else
        errorflag = 1
      endif

      return
      end
