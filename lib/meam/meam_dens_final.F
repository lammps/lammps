c Extern "C" declaration has the form:
c
c  void meam_dens_final_(int *, int *, int *, int *, int *, double *, double *,
c                int *, int *, int *,
c		 double *, double *, double *, double *, double *, double *,
c		 double *, double *, double *, double *, double *, double *, 
c		 double *, double *, double *, double *, double *, int *);
c
c Call from pair_meam.cpp has the form:
c
c  meam_dens_final_(&nlocal,&nmax,&eflag_either,&eflag_global,&eflag_atom,
c            &eng_vdwl,eatom,ntype,type,fmap,
c	     &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],
c	     &arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],gamma,dgamma1,
c	     dgamma2,dgamma3,rho,rho0,rho1,rho2,rho3,frhop,&errorflag);
c

      subroutine meam_dens_final(nlocal, nmax,
     $     eflag_either, eflag_global, eflag_atom, eng_vdwl, eatom,
     $     ntype, type, fmap,
     $     Arho1, Arho2, Arho2b, Arho3, Arho3b, t_ave, tsq_ave,
     $     Gamma, dGamma1, dGamma2, dGamma3,
     $     rho, rho0, rho1, rho2, rho3, fp, errorflag)

      use meam_data
      implicit none

      integer nlocal, nmax, eflag_either, eflag_global, eflag_atom
      integer ntype, type, fmap
      real*8 eng_vdwl, eatom, Arho1, Arho2
      real*8 Arho2b, Arho3, Arho3b
      real*8 t_ave, tsq_ave
      real*8 Gamma, dGamma1, dGamma2, dGamma3
      real*8 rho, rho0, rho1, rho2, rho3
      real*8 fp
      integer errorflag

      dimension eatom(nmax)
      dimension type(nmax), fmap(ntype)
      dimension Arho1(3,nmax), Arho2(6,nmax), Arho2b(nmax)
      dimension Arho3(10,nmax), Arho3b(3,nmax), t_ave(3,nmax)
      dimension tsq_ave(3,nmax)
      dimension Gamma(nmax), dGamma1(nmax), dGamma2(nmax)
      dimension dGamma3(nmax), rho(nmax), rho0(nmax)
      dimension rho1(nmax), rho2(nmax), rho3(nmax)
      dimension fp(nmax)

      integer i, elti
      integer m
      real*8  rhob, G, dG, Gbar, dGbar, gam, shp(3), shpi(3), Z
      real*8  B, denom, rho_bkgd

c     Complete the calculation of density

      do i = 1,nlocal

        elti = fmap(type(i))
        if (elti.gt.0) then
          rho1(i) = 0.d0
          rho2(i) = -1.d0/3.d0*Arho2b(i)*Arho2b(i)
          rho3(i) = 0.d0
          do m = 1,3
            rho1(i) = rho1(i) + Arho1(m,i)*Arho1(m,i)
            rho3(i) = rho3(i) - 3.d0/5.d0*Arho3b(m,i)*Arho3b(m,i)
          enddo
          do m = 1,6
            rho2(i) = rho2(i) + v2D(m)*Arho2(m,i)*Arho2(m,i)
          enddo
          do m = 1,10
            rho3(i) = rho3(i) + v3D(m)*Arho3(m,i)*Arho3(m,i)
          enddo
          
          if( rho0(i) .gt. 0.0 ) then
            if (ialloy.eq.1) then
              t_ave(1,i) = t_ave(1,i)/tsq_ave(1,i)
              t_ave(2,i) = t_ave(2,i)/tsq_ave(2,i)
              t_ave(3,i) = t_ave(3,i)/tsq_ave(3,i)
            else if (ialloy.eq.2) then
              t_ave(1,i) = t1_meam(elti)
              t_ave(2,i) = t2_meam(elti)
              t_ave(3,i) = t3_meam(elti)
            else
              t_ave(1,i) = t_ave(1,i)/rho0(i)
              t_ave(2,i) = t_ave(2,i)/rho0(i)
              t_ave(3,i) = t_ave(3,i)/rho0(i)
            endif
          endif
          
          Gamma(i) = t_ave(1,i)*rho1(i) 
     $         + t_ave(2,i)*rho2(i) + t_ave(3,i)*rho3(i)
          
          if( rho0(i) .gt. 0.0 ) then
            Gamma(i) = Gamma(i)/(rho0(i)*rho0(i))
          end if

          Z = Z_meam(elti)

          call G_gam(Gamma(i),ibar_meam(elti),
     $         gsmooth_factor,G,errorflag)
          if (errorflag.ne.0) return
          if (ibar_meam(elti).le.0) then
            Gbar = 1.d0
          else
            call get_shpfcn(shp,lattce_meam(elti,elti))
            if (mix_ref_t.eq.1) then
              gam = (t_ave(1,i)*shpi(1)+t_ave(2,i)*shpi(2)
     $             +t_ave(3,i)*shpi(3))/(Z*Z)
            else
              gam = (t1_meam(elti)*shp(1)+t2_meam(elti)*shp(2)
     $             +t3_meam(elti)*shp(3))/(Z*Z)
            endif
            call G_gam(gam,ibar_meam(elti),gsmooth_factor,
     $           Gbar,errorflag)
          endif
          rho(i) = rho0(i) * G

          if (mix_ref_t.eq.1) then
            if (ibar_meam(elti).le.0) then
              Gbar = 1.d0
              dGbar = 0.d0
            else
              call get_shpfcn(shpi,lattce_meam(elti,elti))
              gam = (t_ave(1,i)*shpi(1)+t_ave(2,i)*shpi(2)
     $             +t_ave(3,i)*shpi(3))/(Z*Z)
              call dG_gam(gam,ibar_meam(elti),gsmooth_factor,
     $             Gbar,dGbar)
            endif
            rho_bkgd = rho0_meam(elti)*Z*Gbar
          else
            if (bkgd_dyn.eq.1) then
              rho_bkgd = rho0_meam(elti)*Z
            else
              rho_bkgd = rho_ref_meam(elti)
            endif
          endif
          rhob = rho(i)/rho_bkgd
          denom = 1.d0/rho_bkgd

          call dG_gam(Gamma(i),ibar_meam(elti),gsmooth_factor,G,dG)
          
          dGamma1(i) = (G - 2*dG*Gamma(i))*denom
          
          if( rho0(i) .ne. 0.d0 ) then
            dGamma2(i) = (dG/rho0(i))*denom
          else
            dGamma2(i) = 0.d0
          end if

c     dGamma3 is nonzero only if we are using the "mixed" rule for
c     computing t in the reference system (which is not correct, but
c     included for backward compatibility
          if (mix_ref_t.eq.1) then
            dGamma3(i) = rho0(i)*G*dGbar/(Gbar*Z*Z)*denom
          else
            dGamma3(i) = 0.0
          endif

          B = A_meam(elti)*Ec_meam(elti,elti)
          
          if( rhob .ne. 0.d0 ) then
            if (emb_lin_neg.eq.1 .and. rhob.le.0) then
              fp(i) = -B
            else
              fp(i) = B*(log(rhob)+1.d0)
            endif
            if (eflag_either.ne.0) then
              if (eflag_global.ne.0) then
                if (emb_lin_neg.eq.1 .and. rhob.le.0) then
                  eng_vdwl = eng_vdwl - B*rhob
                else
                  eng_vdwl = eng_vdwl + B*rhob*log(rhob)
                endif
              endif
              if (eflag_atom.ne.0) then
                if (emb_lin_neg.eq.1 .and. rhob.le.0) then
                  eatom(i) = eatom(i) - B*rhob
                else
                  eatom(i) = eatom(i) + B*rhob*log(rhob)
                endif
              endif
            endif
          else
            if (emb_lin_neg.eq.1) then
              fp(i) = -B
            else
              fp(i) = B
            endif
          endif
        endif
      enddo
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine G_gam(Gamma,ibar,gsmooth_factor,G,errorflag)
c     Compute G(Gamma) based on selection flag ibar:
c   0 => G = sqrt(1+Gamma)
c   1 => G = exp(Gamma/2)
c   2 => not implemented 
c   3 => G = 2/(1+exp(-Gamma))
c   4 => G = sqrt(1+Gamma)
c  -5 => G = +-sqrt(abs(1+Gamma))
      implicit none
      real*8 Gamma,G
      real*8 gsmooth_factor, gsmooth_switchpoint
      integer ibar, errorflag
      if (ibar.eq.0.or.ibar.eq.4) then
        gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor+1)
        if (Gamma.lt.gsmooth_switchpoint) then
c         e.g. gsmooth_factor is 99, then:
c         gsmooth_switchpoint = -0.99
c         G = 0.01*(-0.99/Gamma)**99
          G = 1/(gsmooth_factor+1)
     $         *(gsmooth_switchpoint/Gamma)**gsmooth_factor
          G = sqrt(G)
        else
          G = sqrt(1.d0+Gamma)
        endif
      else if (ibar.eq.1) then
        G = exp(Gamma/2.d0)
      else if (ibar.eq.3) then
        G = 2.d0/(1.d0+exp(-Gamma))
      else if (ibar.eq.-5) then
        if ((1.d0+Gamma).ge.0) then
           G = sqrt(1.d0+Gamma)
        else
           G = -sqrt(-1.d0-Gamma)
        endif
      else
         errorflag = 1
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dG_gam(Gamma,ibar,gsmooth_factor,G,dG)
c Compute G(Gamma) and dG(gamma) based on selection flag ibar:
c   0 => G = sqrt(1+Gamma)
c   1 => G = exp(Gamma/2)
c   2 => not implemented 
c   3 => G = 2/(1+exp(-Gamma))
c   4 => G = sqrt(1+Gamma)
c  -5 => G = +-sqrt(abs(1+Gamma))
      real*8 Gamma,G,dG
      real*8 gsmooth_factor, gsmooth_switchpoint
      integer ibar
      if (ibar.eq.0.or.ibar.eq.4) then
        gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor+1)
        if (Gamma.lt.gsmooth_switchpoint) then
c         e.g. gsmooth_factor is 99, then:
c         gsmooth_switchpoint = -0.99
c         G = 0.01*(-0.99/Gamma)**99
          G = 1/(gsmooth_factor+1)
     $         *(gsmooth_switchpoint/Gamma)**gsmooth_factor
          G = sqrt(G)
          dG = -gsmooth_factor*G/(2.0*Gamma)
        else
          G = sqrt(1.d0+Gamma)
          dG = 1.d0/(2.d0*G)
        endif
      else if (ibar.eq.1) then
        G = exp(Gamma/2.d0)
        dG = G/2.d0
      else if (ibar.eq.3) then
        G = 2.d0/(1.d0+exp(-Gamma))
        dG = G*(2.d0-G)/2
      else if (ibar.eq.-5) then
        if ((1.d0+Gamma).ge.0) then
           G = sqrt(1.d0+Gamma)
           dG = 1.d0/(2.d0*G)
        else
           G = -sqrt(-1.d0-Gamma)
           dG = -1.d0/(2.d0*G)
        endif
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
