c Declaration in pair_meam.h:
c
c void meam_setup_done(double *)
c
c Call from pair_meam.cpp:
c
c meam_setup_done(&cutmax)
c

      subroutine meam_setup_done(cutmax)
      use meam_data
      implicit none

      real*8 cutmax

      integer nv2, nv3, m, n, p

c     Force cutoff
      cutforce = rc_meam
      cutforcesq = cutforce*cutforce

c     Pass cutoff back to calling program
      cutmax = cutforce

c     Augment t1 term
      t1_meam(:) = t1_meam(:) + augt1 * 3.d0/5.d0 * t3_meam(:)

c     Compute off-diagonal alloy parameters
      call alloyparams

c indices and factors for Voight notation
      nv2 = 1
      nv3 = 1
      do m = 1,3
        do n = m,3
          vind2D(m,n) = nv2
          vind2D(n,m) = nv2
          nv2 = nv2+1
          do p = n,3
            vind3D(m,n,p) = nv3
            vind3D(m,p,n) = nv3
            vind3D(n,m,p) = nv3
            vind3D(n,p,m) = nv3
            vind3D(p,m,n) = nv3
            vind3D(p,n,m) = nv3
            nv3 = nv3+1
          enddo
        enddo
      enddo

      v2D(1) = 1
      v2D(2) = 2
      v2D(3) = 2
      v2D(4) = 1
      v2D(5) = 2
      v2D(6) = 1

      v3D(1) = 1
      v3D(2) = 3
      v3D(3) = 3
      v3D(4) = 3
      v3D(5) = 6
      v3D(6) = 3
      v3D(7) = 1
      v3D(8) = 3
      v3D(9) = 3
      v3D(10) = 1

      nv2 = 1
      do m = 1,neltypes
        do n = m,neltypes
          eltind(m,n) = nv2
          eltind(n,m) = nv2
          nv2 = nv2+1
        enddo
      enddo

c     Compute pair potentials and setup arrays for interpolation
      nr = 1000
      dr = 1.1*rc_meam/nr
      call compute_pair_meam

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fill off-diagonal alloy parameters
      subroutine alloyparams
      use meam_data
      implicit none
      integer i,j,k
      real*8 eb

c Loop over pairs
      do i = 1,neltypes
        do j = 1,neltypes
c Treat off-diagonal pairs
c If i>j, set all equal to i<j case (which has aready been set,
c here or in the input file)
          if (i.gt.j) then
            re_meam(i,j) = re_meam(j,i)
            Ec_meam(i,j) = Ec_meam(j,i)
            alpha_meam(i,j) = alpha_meam(j,i)
            lattce_meam(i,j) = lattce_meam(j,i)
c If i<j and term is unset, use default values (e.g. mean of i-i and j-j)
          else if (j.gt.i) then
            if (Ec_meam(i,j).eq.0.d0) then
              Ec_meam(i,j) = (Ec_meam(i,i)+Ec_meam(j,j))/2.d0 
     $             - delta_meam(i,j)
            endif
            if (alpha_meam(i,j).eq.0.d0) then
              alpha_meam(i,j) = (alpha_meam(i,i)+alpha_meam(j,j))/2.d0
            endif
            if (re_meam(i,j).eq.0.d0) then
              re_meam(i,j) = (re_meam(i,i)+re_meam(j,j))/2.d0
            endif
          endif
        enddo
      enddo

c Cmin(i,k,j) is symmetric in i-j, but not k.  For all triplets
c where i>j, set equal to the i<j element.  Likewise for Cmax.
        do i = 2,neltypes
          do j = 1,i-1
            do k = 1,neltypes
            Cmin_meam(i,j,k) = Cmin_meam(j,i,k)
            Cmax_meam(i,j,k) = Cmax_meam(j,i,k)
          enddo
        enddo
      enddo

c ebound gives the squared distance such that, for rik2 or rjk2>ebound,
c atom k definitely lies outside the screening function ellipse (so 
c there is no need to calculate its effects).  Here, compute it for all
c triplets (i,j,k) so that ebound(i,j) is the maximized over k
      do i = 1,neltypes
        do j = 1,neltypes
          do k = 1,neltypes
            eb = (Cmax_meam(i,j,k)*Cmax_meam(i,j,k))
     $           /(4.d0*(Cmax_meam(i,j,k)-1.d0))
            ebound_meam(i,j) = max(ebound_meam(i,j),eb)
          enddo
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c compute MEAM pair potential for each pair of element types
c

      subroutine compute_pair_meam
      use meam_data
      implicit none

      real*8 r, temp
      integer j,a,b,nv2
      real*8 astar,frac,phizbl
      integer n,nmax,Z1,Z2
      real*8 arat,scrn,scrn2
      real*8 phitmp

      real*8, external :: phi_meam
      real*8, external :: zbl
      real*8, external :: compute_phi
 

c allocate memory for array that defines the potential
      allocate(phir(nr,(neltypes*(neltypes+1))/2))

c allocate coeff memory

      allocate(phirar(nr,(neltypes*(neltypes+1))/2))
      allocate(phirar1(nr,(neltypes*(neltypes+1))/2))
      allocate(phirar2(nr,(neltypes*(neltypes+1))/2))
      allocate(phirar3(nr,(neltypes*(neltypes+1))/2))
      allocate(phirar4(nr,(neltypes*(neltypes+1))/2))
      allocate(phirar5(nr,(neltypes*(neltypes+1))/2))
      allocate(phirar6(nr,(neltypes*(neltypes+1))/2))

c HACK for debug: compute phi_meam for Ni3Si equilibrium spacing
      r = 2.47770216
      phitmp = phi_meam(r,1,2)

c loop over pairs of element types
      nv2 = 0
      do a = 1,neltypes
        do b = a,neltypes
          nv2 = nv2 + 1

c loop over r values and compute
          do j = 1,nr

            r = (j-1)*dr

            phir(j,nv2) = phi_meam(r,a,b)

c if using second-nearest neighbor, solve recursive problem
c (see Lee and Baskes, PRB 62(13):8564 eqn.(21))
            if (nn2_meam(a,b).eq.1) then
c              if (a.ne.b) then
c                write(6,*) 'Second-neighbor currently only valid '
c                write(6,*) 'if element types are the same.'
c                call error(' ')
c              endif
              call get_Zij(Z1,lattce_meam(a,b))
              call get_Zij2(Z2,arat,scrn,lattce_meam(a,b),
     $             Cmin_meam(a,a,b),Cmax_meam(a,a,b))

              if (lattce_meam(a,b).eq.'b1') then
                phir(j,nv2) = phir(j,nv2) -
     $               Z2*scrn/(2*Z1) * phi_meam(r*arat,a,a)
                call get_Zij2(Z2,arat,scrn2,lattce_meam(a,b),
     $               Cmin_meam(b,b,a),Cmax_meam(b,b,a))
                phir(j,nv2) = phir(j,nv2) -
     $               Z2*scrn2/(2*Z1) * phi_meam(r*arat,b,b)
              else
                nmax = 10
                do n = 1,nmax
                  phir(j,nv2) = phir(j,nv2) +
     $                 (-Z2*scrn/Z1)**n * phi_meam(r*arat**n,a,b)
                enddo
              endif

            endif

c if astar <= -3
c   potential is zbl potential
c else if -3 < astar < -1
c   potential is linear combination with zbl potential
c endif
            astar = alpha_meam(a,b) * (r/re_meam(a,b) - 1.d0)
            if (astar.le.-3.d0) then
              phir(j,nv2) = zbl(r,ielt_meam(a),ielt_meam(b))
            else if (astar.gt.-3.d0.and.astar.lt.-1.d0) then
              call fcut(1-(astar+1.d0)/(-3.d0+1.d0),frac)
              phizbl = zbl(r,ielt_meam(a),ielt_meam(b))
              phir(j,nv2) = frac*phir(j,nv2) + (1-frac)*phizbl
            endif

          enddo

c call interpolation
          call interpolate_meam(nv2)
          
        enddo
      enddo

      return
      end


c----------------------------------------------------------------------c      
c Compute MEAM pair potential for distance r, element types a and b
c
      real*8 recursive function phi_meam(r,a,b)result(phi_m)
      use meam_data
      implicit none

      
      integer a,b
      real*8 r
      real*8 a1,a2,a12
      real*8 t11av,t21av,t31av,t12av,t22av,t32av
      real*8 G1,G2,s1(3),s2(3),s12(3),rho0_1,rho0_2
      real*8 Gam1,Gam2,Z1,Z2
      real*8 rhobar1,rhobar2,F1,F2
      real*8 rhoa01,rhoa11,rhoa21,rhoa31
      real*8 rhoa02,rhoa12,rhoa22,rhoa32
      real*8 rho01,rho11,rho21,rho31
      real*8 rho02,rho12,rho22,rho32
      real*8 scalfac,phiaa,phibb
      real*8 Eu
      integer Z12, errorflag
      character*3 latta,lattb

      real*8, external :: erose

c Equation numbers below refer to:
c   I. Huang et.al., Modelling simul. Mater. Sci. Eng. 3:615

c get number of neighbors in the reference structure
c   Nref(i,j) = # of i's neighbors of type j
      call get_Zij(Z12,lattce_meam(a,b))

      call get_densref(r,a,b,rho01,rho11,rho21,rho31,
     $     rho02,rho12,rho22,rho32)

c calculate average weighting factors for the reference structure
      if (lattce_meam(a,b).eq.'c11') then
        scalfac = 1.0/(rho01+rho02) 
        t11av = scalfac*(t1_meam(a)*rho01 + t1_meam(b)*rho02)
        t12av = t11av 
        t21av = scalfac*(t2_meam(a)*rho01 + t2_meam(b)*rho02)
        t22av = t21av 
        t31av = scalfac*(t3_meam(a)*rho01 + t3_meam(b)*rho02)
        t32av = t31av 
      else      
c average weighting factors for the reference structure, eqn. I.8
         call get_tavref(t11av,t21av,t31av,t12av,t22av,t32av,
     $       t1_meam(a),t2_meam(a),t3_meam(a),
     $       t1_meam(b),t2_meam(b),t3_meam(b),
     $       r,a,b,lattce_meam(a,b))
      endif

c for c11b structure, calculate background electron densities
      if (lattce_meam(a,b).eq.'c11') then
         latta = lattce_meam(a,a)
         if (latta.eq.'dia') then
            rhobar1 = ((Z12/2)*(rho02+rho01))**2 +
     $                t11av*(rho12-rho11)**2 +
     $                t21av/6.0*(rho22+rho21)**2 +
     $                121.0/40.*t31av*(rho32-rho31)**2
            rhobar1 = sqrt(rhobar1)
            rhobar2 = (Z12*rho01)**2 + 2.0/3.0*t21av*rho21**2
            rhobar2 = sqrt(rhobar2)
         else
            rhobar2 = ((Z12/2)*(rho01+rho02))**2 +
     $                t12av*(rho11-rho12)**2 +
     $                t22av/6.0*(rho21+rho22)**2 +
     $                121.0/40.*t32av*(rho31-rho32)**2
            rhobar2 = sqrt(rhobar2)
            rhobar1 = (Z12*rho02)**2 + 2.0/3.0*t22av*rho22**2
            rhobar1 = sqrt(rhobar1)
         endif
      else 
c for other structures, use formalism developed in Huang's paper
c
c composition-dependent scaling, equation I.7
c -- note: The shape factors for the individual elements are used, since
c          this function should cancel the F(rho) terms when the atoms are
c          in the reference structure.
c -- GJW (6/12/09): I suspect there may be a problem here... since we should be
c                   using the pure element reference structure, we should also
c                   be using the element "t" values (not the averages compute
c                   above).  It doesn't matter for reference structures with
c                   s(1)=s(2)=s(3)=0, but might cause diffs otherwise.  For now
c                   what's here is consistent with what's done in meam_dens_final.
         Z1 = Z_meam(a)
         Z2 = Z_meam(b)
         if (ibar_meam(a).le.0) then
            G1 = 1.d0
         else
            call get_shpfcn(s1,lattce_meam(a,a))
            Gam1 = (s1(1)*t11av+s1(2)*t21av+s1(3)*t31av)/(Z1*Z1)
            call G_gam(Gam1,ibar_meam(a),gsmooth_factor,G1,errorflag)
         endif
         if (ibar_meam(b).le.0) then
            G2 = 1.d0
         else
            call get_shpfcn(s2,lattce_meam(b,b))
            Gam2 = (s2(1)*t12av+s2(2)*t22av+s2(3)*t32av)/(Z2*Z2)
            call G_gam(Gam2,ibar_meam(b),gsmooth_factor,G2,errorflag)
         endif
         rho0_1 = rho0_meam(a)*Z1*G1
         rho0_2 = rho0_meam(b)*Z2*G2
         
c get number of neighbors in the reference structure
c   Nref(i,j) = # of i's neighbors of type j
c		call get_Zij(Z12,lattce_meam(a,b))
c	
c		call get_densref(r,a,b,rho01,rho11,rho21,rho31,
c		$     rho02,rho12,rho22,rho32)
      
c compute total background density, eqn I.6
         Gam1 = (t11av*rho11+t21av*rho21+t31av*rho31)
         Gam1 = Gam1/(rho01*rho01)
         Gam2 = (t12av*rho12+t22av*rho22+t32av*rho32)
         Gam2 = Gam2/(rho02*rho02)
         call G_gam(Gam1,ibar_meam(a),gsmooth_factor,G1,errorflag)
         call G_gam(Gam2,ibar_meam(b),gsmooth_factor,G2,errorflag)
         rhobar1 = rho01/rho0_1*G1
         rhobar2 = rho02/rho0_2*G2
      endif
c compute embedding functions, eqn I.5
      if (rhobar1.eq.0.d0) then
        F1 = 0.d0
      else
        F1 = A_meam(a)*Ec_meam(a,a)*rhobar1*log(rhobar1)
      endif
      if (rhobar2.eq.0.d0) then
        F2 = 0.d0
      else
        F2 = A_meam(b)*Ec_meam(b,b)*rhobar2*log(rhobar2)
      endif

c compute Rose function, I.16
      Eu = erose(r,re_meam(a,b),alpha_meam(a,b),
     $     Ec_meam(a,b),repuls_meam(a,b),attrac_meam(a,b))

c calculate the pair energy
      if (lattce_meam(a,b).eq.'c11') then
        latta = lattce_meam(a,a)
        if (latta.eq.'dia') then
          phiaa = phi_meam(r,a,a)
          phi_m = (3*Eu - F2 - 2*F1 - 5*phiaa)/Z12
        else
          phibb = phi_meam(r,b,b)
          phi_m = (3*Eu - F1 - 2*F2 - 5*phibb)/Z12
        endif
c     phi_m = 0.d0
      else if (lattce_meam(a,b).eq.'l12') then
        phiaa = phi_meam(r,a,a)
        phi_m = Eu/3. - F1/4. - F2/12. - phiaa
      else      
c     
c potential is computed from Rose function and embedding energy
         phi_m = (2*Eu - F1 - F2)/Z12
c
      endif

c if r = 0, just return 0
      if (r.eq.0.d0) then
        phi_m = 0.d0
      endif
      
      return
      end

c----------------------------------------------------------------------c      
c Shape factors for various configurations
      subroutine get_shpfcn(s,latt)
      implicit none
      real*8 s(3)
      character*3 latt
      if (latt.eq.'fcc'.or.latt.eq.'bcc'.or.latt.eq.'b1') then
        s(1) = 0.d0
        s(2) = 0.d0
        s(3) = 0.d0
      else if (latt.eq.'hcp') then
        s(1) = 0.d0
        s(2) = 0.d0
        s(3) = 1.d0/3.d0
      else if (latt.eq.'dia') then
        s(1) = 0.d0
        s(2) = 0.d0
        s(3) = 32.d0/9.d0
      else if (latt.eq.'dim') then
        s(1) = 1.d0
        s(2) = 2.d0/3.d0
c        s(3) = 1.d0
        s(3) = 0.4d0
      else
        s(1) = 0.0
c        call error('Lattice not defined in get_shpfcn.')
      endif
      return
      end
c------------------------------------------------------------------------------c      
c Average weighting factors for the reference structure
      subroutine get_tavref(t11av,t21av,t31av,t12av,t22av,t32av,
     $     t11,t21,t31,t12,t22,t32,
     $     r,a,b,latt)
      use meam_data
      implicit none
      real*8 t11av,t21av,t31av,t12av,t22av,t32av
      real*8 t11,t21,t31,t12,t22,t32,r
      integer a,b
      character*3 latt
      real*8 rhoa01,rhoa02,a1,a2,rho01,rho02
      if (latt.eq.'fcc'.or.latt.eq.'bcc'.or.latt.eq.'dia'
     $     .or.latt.eq.'hcp'.or.latt.eq.'b1'
     $     .or.latt.eq.'dim') then
c all neighbors are of the opposite type
        t11av = t12
        t21av = t22
        t31av = t32
        t12av = t11
        t22av = t21
        t32av = t31
      else
        a1  = r/re_meam(a,a) - 1.d0
        a2  = r/re_meam(b,b) - 1.d0
        rhoa01 = rho0_meam(a)*exp(-beta0_meam(a)*a1)
        rhoa02 = rho0_meam(b)*exp(-beta0_meam(b)*a2)
        if (latt.eq.'l12') then
          rho01 = 8*rhoa01 + 4*rhoa02
          t11av = (8*t11*rhoa01 + 4*t12*rhoa02)/rho01
          t12av = t11
          t21av = (8*t21*rhoa01 + 4*t22*rhoa02)/rho01
          t22av = t21
          t31av = (8*t31*rhoa01 + 4*t32*rhoa02)/rho01
          t32av = t31
        else
c         call error('Lattice not defined in get_tavref.')
        endif
      endif
      return
      end
c------------------------------------------------------------------------------c      
c Number of neighbors for the reference structure
      subroutine get_Zij(Zij,latt)
      implicit none
      integer Zij
      character*3 latt
      if (latt.eq.'fcc') then
        Zij = 12
      else if (latt.eq.'bcc') then
        Zij = 8
      else if (latt.eq.'hcp') then
        Zij = 12
      else if (latt.eq.'b1') then
        Zij = 6
      else if (latt.eq.'dia') then
        Zij = 4
      else if (latt.eq.'dim') then
        Zij = 1
      else if (latt.eq.'c11') then
        Zij = 10
      else if (latt.eq.'l12') then
        Zij = 12
      else
c        call error('Lattice not defined in get_Zij.')
      endif
      return
      end

c------------------------------------------------------------------------------c      
c Zij2 = number of second neighbors, a = distance ratio R1/R2, and S = second 
c neighbor screening function for lattice type "latt"

      subroutine get_Zij2(Zij2,a,S,latt,cmin,cmax)
      implicit none
      integer Zij2
      real*8 a,S,cmin,cmax
      character*3 latt
      real*8 rratio,C,x,sijk
      integer numscr

      if (latt.eq.'bcc') then
        Zij2 = 6
        a = 2.d0/sqrt(3.d0)
        numscr = 4
      else if (latt.eq.'fcc') then
        Zij2 = 6
        a = sqrt(2.d0)
        numscr = 4
      else if (latt.eq.'dia') then
        Zij2 = 0
        a = sqrt(8.d0/3.d0)
        numscr = 4
        if (cmin.lt.0.500001) then
c          call error('can not do 2NN MEAM for dia')
        endif
      else if (latt.eq.'hcp') then
        Zij2 = 6
        a = sqrt(2.d0)
        numscr = 4
      else if (latt.eq.'b1') then
        Zij2 = 12
        a = sqrt(2.d0)
        numscr = 2
      else
c        call error('Lattice not defined in get_Zij2.')
      endif

c Compute screening for each first neighbor      
      C = 4.d0/(a*a) - 1.d0
      x = (C-cmin)/(cmax-cmin)
      call fcut(x,sijk)
c There are numscr first neighbors screening the second neighbors
      S = sijk**numscr

      return
      end

c------------------------------------------------------------------------------c      
c Calculate density functions, assuming reference configuration
      subroutine get_densref(r,a,b,rho01,rho11,rho21,rho31,
     $     rho02,rho12,rho22,rho32)
      use meam_data
      implicit none
      real*8 r,rho01,rho11,rho21,rho31,rho02,rho12,rho22,rho32
      real*8 a1,a2
      real*8 rhoa01,rhoa11,rhoa21,rhoa31,rhoa02,rhoa12,rhoa22,rhoa32
      real*8 C,Cmin,Cmax,fc,s(3)
      character*3 lat
      integer a,b
      integer Zij1nn,Zij2nn
      real*8 rhoa01nn,rhoa02nn
      real*8 arat,scrn,denom

      a1  = r/re_meam(a,a) - 1.d0
      a2  = r/re_meam(b,b) - 1.d0

      rhoa01 = rho0_meam(a)*exp(-beta0_meam(a)*a1)
      rhoa11 = rho0_meam(a)*exp(-beta1_meam(a)*a1)
      rhoa21 = rho0_meam(a)*exp(-beta2_meam(a)*a1)
      rhoa31 = rho0_meam(a)*exp(-beta3_meam(a)*a1)
      rhoa02 = rho0_meam(b)*exp(-beta0_meam(b)*a2)
      rhoa12 = rho0_meam(b)*exp(-beta1_meam(b)*a2)
      rhoa22 = rho0_meam(b)*exp(-beta2_meam(b)*a2)
      rhoa32 = rho0_meam(b)*exp(-beta3_meam(b)*a2)

      lat = lattce_meam(a,b)
      
      rho11 = 0.d0
      rho21 = 0.d0
      rho31 = 0.d0
      rho12 = 0.d0
      rho22 = 0.d0
      rho32 = 0.d0

      call get_Zij(Zij1nn,lat)

      if (lat.eq.'fcc') then
        rho01 = 12.d0*rhoa02
        rho02 = 12.d0*rhoa01
      else if (lat.eq.'bcc') then
        rho01 = 8.d0*rhoa02
        rho02 = 8.d0*rhoa01
      else if (lat.eq.'b1') then
        rho01 = 6*rhoa02
        rho02 = 6*rhoa01
      else if (lat.eq.'dia') then
        rho01 = 4*rhoa02 
        rho02 = 4*rhoa01 
        rho31 = 32.d0/9.d0*rhoa32*rhoa32
        rho32 = 32.d0/9.d0*rhoa31*rhoa31
      else if (lat.eq.'hcp') then
        rho01 = 12*rhoa02 
        rho02 = 12*rhoa01 
        rho31 = 1.d0/3.d0*rhoa32*rhoa32
        rho32 = 1.d0/3.d0*rhoa31*rhoa31
      else if (lat.eq.'dim') then
        call get_shpfcn(s,'dim')
        rho01 = rhoa02
        rho02 = rhoa01
        rho11 = s(1)*rhoa12*rhoa12
        rho12 = s(1)*rhoa11*rhoa11
        rho21 = s(2)*rhoa22*rhoa22
        rho22 = s(2)*rhoa21*rhoa21
        rho31 = s(3)*rhoa32*rhoa32
        rho32 = s(3)*rhoa31*rhoa31
      else if (lat.eq.'c11') then
        rho01 = rhoa01
        rho02 = rhoa02
        rho11 = rhoa11
        rho12 = rhoa12 
        rho21 = rhoa21
        rho22 = rhoa22
        rho31 = rhoa31
        rho32 = rhoa32
      else if (lat.eq.'l12') then
        rho01 = 8*rhoa01 + 4*rhoa02
        rho02 = 12*rhoa01
        if (ialloy.eq.0) then
          rho21 = 8./3.*(rhoa21-rhoa22)*(rhoa21-rhoa22)
        else
          rho21 = 8./3.*(rhoa21*t2_meam(a)-rhoa22*t2_meam(b))**2
          denom = 8*rhoa01*t2_meam(a)**2 + 4*rhoa02*t2_meam(b)**2
          if (denom.gt.0.) then
            rho21 = rho21/denom * rho01
          endif
        endif
      else
c        call error('Lattice not defined in get_densref.')
      endif

      if (nn2_meam(a,b).eq.1) then
c Assume that second neighbor is of same type, first neighbor
c may be of different type

c        if (lat.ne.'bcc') then
c          call error('Second-neighbor not defined for this lattice')
c        endif

        call get_Zij2(Zij2nn,arat,scrn,lat,
     $       Cmin_meam(a,a,b),Cmax_meam(a,a,b))
        
        a1 = arat*r/re_meam(a,a) - 1.d0
        a2 = arat*r/re_meam(b,b) - 1.d0
        
        rhoa01nn = rho0_meam(a)*exp(-beta0_meam(a)*a1)
        rhoa02nn = rho0_meam(b)*exp(-beta0_meam(b)*a2)

        rho01 = rho01 + Zij2nn*scrn*rhoa01nn

c     Assume Zij2nn and arat don't depend on order, but scrn might
        call get_Zij2(Zij2nn,arat,scrn,lat,
     $       Cmin_meam(b,b,a),Cmax_meam(b,b,a))
        rho02 = rho02 + Zij2nn*scrn*rhoa02nn

      endif

      return
      end

c---------------------------------------------------------------------       
c Compute ZBL potential  
c  
      real*8 function zbl(r,z1,z2)  
      implicit none  
      integer i,z1,z2 
      real*8 r,c,d,a,azero,cc,x  
      dimension c(4),d(4) 
      data c /0.028171,0.28022,0.50986,0.18175/ 
      data d /0.20162,0.40290,0.94229,3.1998/ 
      data azero /0.4685/  
      data cc /14.3997/ 
c azero = (9pi^2/128)^1/3 (0.529) Angstroms 
      a = azero/(z1**0.23+z2**0.23) 
      zbl = 0.0  
      x = r/a  
      do i=1,4  
        zbl = zbl + c(i)*exp(-d(i)*x) 
      enddo
      if (r.gt.0.d0) zbl = zbl*z1*z2/r*cc
      return  
      end  
 
c---------------------------------------------------------------------
c Compute Rose energy function, I.16
c
      real*8 function erose(r,re,alpha,Ec,repuls,attrac)
      implicit none
      real*8 r,re,alpha,Ec,repuls,attrac,astar,a3

      erose = 0.d0
      
      if (r.gt.0.d0) then
        astar = alpha * (r/re - 1.d0)
        a3 = 0.d0
        if (astar.ge.0) then
          a3 = attrac
        else if (astar.lt.0) then
          a3 = repuls
        endif
        erose = -Ec * (1 + astar + a3*(astar**3)/(r/re)) * exp(-astar)
c      erose = -Ec*(1+astar+(-attrac+repuls/r)*(astar**3))*exp(-astar)
      endif

      return
      end
      
c -----------------------------------------------------------------------

      subroutine interpolate_meam(ind)
      use meam_data
      implicit none

      integer j,ind
      real*8 drar

c map to coefficient space

      nrar = nr
      drar = dr
      rdrar = 1.0D0/drar

c phir interp
      do j = 1,nrar
        phirar(j,ind) = phir(j,ind)
      enddo

      phirar1(1,ind) = phirar(2,ind)-phirar(1,ind)
      phirar1(2,ind) = 0.5D0*(phirar(3,ind)-phirar(1,ind))
      phirar1(nrar-1,ind) = 0.5D0*(phirar(nrar,ind)
     $     -phirar(nrar-2,ind))
      phirar1(nrar,ind) = 0.0D0
      do j = 3,nrar-2
        phirar1(j,ind) = ((phirar(j-2,ind)-phirar(j+2,ind)) +
     $       8.0D0*(phirar(j+1,ind)-phirar(j-1,ind)))/12.
      enddo
      
      do j = 1,nrar-1
        phirar2(j,ind) = 3.0D0*(phirar(j+1,ind)-phirar(j,ind)) -
     $       2.0D0*phirar1(j,ind) - phirar1(j+1,ind)
        phirar3(j,ind) = phirar1(j,ind) + phirar1(j+1,ind) -
     $       2.0D0*(phirar(j+1,ind)-phirar(j,ind))
      enddo
      phirar2(nrar,ind) = 0.0D0
      phirar3(nrar,ind) = 0.0D0
      
      do j = 1,nrar
        phirar4(j,ind) = phirar1(j,ind)/drar
        phirar5(j,ind) = 2.0D0*phirar2(j,ind)/drar
        phirar6(j,ind) = 3.0D0*phirar3(j,ind)/drar
      enddo

      end

c---------------------------------------------------------------------
c Compute Rose energy function, I.16
c
      real*8 function compute_phi(rij, elti, eltj)
      use meam_data
      implicit none

      real*8  rij, pp
      integer elti, eltj, ind, kk

      ind = eltind(elti, eltj)
      pp = rij*rdrar + 1.0D0
      kk = pp
      kk = min(kk,nrar-1)
      pp = pp - kk
      pp = min(pp,1.0D0)
      compute_phi = ((phirar3(kk,ind)*pp + phirar2(kk,ind))*pp 
     $     + phirar1(kk,ind))*pp + phirar(kk,ind)

      return
      end
