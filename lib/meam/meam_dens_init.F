c     Extern "C" declaration has the form:
c     
c  void meam_dens_init_(int *, int *, int *, double *, int *, int *, int *, double *,
c		 int *, int *, int *, int *,
c		 double *, double *, double *, double *, double *, double *,
c		 double *, double *, double *, double *, double *, int *);
c     
c     
c     Call from pair_meam.cpp has the form:
c     
c    meam_dens_init_(&i,&nmax,ntype,type,fmap,&x[0][0],
c	       &numneigh[i],firstneigh[i],&numneigh_full[i],firstneigh_full[i],
c	       &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
c	       rho0,&arho1[0][0],&arho2[0][0],arho2b,
c	       &arho3[0][0],&arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],&errorflag);
c     

      subroutine meam_dens_init(i, nmax,
     $     ntype, type, fmap, x,
     $     numneigh, firstneigh, 
     $     numneigh_full, firstneigh_full,
     $     scrfcn, dscrfcn, fcpair, rho0, arho1, arho2, arho2b,
     $     arho3, arho3b, t_ave, tsq_ave, errorflag)

      use meam_data
      implicit none

      integer i, nmax, ntype, type, fmap
      real*8  x
      integer numneigh, firstneigh, numneigh_full, firstneigh_full
      real*8 scrfcn, dscrfcn, fcpair
      real*8 rho0, arho1, arho2
      real*8 arho2b, arho3, arho3b, t_ave, tsq_ave
      integer errorflag
      integer j,jn

      dimension x(3,nmax)
      dimension type(nmax), fmap(ntype)
      dimension firstneigh(numneigh), firstneigh_full(numneigh_full)
      dimension scrfcn(numneigh), dscrfcn(numneigh), fcpair(numneigh)
      dimension rho0(nmax), arho1(3,nmax), arho2(6,nmax)
      dimension arho2b(nmax), arho3(10,nmax), arho3b(3,nmax)
      dimension t_ave(3,nmax), tsq_ave(3,nmax)

      errorflag = 0

c     Compute screening function and derivatives
      call getscreen(i, nmax, scrfcn, dscrfcn, fcpair, x,
     $     numneigh, firstneigh,
     $     numneigh_full, firstneigh_full,
     $     ntype, type, fmap)

c     Calculate intermediate density terms to be communicated
      call calc_rho1(i, nmax, ntype, type, fmap, x, 
     $     numneigh, firstneigh,
     $     scrfcn, fcpair, rho0, arho1, arho2, arho2b,
     $     arho3, arho3b, t_ave, tsq_ave)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getscreen(i, nmax, scrfcn, dscrfcn, fcpair, x, 
     $     numneigh, firstneigh,
     $     numneigh_full, firstneigh_full,
     $     ntype, type, fmap)

      use meam_data
      implicit none

      integer i, nmax
      real*8  scrfcn, dscrfcn, fcpair, x
      integer numneigh, firstneigh, numneigh_full, firstneigh_full
      integer ntype, type, fmap

      dimension scrfcn(numneigh), dscrfcn(numneigh)
      dimension fcpair(numneigh), x(3,nmax)
      dimension firstneigh(numneigh), firstneigh_full(numneigh_full)
      dimension type(nmax), fmap(ntype)

      integer jn,j,kn,k
      integer elti,eltj,eltk
      real*8 xitmp,yitmp,zitmp,delxij,delyij,delzij,rij2,rij
      real*8 xjtmp,yjtmp,zjtmp,delxik,delyik,delzik,rik2,rik
      real*8 xktmp,yktmp,zktmp,delxjk,delyjk,delzjk,rjk2,rjk
      real*8 xik,xjk,sij,fcij,sfcij,dfcij,sikj,dfikj,cikj
      real*8 Cmin,Cmax,delc,ebound,rbound,a,coef1,coef2
      real*8 coef1a,coef1b,coef2a,coef2b
      real*8 dcikj
      real*8 dC1a,dC1b,dC2a,dC2b
      real*8 rnorm,fc,dfc,drinv

      drinv = 1.d0/delr_meam
      elti = fmap(type(i))

      if (elti.gt.0) then

        xitmp = x(1,i)
        yitmp = x(2,i)
        zitmp = x(3,i)
        
        do jn = 1,numneigh
          j = firstneigh(jn)

          eltj = fmap(type(j))
          if (eltj.gt.0) then
          
c     First compute screening function itself, sij
            xjtmp = x(1,j)
            yjtmp = x(2,j)
            zjtmp = x(3,j)
            delxij = xjtmp - xitmp
            delyij = yjtmp - yitmp
            delzij = zjtmp - zitmp
            rij2 = delxij*delxij + delyij*delyij + delzij*delzij
            rij = sqrt(rij2)
            if (rij.gt.rc_meam) then
              fcij = 0.0
              dfcij = 0.d0
              sij = 0.d0
            else
              rnorm = (rc_meam-rij)*drinv
              call screen(i, j, nmax, x, rij2, sij,
     $             numneigh_full, firstneigh_full, ntype, type, fmap)
              call dfcut(rnorm,fc,dfc)
              fcij = fc
              dfcij = dfc*drinv
            endif
            
c     Now compute derivatives
            dscrfcn(jn) = 0.d0
            sfcij = sij*fcij
            if (sfcij.eq.0.d0.or.sfcij.eq.1.d0) goto 100
            rbound = ebound_meam(elti,eltj) * rij2
            do kn = 1,numneigh_full
              k = firstneigh_full(kn)
              if (k.eq.j) goto 10
              eltk = fmap(type(k))
              if (eltk.eq.0) goto 10
              xktmp = x(1,k)
              yktmp = x(2,k)
              zktmp = x(3,k)
              delxjk = xktmp - xjtmp
              delyjk = yktmp - yjtmp
              delzjk = zktmp - zjtmp
              rjk2 = delxjk*delxjk + delyjk*delyjk + delzjk*delzjk
              if (rjk2.gt.rbound) goto 10
              delxik = xktmp - xitmp
              delyik = yktmp - yitmp
              delzik = zktmp - zitmp
              rik2 = delxik*delxik + delyik*delyik + delzik*delzik
              if (rik2.gt.rbound) goto 10
              xik = rik2/rij2
              xjk = rjk2/rij2
              a =  1 - (xik-xjk)*(xik-xjk)
c     if a < 0, then ellipse equation doesn't describe this case and
c     atom k can't possibly screen i-j
              if (a.le.0.d0) goto 10
              cikj = (2.d0*(xik+xjk) + a - 2.d0)/a
              Cmax = Cmax_meam(elti,eltj,eltk)
              Cmin = Cmin_meam(elti,eltj,eltk)
              if (cikj.ge.Cmax) then
                goto 10
c     Note that cikj may be slightly negative (within numerical
c     tolerance) if atoms are colinear, so don't reject that case here
c     (other negative cikj cases were handled by the test on "a" above)
c     Note that we never have 0<cikj<Cmin here, else sij=0 (rejected above)
              else
                delc = Cmax - Cmin
                cikj = (cikj-Cmin)/delc
                call dfcut(cikj,sikj,dfikj)
                coef1 = dfikj/(delc*sikj)
                call dCfunc(rij2,rik2,rjk2,dCikj)
                dscrfcn(jn) = dscrfcn(jn) + coef1*dCikj
              endif
 10           continue
            enddo
            coef1 = sfcij
            coef2 = sij*dfcij/rij
            dscrfcn(jn) = dscrfcn(jn)*coef1 - coef2
 100        continue
            
            scrfcn(jn) = sij
            fcpair(jn) = fcij

          endif
          
        enddo
      
      endif

      return
      end
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine calc_rho1(i, nmax, ntype, type, fmap, x, 
     $     numneigh, firstneigh,
     $     scrfcn, fcpair, rho0, arho1, arho2, arho2b,
     $     arho3, arho3b, t_ave, tsq_ave)

      use meam_data
      implicit none

      integer i, nmax, ntype, type, fmap
      real*8  x
      integer numneigh, firstneigh
      real*8 scrfcn, fcpair, rho0, arho1, arho2
      real*8 arho2b, arho3, arho3b, t_ave, tsq_ave

      dimension type(nmax), fmap(ntype), x(3,nmax)
      dimension firstneigh(numneigh)
      dimension scrfcn(numneigh), fcpair(numneigh)
      dimension rho0(nmax), arho1(3,nmax), arho2(6,nmax)
      dimension arho2b(nmax), arho3(10,nmax), arho3b(3,nmax)
      dimension t_ave(3,nmax), tsq_ave(3,nmax)

      integer jn,j,m,n,p,elti,eltj
      integer nv2,nv3
      real*8 xtmp,ytmp,ztmp,delij(3),rij2,rij,sij
      real*8 ai,aj,rhoa0j,rhoa1j,rhoa2j,rhoa3j,A1j,A2j,A3j
      real*8 G,Gbar,gam,shp(3)
      real*8 ro0i,ro0j
      real*8 rhoa0i,rhoa1i,rhoa2i,rhoa3i,A1i,A2i,A3i

      elti = fmap(type(i))
      xtmp = x(1,i)
      ytmp = x(2,i)
      ztmp = x(3,i)
      do jn = 1,numneigh
        if (scrfcn(jn).ne.0.d0) then
          j = firstneigh(jn)
          sij = scrfcn(jn)*fcpair(jn)
          delij(1) = x(1,j) - xtmp
          delij(2) = x(2,j) - ytmp
          delij(3) = x(3,j) - ztmp
          rij2 = delij(1)*delij(1) + delij(2)*delij(2) 
     $         + delij(3)*delij(3)
          if (rij2.lt.cutforcesq) then
            eltj = fmap(type(j))
            rij = sqrt(rij2)
            ai = rij/re_meam(elti,elti) - 1.d0
            aj = rij/re_meam(eltj,eltj) - 1.d0
            ro0i = rho0_meam(elti)
            ro0j = rho0_meam(eltj)
            rhoa0j = ro0j*exp(-beta0_meam(eltj)*aj)*sij
            rhoa1j = ro0j*exp(-beta1_meam(eltj)*aj)*sij
            rhoa2j = ro0j*exp(-beta2_meam(eltj)*aj)*sij
            rhoa3j = ro0j*exp(-beta3_meam(eltj)*aj)*sij
            rhoa0i = ro0i*exp(-beta0_meam(elti)*ai)*sij
            rhoa1i = ro0i*exp(-beta1_meam(elti)*ai)*sij
            rhoa2i = ro0i*exp(-beta2_meam(elti)*ai)*sij
            rhoa3i = ro0i*exp(-beta3_meam(elti)*ai)*sij
            if (ialloy.eq.1) then
              rhoa1j = rhoa1j * t1_meam(eltj)
              rhoa2j = rhoa2j * t2_meam(eltj)
              rhoa3j = rhoa3j * t3_meam(eltj)
              rhoa1i = rhoa1i * t1_meam(elti)
              rhoa2i = rhoa2i * t2_meam(elti)
              rhoa3i = rhoa3i * t3_meam(elti)
            endif
            rho0(i) = rho0(i) + rhoa0j
            rho0(j) = rho0(j) + rhoa0i
c For ialloy = 2, use single-element value (not average)
            if (ialloy.ne.2) then
              t_ave(1,i) = t_ave(1,i) + t1_meam(eltj)*rhoa0j
              t_ave(2,i) = t_ave(2,i) + t2_meam(eltj)*rhoa0j
              t_ave(3,i) = t_ave(3,i) + t3_meam(eltj)*rhoa0j
              t_ave(1,j) = t_ave(1,j) + t1_meam(elti)*rhoa0i
              t_ave(2,j) = t_ave(2,j) + t2_meam(elti)*rhoa0i
              t_ave(3,j) = t_ave(3,j) + t3_meam(elti)*rhoa0i
            endif
            if (ialloy.eq.1) then
              tsq_ave(1,i) = tsq_ave(1,i) +
     $             t1_meam(eltj)*t1_meam(eltj)*rhoa0j
              tsq_ave(2,i) = tsq_ave(2,i) +
     $             t2_meam(eltj)*t2_meam(eltj)*rhoa0j
              tsq_ave(3,i) = tsq_ave(3,i) +
     $             t3_meam(eltj)*t3_meam(eltj)*rhoa0j
              tsq_ave(1,j) = tsq_ave(1,j) +
     $             t1_meam(elti)*t1_meam(elti)*rhoa0i
              tsq_ave(2,j) = tsq_ave(2,j) +
     $             t2_meam(elti)*t2_meam(elti)*rhoa0i
              tsq_ave(3,j) = tsq_ave(3,j) +
     $             t3_meam(elti)*t3_meam(elti)*rhoa0i
            endif
            Arho2b(i) = Arho2b(i) + rhoa2j
            Arho2b(j) = Arho2b(j) + rhoa2i

            A1j = rhoa1j/rij
            A2j = rhoa2j/rij2
            A3j = rhoa3j/(rij2*rij)
            A1i = rhoa1i/rij
            A2i = rhoa2i/rij2
            A3i = rhoa3i/(rij2*rij)
            nv2 = 1
            nv3 = 1
            do m = 1,3
              Arho1(m,i) = Arho1(m,i) + A1j*delij(m)
              Arho1(m,j) = Arho1(m,j) - A1i*delij(m)
              Arho3b(m,i) = Arho3b(m,i) + rhoa3j*delij(m)/rij
              Arho3b(m,j) = Arho3b(m,j) - rhoa3i*delij(m)/rij
              do n = m,3
                Arho2(nv2,i) = Arho2(nv2,i) + A2j*delij(m)*delij(n)
                Arho2(nv2,j) = Arho2(nv2,j) + A2i*delij(m)*delij(n)
                nv2 = nv2+1
                do p = n,3
                  Arho3(nv3,i) = Arho3(nv3,i) 
     $                 + A3j*delij(m)*delij(n)*delij(p)
                  Arho3(nv3,j) = Arho3(nv3,j)
     $                 - A3i*delij(m)*delij(n)*delij(p)
                  nv3 = nv3+1
                enddo
              enddo
            enddo
          endif
        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine screen(i, j, nmax, x, rijsq, sij,
     $     numneigh_full, firstneigh_full, ntype, type, fmap)
c     Screening function
c     Inputs:  i = atom 1 id (integer)
c     j = atom 2 id (integer)
c     rijsq = squared distance between i and j
c     Outputs: sij = screening function 
      use meam_data
      implicit none

      integer i,j,nmax,k,nk,m
      real*8  x,rijsq,sij
      integer numneigh_full, firstneigh_full
      integer ntype, type, fmap

      dimension x(3,nmax), firstneigh_full(numneigh_full)
      dimension type(nmax), fmap(ntype)

      integer elti,eltj,eltk
      real*8 delxik,delyik,delzik
      real*8 delxjk,delyjk,delzjk
      real*8 riksq,rjksq,xik,xjk,cikj,a,delc,sikj,fcij,rij
      real*8 Cmax,Cmin,rbound

      sij = 1.d0
      elti = fmap(type(i))
      eltj = fmap(type(j))

c     if rjksq > ebound*rijsq, atom k is definitely outside the ellipse
      rbound = ebound_meam(elti,eltj)*rijsq

      do nk = 1,numneigh_full
        k = firstneigh_full(nk)
        eltk = fmap(type(k))
        if (k.eq.j) goto 10
        delxjk = x(1,k) - x(1,j)
        delyjk = x(2,k) - x(2,j)
        delzjk = x(3,k) - x(3,j)
        rjksq = delxjk*delxjk + delyjk*delyjk + delzjk*delzjk
        if (rjksq.gt.rbound) goto 10
        delxik = x(1,k) - x(1,i)
        delyik = x(2,k) - x(2,i)
        delzik = x(3,k) - x(3,i)
        riksq = delxik*delxik + delyik*delyik + delzik*delzik
        if (riksq.gt.rbound) goto 10
        xik = riksq/rijsq
        xjk = rjksq/rijsq
        a = 1 - (xik-xjk)*(xik-xjk)
c     if a < 0, then ellipse equation doesn't describe this case and
c     atom k can't possibly screen i-j
        if (a.le.0.d0) goto 10
        cikj = (2.d0*(xik+xjk) + a - 2.d0)/a
        Cmax = Cmax_meam(elti,eltj,eltk)
        Cmin = Cmin_meam(elti,eltj,eltk)
        if (cikj.ge.Cmax) then
          goto 10
c     note that cikj may be slightly negative (within numerical
c     tolerance) if atoms are colinear, so don't reject that case here
c     (other negative cikj cases were handled by the test on "a" above)
        else if (cikj.le.Cmin) then
          sij = 0.d0
          goto 20
        else
          delc = Cmax - Cmin
          cikj = (cikj-Cmin)/delc
          call fcut(cikj,sikj)
        endif
        sij = sij * sikj
 10     continue
      enddo

 20   continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dsij(i,j,k,jn,nmax,numneigh,rij2,dsij1,dsij2,
     $     ntype,type,fmap,x,scrfcn,fcpair)
c     Inputs: i,j,k = id's of 3 atom triplet
c     jn = id of i-j pair
c     rij2 = squared distance between i and j
c     Outputs: dsij1 = deriv. of sij w.r.t. rik
c     dsij2 = deriv. of sij w.r.t. rjk
      use meam_data
      implicit none
      integer i,j,k,jn,nmax,numneigh
      integer elti,eltj,eltk
      real*8 rij2,rik2,rjk2,dsij1,dsij2
      integer ntype, type, fmap
      real*8 x, scrfcn, fcpair

      dimension type(nmax), fmap(ntype)
      dimension x(3,nmax), scrfcn(numneigh), fcpair(numneigh)

      real*8 dxik,dyik,dzik
      real*8 dxjk,dyjk,dzjk
      real*8 rbound,delc,sij,xik,xjk,cikj,sikj,dfc,a
      real*8 Cmax,Cmin,dCikj1,dCikj2

      sij = scrfcn(jn)*fcpair(jn)
      elti = fmap(type(i))
      eltj = fmap(type(j))
      eltk = fmap(type(k))
      Cmax = Cmax_meam(elti,eltj,eltk)
      Cmin = Cmin_meam(elti,eltj,eltk)

      dsij1 = 0.d0
      dsij2 = 0.d0
      if ((sij.ne.0.d0).and.(sij.ne.1.d0)) then
        rbound = rij2*ebound_meam(elti,eltj)
        delc = Cmax-Cmin
        dxjk = x(1,k) - x(1,j)
        dyjk = x(2,k) - x(2,j)
        dzjk = x(3,k) - x(3,j)
        rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
        if (rjk2.le.rbound) then
          dxik = x(1,k) - x(1,i)
          dyik = x(2,k) - x(2,i)
          dzik = x(3,k) - x(3,i)
          rik2 = dxik*dxik + dyik*dyik + dzik*dzik
          if (rik2.le.rbound) then
            xik = rik2/rij2
            xjk = rjk2/rij2
            a = 1 - (xik-xjk)*(xik-xjk)
            if (a.ne.0.d0) then
              cikj = (2.d0*(xik+xjk) + a - 2.d0)/a
              if (cikj.ge.Cmin.and.cikj.le.Cmax) then
                cikj = (cikj-Cmin)/delc
                call dfcut(cikj,sikj,dfc)
                call dCfunc2(rij2,rik2,rjk2,dCikj1,dCikj2)
                a = sij/delc*dfc/sikj
                dsij1 = a*dCikj1
                dsij2 = a*dCikj2
              endif
            endif
          endif
        endif
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fcut(xi,fc)
c     cutoff function
      implicit none
      real*8 xi,fc
      real*8 a
      if (xi.ge.1.d0) then
        fc = 1.d0
      else if (xi.le.0.d0) then
        fc = 0.d0
      else
        a = 1.d0-xi
        a = a*a
        a = a*a
        a = 1.d0-a
        fc = a*a       
c     fc = xi
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dfcut(xi,fc,dfc)
c     cutoff function and its derivative
      implicit none
      real*8 xi,fc,dfc,a,a3,a4
      if (xi.ge.1.d0) then
        fc = 1.d0
        dfc = 0.d0
      else if (xi.le.0.d0) then
        fc = 0.d0
        dfc = 0.d0
      else
        a = 1.d0-xi
        a3 = a*a*a
        a4 = a*a3
        fc = (1.d0-a4)**2
        dfc = 8*(1.d0-a4)*a3
c     fc = xi
c     dfc = 1.d0
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dCfunc(rij2,rik2,rjk2,dCikj)
c     Inputs: rij,rij2,rik2,rjk2
c     Outputs: dCikj = derivative of Cikj w.r.t. rij
      implicit none
      real*8 rij2,rik2,rjk2,dCikj
      real*8 rij4,a,b,denom

      rij4 = rij2*rij2
      a = rik2-rjk2
      b = rik2+rjk2
      denom = rij4 - a*a
      denom = denom*denom
      dCikj = -4*(-2*rij2*a*a + rij4*b + a*a*b)/denom
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dCfunc2(rij2,rik2,rjk2,dCikj1,dCikj2)
c     Inputs: rij,rij2,rik2,rjk2
c     Outputs: dCikj1 = derivative of Cikj w.r.t. rik
c     dCikj2 = derivative of Cikj w.r.t. rjk
      implicit none
      real*8 rij2,rik2,rjk2,dCikj1,dCikj2
      real*8 rij4,rik4,rjk4,a,b,denom

      rij4 = rij2*rij2
      rik4 = rik2*rik2
      rjk4 = rjk2*rjk2
      a = rik2-rjk2
      b = rik2+rjk2
      denom = rij4 - a*a
      denom = denom*denom
      dCikj1 = 4*rij2*(rij4 + rik4 + 2*rik2*rjk2 - 3*rjk4 - 2*rij2*a)/
     $     denom
      dCikj2 = 4*rij2*(rij4 - 3*rik4 + 2*rik2*rjk2 + rjk4 + 2*rij2*a)/
     $     denom
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      


