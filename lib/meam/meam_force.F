c     Extern "C" declaration has the form:
c     
c  void meam_force_(int *, int *, int *, double *, int *, int *, int *, double *,
c		 int *, int *, int *, int *, double *, double *,
c		 double *, double *, double *, double *, double *, double *,
c		 double *, double *, double *, double *, double *, double *,
c		 double *, double *, double *, double *, double *, int *);
c     
c     Call from pair_meam.cpp has the form:
c     
c    meam_force_(&i,&nmax,&eflag,&eng_vdwl,&ntype,type,fmap,&x[0][0],
c	       &numneigh[i],firstneigh[i],&numneigh_full[i],firstneigh_full[i],
c	       &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
c	       dgamma1,dgamma2,dgamma3,rho0,rho1,rho2,rho3,frhop,
c	       &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],&arho3b[0][0],
c	       &t_ave[0][0],&f[0][0],&strssa[0][0][0],&errorflag);
c     

      subroutine meam_force(i, nmax, eflag, eng_vdwl,
     $     ntype, type, fmap, x,
     $     numneigh, firstneigh, numneigh_full, firstneigh_full,
     $     scrfcn, dscrfcn, fcpair, 
     $     dGamma1, dGamma2, dGamma3, rho0, rho1, rho2, rho3, fp,
     $     Arho1, Arho2, Arho2b, Arho3, Arho3b, t_ave, f, strssa,
     $     errorflag)

      use meam_data
      implicit none

      integer eflag, nmax, ntype, type, fmap
      real*8  eng_vdwl, x
      integer numneigh, firstneigh, numneigh_full, firstneigh_full
      real*8  scrfcn, dscrfcn, fcpair
      real*8  dGamma1, dGamma2, dGamma3
      real*8  rho0, rho1, rho2, rho3, fp
      real*8  Arho1, Arho2, Arho2b
      real*8  Arho3, Arho3b
      real*8  t_ave, f, strssa
      integer errorflag

      dimension type(nmax), fmap(ntype)
      dimension x(3,nmax)
      dimension firstneigh(numneigh), firstneigh_full(numneigh_full)
      dimension scrfcn(numneigh), dscrfcn(numneigh), fcpair(numneigh)
      dimension dGamma1(nmax), dGamma2(nmax), dGamma3(nmax)
      dimension rho0(nmax), rho1(nmax), rho2(nmax), rho3(nmax), fp(nmax)
      dimension Arho1(3,nmax), Arho2(6,nmax), Arho2b(nmax)
      dimension Arho3(10,nmax), Arho3b(3,nmax)
      dimension t_ave(3,nmax), f(3,nmax), strssa(3,3,nmax)

      integer i,j,jn,k,kn,kk,m,n,p,q,strscomp
      integer nv2,nv3,elti,eltj,ind
      real*8 xitmp,yitmp,zitmp,delij(3),delref(3),rij2,rij,rij3
      real*8 delik(3),deljk(3)
      real*8 Eu,astar,astarp
      real*8 pp,phiforce,dUdrij,dUdsij,dUdrijm(3),force,forcem
      real*8 B,r,recip,phi,phip,rhop,a
      real*8 sij,fcij,dfcij,ds(3)
      real*8 a0,a1,a1i,a1j,a2,a2i,a2j
      real*8 a3i,a3j,a3i1,a3i2,a3j1,a3j2
      real*8 G,dG,Gbar,dGbar,gam,shpi(3),shpj(3),Z,denom
      real*8 ai,aj,ro0i,ro0j,invrei,invrej
      real*8 b0,rhoa0j,drhoa0j,rhoa0i,drhoa0i
      real*8 b1,rhoa1j,drhoa1j,rhoa1i,drhoa1i
      real*8 b2,rhoa2j,drhoa2j,rhoa2i,drhoa2i
      real*8 a3,a3a,b3,rhoa3j,drhoa3j,rhoa3i,drhoa3i
      real*8 drho0dr1,drho0dr2,drho0ds1,drho0ds2
      real*8 drho1dr1,drho1dr2,drho1ds1,drho1ds2
      real*8 drho1drm1(3),drho1drm2(3)
      real*8 drho2dr1,drho2dr2,drho2ds1,drho2ds2
      real*8 drho2drm1(3),drho2drm2(3)
      real*8 drho3dr1,drho3dr2,drho3ds1,drho3ds2
      real*8 drho3drm1(3),drho3drm2(3)
      real*8 dt1dr1,dt1dr2,dt1ds1,dt1ds2
      real*8 dt2dr1,dt2dr2,dt2ds1,dt2ds2
      real*8 dt3dr1,dt3dr2,dt3ds1,dt3ds2
      real*8 drhodr1,drhodr2,drhods1,drhods2,drhodrm1(3),drhodrm2(3)
      real*8 arg,arg1,arg2
      real*8 arg1i1,arg1j1,arg1i2,arg1j2,arg2i2,arg2j2
      real*8 arg1i3,arg1j3,arg2i3,arg2j3,arg3i3,arg3j3
      real*8 dsij1,dsij2,force1,force2
      real*8 t1i,t2i,t3i,t1j,t2j,t3j

      errorflag = 0

      strscomp = 0

c     Compute forces atom i

      elti = fmap(type(i))
      xitmp = x(1,i)
      yitmp = x(2,i)
      zitmp = x(3,i)

c     Treat each pair
      do jn = 1,numneigh
        
        if (scrfcn(jn).ne.0.d0) then
          j = firstneigh(jn)

          sij = scrfcn(jn)*fcpair(jn)
          delij(1) = x(1,j) - xitmp
          delij(2) = x(2,j) - yitmp
          delij(3) = x(3,j) - zitmp
          rij2 = delij(1)*delij(1) + delij(2)*delij(2)
     $         + delij(3)*delij(3)
          if (rij2.lt.cutforcesq) then
            rij = sqrt(rij2)
            r = rij
            eltj = fmap(type(j))

c     Compute phi and phip
            ind = eltind(elti,eltj)
            pp = rij*rdrar + 1.0D0
            kk = pp
            kk = min(kk,nrar-1)
            pp = pp - kk
            pp = min(pp,1.0D0)
            phi = ((phirar3(kk,ind)*pp + phirar2(kk,ind))*pp 
     $           + phirar1(kk,ind))*pp + phirar(kk,ind)
            phip = (phirar6(kk,ind)*pp + phirar5(kk,ind))*pp 
     $           + phirar4(kk,ind)
            recip = 1.0d0/r

            if (eflag.eq.1) eng_vdwl = eng_vdwl + phi*sij

c     write(1,*) "force_meamf: phi: ",phi
c     write(1,*) "force_meamf: phip: ",phip

c     Compute pair densities and derivatives
            invrei = 1.d0/re_meam(elti,elti)
            ai = rij*invrei - 1.d0
            ro0i = rho0_meam(elti)
            rhoa0i = ro0i*exp(-beta0_meam(elti)*ai)
            drhoa0i = -beta0_meam(elti)*invrei*rhoa0i
            rhoa1i = ro0i*exp(-beta1_meam(elti)*ai)
            drhoa1i = -beta1_meam(elti)*invrei*rhoa1i
            rhoa2i = ro0i*exp(-beta2_meam(elti)*ai)
            drhoa2i = -beta2_meam(elti)*invrei*rhoa2i
            rhoa3i = ro0i*exp(-beta3_meam(elti)*ai)
            drhoa3i = -beta3_meam(elti)*invrei*rhoa3i

            if (elti.ne.eltj) then
              invrej = 1.d0/re_meam(eltj,eltj)
              aj = rij*invrej - 1.d0
              ro0j = rho0_meam(eltj)
              rhoa0j = ro0j*exp(-beta0_meam(eltj)*aj)
              drhoa0j = -beta0_meam(eltj)*invrej*rhoa0j
              rhoa1j = ro0j*exp(-beta1_meam(eltj)*aj)
              drhoa1j = -beta1_meam(eltj)*invrej*rhoa1j
              rhoa2j = ro0j*exp(-beta2_meam(eltj)*aj)
              drhoa2j = -beta2_meam(eltj)*invrej*rhoa2j
              rhoa3j = ro0j*exp(-beta3_meam(eltj)*aj)
              drhoa3j = -beta3_meam(eltj)*invrej*rhoa3j
            else
              rhoa0j = rhoa0i
              drhoa0j = drhoa0i
              rhoa1j = rhoa1i
              drhoa1j = drhoa1i
              rhoa2j = rhoa2i
              drhoa2j = drhoa2i
              rhoa3j = rhoa3i
              drhoa3j = drhoa3i
            endif

            nv2 = 1
            nv3 = 1
            arg1i1 = 0.d0
            arg1j1 = 0.d0
            arg1i2 = 0.d0
            arg1j2 = 0.d0
            arg1i3 = 0.d0
            arg1j3 = 0.d0
            arg3i3 = 0.d0
            arg3j3 = 0.d0
            do n = 1,3
              do p = n,3
                do q = p,3
                  arg = delij(n)*delij(p)*delij(q)*v3D(nv3)
                  arg1i3 = arg1i3 + Arho3(nv3,i)*arg
                  arg1j3 = arg1j3 - Arho3(nv3,j)*arg
                  nv3 = nv3+1
                enddo
                arg = delij(n)*delij(p)*v2D(nv2)
                arg1i2 = arg1i2 + Arho2(nv2,i)*arg
                arg1j2 = arg1j2 + Arho2(nv2,j)*arg
                nv2 = nv2+1
              enddo
              arg1i1 = arg1i1 + Arho1(n,i)*delij(n)
              arg1j1 = arg1j1 - Arho1(n,j)*delij(n)
              arg3i3 = arg3i3 + Arho3b(n,i)*delij(n)
              arg3j3 = arg3j3 - Arho3b(n,j)*delij(n)             
            enddo

c     rho0 terms
            drho0dr1 = drhoa0j * sij
            drho0dr2 = drhoa0i * sij
            
c     rho1 terms  
            a1 = 2*sij/rij              
            drho1dr1 = a1*(drhoa1j-rhoa1j/rij)*arg1i1
            drho1dr2 = a1*(drhoa1i-rhoa1i/rij)*arg1j1
            a1 = 2.d0*sij/rij
            do m = 1,3
              drho1drm1(m) = a1*rhoa1j*Arho1(m,i)
              drho1drm2(m) = -a1*rhoa1i*Arho1(m,j)
            enddo

c     rho2 terms
            a2 = 2*sij/rij2
            drho2dr1 = a2*(drhoa2j - 2*rhoa2j/rij)*arg1i2 
     $           - 2.d0/3.d0*Arho2b(i)*drhoa2j*sij
            drho2dr2 = a2*(drhoa2i - 2*rhoa2i/rij)*arg1j2
     $           - 2.d0/3.d0*Arho2b(j)*drhoa2i*sij
            a2 = 4*sij/rij2
            do m = 1,3
              drho2drm1(m) = 0.d0
              drho2drm2(m) = 0.d0
              do n = 1,3
                drho2drm1(m) = drho2drm1(m) 
     $               + Arho2(vind2D(m,n),i)*delij(n)
                drho2drm2(m) = drho2drm2(m)
     $               - Arho2(vind2D(m,n),j)*delij(n)
              enddo
              drho2drm1(m) = a2*rhoa2j*drho2drm1(m)
              drho2drm2(m) = -a2*rhoa2i*drho2drm2(m)
            enddo

c     rho3 terms
            rij3 = rij*rij2
            a3 = 2*sij/rij3
            a3a = 6.d0/5.d0*sij/rij
            drho3dr1 = a3*(drhoa3j - 3*rhoa3j/rij)*arg1i3 
     $           - a3a*(drhoa3j - rhoa3j/rij)*arg3i3
            drho3dr2 = a3*(drhoa3i - 3*rhoa3i/rij)*arg1j3 
     $           - a3a*(drhoa3i - rhoa3i/rij)*arg3j3
            a3 = 6*sij/rij3
            a3a = 6*sij/(5*rij)
            do m = 1,3
              drho3drm1(m) = 0.d0
              drho3drm2(m) = 0.d0
              nv2 = 1
              do n = 1,3
                do p = n,3
                  arg = delij(n)*delij(p)*v2D(nv2)
                  drho3drm1(m) = drho3drm1(m)
     $                 + Arho3(vind3D(m,n,p),i)*arg
                  drho3drm2(m) = drho3drm2(m)
     $                 + Arho3(vind3D(m,n,p),j)*arg
                  nv2 = nv2 + 1
                enddo
              enddo
              drho3drm1(m) = (a3*drho3drm1(m) - a3a*Arho3b(m,i))
     $             *rhoa3j
              drho3drm2(m) = (-a3*drho3drm2(m) + a3a*Arho3b(m,j))
     $             *rhoa3i
            enddo

c     Compute derivatives of weighting functions t wrt rij
            t1i = t_ave(1,i)
            t2i = t_ave(2,i)
            t3i = t_ave(3,i)
            t1j = t_ave(1,j)
            t2j = t_ave(2,j)
            t3j = t_ave(3,j)

            ai = 0.d0
            if( rho0(i) .ne. 0.d0 ) then
              ai = drhoa0j*sij/rho0(i)
            end if
            aj = 0.d0
            if( rho0(j) .ne. 0.d0 ) then
              aj = drhoa0i*sij/rho0(j)
            end if

            dt1dr1 = ai*(t1_meam(eltj)-t1i)
            dt1dr2 = aj*(t1_meam(elti)-t1j)
            dt2dr1 = ai*(t2_meam(eltj)-t2i)
            dt2dr2 = aj*(t2_meam(elti)-t2j)
            dt3dr1 = ai*(t3_meam(eltj)-t3i)
            dt3dr2 = aj*(t3_meam(elti)-t3j)

c     Compute derivatives of total density wrt rij, sij and rij(3)
            call get_shpfcn(shpi,lattce_meam(elti,elti))
            call get_shpfcn(shpj,lattce_meam(eltj,eltj))
            drhodr1 = dGamma1(i)*drho0dr1
     $           + dGamma2(i)*
     $           (dt1dr1*rho1(i)+t1i*drho1dr1
     $           + dt2dr1*rho2(i)+t2i*drho2dr1
     $           + dt3dr1*rho3(i)+t3i*drho3dr1)
     $           - dGamma3(i)*
     $           (shpi(1)*dt1dr1+shpi(2)*dt2dr1+shpi(3)*dt3dr1)
            drhodr2 = dGamma1(j)*drho0dr2
     $           + dGamma2(j)*
     $           (dt1dr2*rho1(j)+t1j*drho1dr2
     $           + dt2dr2*rho2(j)+t2j*drho2dr2
     $           + dt3dr2*rho3(j)+t3j*drho3dr2)
     $           - dGamma3(j)*
     $           (shpj(1)*dt1dr2+shpj(2)*dt2dr2+shpj(3)*dt3dr2)
            do m = 1,3
              drhodrm1(m) = 0.d0
              drhodrm2(m) = 0.d0
              drhodrm1(m) = dGamma2(i)*
     $             (t1i*drho1drm1(m)
     $             + t2i*drho2drm1(m)
     $             + t3i*drho3drm1(m))
              drhodrm2(m) = dGamma2(j)*
     $             (t1j*drho1drm2(m)
     $             + t2j*drho2drm2(m)
     $             + t3j*drho3drm2(m))
            enddo

c     Compute derivatives wrt sij, but only if necessary
            if (dscrfcn(jn).ne.0.d0) then
              drho0ds1 = rhoa0j
              drho0ds2 = rhoa0i
              a1 = 2.d0/rij
              drho1ds1 = a1*rhoa1j*arg1i1
              drho1ds2 = a1*rhoa1i*arg1j1
              a2 = 2.d0/rij2
              drho2ds1 = a2*rhoa2j*arg1i2 
     $             - 2.d0/3.d0*Arho2b(i)*rhoa2j
              drho2ds2 = a2*rhoa2i*arg1j2 
     $             - 2.d0/3.d0*Arho2b(j)*rhoa2i
              a3 = 2.d0/rij3
              a3a = 6.d0/(5.d0*rij)
              drho3ds1 = a3*rhoa3j*arg1i3 - a3a*rhoa3j*arg3i3
              drho3ds2 = a3*rhoa3i*arg1j3 - a3a*rhoa3i*arg3j3

              ai = 0.d0
              if( rho0(i) .ne. 0.d0 ) then
                ai = rhoa0j/rho0(i)
              end if
              aj = 0.d0
              if( rho0(j) .ne. 0.d0 ) then
                aj = rhoa0i/rho0(j)
              end if

              dt1ds1 = ai*(t1_meam(eltj)-t1i)
              dt1ds2 = aj*(t1_meam(elti)-t1j)
              dt2ds1 = ai*(t2_meam(eltj)-t2i)
              dt2ds2 = aj*(t2_meam(elti)-t2j)
              dt3ds1 = ai*(t3_meam(eltj)-t3i)
              dt3ds2 = aj*(t3_meam(elti)-t3j)
              drhods1 = dGamma1(i)*drho0ds1
     $             + dGamma2(i)*
     $             (dt1ds1*rho1(i)+t1i*drho1ds1
     $             + dt2ds1*rho2(i)+t2i*drho2ds1
     $             + dt3ds1*rho3(i)+t3i*drho3ds1)
     $             - dGamma3(i)*
     $             (shpi(1)*dt1ds1+shpi(2)*dt2ds1+shpi(3)*dt3ds1)
              drhods2 = dGamma1(j)*drho0ds2
     $             + dGamma2(j)*
     $             (dt1ds2*rho1(j)+t1j*drho1ds2
     $             + dt2ds2*rho2(j)+t2j*drho2ds2
     $             + dt3ds2*rho3(j)+t3j*drho3ds2)
     $             - dGamma3(j)*
     $             (shpj(1)*dt1ds2+shpj(2)*dt2ds2+shpj(3)*dt3ds2)
            endif

c     Compute derivatives of energy wrt rij, sij and rij(3)
            dUdrij = phip*sij
     $           + fp(i)*drhodr1 + fp(j)*drhodr2
            dUdsij = 0.d0
            if (dscrfcn(jn).ne.0.d0) then
              dUdsij = phi
     $             + fp(i)*drhods1 + fp(j)*drhods2
            endif
            do m = 1,3
              dUdrijm(m) = fp(i)*drhodrm1(m) + fp(j)*drhodrm2(m)
            enddo
            
c     Add the part of the force due to dUdrij and dUdsij
            force = dUdrij*recip + dUdsij*dscrfcn(jn)
            do m = 1,3
              forcem = delij(m)*force + dUdrijm(m)
              f(m,i) = f(m,i) + forcem
              f(m,j) = f(m,j) - forcem
              if (strscomp>0) then
                do n = 1,3
                  arg = delij(n)*forcem
                  strssa(n,m,i) = strssa(n,m,i) + arg
                  strssa(n,m,j) = strssa(n,m,j) + arg
                enddo
              endif
            enddo
c     write(1,*)"forcei1: ",f(1,i)," ",f(2,i)," ",f(3,i)
c     write(1,*)"forcej1: ",f(1,j)," ",f(2,j)," ",f(3,j)

c     Now compute forces on other atoms k due to change in sij
            if (sij.eq.0.d0.or.sij.eq.1.d0) goto 100
            do kn = 1,numneigh_full
              k = firstneigh_full(kn)
              if (k.ne.j) then
                call dsij(i,j,k,jn,nmax,numneigh,rij2,dsij1,dsij2,
     $               ntype,type,fmap,x,scrfcn,fcpair)
                if (dsij1.ne.0.d0.or.dsij2.ne.0.d0) then
                  force1 = dUdsij*dsij1
                  force2 = dUdsij*dsij2
                  do m = 1,3
                    delik(m) = x(m,k) - x(m,i)
                    deljk(m) = x(m,k) - x(m,j)
                  enddo
                  do m = 1,3
                    f(m,i) = f(m,i) + force1*delik(m)
                    f(m,j) = f(m,j) + force2*deljk(m)
                    f(m,k) = f(m,k) - force1*delik(m)
     $                   - force2*deljk(m)
                    if (strscomp>0) then
                      do n = m,3
                        arg1 = force1*delik(m)*delik(n)
                        arg2 = force2*deljk(m)*deljk(n)
                        strssa(m,n,i) = strssa(m,n,i) + arg1
                        strssa(m,n,j) = strssa(m,n,j) + arg2
                        strssa(m,n,k) = strssa(m,n,k) + 
     $                       arg1 + arg2
                        if (n.ne.m) then
                          strssa(n,m,i) = strssa(n,m,i) + arg1
                          strssa(n,m,j) = strssa(n,m,j) + arg2
                          strssa(n,m,k) = strssa(n,m,k) 
     $                         + arg1 + arg2
                        endif
                      enddo
                    endif
                  enddo

                endif
              endif
c     end of k loop
            enddo
          endif
 100      continue
        endif
c     end of j loop
      enddo

      return
      end
