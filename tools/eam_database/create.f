C author: X. W. Zhou, xzhou@sandia.gov
c      open(unit=5,file='a.i')
      call inter
c      close(5)
      call writeset
      stop
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c main subroutine.                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inter
      character*80 atomtype,atommatch,outfile,outelem
      namelist /funccard/ atomtype
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16)
      common /pass2/ ielement(16),amass(16),Fr(5000,16),
     *   rhor(5000,16),z2r(5000,16,16),ntypes,blat(16),
     *   nrho,drho,nr,dr,rc,outfile,outelem
      ntypes=0
10    continue
      atomtype='none'
      read(5,funccard)
      if (atomtype .eq. 'none') goto 1200
      open(unit=10,file='EAM_code',form='FORMATTED',status='OLD')
11    read(10,9501,end=1210)atommatch
9501  format(a80)
      if (atomtype .eq. atommatch) then
         ntypes=ntypes+1
         length=len_trim(outfile)
         if (length .eq. len(outfile)) then
            outfile = atomtype
         else
            outfile = outfile(1:length)//atomtype
         endif
         length=len_trim(outelem)
         if (length .eq. len(outelem)) then
            outelem = atomtype
         else
            outelem = outelem(1:length)//' '//atomtype
         endif
         read(10,*) re(ntypes)
         read(10,*) fe(ntypes)
         read(10,*) rhoe(ntypes)
         read(10,*) rhos(ntypes)
         read(10,*) alpha(ntypes)
         read(10,*) beta(ntypes)
         read(10,*) A(ntypes)
         read(10,*) B(ntypes)
         read(10,*) cai(ntypes)
         read(10,*) ramda(ntypes)
         read(10,*) Fi0(ntypes)
         read(10,*) Fi1(ntypes)
         read(10,*) Fi2(ntypes)
         read(10,*) Fi3(ntypes)
         read(10,*) Fm0(ntypes)
         read(10,*) Fm1(ntypes)
         read(10,*) Fm2(ntypes)
         read(10,*) Fm3(ntypes)
         read(10,*) fnn(ntypes)
         read(10,*) Fn(ntypes)
         read(10,*) ielement(ntypes)
         read(10,*) amass(ntypes)
         read(10,*) Fm4(ntypes)
         read(10,*) beta1(ntypes)
         read(10,*) ramda1(ntypes)
         read(10,*) rhol(ntypes)
         read(10,*) rhoh(ntypes)
         blat(ntypes)=sqrt(2.0)*re(ntypes)
         rhoin(ntypes)=rhol(ntypes)*rhoe(ntypes)
         rhoout(ntypes)=rhoh(ntypes)*rhoe(ntypes)
      else
         do 1 i=1,27
1        read(10,*)vtmp
         goto 11
      endif
      close(10)
      goto 10
1210  write(6,*)'error: atom type ',atomtype,' not found'
      stop
1200  continue
      nr=2000
      nrho=2000
      alatmax=blat(1)
      rhoemax=rhoe(1)
      do 2 i=2,ntypes
         if (alatmax .lt. blat(i)) alatmax=blat(i)
         if (rhoemax .lt. rhoe(i)) rhoemax=rhoe(i)
2     continue
      rc=sqrt(10.0)/2.0*alatmax
      rst=0.5
      dr=rc/(nr-1.0)
      fmax=-1.0
      do 3 i1=1,ntypes
      do 3 i2=1,i1
      if ( i1 .eq. i2) then
         do 4 i=1,nr
            r=(i-1.0)*dr
            if (r .lt. rst) r=rst
            call prof(i1,r,fvalue)
            if (fmax .lt. fvalue) fmax=fvalue
            rhor(i,i1)=fvalue
            call pair(i1,i2,r,psi)
            z2r(i,i1,i2)=r*psi
4        continue
      else
         do 5 i=1,nr
            r=(i-1.0)*dr
            if (r .lt. rst) r=rst
            call pair(i1,i2,r,psi)
            z2r(i,i1,i2)=r*psi
            z2r(i,i2,i1)=z2r(i,i1,i2)
5        continue
      endif
3     continue
      rhom=fmax
      if (rhom .lt. 2.0*rhoemax) rhom=2.0*rhoemax
      if (rhom .lt. 100.0) rhom=100.0
      drho=rhom/(nrho-1.0)
      do 6 it=1,ntypes
      do 7 i=1,nrho
         rhoF=(i-1.0)*drho
         call embed(it,rhoF,emb)
         Fr(i,it)=emb
7     continue
6     continue
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the electron density.                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prof(it,r,f)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16)
      f=fe(it)*exp(-beta1(it)*(r/re(it)-1.0))
      f=f/(1.0+(r/re(it)-ramda1(it))**20)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the pair potential.                  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pair(it1,it2,r,psi)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16)
      if (it1 .eq. it2) then
         psi1=A(it1)*exp(-alpha(it1)*(r/re(it1)-1.0))
         psi1=psi1/(1.0+(r/re(it1)-cai(it1))**20)
         psi2=B(it1)*exp(-beta(it1)*(r/re(it1)-1.0))
         psi2=psi2/(1.0+(r/re(it1)-ramda(it1))**20)
         psi=psi1-psi2
      else
         psi1=A(it1)*exp(-alpha(it1)*(r/re(it1)-1.0))
         psi1=psi1/(1.0+(r/re(it1)-cai(it1))**20)
         psi2=B(it1)*exp(-beta(it1)*(r/re(it1)-1.0))
         psi2=psi2/(1.0+(r/re(it1)-ramda(it1))**20)
         psia=psi1-psi2
         psi1=A(it2)*exp(-alpha(it2)*(r/re(it2)-1.0))
         psi1=psi1/(1.0+(r/re(it2)-cai(it2))**20)
         psi2=B(it2)*exp(-beta(it2)*(r/re(it2)-1.0))
         psi2=psi2/(1.0+(r/re(it2)-ramda(it2))**20)
         psib=psi1-psi2
         call prof(it1,r,f1)
         call prof(it2,r,f2)
         psi=0.5*(f2/f1*psia+f1/f2*psib)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the embedding energy.                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine embed(it,rho,emb)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16)
      if (rho .lt. rhoe(it)) then
         Fm33=Fm3(it)
      else 
         Fm33=Fm4(it)
      endif
      if (rho .lt. rhoin(it)) then
         emb=Fi0(it)+
     *       Fi1(it)*(rho/rhoin(it)-1.0)+
     *       Fi2(it)*(rho/rhoin(it)-1.0)**2+
     *       Fi3(it)*(rho/rhoin(it)-1.0)**3
      else if (rho .lt. rhoout(it)) then
         emb=Fm0(it)+
     *       Fm1(it)*(rho/rhoe(it)-1.0)+
     *       Fm2(it)*(rho/rhoe(it)-1.0)**2+
     *       Fm33*(rho/rhoe(it)-1.0)**3
      else
         emb=Fn(it)*(1.0-fnn(it)*log(rho/rhos(it)))*
     *       (rho/rhos(it))**fnn(it)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write out set file.                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeset
      character*80 outfile,outelem
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16)
      common /pass2/ ielement(16),amass(16),Fr(5000,16),
     *   rhor(5000,16),z2r(5000,16,16),ntypes,blat(16),
     *   nrho,drho,nr,dr,rc,outfile,outelem
      character*80 struc
      struc='fcc'
      outfile = outfile(1:index(outfile,' ')-1)//'.set'
      open(unit=1,file=outfile)
      write(1,*)
      write(1,*)
      write(1,*)
      write(1,8)ntypes,outelem
8     format(i5,' ',a24)
      write(1,9)nrho,drho,nr,dr,rc
9     format(i5,e24.16,i5,2e24.16)
      do 10 i=1,ntypes
      write(1,11)ielement(i),amass(i),blat(i),struc
      write(1,12)(Fr(j,i),j=1,nrho)
      write(1,12)(rhor(j,i),j=1,nr)
10    continue
11    format(i5,2g15.5,a8)
12    format(5e24.16)
      do 13 i1=1,ntypes
      do 13 i2=1,i1
      write(1,12)(z2r(i,i1,i2),i=1,nr)
13    continue
      close(1)
      return
      end
