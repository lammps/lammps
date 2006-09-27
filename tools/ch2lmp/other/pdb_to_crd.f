c Reads PDB file, writes out charmm file
c Uses a temp file
c PDB format
c text IATOM  TYPE  RES  IRES      X  Y  Z    W
c  A6   I5  2X A4   A4    I5  4X     3F8.3 6X F6.2
c charmm format
c ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
c   I5    I5  1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5
c
c
        character*80 infile,outfile,line
        character*4 str1,type,res,code,segid,resid,residold,resold
        character*1 chain
        logical loxt(1000)
        write (6,*) 'Give input PDB files, output will be .crd'
1       read (5,'(a)') infile
        i=1
2       i=i+1
        if (infile(i:i).eq.' ') then
         outfile=infile(1:i-1)//'.crd'
        else
         goto 2
        endif
        open (unit=11, file=infile, status='old')
        open (unit=12, file='temppdb', status='unknown')
        open (unit=13, file=outfile, status='new')
        write (13,'(a80)') '* converted from '//infile
        write (13,'(a)') '*'
        do 4 i=1,1000
4       loxt(i)=.false.
        nss=0
        ires=0
        iat=0
        residold='    '
        resold='    '
        do 100 i=1,100000
        read (11,'(a80)',end=1000) line
        read (unit=line,fmt=500) str1
        if (str1.eq.'SSBO') then
          nss=nss+1
          goto 100
        else if (str1.eq.'ATOM') then
        iat= iat+1
        read (unit=line,fmt=500) str1,iatom,type,res,chain,resid,
     &    x,y,z,a,w,code
500     format (a4,2x,i5,1x,a4,1x,a4,a1,a4,4x,3f8.3,2f6.2,6x,a4)
          if ((resid.ne.residold).or.(res.ne.resold)) ires=ires+1
          residold=resid
          resold= res
          if (chain.ne.'    ') then
            segid=chain//code
          elseif (code.ne.'    ') then
            segid=code
          else
            segid='MAIN'
          endif
          if (type.eq.'CD1 ') then
             if (res.eq.'ILE ') type='CD  '
          elseif (type.eq.'OCT1') then
             type='OT1 '
          elseif (type.eq.'OCT2') then
             type='OT2 '
          elseif (type.eq.'OXT ') then
             type='OT2 '
             loxt(ires)=.true.
          endif
c fluch resid left
5         if (resid(1:1).eq.' ') then
            resid=resid(2:4)//' '
            goto 5
          endif
6         if (type(1:1).eq.' ') then
            type=type(2:4)//' '
            goto 6
          endif
          write (12,600) iat,ires,res,type,x,y,z,segid,resid,w
600     format (I5,I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,a4,F10.5)
        else
          goto 100
        endif
100     continue
1000    write (6,*) 'Disulfide bonds', nss
        nres=ires
        write (13,'(i5)') iat
        close (unit=12)
        open (unit=12,file='temppdb',status='old')
        do 200 i=1,100000
        read (12,'(a80)',end=2000) line
        read (unit=line,fmt=600) iatom,ires,res,type,x,y,z,segid,resid,w
        if (loxt(ires).and.(type.eq.'O   ')) type='OT1 '
        write (13,600) iatom,ires,res,type,x,y,z,segid,resid,w
200     continue
2000    close (unit=11)
        close (unit=12)
        close (unit=13)
        goto 1
        end

