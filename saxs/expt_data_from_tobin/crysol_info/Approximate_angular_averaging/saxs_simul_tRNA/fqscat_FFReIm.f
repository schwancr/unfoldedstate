cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This program calculates globbic form-factors from 
c     the input set of protein structures.
c     CNS/XPLOR format for the pdb file is assumed.
c     Otherwise, change res.inp into res_cns.inp 
c     and glob.inp into glob_cns.inp.                                    
c     Alex Grishaev, NIH, 09/19/02 - 2006
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program Fqscat_FFReIm

      parameter(npdbmax  = 110000)
      parameter(nresmax  =   2000)
      parameter(nglobmax =  10000)

      parameter(nq = 32) ! set this to the numbers of experimental data points

      double precision a(6,4),b(6,4),c(6),q,pi
      double precision rij,fact,atvj(6),rm,r0
      double precision rho0,Iqcurr(nq+1),Gq(nq+1),store1,store2
      double precision Fqglob(200,nq+1),Fq2glob(200,nq+1)
      double precision Fqat(6,nq+1),Gqat(6,nq+1),Fqscat(13,nq+1)
      double precision rpdb(nresmax,40,3),dist
      double precision qi(nq+1),Iq0(nq+1),dIq0(nq+1),volglob

      integer i,j,k,l,nfile,ifile,npdb,nres(npdbmax),istore
      integer idtype(nresmax,40),nacc(20),solflag,naverscat(13)
      integer natdat(20),totres,natpdb(nresmax),resbreak(nresmax)
      integer nglob,natglob(200),totglob,idat(nglobmax,20),id0(20)
      integer idres(nglobmax,20),idglob(nglobmax),naver(200)
      integer ida1,ida2,it1,it2,idg,idr1,idr2,scatid(200)

      character fname*4,fnamein*8,atscat(6)*1,atom(npdbmax)*4,flag*1
      character res(npdbmax)*3,temp*4,chain*1
      character resdat(20)*3,atdat(20,40)*4,respdb(nresmax)*3
      character atpdb(nresmax,40)*4,atglob(200,20)*4,resglob(200)*3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi = 3.1415927d0

c     average (rm) and effective (r0) atomic radii (A)...
      rm = 1.616d0
      r0 = 1.616d0

c     solvent density (e/A^3)...
      rho0 = 0.334d0

c     solvent contrast contribution flag (1=yes)
      solflag = 1

      store1 = 0.0d0
      store2 = 0.0d0

c     read database of residues and their atoms...
      open(unit=10,file='res_rna_cns.inp',status='old')
      rewind 10
      do i=1,4
       read(10,'(4x,a3)')resdat(i)
       read(10,*)natdat(i)
       do j=1,natdat(i)
        read(10,'(12x,a4)')atdat(i,j)
       enddo
       read(10,*)
      enddo
      close(unit=10,status='keep')

c     read the atomic form factors...
      open(unit=10,file='formfact.inp',status='old')
      rewind 10
      read(10,*)
      read(10,*)
      do i=1,6
       read(10,10)atscat(i),(a(i,j),j=1,4),(b(i,j),j=1,4),c(i)
10     format(a1,x,9(1x,f8.4))
      enddo
      close(unit=10,status='keep')

c     read the atomic displaced solvent volumes...
      open(unit=10,file='Vj_protna.inp',status='old')
      rewind 10
      read(10,*)
      read(10,*)
      do i=1,6
       read(10,15)atvj(i)
15     format(2x,f7.3)
      enddo
      close(unit=10,status='keep')

c     initalize variables...
      do i=1,200
       naver(i) = 0
       do j=1,nq+1
        Fqglob(i,j) = 0.0d0
        Fq2glob(i,j) = 0.0d0
       enddo
      enddo
           
c     read the definitions of scatterers (globs for our purposes here)...
      open(unit=10,file='scatterers.inp',status='old')
      rewind 10
      read(10,*)nglob
      do i=1,nglob
       read(10,'(12x,a4,1x,a3,4x,i2,6x,i2,1x,20(a4))')
     & atglob(i,1),resglob(i),scatid(i),natglob(i),
     &(atglob(i,j),j=2,natglob(i))
24     format(a3,1x,i2,1x,20(1x,a4))
      enddo
      close(unit=10,status='keep')


      open(unit=10,file='tRNA_zeroC_x3_desm95x5.dat',status='old')
      rewind 10
      read(10,*)
      qi(1) = 0.0d0 ! the first entry is always going to be q=0.0
      do i=2,nq+1
       read(10,*)qi(i)
      enddo
      close(unit=10,status='keep')

c     overall Gaussian expansion factor...
      do i=1,nq+1
       q = qi(i)
       Gq(i) = (r0/rm)**3*
     & dexp(-1.0d0*(4.0d0*pi/3.0d0)**(1.5d0)*pi*q**2*(r0**2-rm**2)/
     & (1.0d0))
      enddo

c     calculate atomic scattering profiles from the 4-Gaussian parametrization...
      do i=1,6
       do j=1,nq+1
        q = qi(j)
        Fqat(i,j) = 0.0d0
        Gqat(i,j) = 0.0d0
        do k=1,4
         Fqat(i,j) = Fqat(i,j)+
     &   a(i,k)*dexp(-1.0d0*b(i,k)*q**2/(4.0d0*pi)**2)
        enddo
        Fqat(i,j) = Fqat(i,j)+c(i)
c       dummy solvent scattering...
        Gqat(i,j) = Gq(j)*atvj(i)*
     &  dexp(-1.0d0*q**2*(atvj(i)**(2.0d0/3.0d0))/(4.0d0*pi))
c       if(solflag.eq.1)then
c        fsat(i,j) = fsat(i,j)-rho0*gsat(i,j)
c       endif
       enddo
      enddo


c     read the pdb files...
      open(unit=10,file='files.inp',status='old')
      rewind 10
      read(10,*)nfile
      do ifile=1,nfile
       read(10,'(a4,a1)')fname,chain
       fnamein = fname//'.pdb'

       call pdbread(npdbmax,nresmax,fnamein,resdat,atdat,natdat,
     &      npdb,totres,respdb,natpdb,atpdb,rpdb,resbreak)
       write(6,*)ifile,'   ',fname,' ',chain,' : ',totres,' residues '

c      establish atom types...
       do i=1,totres
        do j=1,natpdb(i)
         idtype(i,j) = 0
         temp = atpdb(i,j)
         do k=1,6
          if(temp(2:2).eq.atscat(k)) idtype(i,j) = k
         enddo
         if(temp(1:1).eq.'H') idtype(i,j) = 4 ! set type to H if atom starts with H
        enddo
       enddo

c      split the structure into globs...
       totglob = 0
       do i=1,totres
        do idg=1,nglob
         if(respdb(i).eq.resglob(idg))then
          do k=1,natglob(idg)
           id0(k) = 0
          enddo
          do k=1,natpdb(i)
           do l=1,natglob(idg)
            if(atglob(idg,l).eq.atpdb(i,k))then
             id0(l) = k
            endif
           enddo
          enddo
          istore =1
          do k=1,natglob(idg)
           istore = istore*id0(k)
          enddo
          if(istore.ne.0)then
           totglob = totglob+1
           idglob(totglob) = idg
           do k=1,natglob(idg)
            idat (totglob,k) = id0(k)
            idres(totglob,k) = i
           enddo 
          endif
         endif
        enddo
       enddo

       store1 = store1+dble(totglob)/dble(totres)
       store2 = store2+(dble(totglob)/dble(totres))**2

c      calculate globbic form factors...
       do i=1,totglob
        idg = idglob(i)
        do j=1,nq+1
         Iqcurr(j) = 0.0d0
        enddo
        volglob = 0.0d0
        do k=1,natglob(idg)
         idr1 = idres(i,k)
         ida1 = idat (i,k)
         it1  = idtype(idr1,ida1)
         volglob = volglob+atvj(it1)
         do j=1,nq+1
          Iqcurr(j) = Iqcurr(j)+Fqat(it1,j)**2
         enddo
         do l=k+1,natglob(idg)
          idr2 = idres(i,l)
          ida2 = idat (i,l)
          it2 = idtype(idr2,ida2)
          dist = dsqrt((rpdb(idr1,ida1,1)-rpdb(idr2,ida2,1))**2+
     &                 (rpdb(idr1,ida1,2)-rpdb(idr2,ida2,2))**2+
     &                 (rpdb(idr1,ida1,3)-rpdb(idr2,ida2,3))**2)
          do j=1,nq+1
           q = qi(j)
           if(q.ne.0.0d0)then
            fact = dsin(q*dist)/(q*dist)
           else
            fact = 1.0d0
           endif
           Iqcurr(j) = Iqcurr(j)+2.0d0*Fqat(it1,j)*Fqat(it2,j)*fact
          enddo
         enddo
        enddo
        do j=1,nq+1
         q = qi(j)
         Iqcurr(j) = dsqrt(Iqcurr(j))
         Iqcurr(j) = Iqcurr(j)-rho0*Gq(j)*volglob*
     &   dexp(-q**2*(volglob**(2.0d0/3.0d0))/(4.0d0*pi))
         Fqglob(idg,j) = Fqglob(idg,j)+Iqcurr(j)
         Fq2glob(idg,j) = Fq2glob(idg,j)+Iqcurr(j)**2
        enddo
        naver(idg) = naver(idg)+1
       enddo

      enddo
      close(unit=10,status='keep')

      do i=1,nglob
       do j=1,nq+1
        if(naver(i).ne.0)then
         Fqglob(i,j) = Fqglob(i,j)/dble(naver(i))
         Fq2glob(i,j) = Fq2glob(i,j)/dble(naver(i))
        else
         Fqglob(i,j)  = 0.0d0
         Fq2glob(i,j) = 0.0d0
        endif
        Fq2glob(i,j) = Fq2glob(i,j)-Fqglob(i,j)**2
        if(Fq2glob(i,j).lt.0.0d0)then
         Fq2glob(i,j) = 0.0d0
        endif
        Fq2glob(i,j) = dsqrt(Fq2glob(i,j))
       enddo
      enddo

      store1 = store1/dble(nfile)
      store2 = store2/dble(nfile)
      store2 = dsqrt(store2-store1**2)
      write(6,*)store1,' +/- ',store2,' globs per residue'

      do i=1,nglob
       do j=1,nq+1
        Fqscat(scatid(i),j) = 
     &  Fqscat(scatid(i),j) + Fqglob(i,j)*dble(naver(i))
       enddo
       naverscat(scatid(i)) = naverscat(scatid(i))+naver(i)
      enddo
      do i=1,13
       do j=1,nq+1
        if(naverscat(i).ne.0)then
         Fqscat(i,j) = Fqscat(i,j)/dble(naverscat(i))
        else
         Fqscat(i,j) = 0.0d0
        endif
       enddo
      enddo

c     MAIN OUTPUT: average globbic form factors !
c     use these for cns/xplor and as input for isglob program
      open(unit=10,file='Fqglob_FFReIm.out',status='unknown')
      rewind 10
      write(10,1000)(resglob(i),i=1,nglob)
1000  format(18x,170(4x,a3,3x))
      do i=1,nq+1
       write(10,1010)qi(i),(Fqglob(j,i),j=1,nglob)
1010   format(es13.6,2x,170(1x,f9.6))
      enddo
      write(10,'(18x,170(4x,i3,3x))')(naver(i),i=1,nglob)
      write(10,'(18x,170(4x,i3,3x))')(i,i=1,nglob)
      close(unit=10,status='keep')

      open(unit=10,file='Fqscat_FFReIm.out',status='unknown')
      rewind 10
      write(10,'(23x,a)')
     &' C        CH       CH2       CH3         N        NH       NH2'//
     &'       NH3         O        OH         S        SH         P '
      do i=1,nq+1
       write(10,1010)qi(i),(Fqscat(j,i),j=1,13)
      enddo
      write(10,'(18x,170(4x,i3,3x))')(naverscat(i),i=1,13)
      write(10,'(18x,170(4x,i3,3x))')(i,i=1,13)
      close(unit=10,status='keep')


      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine pdbread(npdbmax,nresmax,fnamein,resdat,atdat,natdat,
     & npdb,totres,respdb,natpdb,atpdb,rpdb,resbreak)

	 integer i,j,k,l,m,n,npdb,npdbmax,count,flag,natdat(20)
	 integer numb1,numb2,totres,istore
	 double precision dist
	 character resdat(20)*3,atdat(20,40)*4
	 character fnamein*8,chain*1,text1*4,text2*4,text3*3,fl1*1,fl2*1
	 character text4*1

	 integer nres(npdbmax),resid(npdbmax),irep(npdbmax),break(npdbmax)
c      integer segid(npdbmax)
	 double precision r(npdbmax,3),bfact(npdbmax)
	 character atom(npdbmax)*4,res(npdbmax)*3,flid(npdbmax)*1

	 integer natpdb(nresmax),resbreak(nresmax),maxoff,curres
	 double precision rpdb(nresmax,40,3)
	 character atpdb(nresmax,40)*4,respdb(nresmax)*3,aflag*1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       totres = 0

c      first pass (raw read), form atom-based arrays...
cccccccccccccccccccccccccccccccccccccccccccccccccccc
	 open(unit=11,file=fnamein,status='old')
       npdb = 0
	 do i=1,npdbmax
	  read(11,60,err=65,end=70)
     &  text1,text2,fl1,text3,text4,num1,fl2,x1,x2,x3,x4
60      format(a4,8x,a4,a1,a3,1x,a1,i4,a1,3x,3(f8.3),6x,f6.2)
c       only "ATOM" fields with the correct chain identifier...
c	  if(text4.eq.chain)then
	   if(text1.eq.'ATOM')then
c         only "A" or the primary conformer...
	    if((fl1.eq.' ').or.(fl1.eq.'A').or.(fl1.eq.'a'))then
c          if its name and resid matches database...
           flag = 0
	     do j=1,4
	      if(text3.eq.resdat(j))then
             do k=1,natdat(j)
              if(text2.eq.atdat(j,k)) flag = j
             enddo
            endif
	     enddo
	     if(flag.ne.0)then
	      npdb = npdb+1
            atom(npdb) = text2
	      res(npdb)  = text3
	      nres(npdb) = num1
	      r(npdb,1) = x1
	      r(npdb,2) = x2
	      r(npdb,3) = x3
	      bfact(npdb) = x4
	      resid(npdb) = flag
c	      segid(npdb) = text4
	      flid(npdb) = fl2
	      if((x1.eq.0.0d0).and.(x2.eq.0.0d0).and.(x3.eq.0.0d0))then
	       write(6,*)i,'all coordinates = zero'
	       pause
	      endif
	     endif
	    endif
	   endif
c       endif
65      continue
       enddo
70     continue
       close(unit=11,status='keep')


c      correct for the "32B, 32D, 32C, etc." residue numbering...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 do i=1,npdb
	  irep(i) = 0
       enddo
	 flag = 0
	 maxoff = 0
	 aflag = ' '

	 do i=2,npdb
        
c       from 30...
	  if(flag.eq.0)then
         if(flid(i).ne.' ')then
c         something like 30 -> 32B
	    istore = 0
	    if(flid(i).eq.'A')istore =  1
	    if(flid(i).eq.'B')istore =  2
	    if(flid(i).eq.'C')istore =  3
	    if(flid(i).eq.'D')istore =  4
	    if(flid(i).eq.'E')istore =  5
	    if(flid(i).eq.'F')istore =  6
	    if(flid(i).eq.'G')istore =  7
	    if(flid(i).eq.'H')istore =  8
	    if(flid(i).eq.'I')istore =  9
	    if(flid(i).eq.'J')istore = 10
	    if(flid(i).eq.'K')istore = 11
	    if(istore.eq.0)then
c	     if(flid(i).eq.'P')goto 80
	     write(6,*)'unknown numbering offset identifier'
	     write(6,*)i,res(i),nres(i),flid(i)
	     pause
          endif
	    do j=i,npdb
	     if((nres(j).eq.nres(i)).and.(flid(j).eq.flid(i)))then
	      irep(j) = irep(j)+istore
	     endif
	    enddo
	    flag = 1
	    maxflag = max(istore,maxflag)
	    curres = nres(i)
	    aflag = flid(i)
	    write(6,*)'     ',i,' : corrected numbering for ',
     &    nres(i),flid(i),' by ',istore
c	    PAUSE
	    goto 80
         endif
	  endif
        
c       from 32B,D...
	  if(flag.eq.1)then
	   
         if((nres(i).eq.curres).and.(flid(i).eq.aflag))goto 80
         
         if(nres(i).eq.curres)then
c         something like 32B -> 32D...
          istore = 0
	    if(flid(i).eq.'A')istore =  1
	    if(flid(i).eq.'B')istore =  2
	    if(flid(i).eq.'C')istore =  3
	    if(flid(i).eq.'D')istore =  4
	    if(flid(i).eq.'E')istore =  5
	    if(flid(i).eq.'F')istore =  6
	    if(flid(i).eq.'G')istore =  7
	    if(flid(i).eq.'H')istore =  8
	    if(flid(i).eq.'I')istore =  9
	    if(flid(i).eq.'J')istore = 10
	    if(flid(i).eq.'K')istore = 11
	    if(istore.eq.0)then
	     write(6,*)'unknown numbering offset identifier'
	     write(6,*)i,res(i),nres(i),flid(i)
	     pause
          endif
	    do j=i,npdb
	     if((nres(j).eq.nres(i)).and.(flid(j).eq.flid(i)))then
	      irep(j) = irep(j)+istore
	     endif
	    enddo
	    flag = 1
	    maxflag = max(istore,maxflag)
	    curres = nres(i)
	    aflag = flid(i)
	    write(6,*)'     ',i,' : corrected numbering for ',
     &    nres(i),flid(i),' by ',istore
c	    PAUSE
	    goto 80
         endif
         
	   if(nres(i).ne.curres)then
	    
	    if(flid(i).eq.' ')then
c          something like 32D -> 35...
           do j=i,npdb
	      irep(j) = irep(j)+maxflag
	     enddo
	     write(6,*)'     ',i,' : corrected numbering from',
     &     nres(i),flid(i),' by ',maxflag
C	     PAUSE
	     flag = 0
	     maxflag = 0
	     curres = 0
	     aflag = ' '
	     goto 80
	    endif
          
	    if(flid(i).ne.' ')then
c          something like 32D -> 35C...
           do j=i,npdb
	      irep(j) = irep(j)+maxflag
	     enddo
	     write(6,*)'     ',i,' : corrected numbering from',
     &     nres(i),flid(i),' by ',maxflag
c	     PAUSE
	     flag = 0
	     maxflag = 0
	     curres = 0
	     aflag = ' '
	     istore = 0
	     if(flid(i).eq.'A')istore =  1
	     if(flid(i).eq.'B')istore =  2
	     if(flid(i).eq.'C')istore =  3
	     if(flid(i).eq.'D')istore =  4
	     if(flid(i).eq.'E')istore =  5
	     if(flid(i).eq.'F')istore =  6
	     if(flid(i).eq.'G')istore =  7
	     if(flid(i).eq.'H')istore =  8
	     if(flid(i).eq.'I')istore =  9
	     if(flid(i).eq.'J')istore = 10
	     if(flid(i).eq.'K')istore = 11
	     if(istore.eq.0)then
	      write(6,*)'unknown numbering offset identifier'
	      write(6,*)i,res(i),nres(i),flid(i)
	      pause
           endif
	     do j=i,npdb
	      if((nres(j).eq.nres(i)).and.(flid(j).eq.flid(i)))then
	       irep(j) = irep(j)+istore
	      endif
	     enddo
	     flag = 1
	     maxflag = max(istore,maxflag)
	     curres = nres(i)
	     aflag = flid(i)
	     write(6,*)'     ',i,' : corrected numbering for ',
     &     nres(i),flid(i),' by ',istore
	     goto 80
          endif
          
         endif
         
        endif

80      continue
       enddo
       
	 do i=1,npdb
	  nres(i) = nres(i)+irep(i)
	 enddo

c      correct for the "negative N-terminal" residue numbering...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 if(nres(1).ne.1)then
	  count = nres(1)
c	  write(6,*)'corrected N-terminal numbering : ',nres(1)
	  do i=1,npdb
	   nres(i) = nres(i)-count+1
	  enddo
	 endif

c      form residue-based arrays...
cccccccccccccccccccccccccccccccc
	 do i=1,nresmax
	  natpdb(i) = 0
	 enddo
       do i=1,npdb
	  natpdb(nres(i)) = natpdb(nres(i))+1
	  if(natpdb(nres(i)).gt.natdat(resid(i)))then
	   if((nres(i).ne.1).and.(nres(i).ne.nres(npdb)))then
	    write(6,*)'maximum number of atoms exceeded'
	    write(6,*)natpdb(nres(i)),' ',nres(i),' ',res(i)
	    write(6,*)natdat(resid(i)),resid(i)
 	    do j=1,natpdb(nres(i))
 	     write(6,*)j,atpdb(nres(i),j),rpdb(nres(i),j,1)
 	    enddo
	    pause
	   endif
	  endif
	  atpdb(nres(i),natpdb(nres(i))) = atom(i)
	  respdb(nres(i)) = res(i)
	  resbreak(nres(i)) = break(i)

	  do j=1,3
	   rpdb(nres(i),natpdb(nres(i)),j) = r(i,j)
	  enddo
	 enddo
	 totres = nres(npdb)

	 if(totres.eq.0)then
	  write(6,*)'no residues found'
	  stop
	 endif

       return
	 end