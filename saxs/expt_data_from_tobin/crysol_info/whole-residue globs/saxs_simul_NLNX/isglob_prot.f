cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This program calculates protein globbic correction factors
c     from a set of protein structures (PDB or CNS/XPLOR format)                       
c     Alex Grishaev, NIH, 09/19/02-
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program isglob_prot

      parameter(npdbmax  = 110000)
      parameter(nresmax  =   5000)
      parameter(nglobmax =   3000)

      parameter(ns = 27) ! set to the number of experimentaldata points

      double precision a(5,4),b(5,4),c(5),s,pi,wi,wj,dhist
      double precision ds,bres(5,40),isat(ns),isglob(ns)
      double precision rij,fact,atvj(5),rm,r0,zi(5),store
      double precision rho0,resol,iscurr(ns),gs(ns),store1,store2
      double precision fsglob(50,ns),chi,rfact(2000),time1,time2,time3
      double precision zweight(nresmax,40),fact1,fact2
      double precision blim,fsat(5,ns),gsat(5,ns),pofr(1000)
      double precision rpdb(nresmax,40,3),bpdb(nresmax,40),dist
      double precision rglob(nglobmax,3),bglob(nglobmax,30)
      double precision fsglob2(50,50,ns),fsat2(5,5,ns),israt(ns,2)
      double precision si(ns),is0(ns),dis0(ns)
 
      integer i,j,k,l,nfile,ifile,npdb,nres(npdbmax),istore,count
      integer datres(5,40),idmin,idmax,nskip,idg1,idg2,ihist
      integer idtype(nresmax,40),ic,nave,globflag(nresmax,40)
      integer nrescurr,ndres(5),icurr,nacc(30),solflag,bflag
      integer natdat(20),totres,natpdb(nresmax),resbreak(nresmax)
      integer nglob,natglob(100),totglob,idat(nglobmax,30),id0(30)
      integer idres(nglobmax,30),idglob(nglobmax),naver(50)
      integer ida1,ida2,it1,it2,idg,idr1,idr2

      character fname*4,fnamein*8,atscat(5)*1,atom(npdbmax)*4,flag*1
      character res(npdbmax)*3,atres(5,40)*4,gname*13
      character rescurr*3,temp*4,chain*1,pdbid(100)*4
      character resdat(20)*3,atdat(20,40)*4,respdb(nresmax)*3
      character atpdb(nresmax,40)*4,atglob(100,30)*4,resglob(100)*3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi = 3.1415927d0
      ds = 0.01d0

c     average (rm) and effective (r0) atomic radii (A)...
      rm = 1.616d0
	r0 = 1.0d0*rm

c     solvent density (e/A^3)...
      rho0 = 0.334d0

c     solvent contrast correction flag (1=yes)
      solflag = 1

c     Do you want to use the experimental B-factors (1=yes)?
      bflag = 0
      nave  = 0

c     read database of residues and their atoms...
      open(unit=10,file='res_cns.inp',status='old')
      rewind 10
      do i=1,20
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
      do i=1,5
       read(10,10)atscat(i),(a(i,j),j=1,4),(b(i,j),j=1,4),c(i)
10     format(a1,x,9(1x,f8.4))
      enddo
      close(unit=10,status='keep')

c     read the atomic displaced solvent volumes and atomic charges...
      open(unit=10,file='Vj.inp',status='old')
      rewind 10
      read(10,*)
      read(10,*)
      do i=1,5
       read(10,15)atvj(i),zi(i)
15     format(2x,f7.3,1x,f4.1)
      enddo
      close(unit=10,status='keep')
          
c     read the definitions of globs...
      open(unit=10,file='glob_res_cns.inp',status='old')
      rewind 10
      read(10,*)nglob
      do i=1,nglob
       read(10,24)resglob(i),natglob(i),(atglob(i,j),j=1,natglob(i))
24     format(a3,1x,i2,1x,50(1x,a4))
      enddo
      close(unit=10,status='keep')
 
c     read the experimental scattering vector values...
      open(unit=10,file='nisanxd_x8.dat',status='old')
      rewind 10
      read(10,*)
      do i=1,ns
       read(10,*)si(i)
c      si(i) = 0.1d0*dble(i)
      enddo
      close(unit=10,status='keep')

c     overall Gaussian expansion factor...
      do i=1,ns
       s = si(i)
       gs(i) = (r0/rm)**3*
     & dexp(-1.0d0*(4.0d0*pi/3.0d0)**(1.5d0)*pi*s**2*(r0**2-rm**2)/
     & (1.0d0))
      enddo

c     atomic scattering profiles...
      do i=1,5
       do j=1,ns
        s = si(j)
        fsat(i,j) = 0.0d0
        gsat(i,j) = 0.0d0
        do k=1,4
         fsat(i,j) = fsat(i,j)+
     &   a(i,k)*dexp(-1.0d0*b(i,k)*s**2/(4.0d0*pi)**2)
        enddo
        fsat(i,j) = fsat(i,j)+c(i)
c       dummy solvent scattering...
        gsat(i,j) = gs(j)*atvj(i)*
     &  dexp(-1.0d0*s**2*(atvj(i)**(2.0d0/3.0d0))/(4.0d0*pi))
        if(solflag.eq.1)then
         fsat(i,j) = fsat(i,j)-rho0*gsat(i,j)
        endif
       enddo
      enddo

      do i=1,5
       do j=i,5
        do k=1,ns
         fsat2(i,j,k) = fsat(i,k)*fsat(j,k)
         fsat2(j,i,k) = fsat2(i,j,k)
        enddo
       enddo
      enddo

c      read the globbic scattering profiles...
c      this is the output of the fsglob program !
      open(unit=10,file='fsglob.out',status='old')
      rewind 10
      read(10,*) ! skip the header line
      read(10,*) ! skip the q=0 line
      do i=1,ns
       read(10,26)(fsglob(j,i),j=1,nglob)
26     format(15x,40(1x,f9.6))
      enddo
      close(unit=10,status='keep')

      do i=1,nglob
       do j=i,nglob
        do k=1,ns
         fsglob2(i,j,k) = fsglob(i,k)*fsglob(j,k)
         fsglob2(j,i,k) = fsglob2(i,j,k)
        enddo
       enddo
      enddo

      ic = 0
      do i=1,ns
       israt (i,1)   = 0.0d0
       israt (i,2)   = 0.0d0
      enddo

c     read the pdb files...
      open(unit=10,file='files.inp',status='old')
      rewind 10
      read(10,*)nfile

      do ifile=1,nfile
       read(10,'(a4,a1)')fname,chain
       fnamein = fname//'.pdb'
       pdbid(ifile) = fname

c      initalize scattering intensity profiles and R-factors...
       do i=1,ns
        isat  (i) = 0.0d0
        isglob(i) = 0.0d0
       enddo

       if(mod(ifile,100).eq.0.0) write(6,*)nave,' files processed'

       call pdbread(npdbmax,nresmax,fnamein,chain,resdat,atdat,natdat,
     &             npdb,totres,respdb,natpdb,atpdb,rpdb,bpdb,resbreak)
       write(6,*)ifile,'   ',fname,' ',chain,' : ',totres,' residues '

       nave = nave+1
       ic = ic+1

c      assign nuclear charges...
       do i=1,totres
        do j=1,natpdb(i)
         globflag(i,j) = 0
         idtype  (i,j) = 0
         zweight (i,j) = 0.0d0
         temp = atpdb(i,j)
         do k=1,5
          if(temp(2:2).eq.atscat(k))then
            idtype(i,j) = k
           zweight(i,j) = zi(k)
          endif
         enddo
         if((idtype(i,j).eq.0).and.(temp(1:1).eq.'H'))then
          idtype(i,j) = 4
          zweight(i,j) = zi(4)       
         endif
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
           istore = min(1000,istore*id0(k))
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

       count = 0
       do i=1,totglob
        idq = idglob(i)
        do j=1,natglob(idq)
         idr1 = idres(i,j)
         ida1 = idat (i,j)
         if(globflag(idr1,ida1).eq.0)then
          globflag(idr1,ida1) = 1
          count = count+1
         else
           write(6,*)'multiple glob assignments',i,j,idq,idr1,ida1
           pause
         endif
        enddo
       enddo

       write(6,*)'number of residues = ',totres
       count = 0
       do i=1,totres
        count = count+natpdb(i)
       enddo
       write(6,*)'number of atoms    = ',count


c      Z-weighted globbic centers of mass...
       do i=1,totglob
        idg = idglob(i)
        rglob(i,1) = 0.0d0
        rglob(i,2) = 0.0d0
        rglob(i,3) = 0.0d0
        store = 0.0d0
        do j=1,natglob(idg)
         idr1 = idres(i,j)
         ida1 = idat (i,j)
         store = store+zweight(idr1,ida1)
        enddo
        do j=1,natglob(idg)
         idr1 = idres(i,j)
         ida1 = idat (i,j)
         zweight(idr1,ida1) = zweight(idr1,ida1)/store
         rglob(i,1) = rglob(i,1)+rpdb(idr1,ida1,1)*zweight(idr1,ida1)
         rglob(i,2) = rglob(i,2)+rpdb(idr1,ida1,2)*zweight(idr1,ida1)
         rglob(i,3) = rglob(i,3)+rpdb(idr1,ida1,3)*zweight(idr1,ida1)
        enddo
       enddo

       time1 = secnds(0.0)
c      globbic scattering curve...

       write(6,*)'number of globs    = ',totglob
       count = 0
       do i=1,totglob
        count = count+natglob(idglob(i))
       enddo
       write(6,*)'number of atoms in globs    = ',count

       do i=1,totglob
        idg1 = idglob(i)
        do k=1,ns
         isglob(k) = isglob(k)+fsglob2(idg1,idg1,k)
        enddo
        do j=i+1,totglob
         idg2 = idglob(j)
         dist = dsqrt((rglob(i,1)-rglob(j,1))**2+
     &               (rglob(i,2)-rglob(j,2))**2+
     &               (rglob(i,3)-rglob(j,3))**2)
         do k=1,ns
          s = si(k)
          if(s.ne.0.0d0)then
           fact = dsin(s*dist)/(s*dist)
          else
           fact = 1.0d0
          endif
          isglob(k) = isglob(k)+2.0d0*fsglob2(idg1,idg2,k)*fact
         enddo
        enddo
       enddo

       write(6,*)totglob,'-glob scattering curve done for ',fname
       time2 = secnds(0.0)
       write(6,*)'elapsed time for globbic calculation = ',time2-time1

c      exact atomic scattering curve...
       write(6,*)'calculating atomic scattering curve...'
       count = 0
       do i=1,totres
        write(6,*)i
        do j=1,natpdb(i)
         ida1 = idtype(i,j)
         if(globflag(i,j).ne.1)then
          write(6,*)i,' ',respdb(i),' ',atpdb(i,j)
          goto 94
         endif
         count = count+1
         do m=1,ns
          isat(m) = isat(m)+fsat2(ida1,ida1,m)
         enddo
         do l=j+1,natpdb(i)
          if(globflag(i,l).ne.1) goto 90
          ida2 = idtype(i,l)
          dist = dsqrt((rpdb(i,j,1)-rpdb(i,l,1))**2+
     &                 (rpdb(i,j,2)-rpdb(i,l,2))**2+
     &                 (rpdb(i,j,3)-rpdb(i,l,3))**2)
          do m=1,ns
           s = si(m)
           if(s.ne.0.0d0)then
            fact = dsin(s*dist)/(s*dist)
           else
            fact = 1.0d0
           endif
           isat(m) = isat(m)+2.0d0*fsat2(ida1,ida2,m)*fact
          enddo
90        continue
         enddo
         do k=i+1,totres
          do l=1,natpdb(k)
           if(globflag(k,l).ne.1) goto 92
           ida2 = idtype(k,l)
           dist = dsqrt((rpdb(i,j,1)-rpdb(k,l,1))**2+
     &                  (rpdb(i,j,2)-rpdb(k,l,2))**2+
     &                  (rpdb(i,j,3)-rpdb(k,l,3))**2)
           do m=1,ns
            s = si(m)
            if(s.ne.0.0d0)then
            fact = dsin(s*dist)/(s*dist)
            else
            fact = 1.0d0
            endif
            isat(m) = isat(m)+2.0d0*fsat2(ida1,ida2,m)*fact
           enddo           
92         continue
          enddo
         enddo 
94       continue
        enddo
       enddo
       time3 = secnds(0.0)
       write(6,*)'elapsed time for atomic calculation = ',time3-time2
       write(6,*)istore,count

c      fit globbic curve to the atomic curve...
       fact1 = 0.0d0
       fact2 = 0.0d0
       do i=1,ns
        fact1 = fact1+isat  (i)*isglob(i)
        fact2 = fact2+isglob(i)*isglob(i)
       enddo
       fact = fact1/fact2

       rfact(ifile) = 0.0d0
       fact1 = 0.0d0
       fact2 = 0.0d0
       do i=1,ns
        fact1 = fact1+(isglob(i)-isat(i))**2
        fact2 = fact2+isat(i)**2
       enddo
       rfact(ifile) = dsqrt(fact1/fact2)
       write(6,*)'globic/atomic R-factor = ',rfact(ifile)

c      intermediate output - q, I_all_atom(q), I_glob(q)
       open(unit=11,file='compl_at_glob1.out',status='unknown')
       rewind 11
       do i=1,ns
        write(11,'(3(e13.6,1x))')si(i),isat(i),isglob(i)
       enddo
       close(unit=11,status='keep')
       pause

       do i=1,ns
        fact = isglob(i)/isat(i)
        israt(i,1) = israt(i,1)+fact
        israt(i,2) = israt(i,2)+fact**2
       enddo

      enddo
      close(unit=10,status='keep')

c     THE MAIN OUTPUT - GLOBBIC CORRECTION CURVE!
c     the format is q, average correction, rms correction over the pdb files
      open(unit=11,file='compl_corr1.out',status='unknown')
      rewind 11
      do i=1,ns
       fact1 = israt(i,1)/dble(ic)
       fact2 = israt(i,2)/dble(ic)
       fact2 = fact2-fact1**2
       if(fact2.lt.0.0d0)then
        fact2 = 0.0d0
       endif
       fact2 = dsqrt(fact2)
       israt(i,1) = fact1
       israt(i,2) = fact2
       write(11,'(3(e13.6,1x))')si(i),fact1,fact2
      enddo
      close(unit=11,status='keep')


      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pdbread(npdbmax,nresmax,fnamein,chain,resdat,atdat,
     & natdat,npdb,totres,respdb,natpdb,atpdb,rpdb,bpdb,resbreak)

      integer i,j,k,l,m,n,npdb,npdbmax,count,flag,natdat(20)
      integer numb1,numb2,totres,istore
      double precision dist
      character resdat(20)*3,atdat(20,40)*4
      character fnamein*8,chain*1,text1*4,text2*4,text3*3,fl1*1,fl2*1
      character text4*1

      integer nres(npdbmax),resid(npdbmax),irep(npdbmax),break(npdbmax)
      double precision r(npdbmax,3),bfact(npdbmax)
      character atom(npdbmax)*4,res(npdbmax)*3,flid(npdbmax)*1

      integer natpdb(nresmax),resbreak(nresmax),maxoff,curres
      double precision rpdb(nresmax,40,3),bpdb(nresmax,40)
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
c      only "ATOM" fields with the correct chain identifier...
       if(text4.eq.chain)then
        if(text1.eq.'ATOM')then
c        only "A" or the primary conformer...
         if((fl1.eq.' ').or.(fl1.eq.'A').or.(fl1.eq.'a'))then
c         only if it is one of 20 common residues...
          flag = 0
          do j=1,20
           if(text3.eq.resdat(j))flag = j
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
           flid(npdb) = fl2
           if((x1.eq.0.0d0).and.(x2.eq.0.0d0).and.(x3.eq.0.0d0))then
            write(6,*)i,'all coordinates = zero'
            pause
           endif
          endif
         endif
        endif
c       SELENOMETHIONIES...
        if((text1.eq.'HETA').and.(text3.eq.'MSE'))then
c        only "A" or the primary conformer...
         if((fl1.eq.' ').or.(fl1.eq.'A').or.(fl1.eq.'a'))then
c         only if it is one of 20 common residues...
          flag = 13
          npdb = npdb+1
          atom(npdb) = text2
          res(npdb)  = 'MET'
          nres(npdb) = num1
          r(npdb,1) = x1
          r(npdb,2) = x2
          r(npdb,3) = x3
          bfact(npdb) = x4
          resid(npdb) = flag
          flid(npdb) = fl2
          if((x1.eq.0.0d0).and.(x2.eq.0.0d0).and.(x3.eq.0.0d0))then
           write(6,*)i,'all coordinates = zero'
           pause
          endif
         endif
        endif
       endif
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
       
c      from 30...
       if(flag.eq.0)then
        if(flid(i).ne.' ')then
c        something like 30 -> 32B
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
c          if(flid(i).eq.'P')goto 80
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
c         PAUSE
         goto 80
        endif
       endif
       
c      from 32B,D...
       if(flag.eq.1)then
        
        if((nres(i).eq.curres).and.(flid(i).eq.aflag))goto 80
        
        if(nres(i).eq.curres)then
c        something like 32B -> 32D...
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
c         PAUSE
         goto 80
        endif
        
        if(nres(i).ne.curres)then
         
         if(flid(i).eq.' ')then
c         something like 32D -> 35...
          do j=i,npdb
           irep(j) = irep(j)+maxflag
          enddo
          write(6,*)'     ',i,' : corrected numbering from',
     &     nres(i),flid(i),' by ',maxflag
C          PAUSE
          flag = 0
          maxflag = 0
          curres = 0
          atflag = ' '
          goto 80
         endif
         
         if(flid(i).ne.' ')then
c         something like 32D -> 35C...
          do j=i,npdb
           irep(j) = irep(j)+maxflag
          enddo
          write(6,*)'     ',i,' : corrected numbering from',
     &     nres(i),flid(i),' by ',maxflag
c          PAUSE
          flag = 0
          maxflag = 0
          curres = 0
          atflag = ' '
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

c      correct for the apparent sequence skips that are not real  
c      and check for the hidden chain skips
c      based on the N-C peptide distance ...
cccccccccccccccccccccccccccccccccccccccc
      do i=1,npdb
       irep(i) = 0
       break(i) = 0
      enddo
      do i=2,npdb
       if(nres(i-1).ne.nres(i))then
        
c       get the C_i-1/N_i distance
        numb1 = 0
        numb2 = 0
        do j=max(1,i-40),i-1
         if((atom(j).eq.' C  ').and.(nres(j).eq.nres(i-1)))numb1 = j
        enddo
        do j=i,min(npdb,i+40)
         if((atom(j).eq.' N  ').and.(nres(j).eq.nres(i)))numb2 = j
        enddo
        if(numb1.eq.0)then
         write(6,*)'C not found for residue ',nres(i-1),res(i-1)
         pause
        endif
        if(numb2.eq.0)then
         write(6,*)'N not found for residue ',nres(i),res(i)
         pause
        endif
        dist = dsqrt((r(numb1,1)-r(numb2,1))**2+
     &              (r(numb1,2)-r(numb2,2))**2+
     &              (r(numb1,3)-r(numb2,3))**2)
        
c       real chain break...
        if(dabs(dist-1.3d0).gt.4.0d-1)then
         write(6,100)'unusual C-N distance between ',
     &    res(i-1),nres(i-1),' and ',res(i),nres(i),' : ',dist,' A'
100      format(15x,a29,a3,i4,a5,a3,i4,a3,f6.2,a2)

         if(nres(i)-nres(i-1).eq.1)then
          write(6,*)'ATTN: potentially unknown seq. break @ ',nres(i-1)
          write(6,*)r(numb1,1),r(numb1,2),r(numb1,3)
          write(6,*)r(numb2,1),r(numb2,2),r(numb2,3)
           PAUSE
         endif

         if(nres(i)-nres(i-1).gt.1)then
          write(6,120)'interpreted as a real chain break'
120       format(15x,a33)
          do j=max(1,i-40),i-1
           if(nres(j).eq.nres(i-1))break(j) = 1
          enddo
         endif

         if(nres(i)-nres(i-1).lt.1)then
          write(6,*)'CHECK IT: chain break, sequence down-shift'
            write(6,*)nres(i-1),nres(i)
          do j=max(1,i-40),i-1
           if(nres(j).eq.nres(i-1))break(j) = 1
          enddo
C          PAUSE
         endif

        endif
        
c       no chain break...
        if(dabs(dist-1.3d0).le.1.0d-1)then
         if(nres(i)-nres(i-1).gt.1)then
          write(6,*)'            ','apparent up-skip not real',
     &     nres(i-1),res(i-1),nres(i),res(i)
          do j=i,npdb
           irep(j) = irep(j)-(nres(i)-nres(i-1))+1
          enddo
         endif
         if(nres(i)-nres(i-1).lt.1)then
          write(6,*)'CHECK IT: no chain break, sequence down-shift'
            write(6,*)nres(i-1),nres(i)
          do j=i,npdb
           irep(j) = irep(j)-(nres(i)-nres(i-1))+1
          enddo
C          PAUSE
         endif
        endif
        
       endif
      enddo
      
      do i=1,npdb
       nres(i) = nres(i)+irep(i)
      enddo

c     nisa-specific A-B
      do i=2,npdb
       if(nres(i).lt.nres(i-1))then
        do j=i,npdb
         nres(j) = nres(j)+566
        enddo
        goto 130
       endif
      enddo
130   continue

c     complex-specific B-C
      do i=2,npdb
       if(nres(i).lt.nres(i-1))then
        do j=i,npdb
         nres(j) = nres(j)+566
        enddo
        goto 131
       endif
      enddo
131   continue

c     complex-specific C-D
      do i=2,npdb
       if(nres(i).lt.nres(i-1))then
c       write(6,*)i,nres(i-1); pause
        do j=i,npdb
         nres(j) = nres(j)+261
        enddo
        goto 132
       endif
      enddo
132   continue




c      correct for the "negative N-terminal" residue numbering...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(nres(1).ne.1)then
       count = nres(1)
c       write(6,*)'corrected N-terminal numbering : ',nres(1)
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
         write(6,*)fnamein,i
         write(6,*)i,' ',nres(i),' ',res(i)
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
       bpdb(nres(i),natpdb(nres(i))) = bfact(i)

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