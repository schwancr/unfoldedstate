cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This program calculates globbic form-factors from 
c     the input set of protein structures.
c     CNS/XPLOR format for the pdb file is assumed.
c     Otherwise, change res.inp into res_cns.inp 
c     and glob.inp into glob_cns.inp.                                    
c     Alex Grishaev, NIH, 09/19/02 - 2006
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program fsglob_prot

      parameter(npdbmax  = 110000)
      parameter(nresmax  =   2000)
      parameter(nglobmax =  10000)

      parameter(ns = 50) ! set this to the numbers of experimental data points

      double precision a(5,4),b(5,4),c(5),s,pi
      double precision ds,rres(5,40,3),bres(5,40)
      double precision rij,fact,atvj(5),rm,r0,fsign(200)
      double precision rho0,resol,iscurr(ns+1),gs(ns+1),store1,store2
      double precision fsglob(50,ns+1),fs2glob(50,ns+1)
      double precision blim,fsat(5,ns+1),gsat(5,ns+1)
      double precision rpdb(nresmax,40,3),bpdb(nresmax,40),dist
      double precision rglob(nglobmax,20,3),bglob(nglobmax,20)
      double precision dave(50,20,20),drms(50,20,20),ndave(50,20,20)
      double precision si(ns+1),is0(ns+1),dis0(ns+1)

      integer i,j,k,l,nfile,ifile,npdb,nres(npdbmax),istore
      integer datres(5,40),idmin,idmax,idtype(nresmax,40)
      integer nrescurr,ndres(5),icurr,nacc(20),solflag,bflag
      integer natdat(20),totres,natpdb(nresmax),resbreak(nresmax)
      integer nglob,natglob(100),totglob,idat(nglobmax,20),id0(20)
      integer idres(nglobmax,20),idglob(nglobmax),naver(50)
      integer ida1,ida2,it1,it2,idg,idr1,idr2

      character fname*4,fnamein*8,atscat(5)*1,atom(npdbmax)*4,flag*1
      character res(npdbmax)*3,atres(5,40)*4
      character dumat(20,5,40)*4,rescurr*3,temp*4,chain*1
      character resdat(20)*3,atdat(20,40)*4,respdb(nresmax)*3
      character atpdb(nresmax,40)*4,atglob(100,20)*4,resglob(100)*3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi = 3.1415927d0
      ds = 0.01d0

c     average (rm) and effective (r0) atomic radii (A)...
      rm = 1.616d0
      r0 = 1.616d0

c     solvent density (e/A^3)...
      rho0 = 0.334d0

c     solvent contrast contribution flag (1=yes)
      solflag = 1

c     Do you want to use the experimental B-factors (1=yes)?
      bflag = 0

c     B factor threshold...
      blim = 40.0d0

      store1 = 0.0d0
      store2 = 0.0d0

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

c     read the atomic displaced solvent volumes...
      open(unit=10,file='Vj.inp',status='old')
      rewind 10
      read(10,*)
      read(10,*)
      do i=1,5
       read(10,15)atvj(i)
15     format(2x,f7.3)
      enddo
      close(unit=10,status='keep')

c     initalize variables...
      do i=1,nglob
       naver(i) = 0
       do j=1,ns+1
        fsglob (i,j) = 0.0d0
        fs2glob(i,j) = 0.0d0
       enddo
      enddo
           
c     read the definitions of globs...
      open(unit=10,file='glob_cns.inp',status='old')
      rewind 10
      read(10,*)nglob
      do i=1,nglob
       read(10,24)resglob(i),natglob(i),(atglob(i,j),j=1,natglob(i))
24     format(a3,1x,i2,1x,20(1x,a4))
      enddo
      close(unit=10,status='keep')

c     read the experimental data file...
c     INPUT YOUR DATA HERE!
      open(unit=10,file='TolR_rescale13_x3_desm_x6.dat',status='old')
      rewind 10
      read(10,*)
      si(1) = 0.0d0 ! the first entry is always going to be q=0.0
      do i=2,ns+1
       read(10,*)si(i)
c      si(i) = 0.1d0*dble(i-1)
      enddo
      close(unit=10,status='keep')

c     overall Gaussian expansion factor...
      do i=1,ns+1
       s = si(i)
       gs(i) = (r0/rm)**3*
     & dexp(-1.0d0*(4.0d0*pi/3.0d0)**(1.5d0)*pi*s**2*(r0**2-rm**2)/
     & (1.0d0))
      enddo

c     calculate atomic scattering profiles from the 4-Gaussian parametrization...
      do i=1,5
       do j=1,ns+1
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

      do i=1,nglob
       fsign(i) = 0.0d0
       do j=1,natglob(i)
        temp = atglob(i,j)
        l = 0
        do k=1,5
         if(temp(2:2).eq.atscat(k)) l = k
        enddo
        if((l.eq.0).and.(temp(1:1).eq.'H'))then
         l = 4 ! set type to H if atom starts with H
        endif
        if((l.eq.0).and.(temp(1:1).ne.'H'))then
         write(6,*)i,j,temp,'UNKNOWN SCATTERER'
         pause
        endif
        fsign(i) = fsign(i)+fsat(l,1)
       enddo
       fsign(i) = dsign(1.0d0,fsign(i))
      enddo

c     intialize distance statistics...
      do i=1,50
       do j=1,20
        do k=1,20
         dave (i,j,k) = 0.0d0
         drms (i,j,k) = 0.0d0
         ndave(i,j,k) = 0.0d0
        enddo
       enddo
      enddo

c     read the pdb files...
      open(unit=10,file='files.inp',status='old')
      rewind 10
      read(10,*)nfile
      do ifile=1,nfile
       read(10,'(a4,a1)')fname,chain
       fnamein = fname//'.pdb'

       call pdbread(npdbmax,nresmax,fnamein,chain,resdat,atdat,natdat,
     &      npdb,totres,respdb,natpdb,atpdb,rpdb,bpdb,resbreak)
       write(6,*)ifile,'   ',fname,' ',chain,' : ',totres,' residues '


c      establish atom types...
       do i=1,totres
        do j=1,natpdb(i)
         idtype(i,j) = 0
         temp = atpdb(i,j)
         do k=1,5
          if(temp(2:2).eq.atscat(k)) idtype(i,j) = k
         enddo
         if(temp(1:1).eq.'H') idtype(i,j) = 4 ! set type to H if atom starts with H
        enddo
       enddo

c      split the structure into globs...
       totglob = 0
       do i=1,totres
c       peptide H-N-C=O group (glob #1)...
        if((i.lt.totres).and.(respdb(i+1).ne.'PRO'))then
         idg = 1
         do k=1,natglob(idg)
          id0(k) = 0
         enddo
         do k=1,natpdb(i)
          do l=1,2
           if(atglob(idg,l).eq.atpdb(i,k))then
            if(bpdb(i,k).lt.blim)then
             id0(l) = k
            endif
           endif
          enddo
         enddo
         do k=1,natpdb(i+1)
          do l=3,4
           if(atglob(idg,l).eq.atpdb(i+1,k))then
            if(bpdb(i+1,k).lt.blim)then
             id0(l) = k
            endif
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
          do k=1,2
           idat (totglob,k) = id0(k)
           idres(totglob,k) = i
          enddo 
          do k=3,4
           idat (totglob,k) = id0(k)
           idres(totglob,k) = i+1
          enddo 
         endif
        endif
c       Proline peptide N-C=O group (glob #2)...
        if((i.lt.totres).and.(respdb(i+1).eq.'PRO'))then
         idg = 2
         do k=1,natglob(idg)
          id0(k) = 0
         enddo
         do k=1,natpdb(i)
          do l=1,2
           if(atglob(idg,l).eq.atpdb(i,k))then
            if(bpdb(i,k).lt.blim)then
             id0(l) = k
            endif
           endif
          enddo
         enddo
         do k=1,natpdb(i+1)
          do l=3,3
           if(atglob(idg,l).eq.atpdb(i+1,k))then
            if(bpdb(i+1,k).lt.blim)then
             id0(l) = k
            endif
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
          do k=1,2
           idat (totglob,k) = id0(k)
           idres(totglob,k) = i
          enddo 
          do k=3,3
           idat (totglob,k) = id0(k)
           idres(totglob,k) = i+1
          enddo 
         endif
        endif
c       N-terminal NH3+ group (glob #3)...
        if(i.eq.1)then
         idg = 3
         do k=1,natglob(idg)
          id0(k) = 0
         enddo
         do k=1,natpdb(i)
          do l=1,natglob(idg)
           if(atglob(idg,l).eq.atpdb(i,k))then
            if(bpdb(i,k).lt.blim)then
             id0(l) = k
             endif
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
c       C-terminal COO- group (glob #4)...
        if(i.eq.totres)then
         idg = 4
         do k=1,natglob(idg)
          id0(k) = 0
         enddo
         do k=1,natpdb(i)
          do l=1,natglob(idg)
           if(atglob(idg,l).eq.atpdb(i,k))then
            if(bpdb(i,k).lt.blim)then
             id0(l) = k
            endif
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
c       sidechain globs (#5 and up)...
        do idg=5,nglob
         if(respdb(i).eq.resglob(idg))then
          do k=1,natglob(idg)
           id0(k) = 0
          enddo
          do k=1,natpdb(i)
           do l=1,natglob(idg)
            if(atglob(idg,l).eq.atpdb(i,k))then
             if(bpdb(i,k).lt.blim)then
              id0(l) = k
             endif
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
        do j=1,ns+1
         iscurr(j) = 0.0d0
        enddo
        do k=1,natglob(idg)
         idr1 = idres(i,k)
         ida1 = idat (i,k)
         it1  = idtype(idr1,ida1)
         do j=1,ns+1
          iscurr(j) = iscurr(j)+fsat(it1,j)**2
         enddo
         do l=k+1,natglob(idg)
          idr2 = idres(i,l)
          ida2 = idat (i,l)
          it2 = idtype(idr2,ida2)
          dist = dsqrt((rpdb(idr1,ida1,1)-rpdb(idr2,ida2,1))**2+
     &                 (rpdb(idr1,ida1,2)-rpdb(idr2,ida2,2))**2+
     &                 (rpdb(idr1,ida1,3)-rpdb(idr2,ida2,3))**2)
          dave (idg,k,l) = dave (idg,k,l)+dist
          drms (idg,k,l) = drms (idg,k,l)+dist**2
          ndave(idg,k,l) = ndave(idg,k,l)+1.0d0
          do j=1,ns+1
           s = si(j)
           if(s.ne.0.0d0)then
            fact = dsin(s*dist)/(s*dist)
           else
            fact = 1.0d0
           endif
           iscurr(j) = iscurr(j)+2.0d0*fsat(it1,j)*fsat(it2,j)*fact
          enddo
         enddo
        enddo
        do j=1,ns+1
         iscurr(j) = dsqrt(iscurr(j))
         fsglob (idg,j) = fsglob (idg,j)+iscurr(j)
         fs2glob(idg,j) = fs2glob(idg,j)+iscurr(j)**2
        enddo
        naver(idg) = naver(idg)+1
       enddo

      enddo
      close(unit=10,status='keep')

      do i=1,nglob
       do j=1,natglob(i)
        do k=j+1,natglob(i)
         if(ndave(i,j,k).ne.0.0d0)then
          dave(i,j,k) = dave(i,j,k)/ndave(i,j,k)
          drms(i,j,k) = drms(i,j,k)/ndave(i,j,k)
          drms(i,j,k) = drms(i,j,k)-dave (i,j,k)**2
          if(drms(i,j,k).lt.0.0d0)then
           drms(i,j,k) = 0.0d0
          endif
          drms(i,j,k) = dsqrt(drms(i,j,k))
         endif
        enddo
       enddo
      enddo

      do i=1,nglob
       do j=1,ns+1
        if(naver(i).ne.0)then
         fsglob (i,j) = fsglob (i,j)/dble(naver(i))
         fs2glob(i,j) = fs2glob(i,j)/dble(naver(i))
        else
         fsglob(i,j)  = 0.0d0
         fs2glob(i,j) = 0.0d0
        endif
        fs2glob(i,j) = fs2glob(i,j)-fsglob(i,j)**2
        if(fs2glob(i,j).lt.0.0d0)then
         fs2glob(i,j) = 0.0d0
        endif
        fs2glob(i,j) = dsqrt(fs2glob(i,j))
       enddo
      enddo

      store1 = store1/dble(nfile)
      store2 = store2/dble(nfile)
      store2 = dsqrt(store2-store1**2)
      write(6,*)store1,' +/- ',store2,' globs per residue'

      do i=1,nglob
       do j=1,ns+1
        fsglob(i,j) = fsglob(i,j)*fsign(i)
       enddo
      enddo

c     MAIN OUTPUT: average globbic form factors !
c     use these for cns/xplor and as input for isglob program
      open(unit=10,file='fsglob_tolR_rescale13.dat',status='unknown')
      rewind 10
      write(10,1000)(resglob(i),i=1,nglob)
1000  format(18x,40(4x,a3,3x))
      do i=1,ns+1
       write(10,1010)si(i),(fsglob(j,i),j=1,nglob)
1010   format(e13.6,2x,40(1x,f9.6))
      enddo
      close(unit=10,status='keep')

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
c       only "ATOM" fields with the correct chain identifier...
        if(text4.eq.chain)then
         if(text1.eq.'ATOM')then
c         only "A" or the primary conformer...
          if((fl1.eq.' ').or.(fl1.eq.'A').or.(fl1.eq.'a'))then
c          only if it is one of 20 common residues...
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
c        SELENOMETHIONIES...
         if((text1.eq.'HETA').and.(text3.eq.'MSE'))then
c         only "A" or the primary conformer...
          if((fl1.eq.' ').or.(fl1.eq.'A').or.(fl1.eq.'a'))then
c          only if it is one of 20 common residues...
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
c           if(flid(i).eq.'P')goto 80
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
c          PAUSE
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
c          PAUSE
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
C           PAUSE
           flag = 0
           maxflag = 0
           curres = 0
           atflag = ' '
           goto 80
          endif
          
          if(flid(i).ne.' ')then
c          something like 32D -> 35C...
           do j=i,npdb
            irep(j) = irep(j)+maxflag
           enddo
           write(6,*)'     ',i,' : corrected numbering from',
     &     nres(i),flid(i),' by ',maxflag
c           PAUSE
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
         
c        get the C_i-1/N_i distance
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
     &                (r(numb1,2)-r(numb2,2))**2+
     &                (r(numb1,3)-r(numb2,3))**2)
         
c        real chain break...
         if(dabs(dist-1.3d0).gt.4.0d-1)then
          write(6,100)'unusual C-N distance between ',
     &    res(i-1),nres(i-1),' and ',res(i),nres(i),' : ',dist,' A'
100       format(15x,a29,a3,i4,a5,a3,i4,a3,f6.2,a2)

          if(nres(i)-nres(i-1).eq.1)then
           write(6,*)'ATTN: potentially unknown seq. break @ ',nres(i-1)
           write(6,*)r(numb1,1),r(numb1,2),r(numb1,3)
           write(6,*)r(numb2,1),r(numb2,2),r(numb2,3)
            PAUSE
          endif

          if(nres(i)-nres(i-1).gt.1)then
           write(6,120)'interpreted as a real chain break'
120        format(15x,a33)
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
C           PAUSE
          endif

         endif
         
c        no chain break...
         if(dabs(dist-1.3d0).le.1.0d-1)then
          if(nres(i)-nres(i-1).gt.1)then
           write(6,*)'              ','apparent up-skip not real',
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
C           PAUSE
          endif
         endif
         
        endif
       enddo
       
       do i=1,npdb
        nres(i) = nres(i)+irep(i)
       enddo

c      correct for the "negative N-terminal" residue numbering...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(nres(1).ne.1)then
        count = nres(1)
c        write(6,*)'corrected N-terminal numbering : ',nres(1)
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