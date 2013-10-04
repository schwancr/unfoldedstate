cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SAXS data fitting routines
c      
c     This subroutine initializes arrays needed for SAXS data fitting.
c     It is called by the main CNS program.
c
c     Alexander Grishaev, NIH, 12/01/03 - 7/16/04 - 2/5/05
c     In this version, there is no adjustable scaling for the calc data
c
c     This version is built for the MSG data - no solvent layer + 
c     smaller RAM overhead, no drho is input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsinit

      implicit none
      include 'saxs.inc'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      update = 0
      isaxs  = 0
      nsaxs  = 0

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c
c     SAXS data fitting routines
c
c     This subroutine reads "SAXS" statement from the CNS input script
c     It is called by the main CNS program 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsread

      implicit none
      include 'cns.inc'
      include 'comand.inc'
      include 'timer.inc'
      include 'consta.inc'
      include 'saxs.inc'

      character ftype*7,temp80*80
      integer i,j,k,l,length
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      update = update+1

c     read in the input parameters for the HB potential...
      call pusend('SAXS>')
      done =.false.
      do while(.not.done)
       call nextwd('SAXS>')
       call miscom('SAXS>',used)
       if (.not.used) then

        if (wd(1:4).eq.'HELP')then
         WRITE(DUNIT,'(20X,A)')
     &   'SAXS {<SAXS_data_fit_statement>} END ',
     &   '<SAXS-statement>:== ',
     &   ' FORCeconstant <real> ! force constant for all ',
     &   '      ! assignments in the current class. ',
     &   '      ! {default = 50} '
        
        else if(wd(1:4).eq.'NTYP') then
         call nexti('# of glob types =', nglob)        
        else if(wd(1:4).eq.'NDAT') then
         call nexti('# of SAXS data points to fit =', nsaxs)
        else if(wd(1:4).eq.'KFOR')then
         call nextf('SAXS force constant =', saxsforce)
        else if(wd(1:2).eq.'I0')then
         call nextf('SAXS I(0) =', i0)
                
c       print frequency
        else if(wd(1:4).eq.'PRNF') then
         call nexti('SAXS printing frequency =', prnfrq)
                 
c       print mode...
        else if(wd(1:4).eq.'PRIN') then 
         call saxsprint(punit)
        
c       locations of the variuos SAXS-related files...
        else if(wd(1:4).eq.'FILE') then
         call nextan('file type =',ftype,0,8) ! modified nexta
         do i=1,80
          temp80(i:i)=' '
         enddo
         call nextan('file location =',temp80,1,length)
         if(ftype.eq.'GLOBDEF')then
          lenfile(1) = length
          do i=1,length
           fileloc(1)(i:i) = temp80(i:i)
          enddo
         endif
         if(ftype.eq.'GLOB_FS')then
          lenfile(2) = length
          do i=1,length
           fileloc(2)(i:i) = temp80(i:i)
          enddo
         endif
         if(ftype.eq.'IS_CORR')then
          lenfile(3) = length
          do i=1,length
           fileloc(3)(i:i) = temp80(i:i)
          enddo
         endif
         if(ftype.eq.'SAXSDAT')then
          lenfile(4) = length
          do i=1,length
           fileloc(4)(i:i) = temp80(i:i)
          enddo
         endif
         if(ftype.eq.'H2O_COR')then
          lenfile(5) = length
          do i=1,length
           fileloc(5)(i:i) = temp80(i:i)
          enddo
         endif
         
c       check for the END statement...
        else
         call chkend('SAXS>', done)
        endif
       endif
      enddo  
      done = .false.
      
c     read in the potential files and calculate the forces...
      if(update.eq.1) call saxsread2
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SAXS data fitting term.
c      
c     This subrotine reads in the SAXS-related files.
c     It is called by READSAXS subroutine.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsread2
      
      implicit none
      include 'consta.inc'
      include 'saxs.inc'
      integer i,j,k,l
      double precision si
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     read in glob definitions...
      write(6,5)'reading            glob definitions @',
     &fileloc(1)(1:lenfile(1))
5     format(a,1x,a)
      open(unit=55,file=fileloc(1)(1:lenfile(1)),status='old')
      rewind 55
	read(55,'(i3)')nglob
      do i=1,nglob
	 read(55,10)resglob(i),natglob(i),(atglob(i,j),j=1,natglob(i))
10     format(a3,1x,i2,1x,30(1x,a4))
	enddo
	close(unit=55,status='keep')

c     read in globbic scattering profiles...
      write(6,5)'reading globbic scattering profiles @',
     &fileloc(2)(1:lenfile(2))
      open(unit=55,file=fileloc(2)(1:lenfile(2)),status='old')
      rewind 55
      read(55,*)
	read(55,20)(fs0glob(j),j=1,nglob)   ! globbic f(s) at s=0
	do i=1,nsaxs
	 read(55,20)(fsglob(j,i),j=1,nglob) ! globbic f(s) at experimental s
20     format(15x,40(1x,f9.6))
	enddo
	close(unit=55,status='keep')

c     read the globbic correction factors...
      write(6,5)'reading globbic correction function @',
     &fileloc(3)(1:lenfile(3))
      open(unit=55,file=fileloc(3)(1:lenfile(3)),status='old')
      rewind 55
	do i=1,nsaxs
	 read(55,*)si,iscorr(i)
	 write(6,*)si,iscorr(i)
	enddo
	close(unit=55,status='keep')

c     read the experimental SAXS data file (GNOM format)...
      write(6,5)'reading      experimental SAXS data @',
     &fileloc(4)(1:lenfile(4))
      open(unit=55,file=fileloc(4)(1:lenfile(4)),status='old')
      rewind 55
      read(55,*)
      do i=1,nsaxs
	 read(55,*)s(i),is0(i),dis0(i)
	 write(6,*)s(i),is0(i),dis0(i)
      enddo
	close(unit=55,status='keep')
	do i=1,nsaxs
       is0(i) = is0(i)/i0
       dis0(i) = dis0(i)/i0  ! absolute data uncertainty
	 dis0(i) = dis0(i)**2
	enddo

c     read the water layer correction profile...
      write(6,5)'reading      surface H2O correction @',
     &fileloc(5)(1:lenfile(5))
      open(unit=55,file=fileloc(5)(1:lenfile(5)),status='old')
      rewind 55
	do i=1,nsaxs
	 read(55,*)si,watcorr(i)
	 write(6,*)si,watcorr(i)
       is0(i) = is0(i)*watcorr(i)
	enddo
	close(unit=55,status='keep')

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SAXS data fitting routines
c      
c     This subroutine initializes arrays needed for SAXS data fitting.
c     It is called at the beginning of esaxs during its first pass.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsinit2

      implicit none
      include 'cns.inc'
      include 'mtf.inc'
      include 'saxs.inc'

      integer i,j,k,l,count,totres,natpdb(maxnres),idg,count1
      integer nresid(15000),jdg
 	integer idpdb(maxnres,40),istore,glocc(natom),id0(32)
	character store1*4,store2*4,typei*4,typej*4
      character atpdb(maxnres,40)*4,respdb(maxnres)*3,temp*4
      double precision store
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     assign integer residue numbers...
      count = 1
	store1 = resid(1)
	store2 = segid(1)
	nresid(1) = 1
      do i=2,natom
       if((res(i).eq.'ANI ').or.(res(i).eq.'TIP3')) goto 10
	 if(resid(i).ne.store1)then
	  if(segid(i).eq.store2)then
	   count = count+1
	   store1 = resid(i)
	   store2 = segid(i)
	  endif
	 endif
	 nresid(i) = count
	enddo
10    continue

c     form residue-based arrays within the protein...
	do i=1,maxnres
	 natpdb(i) = 0
	enddo
      natprot = 0
      do i=1,natom
       if((res(i).eq.'ANI ').or.(res(i).eq.'TIP3'))then
        natprot = i-1
        totres = nresid(i-1)
        goto 20
       endif
	 natpdb(nresid(i)) = natpdb(nresid(i))+1
       if(natpdb(nresid(i)).gt.40)then
        write(6,*)'atid ',i,' resid ',nresid(i),' :too many atoms/aa'
        pause
       endif
	 atpdb (nresid(i),natpdb(nresid(i))) = type(i)
	 respdb(nresid(i)) = res(i)(1:3)
	 idpdb (nresid(i),natpdb(nresid(i))) = i
	enddo
      totres = nresid(natom)
20    continue
      write(6,*)natom,' total atoms'
      write(6,*)natprot,' protein atoms'
      write(6,*)totres,' protein residues'

c     split the structure into globs...
      totglob = 0
      do i=1,totres
c      peptide H-N-C=O group (glob #1)...
       if((i.lt.totres).and.(respdb(i+1).ne.'PRO'))then
	  idg = 1
	  do k=1,natglob(idg)
	   id0(k) = 0
	  enddo
	  do k=1,natpdb(i)
	   do l=1,2
	    if(atglob(idg,l).eq.atpdb(i,k))then
	     id0(l) = k
	    endif
         enddo
	  enddo
        do k=1,natpdb(i+1)
	   do l=3,4
	    if(atglob(idg,l).eq.atpdb(i+1,k))then
	     id0(l) = k
	    endif
         enddo
	  enddo
        istore =1; do k=1,natglob(idg); istore = istore*id0(k); enddo
        if(istore.ne.0)then
         totglob = totglob+1
	   idglob(totglob) = idg
         do k=1,2
          idat(totglob,k) = idpdb(i,id0(k))
         enddo 
         do k=3,4
          idat(totglob,k) = idpdb(i+1,id0(k))
         enddo
        endif
	 endif
c      Proline peptide N-C=O group (glob #2)...
       if((i.lt.totres).and.(respdb(i+1).eq.'PRO'))then
	  idg = 2
        do k=1,natglob(idg)
	   id0(k) = 0
	  enddo
	  do k=1,natpdb(i)
	   do l=1,2
	    if(atglob(idg,l).eq.atpdb(i,k))then
	     id0(l) = k
          endif
         enddo
	  enddo
	  do k=1,natpdb(i+1)
	   do l=3,3
	    if(atglob(idg,l).eq.atpdb(i+1,k))then
	     id0(l) = k
	    endif
         enddo
        enddo
        istore =1; do k=1,natglob(idg); istore = istore*id0(k); enddo
        if(istore.ne.0)then
	   totglob = totglob+1
         idglob(totglob) = idg
         do k=1,2
          idat(totglob,k) = idpdb(i,id0(k))
         enddo 
         do k=3,3
          idat(totglob,k) = idpdb(i+1,id0(k))
         enddo
        endif
	 endif
c      N-terminal NH3+ group (glob #3)...
       if(i.eq.1)then
	  idg = 3
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
        istore =1; do k=1,natglob(idg); istore = istore*id0(k); enddo
        if(istore.ne.0)then
	   totglob = totglob+1
	   idglob(totglob) = idg
         do k=1,natglob(idg)
          idat(totglob,k) = idpdb(i,id0(k))
	   enddo
        endif
	 endif
c      C-terminal COO- group (glob #4)...
       if(i.eq.totres)then
	  idg = 4
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
        istore =1; do k=1,natglob(idg); istore = istore*id0(k); enddo
        if(istore.ne.0)then
	   totglob = totglob+1
	   idglob(totglob) = idg
         do k=1,natglob(idg)
          idat(totglob,k) = idpdb(i,id0(k))
	   enddo
        endif
	 endif
c      sidechain globs (#5 and up)...
	 do idg=5,nglob
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
           idat(totglob,k) = idpdb(i,id0(k))
	    enddo
         endif
        endif
       enddo
      enddo
      write(6,*)totglob,' globs    found'

      count = 0
      do i=1,totglob
	 idg = idglob(i)
	 do j=1,natglob(idg)
	  count = count+1
	  glocc(idat(i,j)) = 1
	 enddo
	enddo

      count1 = 0
	do i=1,natom
       if((res(i).eq.'ANI ').or.(res(i).eq.'HOH ')) goto 30
	 if(glocc(i).ne.1)then
	  write(6,*)i,' ',res(i),' ',resid(i),' ',type(i),' not in globs'
	  count1 = count1+1
	 endif
	enddo
30    continue
	write(6,*)count1, ' unused atoms'
	write(6,*)totglob,' globs built'
	write(6,*)count,  ' atoms in globs'

c     assign nuclear charge-based atomic weights within each glob...
      do i=1,totglob
	 idg = idglob(i)
	 store = 0.0d0
	 do j=1,natglob(idg)
        zweight(i,j) = 0.0d0
	  temp = type(idat(i,j))
c       XPLOR/CNS atomic naming convention...
	  if(temp(1:1).eq.'H') zweight(i,j) =  1.0d0
	  if(temp(1:1).eq.'D') zweight(i,j) =  1.0d0
	  if(temp(1:1).eq.'C') zweight(i,j) =  6.0d0
	  if(temp(1:1).eq.'N') zweight(i,j) =  7.0d0
	  if(temp(1:1).eq.'O') zweight(i,j) =  8.0d0
	  if(temp(1:1).eq.'S') zweight(i,j) = 16.0d0
        if(zweight(i,j).eq.0.0d0)then
	   write(6,*)resid(idat(i,j)),type(idat(i,j)),
     &   'unknown nuclear charge'
	   pause
        endif
	  store = store+zweight(i,j)
	 enddo
	 do j=1,natglob(idg)
	  zweight(i,j) = zweight(i,j)/store
	 enddo
      enddo

c     calulate intensity at zero q...
      i0calc = 0.0d0
      do i=1,totglob
	 idg = idglob(i)
       do j=1,totglob
	  jdg = idglob(j)
        i0calc = i0calc+fs0glob(idg)*fs0glob(jdg)
       enddo
      enddo
      write(6,'(3x,es9.3,a)')i0calc,' Icalc(0)'
      afact = 1.0d0/i0calc

1000  continue

      return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SAXS data fitting routine
c       
c     This subroutine does the energies & forces calculation.
c     It is called by the ENERGY subrotine if HBDB flag is active.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsener(totener,key)

      implicit none

	include 'cns.inc'
      include 'mtf.inc'
      include 'coord.inc'
      include 'deriv.inc'
      include 'saxs.inc'
      include 'consta.inc'
      include 'ener.inc' 

      integer i,j,k,l
      integer idg,idg1,idg2,ida,ida1,ida2,jdg,idi,idj,ihist

	double precision rij2,srij,xav,yav,zav,fact,store,r
	double precision dist,dist2
	double precision issum1,issum2,rglob(maxnglob,3)
	double precision gradgi(maxnglob,maxns,3),factx,facty,factz
	double precision gradga(maxnglob,3),t1,t2,dtime,xij,yij,zij
	double precision dis(maxns),gradgc2(3),gradx,grady,gradz
      double precision drx,dry,drz,fsijk,sinc,cosinc
      double precision totener

      character typei*4,typej*4,key*7
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      isaxs = isaxs+1

      if(isaxs.eq.1) call saxsinit2 ! prepare all necessary data

c     accumulate and report inter-atomic pdf's...
      if((mod(isaxs-1,prnfrq).eq.0.0).or.(key.eq.'ANALYZE'))then
       do i=1,natprot
        typei = type(i)
        idi = 0
        if(typei(1:1).eq.'H') idi = 1
        if(typei(1:1).eq.'C') idi = 2
        if(typei(1:1).eq.'N') idi = 3
        if(typei(1:1).eq.'O') idi = 4
        if(typei(1:1).eq.'S') idi = 5
        do j=i+1,natprot
         typej = type(j)
         idj = 0
         if(typej(1:1).eq.'H') idj = 1
         if(typej(1:1).eq.'C') idj = 2
         if(typej(1:1).eq.'N') idj = 3
         if(typej(1:1).eq.'O') idj = 4
         if(typej(1:1).eq.'S') idj = 5
         if(idj.ge.idi)then
          dist = dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
          ihist = 1+int(dist/0.1d0)
          if((idi.ne.0).and.(idj.ne.0))then
           if(ihist.le.1000)then
            gat(idi,idj,ihist) = gat(idi,idj,ihist)+1.0d0
           endif
          endif
         endif
        enddo
       enddo

c      report protein's Rg ...
       xav = 0.0d0
	 yav = 0.0d0
	 zav = 0.0d0
	 do i=1,natprot
	  xav = xav+x(i)
	  yav = yav+y(i)
	  zav = zav+z(i)
 	 enddo
	 xav = xav/dble(natprot)
	 yav = yav/dble(natprot)
	 zav = zav/dble(natprot)
       rgyr = 0.0d0
	 do i=1,natprot
        rgyr = rgyr+(x(i)-xav)**2+(y(i)-yav)**2+(z(i)-zav)**2
	 enddo
       rgyr = sqrt(rgyr/dble(natprot))
       
       open(unit=60,file='gat.dat',status='unknown'); rewind 60
       write(60,'(a,1x,i6)')   'Isaxs= ',isaxs
       write(60,'(a,1x,f5.2)') '   Rg= ',rgyr
       do k=1,1000
        write(60,10) (dble(k)-0.5d0)*0.1d0,((gat(i,j,k),j=i,5),i=1,5)
10      format(f6.2,50(1x,i7))
       enddo
       close(unit=60,status='keep')
       do i=1,5
        do j=1,5
         do k=1,1000
          gat(i,j,k) = 0
         enddo
        enddo
       enddo

      endif

c     Z-weighted globbic centers of mass (protein)...
      t1 = secnds(0.0)
      do i=1,totglob
       idg = idglob(i)
	 rglob(i,1) = 0.0d0
	 rglob(i,2) = 0.0d0
	 rglob(i,3) = 0.0d0
	 do j=1,natglob(idg)
	  rglob(i,1) = rglob(i,1)+x(idat(i,j))*zweight(i,j)
	  rglob(i,2) = rglob(i,2)+y(idat(i,j))*zweight(i,j)
	  rglob(i,3) = rglob(i,3)+z(idat(i,j))*zweight(i,j)
       enddo
      enddo

c     calculate globbic scattering curve...
      do k=1,nsaxs
	 is(k) = 0.0d0
	enddo

	do i=1,totglob
	 do k=1,nsaxs
        do l=1,3
 	   gradgi(i,k,l) = 0.0d0
        enddo
	 enddo
      enddo

	do i=1,totglob
 	 idg1 = idglob(i)
	 do k=1,nsaxs
	  is(k) = is(k)+fsglob(idg1,k)*fsglob(idg1,k)
	 enddo
	 do j=i+1,totglob
	  idg2 = idglob(j)
	  drx = rglob(i,1)-rglob(j,1)
	  dry = rglob(i,2)-rglob(j,2)
	  drz = rglob(i,3)-rglob(j,3)
	  dist2 = drx*drx+dry*dry+drz*drz
        dist  = dsqrt(dist2)
	  do k=1,nsaxs
	   srij   = s(k)*dist
         sinc   = dsin(srij)/srij
         cosinc = (dcos(srij)-sinc)/dist2
	   is(k) = is(k) + 2.0d0*fsglob(idg1,k)*fsglob(idg2,k)*sinc
 	   gradgi(i,k,1) = 
     &   gradgi(i,k,1) + fsglob(idg1,k)*fsglob(idg2,k)*cosinc*drx
 	   gradgi(i,k,2) = 
     &   gradgi(i,k,2) + fsglob(idg1,k)*fsglob(idg2,k)*cosinc*dry
         gradgi(i,k,3) = 
     &   gradgi(i,k,3) + fsglob(idg1,k)*fsglob(idg2,k)*cosinc*drz
 	   gradgi(j,k,1) = 
     &   gradgi(j,k,1) - fsglob(idg1,k)*fsglob(idg2,k)*cosinc*drx
 	   gradgi(j,k,2) = 
     &   gradgi(j,k,2) - fsglob(idg1,k)*fsglob(idg2,k)*cosinc*dry
         gradgi(j,k,3) = 
     &   gradgi(j,k,3) - fsglob(idg1,k)*fsglob(idg2,k)*cosinc*drz
	  enddo
       enddo
	enddo

c     correct the raw curve for the systematic globbic errors...
c     calculate chi...
c     calculate expt/calc difference...
      chi = 0.0d0
      do i=1,nsaxs
	 is(i)  = afact*is(i)/iscorr(i)
	 chi    = chi+(is0(i)-is(i))*(is0(i)-is(i))/dis0(i)
	 dis(i) = 2.0d0*(is0(i)-is(i))/dis0(i)
	enddo
      totener = chi
	chi = dsqrt(chi/dble(nsaxs))

 	write(6,'(a,i6,a,f10.5)')
     &'SAXS:  Isaxs = ',isaxs,'      Chi = ',chi
      if(mod(isaxs-1,prnfrq).eq.0.0)then
	 open(unit=11,file='is_exp_calc.out',status='unknown')
	 rewind 11
	 do i=1,nsaxs
	  write(11,1010)s(i),is0(i),is(i),dsqrt(dis0(i)),chi
1010    format(es13.6,2x,100(1x,es12.5))
	 enddo
	 close(unit=11,status='keep')
      endif

      if(key.eq.'ANALYZE') goto 2000

c     calculation of the atomic forces...
      do i=1,totglob
c      full globbic gradient...
	 gradgc2(1) = 0.0d0
	 gradgc2(2) = 0.0d0
	 gradgc2(3) = 0.0d0
	 do k=1,nsaxs
	  gradgc2(1) = gradgc2(1) + dis(k)*afact*gradgi(i,k,1)
	  gradgc2(2) = gradgc2(2) + dis(k)*afact*gradgi(i,k,2)
	  gradgc2(3) = gradgc2(3) + dis(k)*afact*gradgi(i,k,3)
	 enddo
c      atomic gradients...
       idg = idglob(i)
	 do j=1,natglob(idg)
	  ida   = idat(i,j)
	  fact  = saxsforce*zweight(i,j)
	  gradx = fact*gradgc2(1)
	  grady = fact*gradgc2(2)
	  gradz = fact*gradgc2(3)
        dx(ida) = dx(ida)-gradx
        dy(ida) = dy(ida)-grady
        dz(ida) = dz(ida)-gradz
	 enddo
      enddo
      t2 = secnds(0.0)
      dtime = t2-t1
 	write(6,'(a,1x,f6.2,1x,a)')'SAXS: ',dtime,
     &'seconds for energy/forces calculation'

c     write to the console...
      if(mod(isaxs-1,prnfrq).eq.0.0)then
       call saxsprint(6)
      endif
	 
2000  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsprint(nunit)
      
      include 'cns.inc'
      include 'mtf.inc'
      include 'coord.inc'
      include 'saxs.inc'
      integer nunit,i,j,k,l,ihist,istore
      double precision dist,hist(1000),dr,wi,wj,store,dmax,esaxs
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(nunit.ne.6)then
       call saxsener(esaxs,'ANALYZE')
      endif
      
      dr = 0.2d0

      do i=1,1000
       hist(i) = 0.0d0
      enddo

      do i=1,natom-1
       wi = 0.0d0
       if(type(i).eq.'OO  ') goto 20
	 if(type(i)(1:1).eq.'H') wi =  1.0d0
	 if(type(i)(1:1).eq.'D') wi =  1.0d0
	 if(type(i)(1:1).eq.'C') wi =  6.0d0
	 if(type(i)(1:1).eq.'N') wi =  7.0d0
	 if(type(i)(1:1).eq.'O') wi =  8.0d0
	 if(type(i)(1:1).eq.'S') wi = 16.0d0
       do j=i+1,natom
        wj = 0.0d0
       if(type(j).eq.'OO  ') goto 10
        if(type(j)(1:1).eq.'H') wj =  1.0d0
	  if(type(j)(1:1).eq.'D') wj =  1.0d0
	  if(type(j)(1:1).eq.'C') wj =  6.0d0
	  if(type(j)(1:1).eq.'N') wj =  7.0d0
	  if(type(j)(1:1).eq.'O') wj =  8.0d0
	  if(type(j)(1:1).eq.'S') wj = 16.0d0
        dist = sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
        ihist = 1+int(dist/dr)
        if((ihist.gt.0).and.(ihist.le.1000))then
         hist(ihist) = hist(ihist)+wi*wj
        endif
10     continue
       enddo
20     continue
      enddo

      istore = 0
      store = 0.0d0
      do i=1,1000
       store = store+hist(i)
       if(hist(i).ne.0.0d0) istore = i
      enddo
      dmax = dr*(dble(istore)-0.5d0)
      do i=1,1000
       hist(i) = hist(i)/store
      enddo

      write(nunit,'(a)')'REMARK SAXS  =============================='//
     &'================================='
      write(nunit,'(a,1x,f7.5)')'REMARK SAXS  Chi    = ',chi
      write(nunit,'(a,1x,f7.3)')'REMARK SAXS  Force  = ',saxsforce
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  Rgyr   = ',rgyr
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  Dmax   = ',dmax
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  '
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  Fit                  s'//
     &'       Icalc(s)     Iexpt(s)     Sigma(s)'
      do i=1,nsaxs
       write(nunit,100)'REMARK SAXS  Fit     ',
     & s(i),is(i),is0(i),dsqrt(dis0(i))
100    format(a,1x,es13.6,2x,3(1x,es12.5))
      enddo
      if(nunit.eq.6) goto 200
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  '
      write(nunit,'(a,1x,f7.5)')'REMARK SAXS            Weighted P(r)'
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  '
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS            r,A'//
     &'             P(r)'
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  '

      do i=1,1000
       write(nunit,'(a,1x,f6.2,2x,a,es12.5)')'REMARK SAXS          ',
     & dr*(dble(i)-0.5d0),'        ',hist(i)
       if((hist(i).ne.0.0d0).and.(hist(i+1).eq.0.0d0)) goto 200
      enddo

200   continue
      write(nunit,'(a)')'REMARK SAXS  =============================='//
     &'================================='


      return
      end