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
c
c     exptl data and their uncertainties are now rescaled by the 
c     water & globbic corrections, 10/05/07
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine saxsinit

      implicit none
      include 'saxs.inc'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      update = 0
      nsegsaxs = 0
      isaxs  = 0
      nq     = 0

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
      integer i,j,k,l,length,isegm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      update = update+1   
      isegm  = 0   

c     read in the input parameters for the SAXS potential...
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
        
        else if(wd(1:4).eq.'NDAT') then
         call nexti('# of SAXS data points to fit =', nq)

        else if(wd(1:4).eq.'KFOR')then
         call nextf('SAXS force constant =', saxsforce)

        else if(wd(1:2).eq.'I0')then
         call nextf('SAXS I(0) =', I0expt)
                
c       update frequency
        else if(wd(1:4).eq.'UPDF') then
         call nexti('SAXS update frequency =', updfrq)

c       grid dimension
        else if(wd(1:4).eq.'NGRD') then
         call nexti('grid order =', ngrid)

c       grid type
        else if(wd(1:4).eq.'GRID') then
         call nextan('grid type ',gridtype ,0,4) 

c       number of active segments
        else if(wd(1:4).eq.'SEGM') then
         isegm = isegm+1
         nsegsaxs = nsegsaxs+1
         call nextan('segment name ',segmsaxs(isegm) ,0,4) 

c       min residue # 
        else if(wd(1:4).eq.'NMIN') then
         call nexti('min residue =', nminsaxs(isegm))

c       max residue # 
        else if(wd(1:4).eq.'NMAX') then
         call nexti('max residue =', nmaxsaxs(isegm))
              
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

         if(ftype.eq.'SCATDEF')then
          lenfile(1) = length
          do i=1,length
           fileloc(1)(i:i) = temp80(i:i)
          enddo
         endif
         if(ftype.eq.'SCAT_FQ')then
          lenfile(2) = length
          do i=1,length
           fileloc(2)(i:i) = temp80(i:i)
          enddo
         endif
         if(ftype.eq.'GRD_UNC')then
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
      double precision qi,fact1,fact2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     read in scatterer type definitions...
      write(6,5)'reading       scatterer definitions @',
     &fileloc(1)(1:lenfile(1))
5     format(a,1x,a)
      open(unit=55,file=fileloc(1)(1:lenfile(1)),status='old')
      rewind 55
	read(55,*)nscatname
      do i=1,nscatname
	 read(55,10)scatatom(i),scatres(i),scatid(i)
10     format(13x,a4,a4,3x,i2)
	enddo
	close(unit=55,status='keep')

c     read in atomic scattering profiles...
      write(6,5)'reading         scattering profiles @',
     &fileloc(2)(1:lenfile(2))
      open(unit=55,file=fileloc(2)(1:lenfile(2)),status='old')
      rewind 55
      read(55,*)
	read(55,20)(Fq0(j),j=1,14)   !  f(q) at q=0
	do i=1,nq
	 read(55,20)(Fq(j,i),j=1,14) !  f(q) at experimental q
20     format(15x,40(1x,f9.6))
	enddo
	close(unit=55,status='keep')

c     read the grid correction factors...
      write(6,5)'reading grid-based uncertainty      @',
     &fileloc(3)(1:lenfile(3))
      open(unit=55,file=fileloc(3)(1:lenfile(3)),status='old')
      rewind 55
	do i=1,nq
	 read(55,*)qi,gridunc(i)
	 write(6,*)qi,gridunc(i)
	enddo
	close(unit=55,status='keep')

c     read the experimental SAXS data file (GNOM format)...
      write(6,5)'reading      experimental SAXS data @',
     &fileloc(4)(1:lenfile(4))
      open(unit=55,file=fileloc(4)(1:lenfile(4)),status='old')
      rewind 55
      read(55,*)
      do i=1,nq
	 read(55,*)q(i),Iqexpt(i),dIqexpt(i)
	 write(6,*)q(i),Iqexpt(i),dIqexpt(i)
      enddo
	close(unit=55,status='keep')

c     read the water layer correction profile...
      write(6,5)'reading      surface H2O correction @',
     &fileloc(5)(1:lenfile(5))
      open(unit=55,file=fileloc(5)(1:lenfile(5)),status='old')
      rewind 55
	do i=1,nq
	 read(55,*)qi,watcorr(i)
	 write(6,*)qi,watcorr(i)
	enddo
	close(unit=55,status='keep')

c     rescale the experimental data
	do i=1,nq
       Iqexpt (i) = Iqexpt (i)/I0expt  ! rescale expt to I(0) = 1
       dIqexpt(i) = dIqexpt(i)/I0expt  ! rescale expt to I(0) = 1
       fact1 = dIqexpt(i)/Iqexpt(i)           ! include grid uncertainty into the total uncertainty
       fact2 = dsqrt(fact1**2+gridunc(i)**2)  ! include grid uncertainty into the total uncertainty
       dIqexpt(i) = Iqexpt(i)*fact2           ! include grid uncertainty into the total uncertainty
	 dIqexpt(i) = dIqexpt(i)**2   ! square it
       Iqcorr(i) = 1.0d0/watcorr(i) ! apply to Icalc
c      Iqcorr(i) = 1.0d0
	enddo

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

      integer i,j,k,l,id,jd,flag,nfib,count
      double precision store,fib(15),tp(512,3),h,pi
      character store1*4,store2*4
      integer nresid(10000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi = 3.1415926535898d0

c     assign integer residue numbers...
      count = 1
	store1 = resid(1)
	nresid(1) = 1
      do i=2,natom
       if((res(i).eq.'ANI ').or.(res(i).eq.'TIP3')) goto 5
       flag = 0
       do j=1,nsegsaxs
        if(segid(i).eq.segmsaxs(j)) flag = 1
       enddo
       if(flag.eq.1)then
	  if(resid(i).ne.store1)then
	   count = count+1
	   store1 = resid(i)
	  endif
	  nresid(i) = count
       endif
	enddo
5     continue

c     build the list of scatterers and their atom types...
      nscat = 0
      do i=1,natom
       flag = 0
       do j=1,nsegsaxs
        if(segid(i).eq.segmsaxs(j))then
         if((nresid(i).ge.nminsaxs(j)).and.
     &      (nresid(i).le.nmaxsaxs(j)))then
          flag = 1
         endif
        endif
       enddo
       if(flag.eq.0)goto 10
       flag = 0
       do j=1,nscatname
        if((res(i).eq.scatres(j)).and.(type(i).eq.scatatom(j)))then
         nscat = nscat+1
         scat(nscat) = i
         scatype(nscat) = scatid(j)
         goto 10
        endif
       enddo
10     continue
      enddo
	write(6,*)nscat, ' total scatterers'


c     calulate intensity at zero q...
      I0calc = 0.0d0
      do i=1,nscat
	 id = scatype(i)
       do j=1,nscat
	  jd = scatype(j)
        I0calc = I0calc+fq0(id)*fq0(jd)
       enddo
      enddo
      write(6,'(a,3x,es12.5)')'Icalc(0) = ',I0calc
      Afact = 1.0d0/I0calc


c     define the directional grid for averaging of the scattering amplitude...
      if(gridtype.eq.'FIBO')then
       if(ngrid.gt.15)then
        write(6,*)' too many equi-spaced vectors, decrease ngrid to <16'
       endif

       nfib = ngrid
       fib(1) = 1.0d0
	 fib(2) = 2.0d0
	 do i=3,nfib
	  fib(i) = fib(i-1)+fib(i-2)
	 enddo
       ngrid = fib(nfib)+1

	 do i=1,ngrid
	  tp(i,1) = dacos(1.0d0-2.0d0*dble(i-1)/fib(nfib))
	  tp(i,2) = 2.0d0*pi*dmod(i*fib(nfib-1)/fib(nfib),1.0d0)
        grid0(i,1) = dsin(tp(i,1))*dcos(tp(i,2))
        grid0(i,2) = dsin(tp(i,1))*dsin(tp(i,2))
        grid0(i,3) = dcos(tp(i,1))
	 enddo
      endif

      if(gridtype.eq.'SPIR')then
       if(ngrid.gt.1000)then
        write(6,*)' too many equi-spaced vectors, decrease ngrid to <1K'
       endif
       do i=1,ngrid
        h = -1.0d0+2.0d0*dble(i-1)/dble(ngrid-1)
        tp(i,1) = dacos(h)
        if((i.eq.1).or.(i.eq.ngrid))then
         tp(i,2) = 0.0d0
        else
         tp(i,2) = (tp(i-1,2) + 3.6d0/dsqrt(dble(ngrid)*(1.0d0-h**2)))
        endif
        grid0(i,1) = dsin(tp(i,1))*dcos(tp(i,2))
        grid0(i,2) = dsin(tp(i,1))*dsin(tp(i,2))
        grid0(i,3) = dcos(tp(i,1))
       enddo
      endif

c     intialize axes for random grid rotations
      nseed = 2*int(secnds(0.0))+1299709
      nfib = 12
      fib(1) = 1.0d0
	fib(2) = 2.0d0
	do i=3,nfib
	 fib(i) = fib(i-1)+fib(i-2)
	enddo
      ngridmc = fib(nfib)+1
	do i=1,ngridmc
	 tp(i,1) = dacos(1.0d0-2.0d0*dble(i-1)/fib(nfib))
	 tp(i,2) = 2.0d0*pi*dmod(i*fib(nfib-1)/fib(nfib),1.0d0)
       gridmc(i,1) = dsin(tp(i,1))*dcos(tp(i,2))
       gridmc(i,2) = dsin(tp(i,1))*dsin(tp(i,2))
       gridmc(i,3) = dcos(tp(i,1))
	enddo
      igridmc = 0


      return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SAXS data fitting routine
c       
c     This subroutine does the energies & forces calculation.
c     It is called by the ENERGY subrotine if SAXS flag is active.
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

      integer i,j,k,l,iq,igrid,iscat,id,imod

	double precision xav,yav,zav
	double precision gradscat(maxnq,maxnscat,3)
	double precision T1,T2,dtime
	double precision diq(maxnq),gradx,grady,gradz
      double precision totener
      double precision junk,theta,cc,ss,tt,rotm(3,3),rotv(3),axx,axy,axz
      double precision qr,mult
      double precision cosqr(maxnscat),sinqr(maxnscat)
      double precision ReAq,ImAq,fact

      character key*7
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      isaxs = isaxs+1

      if(isaxs.eq.1) call saxsinit2 ! prepare all necessary data

      if(mod(isaxs-1,3).ne.0.0d0) goto 1999

c     re-initialize grid and re-center every 50 steps

      if(mod(isaxs-1,updfrq).eq.0.0)then
       igridmc = igridmc+1
       if(igridmc.gt.ngridmc) igridmc = igridmc-ngridmc
       axx = gridmc(igridmc,1)
       axy = gridmc(igridmc,2)
       axz = gridmc(igridmc,3)
       junk = ran(nseed)
       theta = 360.0d0*ran(nseed)
       write(6,*)'random grid rotation through vector # ',igridmc
       write(6,'(a,3(f9.6))')'axis  = ',axx,axy,axz
       write(6,*)'angle = ',theta
       cc = dcosd(theta)
       ss = dsind(theta)
       tt = 1.0d0-dcosd(theta)
       rotm(1,1) = tt*axx*axx+cc
       rotm(1,2) = tt*axx*axy-ss*axz
       rotm(1,3) = tt*axx*axz+ss*axy
       rotm(2,1) = tt*axx*axy+ss*axz
       rotm(2,2) = tt*axy*axy+cc
       rotm(2,3) = tt*axy*axz-ss*axx
       rotm(3,1) = tt*axx*axz-ss*axy
       rotm(3,2) = tt*axy*axz+ss*axx
       rotm(3,3) = tt*axz*axz+cc
       do i=1,ngrid
        rotv(1) = 
     &  grid0(i,1)*rotm(1,1)+grid0(i,2)*rotm(1,2)+grid0(i,3)*rotm(1,3)
        rotv(2) = 
     &  grid0(i,1)*rotm(2,1)+grid0(i,2)*rotm(2,2)+grid0(i,3)*rotm(2,3)
        rotv(3) = 
     &  grid0(i,1)*rotm(3,1)+grid0(i,2)*rotm(3,2)+grid0(i,3)*rotm(3,3)
        grid(i,1) = rotv(1)
        grid(i,2) = rotv(2)
        grid(i,3) = rotv(3)

        do j=1,nq
         qx(j,i) = q(j)*grid(i,1)
         qy(j,i) = q(j)*grid(i,2)
         qz(j,i) = q(j)*grid(i,3)
        enddo
       enddo

       xav = 0.0d0
	 yav = 0.0d0
	 zav = 0.0d0
	 do i=1,nscat
	  xav = xav+x(scat(i))
	  yav = yav+y(scat(i))
	  zav = zav+z(scat(i))
 	 enddo
	 xav = xav/dble(nscat)
	 yav = yav/dble(nscat)
	 zav = zav/dble(nscat)
       write(6,'(a,3(f8.3))')'scattering center = ',xav,yav,zav
       do i=1,natom
        x(i) = x(i)-xav
        y(i) = y(i)-yav
        z(i) = z(i)-zav
       enddo

      endif

c     report Rg...
      if((mod(isaxs-1,updfrq).eq.0.0).or.(key.eq.'ANALYZE'))then
       Rgyr = 0.0d0
	 do i=1,nscat
        Rgyr = Rgyr+x(scat(i))**2+y(scat(i))**2+z(scat(i))**2
	 enddo
       Rgyr = sqrt(Rgyr/dble(nscat))
       write(6,'(a,f7.3)')'Rg = ',Rgyr
      endif

      T1 = secnds(0.0)

c     Calculate predicted scattering curve (update every step)
c     do iq=1,nq
c     Iqcalc(iq) = 0.0d0
c      do igrid=1,ngrid
c       ReAq(iq,igrid) = 0.0d0
c       ImAq(iq,igrid) = 0.0d0
c       do iscat=1,nscat
c        qr = qx(iq,igrid)*x(scat(iscat))+
c    &        qy(iq,igrid)*y(scat(iscat))+
c    &        qz(iq,igrid)*z(scat(iscat))
c        cosqr(iq,igrid,iscat) = dcos(qr)
c        sinqr(iq,igrid,iscat) = dsin(qr)
c        ReAq(iq,igrid) = ReAq(iq,igrid)+
c    &                    Fq(scatype(iscat),iq)*cosqr(iq,igrid,iscat)
c        ImAq(iq,igrid) = ImAq(iq,igrid)+
c    &                    Fq(scatype(iscat),iq)*sinqr(iq,igrid,iscat)
c       enddo
c       Iqcalc(iq) = Iqcalc(iq)+ReAq(iq,igrid)**2+ImAq(iq,igrid)**2
c      enddo
c      Iqcalc(iq) = Iqcalc(iq)/dble(ngrid)
c     enddo
c     T2 = secnds(0.0)


      do iq=1,nq
       Iqcalc(iq) = 0.0d0
       do iscat=1,nscat
        gradscat(iq,iscat,1) = 0.0d0
        gradscat(iq,iscat,3) = 0.0d0
        gradscat(iq,iscat,3) = 0.0d0
       enddo
      enddo

      do iq=1,nq
       do igrid=1,ngrid
        ReAq = 0.0d0
        ImAq = 0.0d0
        do iscat=1,nscat
         qr = qx(iq,igrid)*x(scat(iscat))+
     &        qy(iq,igrid)*y(scat(iscat))+
     &        qz(iq,igrid)*z(scat(iscat))
         cosqr(iscat) = dcos(qr)
         sinqr(iscat) = dsin(qr)
         ReAq = ReAq+Fq(scatype(iscat),iq)*cosqr(iscat)
         ImAq = ImAq+Fq(scatype(iscat),iq)*sinqr(iscat)
        enddo
        Iqcalc(iq) = Iqcalc(iq)+ReAq**2+ImAq**2
        do iscat=1,nscat
         fact = (cosqr(iscat)*ImAq-sinqr(iscat)*ReAq)
         gradscat(iq,iscat,1) = gradscat(iq,iscat,1)+qx(iq,igrid)*fact
         gradscat(iq,iscat,2) = gradscat(iq,iscat,2)+qy(iq,igrid)*fact
         gradscat(iq,iscat,3) = gradscat(iq,iscat,3)+qz(iq,igrid)*fact
        enddo
       enddo
       Iqcalc(iq) = Iqcalc(iq)/dble(ngrid)
      enddo

      Chi = 0.0d0
      do i=1,Nq
	 Iqcalc(i)  = Afact*Iqcorr(i)*Iqcalc(i) ! normalize by I(0), apply H2O correction
	 Chi    = Chi+(Iqexpt(i)-Iqcalc(i))**2/dIqexpt(i)
	 dIq(i) = 2.0d0*(Iqcalc(i)-Iqexpt(i))*Afact*Iqcorr(i)/
     &         (dIqexpt(i)*dble(Nq))
	enddo
	Chi = dsqrt(Chi/dble(Nq))
      totener = Chi
      if(mod(isaxs-1,1000).eq.0.0)then
	 open(unit=11,file='Iq_exp_calc.out',status='unknown')
	 rewind 11
       write(11,'(a,i5)')'isaxs = ',isaxs
	 do i=1,nq
	  write(11,1010)q(i),Iqexpt(i),Iqcalc(i),dsqrt(dIqexpt(i)),chi
1010    format(es13.6,2x,100(1x,es12.5))
	 enddo
	 close(unit=11,status='keep')
      endif

c     if(key.eq.'ANALYZE') goto 2000

c     calculation of the atomic forces (update every step)...
c	do iscat=1,nscat
c      id = scat(iscat)
c      gradx = 0.0d0
c      grady = 0.0d0
c      gradz = 0.0d0
c	 do iq=1,nq
c       fact = 2.0d0*Fq(scatype(iscat),iq)/dble(ngrid)
c       gradscat(iq,iscat,1) = gradscat(iq,iscat,1)*fact
c       gradscat(iq,iscat,2) = gradscat(iq,iscat,2)*fact
c       gradscat(iq,iscat,3) = gradscat(iq,iscat,3)*fact
c       gradx = gradx+gradscat(iq,iscat,1)*dIq(iq)
c       grady = grady+gradscat(iq,iscat,2)*dIq(iq)
c       gradz = gradz+gradscat(iq,iscat,3)*dIq(iq)
c	 enddo
c      dx(id) = dx(id)+gradx*saxsforce
c      dy(id) = dy(id)+grady*saxsforce
c      dz(id) = dz(id)+gradz*saxsforce
c     enddo


c     calculation of the atomic forces (do not update every step)...
 	do iscat=1,nscat
       id = scat(iscat)
       gradscat0(iscat,1) = 0.0d0
       gradscat0(iscat,2) = 0.0d0
       gradscat0(iscat,3) = 0.0d0
 	 do iq=1,nq
        fact = 2.0d0*Fq(scatype(iscat),iq)/dble(ngrid)
        gradscat0(iscat,1) = 
     &  gradscat0(iscat,1)+gradscat(iq,iscat,1)*dIq(iq)
        gradscat0(iscat,2) = 
     &  gradscat0(iscat,2)+gradscat(iq,iscat,2)*dIq(iq)
        gradscat0(iscat,3) = 
     &  gradscat0(iscat,3)+gradscat(iq,iscat,3)*dIq(iq)
       enddo
      enddo

1999  continue

c     do not update every step...
 	do iscat=1,nscat
       id = scat(iscat)
       dx(id) = dx(id)+gradscat0(iscat,1)*saxsforce
       dy(id) = dy(id)+gradscat0(iscat,2)*saxsforce
       dz(id) = dz(id)+gradscat0(iscat,3)*saxsforce
      enddo

      T2 = secnds(0.0)
      dtime = T2-T1

c     if(mod(isaxs,10).eq.0.0)then
 	 write(6,*)
 	 write(6,'(a,i6,a,f10.5,a,f8.4,x,a)')
     & ' SAXS:  step = ',isaxs,'  Chi = ',chi,
     &  '        ',dtime,' s calculation time'
c     endif

c     write to the console...
      if(mod(isaxs-1,updfrq).eq.0.0)then
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
      
      dr = 0.5d0

      do i=1,1000
       hist(i) = 0.0d0
      enddo

      do i=1,nscat
	 wi = Fq0(scatype(i))
       do j=i+1,nscat
	  wj = Fq0(scatype(j))
        dist = dsqrt((x(scat(i))-x(scat(j)))**2+
     &               (y(scat(i))-y(scat(j)))**2+
     &               (z(scat(i))-z(scat(j)))**2)
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
      open(unit=10,file='Pr.dat',status='unknown'); rewind 10
      write(10,'(a,i5)')'isaxs = ',isaxs
      do i=1,1000
       write(10,'(f6.2,2x,es12.5)')dr*(dble(i)-0.5d0),hist(i)
       if((hist(i).ne.0.0d0).and.(hist(i+1).eq.0.0d0)) goto 21
      enddo
21    continue
      close(unit=10,status='keep')


      write(nunit,'(a)')'REMARK SAXS  =============================='//
     &'================================='
      write(nunit,'(a,1x,f7.5)')'REMARK SAXS  Chi    = ',chi
      write(nunit,'(a,1x,f7.3)')'REMARK SAXS  Force  = ',saxsforce
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  Rgyr   = ',rgyr
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  Dmax   = ',dmax
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  '
      write(nunit,'(a,1x,f7.4)')'REMARK SAXS  Fit                  q'//
     &'       Icalc(q)     Iexpt(q)     Sigma(q)'
      do i=1,nq
       write(nunit,100)'REMARK SAXS  Fit     ',
     & q(i),Iqcalc(i),Iqexpt(i),dsqrt(dIqexpt(i))
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