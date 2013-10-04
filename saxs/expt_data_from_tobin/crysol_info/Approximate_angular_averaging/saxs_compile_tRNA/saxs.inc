cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This is the data block for the 
c     SAXS data fitting pseudo-potential

      integer maxnq,maxnscat,maxscatname
      parameter (maxnq       =  100) ! max number of data points
      parameter (maxnscat    = 6000) ! max # of scatterers
      parameter (maxscatname =  200) ! 85 for all heavy atoms in RNA 

      integer update,nsegsaxs,isaxs,nq
      integer scatid(maxscatname),nscatname,updfrq,saxsflag,lenfile(5)
      integer ngrid,ngridmc,igridmc,nscat,scat(maxnscat)
      integer scatype(maxnscat),nseed,nminsaxs(32),nmaxsaxs(32)

      double precision saxsforce,I0expt,Fq0(20),Fq(20,maxnq)
	  
	  double precision gridunc(maxnq),q(maxnq),Iqexpt(maxnq)
      double precision dIqexpt(maxnq),watcorr(maxnq),Iqcorr(maxnq)
      double precision I0calc,Iqcalc(maxnq),afact,rgyr,chi
      double precision grid(512,3),grid0(512,3),gridmc(512,3)
	  double precision qx(maxnq,512),qy(maxnq,512),qz(maxnq,512)
	  double precision gradscat0(maxnscat,3)

      character segmsaxs(32)*4,fileloc(5)*80       
	  character scatatom(maxscatname)*4,scatres(maxscatname)*4
	  character gridtype*4

      common /saxsintvar/  update,nsegsaxs,isaxs,nq,scatid,nscatname,
     &                     updfrq,saxsflag,lenfile,ngrid,ngridmc,
     &                     igridmc,nscat,scat,scatype,nseed,
     &                     nminsaxs,nmaxsaxs

      common /saxsfpvar /  saxsforce,I0expt,Fq,Fq0,gridunc,q,Iqexpt,
     &                     dIqexpt,watcorr,Iqcorr,I0calc,Iqcalc,
     &                     afact,rgyr,chi,grid,grid0,gridmc,qx,qy,qz,
     &                     gradscat0

      common /saxscharvar/ segmsaxs,fileloc,scatatom,scatres,gridtype