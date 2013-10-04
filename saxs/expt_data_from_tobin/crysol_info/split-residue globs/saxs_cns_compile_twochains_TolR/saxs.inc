cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This is the data block for the 
c     SAXS data fitting pseudo-potential

      integer maxnq,maxnglob,maxglobtype,maxglobnat,maxnres
      parameter (maxnq       =  200)
      parameter (maxnglob    = 2500)
      parameter (maxglobtype =   37) ! 37 for split-res, 20 for whole-res
      parameter (maxglobnat  =   16) ! 16 for split-res, 25 for whole-res
      parameter (maxnres     = 2000)

      integer          nsaxs,nglob,totglob,natglob(maxglobtype)
      integer          idat(maxnglob,maxglobnat),idglob(maxnglob)
      integer          isaxs,saxsflag,gat(5,5,1000),lenfile(5)
      integer          prnfrq,update,natsaxs,nsegsaxs

      double precision Iqexpt(maxnq),dIqexpt(maxnq),q(maxnq),saxsforce
      double precision Iqcalc(maxnq),Iqcorr(maxnq)
      double precision fqglob(maxglobtype,maxnq),I0expt,I0calc,chi
      double precision fq0glob(maxglobtype)
      double precision zweight(maxnglob,maxglobnat),rgyr,afact
      double precision watcorr(maxnq),globcorr(maxnq)
      character        atglob(maxglobtype,maxglobnat)*4
      character        resglob(maxglobtype)*3
      character        fileloc(5)*80
      character        segmsaxs(32)*4


      common /saxsintvar/  nsaxs,nglob,totglob,natglob,idat,idglob,
     &                     isaxs,saxsflag,gat,lenfile,prnfrq,update,
     &                     natsaxs,nsegsaxs
      common /saxsfpvar /  Iqexpt,dIqexpt,q,saxsforce,Iqcalc,Iqcorr,
     &                     fqglob,zweight,I0expt,I0calc,
     &                     fq0glob,rgyr,chi,afact,watcorr,globcorr
      common /saxscharvar/ atglob,resglob,fileloc,segmsaxs