cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This is the data block for the 
c     SAXS data fitting pseudo-potential

      integer maxns,maxnglob,maxglobtype,maxglobnat,maxnres
      parameter (maxns       =  200)
      parameter (maxnglob    = 2500)
      parameter (maxglobtype =   37)
      parameter (maxglobnat  =   16)
      parameter (maxnres     = 1000)

      integer          nsaxs,nglob,totglob,natglob(maxglobtype)
      integer          idat(maxnglob,maxglobnat),idglob(maxnglob)
      integer          isaxs,saxsflag,gat(5,5,1000),lenfile(5)
      integer          prnfrq,update,natprot

      double precision is0(maxns),dis0(maxns),s(maxns),saxsforce
      double precision is(maxns),iscorr(maxns)
      double precision fsglob(maxglobtype,maxns),i0,i0calc,chi
      double precision fs0glob(maxglobtype)
      double precision zweight(maxnglob,maxglobnat),rgyr,afact
      double precision watcorr(maxns)
      character        atglob(maxglobtype,maxglobnat)*4
      character        resglob(maxglobtype)*3
      character        fileloc(5)*80


      common /saxsintvar/  nsaxs,nglob,totglob,natglob,idat,idglob,
     &                     isaxs,saxsflag,gat,lenfile,prnfrq,update,
     &                     natprot
      common /saxsfpvar /  is0,dis0,s,saxsforce,is,iscorr,
     &                     fsglob,zweight,i0,i0calc,
     &                     fs0glob,rgyr,chi,afact,watcorr
      common /saxscharvar/ atglob,resglob,fileloc