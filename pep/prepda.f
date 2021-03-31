      subroutine PREPDA(buff)

      implicit none


c*** start of declarations inserted by spag
      integer   i, lrec, mrec, ndim
c*** end of declarations inserted by spag


c
c d. white  april 1974  subroutine prepda
c
c z. goldberg   march 1980 -- expanded common filtim
c
c
c initialize direct access data set
c use dfile paul macneil june, 1978
c
c parameter
      real*10 buff(1000)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'filtda.inc'
c
c local
c
      integer  firetc


c
c maxp = maximum number of parameters
c Maxe = maximum number of epochs
c maxp, Maxe used for allocation purposes only
      if(Nparam.gt.Maxp)
     .     call SUICID('MAXP TOO SMALL, STOP IN PREPDA  ', 8)
      do i = 1, 1000
         buff(i) = 0._10
      end do
c
c mfile = norm eqn
c max number of records = Maxe*(maxp + 1)
c max record lgth = maxp "double" words
      Mfile = 80
      mrec  = Maxe*(Maxp + 1)
      ndim  = Maxp
      lrec  = 4*u_stdsz*ndim
c
c open a direct access file
c
c Unbelievably, IBM's VS fortran requires knowledge of the file
c size BEFORE the open statement.  Thus the non-portable fileinf call.
c

c     nonportable IBMism

      call fileinf(firetc, 'MAXREC', mrec)
      open(Mfile, access='DIRECT', recl=lrec, status='SCRATCH')
      call PREPDB(buff, ndim, mrec, Mfile)
c
c lfile = smear matrix
c max number of records = Maxe
c max record length = 400 items (20 * 20)
      Lfile = 81
      mrec  = Maxe
      ndim  = 400
      lrec  = 4*u_stdsz*ndim

c     nonportable IBMism

      call fileinf(firetc, 'MAXREC', mrec)
      open(Lfile, access='DIRECT', recl=lrec, status='SCRATCH')
      call PREPDB(buff, ndim, mrec, Lfile)
c
c define direct access for filter-predct link
c max number of records=Maxe
c max record length=(maxp+2) "double" words
      Kfile = 82
      mrec  = Maxe
      ndim  = Maxp+2
      lrec  = 4*u_stdsz*ndim

c     nonportable IBMism

      call fileinf(firetc, 'MAXREC', mrec)
      open(Kfile, access='DIRECT', recl=lrec, status='SCRATCH')

      call PREPDB(buff, ndim, mrec, Kfile)

      return
      end
