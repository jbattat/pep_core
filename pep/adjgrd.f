      subroutine ADJGRD(ind,ll,aplnt,npl)
 
      implicit none
c
c        r.b. goldstein sept. 1978
c        called from adjust for grid topography model parameters
c
c arguments
      integer ind,ll
      character*8 aplnt
      integer*2 npl
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'inodta.inc'
      include 'scoef4.inc'
      real*4 Grid(u_stdsz,4,1000/u_stdsz)
      integer*2 Lgrid(4,1000)
      equivalence (Pzhar,Grid),(Lpzhar,Lgrid)
 
c local
      real*10 fact, old
      integer   i,ishp,k
      character*1 astrik(3)/'*','&',' '/
 
      fact = 1._10
      ishp = 0
      call LINCHK
      do i = 1,1000/u_stdsz
         do k = 1,u_stdsz
            ishp = ishp + 1
            if(ishp.gt.ll) goto 100
            if(Lgrid(ind,ishp).gt.0) then
               old = Grid(k,ind,i)
               call ADJAST(old,fact)
               write(Iout,10) astrik(Ntype),N,ishp,aplnt,npl,old,
     .                         Adj,Nwv,Sig,Fract
   10          format(1x,a1,i4,'. LGRID(',i4,')= 1 GRID(KM) OF ',
     .                a8,i3,8x,1pd22.15,1pd16.8,1pd22.15,1pd10.3,
     .                0pf8.3)
               if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,
     .            ishp,aplnt,npl,old,Adj,Nwv,Sig,Fract
               if(Keepit) Grid(k,ind,i) = Nwv
            endif
         end do
      end do
  100 call FRSTAT(1,2,'GRD HGTS')
 
      return
      end
