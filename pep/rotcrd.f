      subroutine ROTCRD(jd,fract,x,nvel,nplntr,n)
 
      implicit none
c
c m.e.ash,r.w.king sept 1972  subroutine rotcrd
c perform interpolation for moon,planet or earth rotation coordinate
c
c for earth rotation, interpolation indices are determined for
c receive (n=1) and send (n=2) times

c arguments
      integer*4 jd,nvel,nplntr,n
      real*10 fract,x(6)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'pqind.inc'
      include 'rotdta.inc'

c local
      real*10 ypr(7,3,6),yer(7,3,6)
      integer*4 imrcl/6/,iercl/7/,lsw
 
      if(nplntr.eq.3) then
c
c read earth rotation tape
         lsw = n - 2
         call EVTRP(jd,fract,0,lsw,iercl,yer,x,Earth,i_mxplprt+1,
     .              Jder,Fer,Per(1,n))
         if(jd.gt.0 .and. nvel.gt.0)
     .       call EVTRP(jd,fract,1,0,iercl,yer,x,Earth,i_mxplprt+1,
     .       Jder,Fer,Per(1,n))
      else
 
         call EVTRP(jd,fract,0,-1,imrcl,ypr,x,Plnmon,i_mxplprt+1,
     .              Jdpr,Fpr,Ppr)
         if(jd.gt.0 .and. nvel.gt.0)
     .       call EVTRP(jd,fract,1,0,imrcl,ypr,x,Plnmon,i_mxplprt+1,
     .       Jdpr,Fpr,Ppr)
      endif
 
      return
      end
