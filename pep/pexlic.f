      subroutine PEXLIC
c set up partial derivative controls for pulsar companions
c Note: any exo-planets from the obslib must also be input and must
c be in the same order.
      implicit none

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdat.inc'
      include 'empcnd.inc'
      include 'ennips.inc'
      include 'funcon.inc'
      include 'lcntrl.inc'
      include 'ltrapx.inc'
      include 'maxpxdat.inc'
      include 'mtrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'zeroes.inc'

c local
      integer*4 l,m,n

      Nmpex = 0
      if(Nplsr.gt.0) then
         m=1
         n=0
         do l=1,Numpln
            if(Nplnt(l).gt.0 .and. Npcent(l).eq.-4 .and. 
     .       Aplnt(l)(1:4).eq.Spotf) then
               n=n+1
               if(n.gt.maxpex) call SUICID(
     . 'TOO MANY EXO-PLANETS, STOP IN PEXLIC',9)
               Klanex(n)=l
               Nplex(n)=Nplnt(l)
               if(Iabs1.gt.0 .and. 
     .          m.le.Mmpex .and. Mplex(m).eq.Nplex(n)) then
                  call LVTBDY(Lpex(1,n),Lpl(1,l),Mpex(1,m),u_nmbod)
                  m=m+1
               else
                  call LVTBDY(Lpex(1,n),Lpl(1,l),izero2,u_nmbod)
               endif
               if(Pcond(3,l).ne.90._10 .or. Pcond(4,l).ne.0._10 .or. 
     .          Pcond(9,l).le.0._10) call SUICID(
     . 'UNCONVENTIONAL PULSAR COMPANION ELEMENTS, STOP IN PEXLIC',14)
               call JNITL(SQRT(Pcond(1,l)**3*(Twopi/Pcond(9,l))**2),
     .          Pcond(1,l), Elptn(1,n), 1, Dynpt(1,1,n))
            endif
         end do
         if(Iabs1.gt.0 .and. m.le.Mmpex)
     .    call SUICID('EXO-PLANET NOT INPUT',5)
         Nmpex=n
      endif
      return
      end
