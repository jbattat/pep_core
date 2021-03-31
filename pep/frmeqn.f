      subroutine FRMEQN(rst,ntop)
 
      implicit none
c
c m.e. ash - february 1970 - subroutine frmeqn
c equinox-equator-latitude logic routine for forming
c normal equations from saved normal equations
c
c arguments
      integer ntop
      real*10 rst(1000)

c array dimensions
      include 'globdefs.inc'
c commons
      include 'anctrl.inc'
      include 'eqenox.inc'
      include 'eqenxm.inc'
      include 'restor.inc'
      include 'wrkcomrs.inc'

c local
      integer   i,j,k

      if(Mumeqn.gt.0) then
         do k = 1,Mumeqn
 
            if(Numeqn.gt.0) then
               do j = 1,Numeqn
                  if(Eqnsit(j).eq.Eqnst1(k)) then
                     if(Eqnser(j).eq.Eqnsr1(k)) then
                        Nrst = Leqn1(j) - 1
                        call FRMBDY(rst,Leqn(1,j),Meqn(1,k),-3,ntop)
                        goto 50
                     endif
                  endif
               end do
            endif
 
            do i = 1,3
               if(Meqn(i,k).gt.0) Nsav = Nsav + 1
            end do
 
   50    end do
      endif
 
      return
      end
