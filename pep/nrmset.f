      subroutine NRMSET(kall)
 
      implicit none
c
c m.e.ash    feb 1967    subroutine nrmset
c the right side and the lower diagonal half of the symmetric
c coefficient matrix of the normal equations are initialized
c
c parameters
      integer*4 kall
c kall - if 0 then clear statistics and equations
c if 1 then clear only equations
c
c commons
      include 'fcntrl.inc'
      include 'nrmmat.inc'
      include 'numnum.inc'
      include 'rtside.inc'
c local variables
      integer*4 i,length
c
c initialize to zero
      if(kall.le.0) then
         Measmt = 0
         do i = 1, 3
            Ermeas(i) = 0.0_10
         end do
         Sumaps = 0._10
         Sumzns = 0._10
         Sumzsm = 0._10
      endif
      length = (Nparam*(Nparam+1))/2
      if(length.gt.Nrmsiz) call SUICID(
     .' SIZE OF NORMAL EQUATIONS EXCEEDS STORAGE IN /NRMMAT/ LABELED COM
     .MON, STOP IN NRMSET',21)
 
      call ZFILL(Side,16*Nparam)
      call ZFILL(B,16*length)
      call ZFILL(Vectk,16*Nparam)
 
      return
      end
