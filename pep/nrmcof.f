      subroutine NRMCOF(deriv, dermd, numpar, nmp2, ict19, lhsflg,
     .                  rhsflg, b)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   l
 
c*** end of declarations inserted by spag
 
 
c
c m.ash   nov 1970    subroutine nrmcof
c increment right side and coefficient matrix of normal equations
c
 
      real*10 deriv(200), dermd(200)
      real*10 b(*)
      integer*2 nmp2, ict19, numpar
      logical   lhsflg, rhsflg
 
      include 'bernum.inc'
      include 'rtside.inc'
 
c local
      integer*4 i, k, ip, irow, kp, idimen
      real*10 wxyz,dermdi
 
c
c if ict(19) = 0, dermd has dimension numpar (no correlated obs.)
c = 1, dermd has dimension nmp2 (correlates obs.)
 
      idimen = numpar
      if(ict19.gt.0 .and. nmp2.gt.numpar) idimen = nmp2
      do i = 3, idimen
         ip = Iptr(i)
         if(ip.gt.0) then
            dermdi=dermd(i)
            wxyz = dermdi*Wobsth
            if(rhsflg) Side(ip) = Side(ip) + deriv(2)*wxyz
            if(lhsflg) then
               Vectk(ip) = Vectk(ip) + dermdi*Wobst1
               irow = (ip*(ip-1))/2
               do k = 3, numpar
                  kp = Iptr(k)
                  if(kp.gt.0 .and. kp.le.ip) then
                     l    = irow + kp
                     b(l) = b(l) + deriv(k)*wxyz
                  endif
                  end do
            endif
         endif
         end do

      return
      end
