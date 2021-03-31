      subroutine SNORM5(sidfi,pvcol,iersol)
 
      implicit none
c
c m.ash   jan 1972    subroutine snorm5
c clean up solution, calculate error in solution
c
c arguments
      real*10 sidfi(1000),pvcol(1000)
      integer*4 iersol
c
c common
      include 'fcntrl.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'nrmmat.inc'
      include 'rtside.inc'

c local
      real*10 err,quanti,scl,sum(2)
      integer*4 i,i0,ii,ik,k,k1,l,linexp
      character*80 mesg/'   GAUSS-SEIDEL CLEANUP OF SOLUTION AND WRITING
     . OUT SOLUTION & ERROR QUANTITIES '/
c
c avoid printing on jout
      call PAGSET(-1,0)
c
c clean up of the solution
c
      iersol = 0
      if(Ict(48).gt.0) then
         call GAUSEI(B,sidfi,Side,Nparam,Ict(48),Eps(14),iersol)
         Line = 60
      endif
c
c write out solution and error quantities
c
      if(Fict(6).ne.1) then
         linexp = 4 + (Nparam - 1)/8
         call PAGCHK(60,linexp,0)
         write(Iout,50) Nparam,(sidfi(i),i = 1,Nparam)
   50    format('-THE',i4,
     .          ' RIGHT HAND SIDES OF THE NORMAL EQUATIONS ARE'/
     .          (4x,1p,8D16.8))
         call PAGCHK(60,linexp,0)
         write(Iout,100) Nparam,(Side(i),i = 1,Nparam)
  100    format('-THE',i4,' SOLUTIONS OF THE NORMAL EQUATIONS ARE'/
     .          (4x,1p,8D16.8))
c
c substitute solution into normal equations and subtract right
c hand side
c also sum the relative error
         quanti = 0._10
         err    = 0._10
         ii     = 0
         do i = 1,Nparam
            i0     = ii
            ii     = ii + i
            sum(1) = 0._10
            sum(2) = 0._10
            k1     = 1
            ik     = i0
            do k = 1,Nparam
               ik = ik + k1
               if(k.ge.i) k1 = k
               call XLOAD8(Side(k))
               call XMUL8(B(ik))
               call XADD(sum)
               call XSTORE(sum)
               end do
            call XSUB8(sidfi(i))
            call STORND(sum)
            scl = B(ii)
            if(scl.eq.0._10) scl = 1E-50_10
            err    = err + ABS(sum(1)/scl)
            quanti = quanti + ABS(Side(i))
            pvcol(i) = sum(1)
            end do
         linexp = 4 + (Nparam - 1)/16
         call PAGCHK(60,linexp,0)
         write(Iout,150) (pvcol(i),i = 1,Nparam)
  150    format('-THE RESULT OF SUSTITUTING THE SOLUTIONS IN THE NOR',
     .   'MAL EQUATIONS AND SUBTRACTING THE RIGHT HAND SIDES IS'/
     .    (4x,1p,16D8.1))
c
c divide above differences by non-zero right sides
         do i = 1,Nparam
            if(sidfi(i).ne.0._10) then
               pvcol(i) = pvcol(i)/sidfi(i)
            else
               pvcol(i) = 0._10
            endif
            end do
         call PAGCHK(58,linexp,0)
         write(Iout,200) (pvcol(i),i = 1,Nparam)
  200    format(
     .'-THE RESULT OF DIVIDING THESE DIFFERENCES BY THE RIGHT SIDES OF T
     .HE NORMAL EQUATIONS (FOR THE NON-ZERO RIGHT SIDES) IS'/
     .  (4x,1p,16D8.1))
c
c write out relative error of the solution
         write(Iout,250) err,quanti
  250    format('0 ERROR',1pd12.2,'  RELATIVE TO',1pd12.2)
         Line = Line + 2
c
c write out timer information
         k = 41
 
c decide whether to print first half of message
         if(Ict(48).gt.0) k = 1
         l = (81 - k)/4
         call TIMRIT(mesg(k:80),l)
      endif
 
      return
      end
