      subroutine SMLROT(prec,dpranh)
 
      implicit none
c
c m.e.ash    may 1971    subroutine smlrot
c multiply precession matrix by small rotations
c
c parameters
      real*10 prec(3,3),dpranh(3),dprmat(3,3)
c local
      real*10 temp(3,3),dum(3)
      integer*4 i,j

      do i=1,3
         dum(1)=prec(i,1)
     .            + (-prec(i,2)*dpranh(3) + prec(i,3)*dpranh(2))
         dum(2)=prec(i,2)
     .            + (prec(i,1)*dpranh(3) - prec(i,3)*dpranh(1))
         dum(3)=prec(i,3)
     .            + (-prec(i,1)*dpranh(2) + prec(i,2)*dpranh(1))
         do j=1,3
            prec(i,j)=dum(j)
         end do
      end do
      return

      entry FULROT(prec,dprmat)
c multiply precession matrix by possibly small rotations expressed
c in the form of a full matrix
c alternate entry point added 2013 Feb 22 - J.F.Chandler
      do i=1,3
         do j=1,3
            temp(i,j)=prec(i,j)
         end do
      end do
      call PRODCT(temp,dprmat,prec,3,3,3)
      return
      end
