      subroutine LIBCHK(ntype,cndx,idd)
 
      implicit none

c
c m.e.ash    dec 1968      subroutine libchk
c check consistency of control vectors on observation library tape
c insert initial conditions from observation library tape

c array dimensions
      include 'globdefs.inc'

c arguments
      integer*4 ntype,idd
      real*10 cndx(u_nmbod)

c commons
      include 'lcntrl.inc'
      include 'mtrapx.inc'
      integer*2 mbdx(u_nmbod,5)
      equivalence (Memx(1),mbdx(1,1))
      include 'namtim.inc'
      include 'empcnd.inc'

c local variables
      integer*4 i,j,ngo,ntype1,ntype2
 
      ntype1 = iabs(ntype) - 4
      ntype2 = min0(5,ntype1 + 4)
      ngo    = 0
 
      do i = 1, 6
         if(Lpl(i,ntype1).gt.0) then
 
            if(ntype.le.0) then
               if(Msbx(i).le.0) goto 100
            else if(mbdx(i,ntype2).le.0) then
               goto 100
            endif
 
            if(ngo.le.0) then
               ngo = 1
               Jdpl0(ntype1) = idd
               do j = 1, 6
                  Pcond(j,ntype1) = cndx(j)
               end do
            endif
         endif
 
      end do
 
      return
 
100   continue
      call SUICID('INCONSISTENCY BETWEEN INPUT L VECTOR AND M VECTOR'//
     .            ' ON OBSERVATION LIBRARY TAPE, STOP LIBCHK  ', 23)
 
      end
