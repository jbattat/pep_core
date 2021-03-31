      subroutine DTPUN(ncard)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, ncard, ndt, numdt1, numdt2
 
c*** end of declarations inserted by spag
 
 
c
c f.a. kreimendahl  december 1979  subroutine dtpun
c punch whole dt table if any are adjusted
c
 
      include 'inodta.inc'
      include 'dtparm.inc'
 
      character*4 blank/'    '/,qdt/'DT('/
 
      if( Numdt .gt. 0 ) then
 
c test if any dt's are adjusted
         do i = 1, Numdt
            if( Ldt(i) .gt. 0 ) go to 100
 
c test if dt table includes wobble
            if( Jddt0 .le. 0 ) then
               if( Ldt(i+200) .gt. 0 .or. Ldt(i+400) .gt. 0 ) go to 100
            endif
         end do
      endif
      return
 
c punch l-vector for ut1 dt's
  100 write(Ipunch, 200) Numdt, Jddt0, (blank, i, Ldt(i), i = 1, Numdt)
  200 format(' NUMDT=', i3, ', JDDT0=',
     .       i7/(6(1A1,'LDT(',i3,')=',i1,',')))
      ncard  = ncard + 2 + (Numdt - 1)/6
      numdt1 = Numdt + 200
      numdt2 = Numdt + 400
 
c punch l-vector for xwob & ywob dt's
      write(Ipunch, 300) (blank, i, Ldt(i), i = 201, numdt1)
      write(Ipunch, 300) (blank, i, Ldt(i), i = 401, numdt2)
  300 format(6(1A1,'LDT(',i3,')=',i1,','))
      ncard = ncard + 2*((Numdt-1)/6 + 1)
      ndt   = Numdt
      if( Jddt0 .le. 0 ) ndt = 400 + Numdt
      do i = 1, Numdt
 
c punch ut1, xwob & ywob dt's
         write(Ipunch, 350) i, Jddt(i),
     .                      (qdt, j, Dt(j), j = i, ndt, 200)
  350    format(' JDDT(', i3, ')=', i7, ',',
     .          (t22,2(a3,i3,')=',1pe13.6,',')))
      end do
      ncard = ncard + Numdt*(1 + ndt/400)
 
      return
      end
