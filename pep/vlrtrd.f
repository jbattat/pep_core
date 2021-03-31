      subroutine VLRTRD(v1, v2, v3, scal, n, nvlcn)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 DOT, scal
      integer   i, ivel, ixst, n, ncor, nvlcn
 
c*** end of declarations inserted by spag
 
 
c        j.f.chandler   1978 may
c        calculate apparent velocities and time rates due to light time
c        also calculate light time iteration correction factors for ptls
c        input:
c          v1,v2,v3 - true velocities of receive, reflect, send points
c          scal - 1/c in appropriate units
c          n - control integer with packed bits ---
c                1: do receive-reflect (must be on)
c                2: do reflect-send
c                4: do velocity setup in xsitep
c                8: correct velocity vectors that were passed
c          nvlcn - assumed zero for now.  if greater than zero, then
c                  velocities must be offset by corresponding 'center
c                  of mass of solar system' velocities, which should be
c                  available in common somewhere
c
      real*10 v1(3), v2(3), v3(3), vs(3), vt(3)
 
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'rtrdvl.inc'
 
      ncor = mod(n, 4)
      ixst = mod(n, 8)/4
      ivel = mod(n, 16)/8
      if( Jct(67) .gt. 0 ) then
 
c get relative apparent time rates
         Omv1 = 1._10 - DOT(Xsitp0(1,1), v1)*scal
         Omv2 = 1._10 - DOT(Xsitp0(1,1), v2)*scal
         do i = 1, 3
            V2fac(i) = v2(i)*(scal/Omv2)
         end do
         Dt2dt1 = Omv1/Omv2
         if( ncor .ge. 2 ) then
            Opv2   = 1._10 + DOT(Xsitp0(1,2), v2)*scal
            Opv3   = 1._10 + DOT(Xsitp0(1,2), v3)*scal
            Dt3dt2 = Opv2/Opv3
            Dt3dt1 = Dt3dt2*Dt2dt1
            do i = 1, 3
               V3fac(i) = v3(i)*(scal/Opv3)
               Vxfac(i) = scal*(v2(i) - v3(i)*Dt3dt2)
 
            end do
         endif
c
c correction factors set to unity in cmpar3
      else if( ixst .ne. 1 ) then
         return
      endif
c
c get retarded velocities
      if( ixst + ivel .ne. 0 ) then
 
c do first retarded velocity
         do i = 1, 3
            vs(i) = v2(i)*Dt2dt1
            if( ixst .eq. 1 ) Xsitep(i + 3, 1) = v1(i) - vs(i)
            if( ivel .eq. 1 ) v2(i) = vs(i)
         end do
         if( ncor .ge. 2 ) then
 
c do second retarded velocity
            do i = 1, 3
               vt(i) = v3(i)*Dt3dt1
               if( ixst .eq. 1 ) Xsitep(i + 3, 2) = vt(i) - vs(i)
               if( ivel .eq. 1 ) v3(i) = vt(i)
            end do
         endif
      endif
      return
      end
