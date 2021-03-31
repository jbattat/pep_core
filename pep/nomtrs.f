      subroutine NOMTRS(side, temp, nptr, npi, wtrans)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, npi, npnpx2
 
c*** end of declarations inserted by spag
 
 
c
c paul macneil january, 1978 subroutine nomtrs
c
c read and transform initial condition offsets,
c correct the propagated residual
 
      include 'fcntrl.inc'
      include 'filtds.inc'
      include 'filtim.inc'
 
      real*10   wtrans(npi, Npnp), temp(npi)
      real*10 side(npi), w(40)
      integer*2 nptr(Npnp)
 
      if(Iterat.le.1 .and. Jct(56).lt.2 .and.
     .      .not.Filflg(1) .and. Fict(8).le.0) return
 
c read in initial condition offsets
      npnpx2 = Npnp*2
      read(Iconof) (w(i), i = 1, npnpx2)
c
c form i(k/k-1)*w
      do i = 1, npi
         temp(i) = 0.0
         do j = 1, Npnp
 
c ic deltas follow ics; deltas used
            temp(i) = temp(i) + wtrans(i, j)*w(j + Npnp)
            end do
         end do
c
c correct u for transformed nominals
      do i = 1, npi
         side(i) = side(i) - temp(i)
         end do
c
      return
      end
