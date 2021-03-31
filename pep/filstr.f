      subroutine FILSTR(b,side,buff,np,ipoch,jd1,jd2,kode)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, irow, j, jread, k
 
c*** end of declarations inserted by spag
 
 
c
c d. white  subroutine filstr  september 1974
c
c z. goldberg   march 1980 -- expanded common filtim
c
c
c         restore filtered normal eqn for one epoch
c         kode=1 => restore forward only (last epoch)
c         kode=2 => restore backward only (first epoch)
c         kode=3 => add forward + backward (intermediate epoch)
c         kode<0 => read only (for skipping)
c                   -1 forward, -2 backward, -3 both
c
c         parameters
      integer*4 np, ipoch, kode
      real*10 b(1),side(np),buff(np)
      real*10 jd1, jd2
c
c common
      include 'cureph.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
c
c local
      real*10 j1, j2
      character*8 msg(2) /'CHECK', 'FILSTR'/
      character*1 msgl(8),nums(4) /'1', '2', '3', '4'/
      equivalence (msg,msgl)
c
c formats
  100 format(/' EPOCH NUMBER', i4, 5x, 'FROM', f13.4, ' TO', f13.4)
c
c see if just reading
      jread = 0
      if(kode.le.0) then
         jread = 1
         kode  = -kode
c
c page heading
      else if(Fict(6).ne.1) then
         call PAGSET('RESTORE FILTERED SNE', 5)
         call NEWPG
      endif
c
c id record
      read(Outsne) k,j1,j2,i
      if(Fict(6).ne.1) then
         call PAGCHK(60,2,0)
         write(Iout,100) i,j1,j2
      endif
 
      Ephnum = i
      Datea  = j1
      Dateb  = j2
 
c checks
      j = 1
      if(k.eq.kode) then
         j = j + 1
         if(ABS(j1-jd1).le.1.0E-6_10) then
            j = j + 1
            if(ABS(j2-jd2).le.1.0E-6_10) then
               j = j + 1
               if(i.eq.ipoch) goto 200
            endif
         endif
      endif
      msgl(7) = nums(j)
      call SUICID(msg,4)
c
c read in first set (forward or backward)
  200 read(Outsne) side
      do i = 1, np
         read(Outsne) buff
         if(jread.ne.1) then
            irow = i*(i - 1)/2
            do j = 1, i
               b(irow + j) = buff(j)
            end do
         endif
      end do
c
c add second set (forward)
      if(kode.eq.3) then
         read(Outsne) buff
         if(jread.ne.1) then
            do i = 1, np
               side(i) = side(i) + buff(i)
            end do
         endif
         do i = 1, np
            read(Outsne) buff
            if(jread.ne.1) then
               irow = i*(i - 1)/2
               do j = 1, i
                  k    = irow + j
                  b(k) = b(k) + buff(j)
               end do
            endif
         end do
      endif
      return
      end
