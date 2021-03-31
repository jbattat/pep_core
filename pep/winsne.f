      subroutine WINSNE(b,nms,side,edriv,mparam,lhsflg,rhsflg,measmt,
     .                  ermeas)
 
      implicit none
c
c d. white  april 1974  subroutine winsne
c
c z. goldberg   march 1980 -- expanded common filtim
c
c
c read direct access and rewrite into filter input sne
c
c parameters
      integer*4 mparam,measmt
      real*10 b(1),side(mparam),edriv(mparam),ermeas(3)
      character*8 nms(1)
      logical*4 lhsflg,rhsflg
c note that b and nms are work areas that may share storage
c
c common
      include 'aprtbf.inc'
      include 'fcntrl.inc'
      include 'filmis.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
c
c local
      real*10 jd1,jd2
      integer   i,ipoch,j,jrow,k,kode,nsave
      character*4 wside(3)/'RHS','LHS','R&L'/
      integer*4 null/0/
c
c write heading on iout
      write(Iout,100) Iterat,Heding,Date,Npage
  100 format('1CREATE FILTER INPUT SNE (ITERAT=',i2,')',5x,18A4,
     .       1x,2A4,' PAGE',i5)
      Npage = Npage+1
      write(Iout,200)
  200 format(/' ','SIDE SAVED',3x,'EPOCH START',4x,'EPOCH END',
     .       6x,'NUMBER OF OBS')
c
c write title to insne (heading and date from this run)
      rewind Insne
      write(Insne) Heding,Date,Nparam,Nepoch
c
c write residual summary to insne
      write(Insne) measmt,ermeas
c
c write names
      rewind Ibuf1
      j = 2*Nparam
      read(Ibuf1) (nms(i),i = 1,j)
      write(Insne) (nms(i),i = 1,j)
      rewind Ibuf1
c
c set up for loop
      kode = 0
      if(lhsflg) kode = kode+2
      if(rhsflg) kode = kode+1
c
c loop through epochs
      do i = 1,Nepoch
c
c read in from direct access
         call RFILDA(b,side,i-1,Nparam,edriv,lhsflg,rhsflg)
c
c write out data for this epoch
         jd1 = Fep(i)
         jd2 = Fep(i+1)
         write(Iout,300) wside(kode),jd1,jd2,Nwrite(i)
  300    format(1x,3x,a4,3x,2F15.4,10x,i6)
         write(Insne) kode,jd1,jd2,i,Nwrite(i)
         if(Nwrite(i).gt.0) then
            nsave = mparam
            if(rhsflg) write(Insne) side
            if(lhsflg) then
               do j = 1,nsave
                  edriv(j) = 0.0_10
               end do
               do j = 1,nsave
                  jrow = j*(j-1)/2
                  do k = 1,j
                     edriv(k) = b(jrow+k)
                  end do
                  write(Insne) edriv
               end do
            endif
         endif
      end do
c
c write reverse data if needed
      if(Fict(2).ne.0) then
         write(Insne) null
         write(Iout,500)
  500    format(/' BACKWARDS DATA IS BEING WRITTEN')
         jd1 = jd2
         do i = 1,Nepoch
            ipoch = Nepoch+1-i
c
c read direct access
            call RFILDA(b,side,ipoch-1,mparam,edriv,lhsflg,
     .                  rhsflg)
c
c data
            jd1 = Fep(ipoch)
            jd2 = Fep(ipoch+1)
            write(Insne) kode,jd1,jd2,ipoch,Nwrite(ipoch)
            if(Nwrite(ipoch).gt.0) then
               if(rhsflg) write(Insne) side
               if(lhsflg) then
                  do j = 1,Nparam
                     edriv(j) = 0.0_10
                  end do
                  do j = 1,Nparam
                     jrow = j*(j-1)/2
                     do k = 1,j
                        edriv(k) = b(jrow+k)
                     end do
                     write(Insne) edriv
                  end do
               endif
            endif
         end do
      endif
c
c write totals
      write(Iout,400) Twrite,Nor
  400 format(/' ','TOTAL OBS INCLUDED',i6,' OUT OF RANGE',i5)
      return
      end
