      subroutine RSMAT(s, snames, npnp1)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, mpnp, ne, nerd, nes, npnp1
 
c*** end of declarations inserted by spag
 
 
c
c d. white  subroutine rsmat  september 1974
c sweigh scaling added  paul macneil june, 1978
c skip irrelevant epochs on smat - zmg 2/85
c
c read smear matrix into direct access
c
c   n o t e: s and snames are presumed to share storage in the caller
c
c parameters
      character*8 snames(npnp1, npnp1)
      real*10 s(npnp1, npnp1)
c
c common
      include 'filtda.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
c
c
c local
      real*10 jd1, jd2, jds1, jds2
      character*80 title,planet
c
c read heading records
      read(Smat,50) title, mpnp, nes
   50 format(a80/ 2i5)
      write(Iout, 100) title, mpnp, nes
  100 format('-TITLE OF SMEARING MATRIX DATA SET'/ 1x, a80/
     .       '  NUMBER OF PROCESS NOISE PARAMETERS =', i4,
     .       '  NUMBER OF EPOCHS =', i5)
      read(Smat,150) title, planet
  150 format(a80/ a8)
      write(Iout, 200) title, planet
  200 format('-TITLE OF INTEGRATION USED IN SMEAR'/1x, A80, 10x,
     .       'PLANET = ', a8)
      read(Smat,250) ((snames(i,j),i=1,2), j = 1, mpnp)
  250 format(1x,a8,1x,a8)
      write(Iout, 300) ((snames(i,j),i=1,2), j = 1, mpnp)
  300 format('-NAMES OF PROCESS NOISE PARAMETERS'/(1x,a8,1x,a8))
c
c checks on process params
      if(mpnp.ne.Npnp) goto 350
      do i = 1, Npnp
         if(pnames(1,i).ne.snames(1,i) .or.
     .      pnames(2,i).ne.snames(2,i)) goto 350
         end do
      if(nes.lt.Nepoch-1) goto 360
      goto 400
  350 call SUICID('PARAMETER CHECK IN RSMAT', 6)
  360 call SUICID('EPOCH CHECK IN RSMAT',5)
 
c read into direct access
c loop on number of smear epochs
  400 ne   = Nepoch - 1
      nerd = 0
      do i = 1, ne
         jd1 = (Fep(i) + Fep(i+1))/2._10
         jd2 = (Fep(i+1) + Fep(i+2))/2._10
         do while( .true. )
c read smear interval and S matrix
            read(Smat,405) jds1, jds2, s
  405       format(2f25.16/ (1p,3d25.17))
            nerd = nerd + 1
c
c check epoch alignment
            if(ABS(jd1-jds1).le.Feps(3) .and.
     .         ABS(jd2-jds2).le.Feps(3) ) then
c
c scale by feps(2)
               do j = 1, Npnp
                  do k = 1, Npnp
                     s(k, j) = s(k, j)*Feps(2)*Sweigh(i)
                     end do
                  end do
c
c write to direct access
               write(Lfile, rec = i) s
               go to 450
            else if(jds1.le.jd1-Feps(3) .and.
     .              jds2.le.jd2-Feps(3) ) then
c
c skip last epoch read, & check number remaining
               write(Iout, 410) nerd, Smat, jds1, jds2, jd1, jd2
  410          format('-SKIPPING INPUT EPOCH #', i2, ' ON UNIT',i3/
     .                6x, 'TIMES READ:  ', 2F20.9/
     .                6x, 'TIMES SOUGHT:', 2F20.9)
               if(nes-nerd.le.ne-i) then
                  write(Iout, 420)
  420             format(
     .         '-*** TOO FEW EPOCHS REMAIN ON SMAT INPUT DATA SET ***')
                  goto 360
               end if
            else
c
c fatal epoch mismatch
               write(Iout, 430) Smat, i, jd1, jd2, jds1, jds2
  430          format('-*** NO EPOCH FOUND ON UNIT', i3,
     .                ' FOR SMEAR INTERVAL #', i2, ', FROM', f18.9,
     .                ' TO', f18.9/ 5x, 'LAST EPOCH READ:', 34x,
     .                'FROM', f18.9, ' TO', f18.9)
               goto 360
            end if
            end do
  450    end do
      return
      end
