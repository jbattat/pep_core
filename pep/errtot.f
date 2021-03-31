      subroutine ERRTOT(ncodg)
 
      implicit none

c m.e.ash   april 1967    subroutine errtot
c printout total error analysis
      integer*2 ncodg
c ncodg=0, error analysis for all observations. otherwise,
c error analysis of observations of type ncodg made of body

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdat.inc'
      character*8 pname
      equivalence (Comcon(127),pname)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      real*10 erstf(3,3)
      equivalence (erstf(1,1),Erquan(17))
      include 'obstuf.inc'
      include 'statsrad.inc'
      include 'timstf.inc'
c
c           nit(1) = page number at start of processing all observations
c           nit(2) = page number at start of processing of observations
c                    of type ncodg made of body aplnt(klam)
c           nit(3) = number of processed records for all observations
c           nit(4) = number of non-processed records for all
c                    observations
c           nit(5) = number of processed records for observations
c                    of type ncodg made of body aplnt(klam)
c           nit(6) = number of non-processed records for observations
c                    of type ncodg made of body aplnt(klam)
c           nit(7) = number of measurements deleted in least squares
c                    and error analysis for all observations
c           nit(8) = number of measurements included in least squares
c                    and error analysis for all observations
c           nit(9) = number of 1st measurements deleted in least squares
c                    and error analysis for observations
c                    of type ncodg made of body aplnt(klam)
c           nit(10)= number of 2nd measurements deleted in least squares
c                    and error analysis for observations
c                    of type ncodg made of body aplnt(klam)
c           nit(11)=number of 1st measurements included in least squares
c                    and error analysis for observations
c                    of type ncodg made of body aplnt(klam)
c           nit(12)=number of 2nd measurements included in least squares
c                    and error analysis for observations
c                    of type ncodg made of body aplnt(klam)
c           nit(13)= real timer at start of processing all observations
c           nit(14)= task timer at start of processing all observations
c           nit(15)= real timer at start of processing observations
c                    of type ncodg made of body aplnt(klam)
c           nit(16)= task timer at start of processing observations
c                    of type ncodg made of body aplnt(klam)
c
c local
      real*10 fnobs,ww(2)
      integer*4 i,j,k,king,kong,ler1,ler2,ling,mpage,nn2,nnnn(2,2)
      real*4    a(2,5)
      character*4 form(9)/'(2(A','32,','I10,','I14/','),(A',
     .    '32,1','P2E1','4.5)',')   '/
      character*4 form1(4)/'   /','I14/','P1E1','P2E1'/
      character*32 words(2)/
     .  ' NUMBER OF MEASUREMENTS DELETED ',
     .  ' NUMBER OF MEASUREMENTS INCLUDED'/
      character*32 wores(5)/
     .  '          AVERAGE (OBS-TH)/ERROR',
     .  '       AVERAGE ABS(OBS-TH)/ERROR',
     .  ' ROOT MEAN SQUARE (OBS-TH)/ERROR',
     .  '     AVERAGE ((OBS-TH)/ERROR)**2',
     .  '         SUM ((OBS-TH)/ERROR)**2'/
c
c new page if necessary
c increment line count in advance
      call PAGCHK(58,14,0)
      if(ncodg.gt.0) then
c
c for observations of type ncodg of body
         king    = 5
         ling    = 10
         ler1    = 7
         ler2    = 9
         kong    = 0
         nn2     = 2
         form(4) = form1(2)
         form(7) = form1(4)
         mpage   = Npage - Nit(2)
         write(Iout,50) Svtpbs,pname,mpage,Nit(2),Labsv
   50    format('-ERROR ANALYSIS FOR ',4A4,' OBSERVATIONS OF ',a8,
     .          ' PROCESSED DURING THE LAST',i5,
     .          ' PAGES STARTING ON PAGE',i5/30x,10A4)
      else
c
c for all observations
         king    = 3
         ling    = 7
         ler1    = 10
         ler2    = 12
         kong    = 2
         nn2     = 1
         form(4) = form1(1)
         form(7) = form1(3)
         mpage   = Npage - Nit(1)
         write(Iout,100) mpage,Nit(1)
  100    format('-ERROR ANALYSIS FOR ALL OBSERVATIONS PROCESSED IN THE',
     .   ' COMPAR LINK DURING THE PAST',i5,' PAGES STARTING ON PAGE',
     .    i5/)
      endif
c
c error analysis
      do i = 1,nn2
         kong = kong + 1
         nnnn(i,1) = Nit(king + 3 + i)
         nnnn(i,2) = Nit(ling + i)
         fnobs = nnnn(i,2)
         if(fnobs.gt.0._10) then
            a(i,1) = erstf(1,kong)/fnobs
            a(i,2) = erstf(2,kong)/fnobs
            ww(1)  = erstf(3,kong)/fnobs
            a(i,3) = SQRT(ww(1))
            a(i,4) = ww(1)
            a(i,5) = erstf(3,kong)
         else
            do j = 1,5
               a(i,j) = 0.0E0
               end do
         endif
         end do
      write(Iout,form) (words(j),(nnnn(i,j),i=1,nn2),j=1,2),
     .                  (wores(j),(a(i,j),i=1,nn2),j=1,5)
c
c statistics for error reads
      write(Iout,200) (stat(i),i = ler1,ler2)
  200 format(
     .    ' NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (',
     .    i5,',',i5,',',i5,')')
      do i = ler1,ler2
         stat(i) = 0
         end do
c
c printout time elapsed
      write(Iout,300) Nit(king),Nit(king + 1),Nparam
  300 format(/i8,' +',i5,' OBSERVATION RECORDS PROCESSED (NPARAM=',
     .       i4,')')
      call TIMRTC('   PROCESSING OBSERVATIONS  ',7,Nit(king+10))
c
c re-initialize
      if(ncodg.gt.0) then
         Nit(2)  = Npage
         Nit(15) = Ireal0
         Nit(16) = Itotsk
         do i = 3,4
            Nit(i)   = Nit(i) + Nit(i + 2)
            Nit(i+2) = 0
            end do
         k = 7
         do i = 7,8
            k = k + 2
            Nit(i)   = Nit(i) + Nit(k) + Nit(k + 1)
            Nit(k)   = 0
            Nit(k+1) = 0
            end do
         do i = 1,3
            erstf(i,3) = erstf(i,1) + erstf(i,2) + erstf(i,3)
            erstf(i,1) = 0._10
            erstf(i,2) = 0._10
            end do
      endif
      return
      end
