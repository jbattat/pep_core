      subroutine SCRAAD(jd,fract)
 
      implicit none
 
c arguments
      integer*4 jd
      real*10 fract

c     m.e.ash, l.friedman-  may 1970  --subroutine scraad
c     observing body  tape is read either forward or backward in time
c     in the variable tabular interval mode
c           jd   = julian day no.
c           fract= fraction of day from midnight
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'scdta.inc'
      include 'scdtavtp.inc'
      include 'trpcom.inc'
c
c local variables
      real*10 xx,yy,zz,zz1
      integer*4 i,iord,j,jord,k,l,m,m1,m2,mm,n,numpt,numpt1

c see if the correct records are in storage
  100 xx = jd - Jdsc(2)
      xx = xx + (fract - Fsc(2))
 
c xx= given time - time of record 2
      yy = jd - Jdsc(5)
      yy = yy + (fract - Fsc(5))
 
c yy= given time - time of record 5
      if(xx*yy.le.0._10) return
c
c see if correct records are ahead on tape
      if(ABS(xx).lt.ABS(yy)) then
c
c correct records are behind on tape
         mm = 1
         zz = xx + Gc2(1)
 
c zz= given time - time of record 1
         n = 8
         if(zz*xx.lt.0.0_10) then
            zz1 = zz + Gc1(1)
 
c zz1= given time - time of record 0
            n = 9
            if(zz*zz1.lt.0.0_10) then
c
c correct records are more than two behind on tape
               rewind Jsc
               write(Iout,110) jd,fract,Jsc
  110          format(i17,f17.13,' TOO FAR BEHIND ON TAPE', i3,
     .                ' FOR SATELLITE-PROBE, REWIND TAPE')
               do i = 1, 2
                  read(Jsc,err=120)
                  goto 130
  120             read(Jsc)
  130          end do
               goto 200
            endif
         endif
c
c backspace records
         do i = 1, n
            backspace Jsc
         end do
      else
c
c see if we move 5 tabular points in storage
         zz = yy - Gc1(6)
 
c zz= given time - time of record 6
         if(zz*yy.gt.0._10) then
c
c see if we move 4 tabular points in storage
            zz1 = zz - Gc2(6)
 
c zz1= given time - time of record 7
            if(zz*zz1.gt.0._10) then
c
c move  3  tabular point in storage preparatory to read rec.
               m1 = 3
               m2 = 3
               mm = 5
            else
               m1 = 4
               m2 = 2
               mm = 5
            endif
         else
            m1 = 5
            m2 = 1
            mm = 6
         endif
         do while( .true. )
c
c move tabular points in storage
            do l = 1, m1
               m = l + m2
               Jdsc(l)   = Jdsc(m)
               Fsc(l)    = Fsc(m)
               Iorcer(l) = Iorcer(m)
               Jorcer(l) = Jorcer(m)
               Numprc(l) = Numprc(m)
               Gc1(l)    = Gc1(m)
               Gc2(l)    = Gc2(m)
               iord      = Iorcer(l)
               do j = 1, iord
                  do i = 1, 3
                     Sprc(l,j,i) = Sprc(m,j,i)
                  end do
               end do
               if(Numprc(l).gt.0) then
                  numpt = Numprc(l)
                  jord  = Jorcer(l)
                  do k = 1, numpt
                     do j = 1, jord
                        do i = 1, 3
                           Dsprc(l,j,i,k) = Dsprc(m,j,i,k)
                        end do
                     end do
                  end do
               endif
            end do
            if(m1.gt.3) goto 300
            m2 = 1
c
c search ahead on tape
            read(Jsc) Jdsc(4),Fsc(4),iord,jord,Numprc(4),numpt1,
     .                ((Sprc(4,j,i),i=1,3),j = 1,iord),
     .                (((Dsprc(4,j,i,k),i=1,3),j=1,jord),k = 1,
     .                numpt1), Gc1(4),Gc2(4)
            Iorcer(4) = iord
            Jorcer(4) = jord
 
            zz  = zz1
            zz1 = zz - Gc2(4)
            if(zz*zz1.le.0) goto 300
         end do
      endif
c
c read variable tabular points into storage
      entry SCRAD1(jd)
  200 mm   = 1
  300 do l = mm, 6
         read(Jsc) Jdsc(l),Fsc(l),iord,jord,Numprc(l),numpt1,
     .             ((Sprc(l,j,i),i=1,3),j = 1,iord),
     .             (((Dsprc(l,j,i,k),i=1,3),j=1,jord),k = 1,numpt1),
     .             Gc1(l),Gc2(l)
         Iorcer(l) = iord
         Jorcer(l) = jord
      end do
      if(jd.gt.0) goto 100
c
c get a sharper criterion for saying point is before start
      if(Idirsc.le.0) then
         Frsc2 = Fsc(2)
         Jdsc2 = Jdsc(2)
      else
         Frsc1 = Fsc(2)
         Jdsc1 = Jdsc(2)
      endif
 
      return
      end
