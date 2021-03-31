      subroutine SNORM6(tempa,tempb,b,scale)

      implicit none
c
c r.reasenberg  feb 73  subroutine snorm6
c version of 17 april 1975  r.reasenberg
c
c
c        this subroutine provides the services required to calculate
c        and display the series-by-series adjustments.  it is called
c        at different entry points by - - nrmadd, snorml, solctl,
c        and soluse.
c
      include 'fcntrl.inc'

c parameters to entry SNORM6
      real*10 tempa(Nparam),tempb(Nparam),b(10),scale(Nparam)
c parameters to entry SETSD
      integer*4 isdset,isolset,inmset
c parameters to entry SVSD
      integer*4 ntp,nsq
      real*10 side(Nparam)
c parameters to entry ADDSD (also, side as above)
      real*10 save(Nparam)
c parameters to entry ADJSER
      real*10 bb(Nparam,23),solut(Nparam),sigma(Nparam)
      integer*4 iout
      character*8 nms(2,Nparam)
c typically, b and bb share storage

c local
      integer*4 sdset,solset,nmset
      character*8 blnk/'        '/
      real*10 solutn,sum
      integer   j,k,k1,kl,klz,l,m,n,nl,nneg,nparax,nums
      integer*2 ntpv(20),nsqv(20),numv(20)
      integer*4 nlist/8/, mline/55/, nid/3/, num/0/, nv(3)
      real*10 rms,rmsn
      real*10 fnum
c
c - - - - - - - - - - - - - - - - - - -
      if(num.gt.0) then
         rewind sdset
         rewind solset
c find ser-by-ser solution - - 1) read rhs from data set sdset,
c 2) multiply by the covariance matrix, and 3) write onto data
c set solset
c data sets defined at entry setsd
         do j = 1,num

            read(sdset,end=20,err=20) nid,nv,tempa

            if(nv(3).eq.j) then

               klz = 0
               do k = 1,Nparam
                  sum = 0.0
                  do l = 1,k
                     kl  = klz + l
                     sum = sum + b(kl)*tempa(l)*scale(l)
                  end do
                  if(k.lt.Nparam) then
                     klz = klz + k
                     k1  = k + 1
                     do l = k1,Nparam
                        kl  = kl + l - 1
                        sum = sum + b(kl)*tempa(l)*scale(l)
                     end do
                  endif
                  tempb(k) = sum*scale(k)
               end do
               write(solset) nid,nv,tempb
               goto 50
            endif

   20       write(6,40)
   40       format('0ERROR IN SNORM6')
            num = 0
            return
   50    end do
         endfile solset
         rewind sdset
         rewind solset
      endif
      return



c - - - - - - - - - - - - - - - - - - -
      entry SETSD(isdset,isolset,inmset)
      sdset  = isdset
      solset = isolset
      nmset  = inmset
      num    = 0
      rewind sdset
      return



c - - - - - - - - - - - - - - - - - - -
      entry SVSD(ntp,nsq,side)
      num = num + 1
      write(sdset) nid,ntp,nsq,num,side
      do j = 1,Nparam
         side(j) = 0.0
      end do
      return



c - - - - - - - - - - - - - - - - - - -
      entry ADDSD(side,save)
      if(num.gt.0) then
         endfile sdset
         rewind sdset
         do j = 1,num
            read(sdset,end=200,err=200) nid,nv,save
            if(nv(3).ne.j) goto 200
            do k = 1,Nparam
               side(k) = side(k) + save(k)
            end do
         end do
         rewind sdset
      else
         write(6,100)
  100    format('0ADDSD E.P. OF SNORM6 CALLED WITH NUM.LE.0')
      endif
      return

  200 write(6,300)
  300 format('0ERROR IN ADDSD  - -  STOP')
      stop



c - - - - - - - - - - - - - - - - - - -
      entry ADJSER(solut,iout,bb,nms,sigma)
      if(num.le.0) return
      if(nmset.le.0) goto 600

c get the list of names from data set nmset=ibuf1
      rewind nmset
      read(nmset,end=400,err=400) nms

c read scale factors, using trick to avoid nominals
      read(nmset,end=400,err=400)
     .     (bb(m,23),bb(m,23),m = 1,Nparam)
      rewind nmset
      goto 700

  400 write(iout,500) nmset
  500 format('0ERROR READING LIST OF NAMES IN SUBROUTINE SNORM6',
     .       ' AT ENTRY ADJSER FORM DATASET',i5)
      rewind nmset
  600 do m = 1,Nparam
         bb(m,23) = 1._10
         nms(1,m) = blnk
         nms(2,m) = blnk
      end do

  700 do m = 1,Nparam
         bb(m,22) = 0._10
         bb(m,21) = 0._10
      end do
      nums = num

c find maximum number of adjustments to be printed  nparax
      nparax = Ict(18)/10
      nparax = iabs(nparax)
      if(nparax.eq.0) nparax = Nparam
      if(nparax.gt.Nparam) nparax = Nparam
      do while(.true.)
c start of loop for reading solutions, a maximun of nlist
c per batch
         nl   = min0(nlist,nums)
         Line = 1000
         do n = 1,nl
            read(solset,end=800,err=800) nid,nv,(bb(m,n),m=1,nparax)
            ntpv(n) = nv(1)
            nsqv(n) = nv(2)
            numv(n) = nv(3)
            do m = 1,nparax
               bb(m,n)  = bb(m,n)*bb(m,23)
               bb(m,21) = bb(m,21) + bb(m,n)
               bb(m,22) = bb(m,22) + bb(m,n)**2
            end do
         end do

         do m = 1,nparax
            if(Line.ge.mline) then
               Line = 3

c write out page heading
               write(iout,710) Iterat,Heding,Date,Npage
  710          format('1LEAST SQUARES ANALYSIS  (ITERAT=',i3,')',4x,
     .                18A4,1x,2A4,' PAGE',i5)
  720          format(20x,8(5x,i3,i4))
               write(iout,720) (numv(n),ntpv(n),n = 1,nl)
               write(iout,730) (nsqv(n),n = 1,nl)
  730          format(20x,8(5x,i7))
               Npage = Npage + 1
            endif
            Line = Line + 1
  740       format(i5,2(1x,a8),t127,i5,t25,1p,8D12.4)
            write(iout,740) m,nms(1,m),nms(2,m),m,
     .                       (bb(m,n),n = 1,nl)
         end do
         nums = nums - nl
         if(nums.le.0) then
c end of loop reading and writing solutions and forming
c solution statistics
            Line = 1000
            fnum = num
            nneg = 0

            do m = 1,nparax
               if(Line.ge.mline) then

c write out page heading
                  write(iout,710) Iterat,Heding,Date,Npage
                  Npage = Npage + 1
                  write(iout,750)
  750             format(4x,'M',25x,'SOLUTION',7x,
     .                   'SUM PARTS - SOLUTION',8x,'R M S',11x,
     .                   'R M S /SIGMA',/)
                  Line = 4
               endif
               Line = Line + 1
               bb(m,20) = bb(m,21) - solut(m)
               rms = bb(m,22) - bb(m,21)**2/fnum
               if(rms.lt.0) then
                  nneg = nneg + 1
               else if(rms.ne.0) then
                  rms = SQRT(rms)
               endif
               rmsn   = rms/sigma(m)/bb(m,23)
               solutn = solut(m)*bb(m,23)
               write(iout,760) m,nms(1,m),nms(2,m),solutn,
     .                          bb(m,20),rms,rmsn,m
  760          format(i5,2(1x,a8),1p,4D20.10,i5)
            end do
            if(nneg.gt.0) write(iout,780) nneg
  780       format('0THERE ARE',i5,' ERRORS IN RMS')
            rewind solset
            return
         endif
      end do

  800 write(6,900)
  900 format('0ERROR IN ADJSER')
      return

      end
