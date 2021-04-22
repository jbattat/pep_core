      subroutine SSREED(jd,npl,itp)

      implicit none

c
c j.f.chandler   jan 2019   subroutine ssreed
c planet tape is read either forward or backward in time

c parameters
      integer*4 jd,itp
      integer*2 npl

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'tapdta.inc'
      include 'tapdtp.inc'
      include 'trpcom.inc'

c local variables
      integer   i,ivl,j,k,l,l1,l2,mm,n
c
      if(jd.le.Jdss1(npl)) goto 900
      if(Jdss2(npl).le.jd) goto 900
  100 n = Jdss(2,npl) - jd
      if(Idirz(npl+3).lt.0) then

         if(n.lt.0) then
c
c backspace logic
            n = -n
            goto 200
         else if(n.eq.0) then
            return
         else
            n = Jdss(3,npl) - jd
         endif

      else if(n.lt.0) then
         n = jd - Jdss(3,npl)
      else if(n.eq.0) then
         return
      else
         goto 200
      endif
      if(n.lt.0) return
      if(n.ne.0) then

         n = n/Intss5(npl) - 1
         if(n.lt.0) then
         else if(n.eq.0) then

            Jdss(1,npl)  = Jdss(3,npl)
            Issvel(1,npl)= Issvel(3,npl)
            Fss(1,npl)   = Fss(3,npl)
            ivl = Issvel(1,npl)
            do k = 1,5
               n = k + 10
               do j = 1,Iparss(npl)
                  do i = 1,ivl
                     Planss(i,j,k,npl) = Planss(i,j,n,npl)
                  end do
               end do
            end do
            mm = 2
            goto 400
         else
            n = n - 1
            if(n.gt.0) then
               do i = 1,n
                  read(itp)
               end do
            endif
            goto 300
         endif
      endif

      Jdss(1,npl)   = Jdss(2,npl)
      Issvel(1,npl) = Issvel(2,npl)
      Fss(1,npl)    = Fss(2,npl)
      Jdss(2,npl)   = Jdss(3,npl)
      Issvel(2,npl) = Issvel(3,npl)
      Fss(2,npl)    = Fss(3,npl)
      ivl = max0(Issvel(1,npl),Issvel(2,npl))
      do k = 1,10
         n = k + 5
         do j = 1,Iparss(npl)
            do i = 1,ivl
               Planss(i,j,k,npl) = Planss(i,j,n,npl)
            end do
         end do
      end do
      mm = 3
      goto 400
  200 n = (n - 1)/Intss5(npl) + 4
      do i = 1,n
         backspace itp
      end do
      goto 300
c
c
c read planet data into storage
      entry SSRED1(jd,npl,itp)
  300 mm     = 1
  400 do l   = mm,3
         l2  = l*5
         l1  = l2 - 4
         read(itp) Jdss(l,npl),Fss(l,npl),ivl,
     .    (((Planss(i,j,k,npl),i=1,ivl),j=1,Iparss(npl)),k = l1,l2)
         Issvel(l,npl) = ivl
      end do
      if(jd.gt.0) goto 100
      return
c
  900 write(Iout,1000) jd,itp,npl,Jdss1(npl),Jdss2(npl)
 1000 format(i17,' IS NOT ON DATA SET',i3,' FOR PLANET',i3,
     . ' (',i7,'-',i7,')')
      jd = 0
      return
      end
