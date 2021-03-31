      subroutine PRDSID
 
      implicit none

c m.ash   nov 1970    subroutine prdsid
c create derivative vector for prediction

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bernum.inc'
      include 'fcntrl.inc'
      include 'ktrap.inc'
      include 'nrmmat.inc'
      include 'prdmov.inc'
      include 'rtsidesl.inc'

c local
      real*10 dws1,weight,ws1
      integer   i,ii,ip,j,jp,k,ll
 
      do j = Num1,Num2
         k = j
         if(mod(Jacc,2).ne.0) k = 2
c
c error statistics for old obs-th residuals
         Neas22(k) = Neas22(k) + 1
         Istrik(j) = 2
         if(Deriv(1,j).gt.0.0_10) then
            if(ABS(Deriv(2,j)).lt.Eps(Jacc+j)*Deriv(1,j)) then
               Istrik(j) = 1
               Meas22(k) = Meas22(k) + 1
               weight    = Deriv(2,j)/Deriv(1,j)
               Wobsth    = 1._10/Deriv(1,j)**2
               Erm22(1,1,k) = Erm22(1,1,k) + weight
               Erm22(2,1,k) = Erm22(2,1,k) + ABS(weight)
               Erm22(3,1,k) = Erm22(3,1,k) + weight**2
               Erm22(6,1,k) = Erm22(6,1,k) + Deriv(2,j)**2
               Erm22(7,1,k) = Erm22(7,1,k) + Wobsth
            endif
         endif
c
c calculate predicted observed minus theory
         ws1 = 0._10
         do i = 3,Numpar
            ip = Iptr(i)
            if(ip.gt.0) ws1 = ws1 + Deriv(i,j)*Solut(ip)
         end do
         Obsth(j) = Deriv(2,j) - ws1
c
c calculate uncertainty of prediction
         if(Ict(14).gt.0) then
            ws1 = 0.0_10
            do i = 3,Numpar
               ip = Iptr(i)
               if(ip.gt.0) then
                  do ii = i,Numpar
                     jp = Iptr(ii)
                     if(jp.gt.0) then
                        if(jp.lt.ip) then
                           ll = (ip*(ip-1))/2 + jp
                        else
                           ll = (jp*(jp-1))/2 + ip
                        endif
                        dws1 = B(ll)*Deriv(i,j)*Deriv(ii,j)
                        ws1  = ws1 + dws1
                        if(i.ne.ii) ws1 = ws1 + dws1
                     endif
                  end do
               endif
            end do
            Ucert(j) = SQRT(ABS(ws1))
            if(ws1.lt.0._10) then
               Ucert(j) = -Ucert(j)
               if(Negu.eq.0)
     .              call SUICID('NEGATIVE VARIANCE OF PREDICT',-7)
               Negu = 1
            endif
         endif
c
c
c error statistics for new obs-theory residuals
         if(Istrik(j).le.1) then
            weight = Obsth(j)/Deriv(1,j)
            Erm22(1,2,k) = Erm22(1,2,k) + weight
            Erm22(2,2,k) = Erm22(2,2,k) + ABS(weight)
            Erm22(3,2,k) = Erm22(3,2,k) + weight**2
            Erm22(6,2,k) = Erm22(6,2,k) + Obsth(j)**2
            Erm22(7,2,k) = Erm22(7,2,k) + Wobsth
         endif
 
      end do
 
      return
      end
