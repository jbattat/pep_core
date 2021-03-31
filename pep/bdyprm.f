      subroutine BDYPRM(lvec,cond,prefix,nplic,k,names,iskale,xnom)
 
      implicit none
c
c d. white  april 1973  subroutine bdyprm
c
c generates names for ic's and con's

c array dimensions
      include 'globdefs.inc'
c
c parameters
      character*8 prefix,names(2,1)
      real*10 cond(u_nmbod),iskale(1),xnom(1)
      integer*2 lvec(u_nmbod)

c external function
      real*10 SKALE

c locals
      integer*4 i,j,k,l,nplic
      character*8 orbels(6) /'A', 'E', 'INC', 'ASC', 'PER', 'ANOM'/
      character*8 ardbrd(6) /'R', 'RA', 'DECL', 'V*', 'FLAZ', 'FLANG'/
      character*8 cartes(6) /'X', 'Y', 'Z' ,'X''', 'Y''', 'Z'''/
      integer*2 orbic(6)/0, 0, 1, 1, 1, 1/
      integer*2 ardic(6)/0, 1, 1, 0, 1, 1/
      integer*2 plnprm(u_nmbod-6)/3, 4, 4, 4, 5, 5, 4, 6, 5, 14*0, 2/
      integer*2 prbprm(u_nmbod-6)/22*0, 4, 2/
      integer*2 mnprm(u_nmbod-6)/3, 22*0, 2/
      integer*2 embprm(u_nmbod-6)/23*0, 2/
      character*8 ccond /'COND( )'/, CON1 /'CON( )'/, CON2 /'CON(  )'/
      character*1 nums(10)/'1','2','3','4','5','6','7','8','9','0'/
      character*8 temp
      character*1 tem(8)
      equivalence(tem, temp)
      integer*4 nplnt, icnd
c
c initial conditions
      icnd  = nplic/100
      nplnt = nplic - icnd*100
      icnd  = iabs(icnd) - 2
      do i = 1, 6
         if(lvec(i).ne.0) then
            k = k + 1
            names(1, k) = prefix
            xnom(k)     = cond(i)
            if(nplnt.lt.0) then
               temp   = ccond
               tem(6) = nums(i)
               names(2, k) = temp
            else if(icnd.lt.0) then
               names(2, k) = cartes(i)
            else if(icnd.eq.2) then
               names(2, k) = ardbrd(i)
               j = ardic(i)
               iskale(k) = SKALE(80 + j)
               xnom(k)   = xnom(k)/iskale(k)
            else
               names(2, k) = orbels(i)
               j = orbic(i)
               iskale(k) = SKALE(80 + j)
               xnom(k)   = xnom(k)/iskale(k)
            endif
         endif
      end do
c
c other body parameters
      do i = 7, u_nmbod
         if(lvec(i).eq.0) return
         k = k + 1
         l = lvec(i)
         names(1, k) = prefix
         j = plnprm(l)
         if(nplnt.le.0 .or. nplnt.gt.30) j   = prbprm(l)
         if(nplnt.eq.-3 .or. nplnt.eq.-10) j = 0
         if(nplnt.eq.3) j  = embprm(l)
         if(nplnt.eq.10) j = mnprm(l)
         iskale(k) = SKALE(80 + j)
         xnom(k)   = cond(l + 6)/iskale(k)
         if(l.ge.10) then
            temp = con2
            call EBCDI(l, tem(5), 2)
         else
            temp   = con1
            tem(5) = nums(l)
         endif
         names(2, k) = temp
      end do
      return
      end
