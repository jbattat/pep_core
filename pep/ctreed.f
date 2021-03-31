      subroutine CTREED(jdt,frt,tctat,ctmon)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 ctint, frt, tctat, TERPF
      integer   i, inta, intmax, jd1, jd2, jda, jdb, jdt, k, m, n, nb,
     .          nn, npgsv, nplx, npr, nr, ntbct
 
c*** end of declarations inserted by spag
 
 
c
c           j.f.chandler - 1984 july - subroutine ctreed
c           interpolate ct-at from external integration
c
c     jdt   - julian day number of desired epoch
c             (returned as 0 if date not on data set)
c     frt   - fraction of day since midnight (ct) on 'jdt'
c     tctat - output value of ct-at in sec (iteration is pointless)
c     ctmon - if true, tctat includes monthly term
c
c     enter at ctrd1 to initialize
c              (must have jdt=0)
c
 
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      include 'tabval.inc'
 
c temporary storage
      common/YVECT/ Title(16),Prmtct(100)
      character*8 Title
      real*10 Prmtct
      real*10 mass9(9)
      equivalence (mass9,Prmtct)
 
      real*10 ctatvl(200),yvct(7,2),pct(4)
      logical*4 warned/.false./
      integer*4 jdct(3),nprct(2),nprv(3)
      character*8 montrm/'CT - AT '/
      logical*4 monthl, ctmon
 
 
      ctmon = .false.
 
c check requested time
      if(Nk1.lt.0 .and. Npage.ne.npgsv) warned = .false.
      npgsv = Npage
      if(Ictat.le.0 .or. jdt.lt.jd1 .or. jdt.ge.jd2) then
         if(.not. (warned)) then
            warned = .true.
            call PAGCHK(60,2,0)
            write(Iout,20) jdt,Ictat,jd1,jd2
   20       format('0******* WARNING:', i8, ' NOT ON CT-AT DATA SET',
     .             i3,'  (',i7,'-',i7,')')
            if(Mout.gt.0) write(Mout,20) jdt,Ictat,jd1,jd2
         endif
         jdt = 0
         return
      endif
 
 
c check position on tape
 
  100 if(jdt.lt.jda) then
 
c must back up
         n = 3 + (jdct(1)-jdt-1)/intmax
         do i = 1, n
            backspace Ictat
         end do
         m = 1
 
      else if(jdt.lt.jdb) then
 
c set up y-vectors
         nb = (jdt - jdct(1))/inta
         if(nb.ne.ntbct) then
            ntbct = nb
            nr    = nb - 4
            do i = 1, 10
               Tabvl(i) = ctatvl(nr + i)
            end do
            Ntab1 = 1
            Ntab2 = 2
            call YCOFF(yvct)
         endif
 
c interpolate value and return
         pct(1) = ((jdt-jdct(1)-ntbct*inta) + frt)/ctint
         pct(2) = pct(1)**2
         pct(3) = 1._10 - pct(1)
         pct(4) = pct(3)**2
         tctat  = TERPF(pct,yvct)
         ctmon  = monthl
         return
 
      else
 
c must skip forward
         n = (jdt - jdct(3))/intmax
         if(n.le.0) then
 
c preserve one record in core, then read one more
            jdct(1)  = jdct(2)
            nprct(1) = nprct(2)
            n  = nprct(1)
            nn = nprv(2)
            do i = 1, n
               ctatvl(i) = ctatvl(i + nn)
            end do
            m = 2
         else
            if(jdt.gt.jd2-intmax) n=n-1
            do i = 1, n
               read(Ictat)
            end do
            m = 1
         endif
 
      endif
 
 
c read one or two more records
  200 do k  = m, 2
         nn = nprv(k)
         read(Ictat) jdct(k),npr,(ctatvl(i+nn),i = 1,npr)
         nprv(k + 1) = nn + npr
         nprct(k)    = npr
      end do
      jdct(3) = jdct(2) + inta*npr
      jda     = jdct(1) + 4*inta
      jdb     = jdct(3) - 5*inta
      ntbct   = 9999
      if(jdt.eq.0) return
 
      goto 100

c
c-----------------------------------------------------------------------
c read headers
c
      entry CTRD1(jdt)
 
      warned  = .false.
      nprv(1) = 0
      if(Ictat.le.0) return
      read(Ictat) Title
      read(Ictat) jd1,jd2,inta,nplx,Prmtct
      monthl = (Title(1).eq.montrm)
 
c print constants
      call PAGCHK(60,7,0)
 
      write(Iout,300) Ictat,Title,jd1,jd2,inta,nplx,mass9
  300 format(
     . '0INFORMATION FROM HEADER RECORDS OF CT-AT INTEGRATION DATA SET'
     . , i3, '  TITLE='/1x, 16A8/' JD1=', i8, ' JD2=', i8, '  INT=',
     . i4, '  CT-AT FOR PLANET', i3/' MASSES=', (t10,1p,5D22.15))
 
      if(monthl) then
         write(Iout,500)
  500    format(' MONTHLY TERM INCLUDED')
      else
         write(Iout,400)
  400    format(' MONTHLY TERM NOT INCLUDED')
      endif

      ctint  = inta
      jd1    = jd1 + 4*inta
      jd2    = jd2 - 5*inta
      intmax = 100*inta
      m      = 1
 
      goto 200
 
      end
