      subroutine BODPUN(nplnt,ncentr,name,itape,jd0,l,cond,nzone,
     .                  lzhar,zhar,ntess,lchar,char,lshar,shar,
     .                  nsize,klam,icnd,nshape,scontl,ncard)

      implicit none


c*** start of declarations inserted by spag
      integer   itape,jd0,jf,klam,n12,ncard,nf,nshape,
     .          nsize,nt1,nz1

c*** end of declarations inserted by spag


c
c ash / amuchastegui - october 1969 - subroutine bodpun
c punch adjusted body parameters
c
c
c array dimensions
      include 'globdefs.inc'
c
c parameters
      character*8 name
      integer*2 nplnt,ncentr,l(u_nmbod),nzone,lzhar(nsize,20),ntess,
     .          lshar(nsize,100),lchar(nsize,100),icnd
      real*10 cond(u_nmbod),zhar(nsize,20),char(nsize,100),
     .          shar(nsize,100)
      real*4    scontl(nsize,9)
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'

c subscript stuff + temporaries
      common/WRKCOM/ Bnum((u_mxtes*(u_mxtes+1))/2-1),
     . Con(u_nmbod-6),Icon(u_nmbod-6)
      character*8 Bnum
      real*10 Con
      integer*4 Icon
c
c local
      character*4 blank/'    '/
      integer*4 n1(5)
      integer   i, k, m
      character*4 fortag(8)/'FA( ', 'FB( ', 'FC( ', 'FD( ', 'FAP(',
     .          'FCP(', 'FAPP', 'FBPP'/
c
c initialize flags
      do i = 1, 5
         n1(i) = 0
         end do
c
c see if initial conditions were adjusted
      do i = 1, 6
         if(l(i).gt.0) then
            n1(1) = 1
            goto 100
         endif
         end do
c
c see if parameters were adjusted
  100 do i = 7, u_nmbod
         if(l(i).le.0) goto 200
         n1(2) = n1(2) + 1
         end do
c
c see if  harmonics were adjusted
  200 if(nshape.le.0) then
         nz1 = nzone - 1
         nt1 = (ntess*(ntess+1))/2 - 1

c fourier shape
      else if(nshape.gt.1) then

c grid shape
         nz1 = 1000
         nt1 = 0
      else
         nz1 = 122
         nt1 = 0
      endif
      do i = 1, nz1
         if(lzhar(klam,i).gt.0) then
            n1(3) = 1
            goto 300
         endif
         end do
  300 if(nt1.gt.0) then
         do i = 1, nt1
            if(lchar(klam,i).gt.0) then
               n1(4) = 1
               goto 350
            endif
            end do
  350    do i = 1, nt1
            if(lshar(klam,i).gt.0) then
               n1(5) = 1
               goto 400
            endif
            end do
      endif
c
c see if there is to be punched output
  400 do i = 1, 5
         if(n1(i).gt.0) goto 410
         end do
      return

  410 write(Ipunch, 420) nplnt, ncentr, name, icnd, jd0, itape
  420 format('*OBJECT', i5/' NCENTR=', i3, ', NAME=''', a8, '''',
     .       ', ICND=', i2, ', JD0=', i8, ',ITAPE=', i2, ',')
      ncard = ncard + 2
c
c punch adjusted initial conditions
      if(n1(1).gt.0) then
         write(Ipunch, 430) (l(k), k = 1, 6)
  430    format('  L=', 6(i2,','))
         ncard = ncard + 1
         if(Jct(13).gt.0) then
            write(Ipunch, 435)
  435       format('  JTYPE=6')
            ncard = ncard + 1
         endif
         write(Ipunch, 440) (k, cond(k), k = 1, 6)
  440    format(2('  COND(',i1,')=',1pd26.19,','))
         ncard = ncard + 3
      endif
c
c punch adjusted parameters
      if(n1(2).gt.0) then
         n12 = n1(2)
         do m = 1, n12
            k = l(m+6)
            Con(m)  = cond(k+6)
            Icon(m) = k
            end do
         write(Ipunch, 450) (l(k+6), k = 1, n12)
  450    format('  L(7)=', (3x,20(i2,',')))
         ncard = ncard + (n12 - 1)/20 + 1
         write(Ipunch, 460) (blank, Icon(k), Con(k), k = 1, n12)
  460    format(2(a2,'CON(',i2,')=',1pd26.19,','))
         ncard = ncard + (n12 - 1)/2 + 1
      endif
c
      if(n1(3)+n1(4)+n1(5) .le. 0) return
      if(nshape.gt.0) then
c
c extra shape models
         write(Ipunch, 470) nshape
  470    format('  NSHP=', i2, ',')
         ncard = ncard + 1
      endif

      if(nshape.eq.2) then
c
c grid model
         call GRDPUN(lzhar,zhar,nsize,klam,nshape,scontl,ncard)
      else if(nshape.eq.1) then

c fourier model
         jf = 0
         do nf = 1, 6
            k = 0
            do m = 1, 20
               if(lzhar(klam,jf+m).gt.0) then
                  k = k + 1
                  Icon(k) = m
               endif
               end do
            if(k.gt.0) then
               write(Ipunch,480) (blank,fortag(nf),Icon(m),
     .              lzhar(klam,Icon(m)+jf), m = 1,k)
  480          format(5(a2,'L',a4,i2,')=',i2,','))
               ncard = ncard + (k - 1)/5 + 1
               write(Ipunch,490) (blank,fortag(nf),Icon(m),
     .              zhar(klam,Icon(m)+jf), m = 1,k)
  490          format(2(a2,a4,i2,')=',1pd26.19,','))
               ncard = ncard + (k - 1)/2 + 1
            endif
            jf = jf + 20
            end do
         do nf = 7, 8
            if(lzhar(klam,nf+114).gt.0) then
               write(Ipunch,500) fortag(nf),lzhar(klam,nf+114),
     .                           fortag(nf),zhar(klam,nf+114)
  500          format('  L',a4,'=',i2,', ',a4,'=',1pd26.19,',')
               ncard = ncard + 1
            endif
            end do
         write(Ipunch,510) (m,scontl(klam,m), m = 1,2)
  510    format(2('  TLAT(',i1,')=',f8.3,','))
         ncard = ncard + 1
      else
c
c punch adjusted zonal harmonic coefficients
         write(Ipunch,520) nzone,ntess
  520    format(' NZONE=', i3, ',', 3x, 'NTESS=', i3, ',')
         ncard = ncard + 1
         if(n1(3).gt.0) then
            write(Ipunch,530) (blank,m+1,lzhar(klam,m),m=1,nz1)
  530       format(6(a2,'LJ(',i2,')=',i2,','))
            ncard = ncard + (nz1 - 1)/6 + 1
            write(Ipunch,540) (blank,m+1,zhar(klam,m),m=1,nz1)
  540       format(2(1x,a2,'J(',i2,')=',1pd26.19,','))

            ncard = ncard + (nz1 - 1)/2 + 1
         endif
c
c punch adjusted tess. cosine harmonic coefficients
         if(n1(4).gt.0) then
            write(Ipunch,550) (blank,Bnum(m),lchar(klam,m),m=1,nt1)
  550       format(5(a1,'LC',a7,'=',i2,','))
            ncard = ncard + (nt1 - 1)/5 + 1
            write(Ipunch,560) (blank,Bnum(m),char(klam,m),m=1,nt1)
  560       format(2(a1,'C',a7,'=',1pd26.19,','))
            ncard = ncard + (nt1 - 1)/2 + 1
         endif
c
c punch adjusted tess.   sine harmonic coefficients
         if(n1(5).gt.0) then
            write(Ipunch,570) (blank,Bnum(m),lshar(klam,m),m=1,nt1)
  570       format(5(a1,'LS',a7,'=',i2,','))
            ncard = ncard + (nt1 - 1)/5 + 1
            write(Ipunch,580) (blank,Bnum(m),shar(klam,m),m=1,nt1)
  580       format(2(a1,'S',a7,'=',1pd26.19,','))
            ncard = ncard + (nt1 - 1)/2 + 1
         endif
      endif

      return
      end
