      subroutine ASTCRD(jd, fract)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 fk, fl, fract
      integer   i, ivl, ivl11, ivl41, ivl81, j, jd, k, k1, k2, k3, k4,
     .          l, l1, l2, l3, l4, l5, l6, m1
      integer   m2, m3, mm, mtab, n, n3, nasu, ni, nj, nk, nn, ntab,
     .          ntqn
 
c*** end of declarations inserted by spag
 
 
c  j.f.chandler - 1976 july
c  override elliptic orbit calculations for bodies in fmmips
c           y-vectors are calculated as they are needed
c           jd   = input julian day number
c           fract= input fraction of day
c            yast(.,i)= output position for perturbing body i
c           ler     =  0 positions only needed for perturbing bodies
c           ler     =  1 velocities needed in relativity calculation
c           ler     =  2 velocities and accelerations needed in rel.cal.
c  note - present configuration assumes ler = 0
c
c common
      include 'b2dtaint.inc'
      include 'b2ydta.inc'
      include 'inodta.inc'
      include 'prtpin.inc'
      include 'scdtaint.inc'
 
      character*76 gmess/' JD=*******   READ ERROR ON PERTURBING ASTEROI
     .D DATASET **, STOP IN ASTCRD  '/
      real*10 dd(2, 3)
      equivalence (dd(1,1), D1)
      integer*2 itab(3)/36, 6, 1/, ntq(12)
 
      if(Nast.le.0) return
 
c see if jd is on tape
      if(jd.gt.Jdb21) go to 300
 
c not on tape
  100 call MVC('   NOT FOUND',1,12, gmess,13)
 
c encode jd and jpert into message
  200 call EBCDIX(jd, gmess, 5, 7)
      call EBCDIX(Jpert, gmess, 57, 2)
      call SUICID(gmess, 19)
  300 if(jd.ge.Jdb22) go to 100
c read records from tape ( if necessary)
c get correct three records of perturbing sat. tape
      fk = ((jd-Jdb2(2)) + fract - Frb2(2))/dd(1, 1)
      if(fk.lt.0) then
c
c correct records are behind on the tape
         nn = fk + 1.E-4_10
         nn = 4 - nn/40
         do i = 1, nn
            backspace Jpert
         end do
      else
         nn = fk/40._10
         if(nn.eq.0) go to 500
         mm = nn - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on the tape
c storage must be shifted
            Jdb2(1) = Jdb2(nn + 1)
            Frb2(1) = Frb2(nn + 1)
            Ivl1(1) = Ivl1(nn + 1)
            Ivl4(1) = Ivl4(nn + 1)
            Ivl8(1) = Ivl8(nn + 1)
            l1  = 1
            l2  = 40
            l3  = 1
            l4  = 10
            l5  = 1
            l6  = 5
            m1  = nn*40
            m2  = nn*10
            m3  = nn*5
            n   = 2
            ivl = max0(Ivl1(1), Ivl4(1), Ivl8(1))
            do while( .true. )
               do i = 1, ivl
                  if(Na1.gt.0) then
                     do j = l1, l2
                        mm = j + m1
                        do k = 1, Na1
                           Bod1(i, j, k) = Bod1(i, mm, k)
                        end do
                     end do
                  endif
                  if(Na4.gt.0) then
                     do j = l3, l4
                        mm = j + m2
                        do k = 1, Na4
                           Bod4(i, j, k) = Bod4(i, mm, k)
                        end do
                     end do
                  endif
                  if(Na8.gt.0) then
                     do j = l5, l6
                        mm = j + m3
                        do k = 1, Na8
                           Bod8(i, j, k) = Bod8(i, mm, k)
                        end do
                     end do
                  endif
               end do
               if(nn.gt.1) go to 400
               nn = 2
               Jdb2(2) = Jdb2(3)
               Frb2(2) = Frb2(3)
               Ivl1(2) = Ivl1(3)
               Ivl4(2) = Ivl4(3)
               Ivl8(2) = Ivl8(3)
               l1  = 41
               l2  = 80
               l3  = 11
               l4  = 20
               l5  = 6
               l6  = 10
               n   = 3
               ivl = max0(Ivl1(2), Ivl4(2), Ivl8(2))
            end do
         else if(mm.ne.0) then
 
c correct records are ahead on tape
            do i = 1, mm
               read(Jpert, end=100, err=310)
               go to 320
  310          read(Jpert)
 
c second read of error record might not be needed if system changes
  320       end do
         endif
      endif
c
c read perturbing satellite tape
      entry RDATAP(jd)
 
c where to start in array
      n     = 1
  400 do l  = n, 3
         l2 = 40*l
         l1 = l2 - 39
         l4 = 10*l
         l3 = l4 - 9
         l6 = 5*l
         l5 = l6 - 4
         if(Na1.le.0) l2 = l1
         if(Na4.le.0) l4 = l3
         if(Na8.le.0) l6 = l5
         read(Jpert,end=100,err=200) Jdb2(l),Frb2(l),ivl11,ivl41,ivl81,
     .        (((Bod1(i,j,k),i=1,ivl11),j=l1,l2), k=1,max0(1,0+Na1)),
     .        (((Bod4(i,j,k),i=1,ivl41),j=l3,l4), k=1,max0(1,0+Na4)),
     .        (((Bod8(i,j,k),i=1,ivl81),j=l5,l6), k=1,max0(1,0+Na8))
         Ivl1(l) = ivl11
         Ivl4(l) = ivl41
         Ivl8(l) = ivl81
      end do
      do i = 1, Nast
         ntq(i) = -99
      end do
      if(jd.eq.0) return
  500 fk = ((jd-Jdb2(2)) + fract - Frb2(2))
      ni = 1
      do n = 1, Nast
         nasu = Nas(n)
         if(nasu.gt.0 .and. Kpa(n).ge.0) then
 
c set up indexing for 'equivalence'  of bod48 bod1
            nj = n - Na1
            if(nj.le.0) then
               n3 = (n - 1)*120
            else
               nk = nj - Na4
               if(nk.le.0) then
                  n3 = (nj - 1)*30 + 480
                  ni = 2
               else
                  n3 = (nk - 1)*15 + 600
                  ni = 3
               endif
            endif
            fl = fk/dd(1, ni)
 
c compute p vector and tabular indices
            ntab = fl
            P(1) = fl - ntab
            ntab = ntab + 1
            mtab = ntab + 1
            P(3) = 1.0_10 - P(1)
            P(2) = P(1)**2
            P(4) = P(3)**2
            ntqn = ntq(n)
            if(ntab.ne.ntqn) then
               k1 = ntab + itab(ni)
               k3 = 1
               if(mtab.ne.ntqn) then
                  if(ntab.eq.(ntqn+1)) then
                     k1 = k1 + 1
                     k3 = 2
                  else
                     k2 = 2
                     go to 510
                  endif
               endif
               k2 = 1
               k4 = 3 - k3
c
c shift y-vectors
c this would be  limvel if there were room for vel or acc.
               ivl = 3
               do i = 1, 5
                  do j = 1, ivl
                     Ybdast(i, j, k4, n) = Ybdast(i, j, k3, n)
                  end do
               end do
  510          ntq(n) = ntab
               k1     = k1 + n3
               call YPRTCD(Bod1(1,k1,1), Ybdast(1,1,k3,n), k2)
            endif
            call ASTERP(Ybdast(1,1,1,n), dd(1,ni), nasu)
         endif
      end do
      return
      end
