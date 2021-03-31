      subroutine PRTCRD(jd,fract)
 
      implicit none
c
c m.e.ash   aug 1968   subroutine prtcrd (n-body format)
c determination of perturbing planet coordinates
c revised 5/26/66   to reduce storage requirements
c y-vectors are calculated as they are needed
c
c arguments
      integer*4 jd
      real*10 fract
c           jd   = input julian day number
c           fract= input fraction of day
c           xpert(.,i)= output position,velocity,acceleration for
c                       perturbing body i (vel,acc needed for rel)
c           rpert(i),rpert2(i),rpert3(i)= output distance, square, cube
c                       for perturbing body i
c           xpert3(i,j) =xpert(i,j)/rpert3(j)   (i=1,3) (j=1,10)
c           ler     =  0 positions only needed for perturbing bodies
c           ler     =  1 velocities needed in relativity calculation
c           ler     =  2 velocities and accelerations needed in rel.cal.
c           merc =position of mercury from perturbing planet tape
c           body =position of venus through pluto from tape
c           moon =position of moon from tape
c     mdstau =moon distance unit on perturbing planet tape in au

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bddtaint.inc'
      include 'bdydta.inc'
      integer*4 zbdydt/29916/   !r8=16188,r10=29916
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'orblun.inc'
      include 'param.inc'
      include 'prtcod.inc'
      include 'prtpin.inc'
      include 'scdtaint.inc'

c local
      real*10 fract1,frsim,tsim
      integer   i,ivl,j,jdsim,k,k1,k2,k3,k4,k5,l,l1,l2,
     .          l3,l4,l5,l6,limvl,m1,m2,m3,mm,mtab,mvl
      integer   n,ngo,nn,ntab,ntbody,ntmerc,ntmoon
      character*76 gmess/  ' JD=*******   READ ERROR ON PERTURBING PLANE
     1T DATA SET  **, STOP IN PRTCRD  '/
c
c see if jd is on perturbing planet tape
      if(jd.gt.Jdbd1. and. jd.lt.Jdbd2) goto 300
 
c not on tape
  100 call MVC('   NOT FOUND',1,12,gmess,13)
c
c error message for read redundency on perturbing planet tape
  200 call EBCDIX(jd,gmess,5,7)
      call EBCDIX(Ipert,gmess,57,2)
      call SUICID(gmess,19)
c
c get correct three records of perturbing planet tape into
c storage
  300 nn = Ibdsgn*(jd - Jdbd(2))
      if(nn.lt.0) then
c
c correct records are behind on the tape
         nn = 4 - (nn + 1)/20
         if(Ipert.le.0) then
 
            nn = 1
            goto 400
         else
            do i = 1,nn
               backspace Ipert
               end do
         endif
      else if(nn.eq.0) then
         goto 500
      else
c
c correct records are ahead on tape
         n = nn/20
         if(n.eq.0) goto 500
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on the tape
c storage must be shifted
            l1 = 1
            l2 = 10
            l3 = 1
            l4 = 5
            l5 = 1
            l6 = 40
            m1 = n*10
            m2 = n*5
            m3 = n*40
            nn = 2
            Jdbd(1) = Jdbd(n + 1)
            Frct(1) = Frct(n + 1)
            Ivel(1) = Ivel(n + 1)
            Mvel(1) = Mvel(n + 1)
            ivl     = Ivel(1)
            mvl     = Mvel(1)
            do while( .true. )
               do i = 1,ivl
                  do j = l1,l2
                     mm = j + m1
                     Merc(i,j) = Merc(i,mm)
                     end do
                  do j = l3,l4
                     mm = j + m2
                     do k = 1,8
                        Body(i,j,k) = Body(i,mm,k)
                        end do
                     end do
                  end do
               do j = l5,l6
                  mm = j + m3
                  do i = 1,mvl
                     Mon(i,j) = Mon(i,mm)
                     end do
                  end do
               do j = l5,l6
                  mm = j + m3
                  Psid(j) = Psid(mm)
                  Epsd(j) = Epsd(mm)
                  do i = 1,3
                     Librt(j,i) = Librt(mm,i)
                     end do
                  end do
               if(n.eq.2) goto 400
               Jdbd(2) = Jdbd(3)
               Frct(2) = Frct(3)
               Ivel(2) = Ivel(3)
               Mvel(2) = Mvel(3)
               ivl     = Ivel(2)
               mvl     = Mvel(2)
               n  = 2
               nn = 3
               l1 = 11
               l2 = 20
               l3 = 6
               l4 = 10
               l5 = 41
               l6 = 80
               end do
         else if(mm.ne.0) then
            if(Ipert.gt.0) then
               do i = 1,mm
                  read(Ipert,end=100,err=310)
                  goto 320
  310             read(Ipert)
 
c second read of error record might not be needed if system changes
  320             end do
            endif
         endif
         nn = 1
         goto 400
      endif
c
c read perturbing planet tape
      entry RDYCAL(jd)
      if(jd.gt.0 .or. Ipert.gt.0) then
         nn = 1
      else
 
c initial setup for pretending to read tape
         call ZFILL(Merc,zbdydt)
         if(Kiss(10).ge.0 .and. Ict(50).gt.0)
     .       call LUNSET(Gauss,Mass(3),2440001,0._10,0)
         return
      endif
  400 do l  = nn,3
         l2 = l*10
         l1 = l2 - 9
         l4 = l*5
         l3 = l4 - 4
         l6 = l*40
         l5 = l6 - 39
         if(Ipert.gt.0) then
 
c read a record from tape
            read(Ipert,end=100,err=200) Jdbd(l),Frct(l),ivl,
     .           mvl,((Merc(i,j),i=1,ivl),j = l1,l2),
     .           (((Body(i,j,k),i=1,ivl),j=l3,l4),k = 1,8),
     .           ((Mon(i,j),i=1,mvl),j = l5,l6),
     .           (Psid(i),Epsd(i),i = l5,l6),
     .           ((Librt(j,i),i=1,3),j = l5,l6)
            Ivel(l) = ivl
            Mvel(l) = mvl
            if(Kiss(10).ge.0 .and. Nmoon.gt.0) then
               do j = l5,l6
                  do i = 1,mvl
                     Mon(i,j) = Mon(i,j)*Mdstau
                     end do
                  end do
            endif
         else
 
c simulate reading tape
            Ivel(l) = 1
            Mvel(l) = 6
            Frct(l) = 0._10
            if(l.gt.1) Jdbd(l) = Jdbd(l - 1) + 20*Ibdsgn
            if(l.eq.1) Jdbd(l) = (jd - 1)/20*20 - 19
            if(l.eq.1 .and. Ibdsgn.lt.0) Jdbd(l) = Jdbd(l) + 60
            if(Ict(50).gt.0) then
               tsim = Jdbd(l)
               do j = l5,l6
 
c simulate moon coordinates
                  if(Kiss(10).ge.0) then
                     jdsim = tsim
                     frsim = tsim - jdsim
                     call LUNORB(jdsim,frsim,0)
                     do i = 1,6
                        Mon(i,j) = Ylun(i)
                        end do
                  endif
 
c simulate earth nutation
                  if(Jct(21).le.0) then
                     call PEPNUT((tsim-2451545.5_10)/36525._10,Psid(j),
     .                        Epsd(j))
                  else
                     call IAU2000A(tsim-0.5_10,Psid(j),Epsd(j))
                  endif
                  tsim = tsim + Dmoon(1)
                  end do
            endif
         endif
         end do
      ntmoon = -99
      ntbody = -99
      ntmerc = -99
      if(jd.eq.0) return
c
c
c some setup for moon, nutation, libration interpolation
      nn     = Ibdsgn*(jd - Jdbd(2))
  500 fract1 = fract*Dir
      if(Kiss(10).lt.0 .and. Kiss(11).le.0 .and. Kiss(12).le.0) goto 900
      ntab = 2*nn + 1
      P(1) = fract - 0.5_10
      if(P(1).lt.0._10) then
         P(1) = fract
      else
         ntab = ntab + Ibdsgn
      endif
      P(1)  = P(1)/Dmoon(1)
      Nptrp = 3
  600 do while( .true. )
c
c compute p vector and tabular indices
         if(P(1).lt.0._10) then
            P(1) = P(1) + 1._10
            ntab = ntab - 1
         endif
         mtab = ntab + 1
         P(3) = 1._10 - P(1)
         P(2) = P(1)**2
         P(4) = P(3)**2
         if(Nptrp.eq.1) then
c
c perform mercury interpolation
            if(ntab.ne.ntmerc) then
               if(mtab.eq.ntmerc) then
                  k1 = ntab + 6
                  k3 = 1
               else if(ntab.eq.ntmerc+1) then
                  k1 = mtab + 6
                  k3 = 2
               else
                  k1 = ntab + 6
                  k2 = 2
                  k3 = 1
                  goto 610
               endif
               k2    = 1
               k4    = 3 - k3
               limvl = Limvel
               if(Ler.lt.0 .and. Nptrp.eq.Lps) limvl = 6
               do i = 1,5
                  do j = 1,limvl
                     Ymerc(i,j,k4) = Ymerc(i,j,k3)
                     end do
                  end do
  610          call YPRTCD(Merc(1,k1),Ymerc(1,1,k3),k2)
               ntmerc = ntab
            endif
            call PRTERP(Ymerc,Dmerc)
            return
         else if(Nptrp.eq.2) then
            ngo = 3
            if(ntab.ne.ntbody) then
               if(mtab.eq.ntbody) then
                  k1 = ntab + 1
                  k3 = 1
               else if(ntab.eq.ntbody+1) then
                  k1 = mtab + 1
                  k3 = 2
               else
                  ngo = 1
                  k1  = ntab + 1
                  k2  = 2
                  k3  = 1
                  goto 620
               endif
               k2     = 1
               ngo    = 2
  620          ntbody = ntab
            endif
            do while( .true. )
               if(Kiss(Nptrp).ge.0) then
                  k5 = Nptrp - 1
                  if(ngo.eq.2) then
                     k4    = 3 - k3
                     limvl = Limvel
                     if(Ler.lt.0 .and. Nptrp.eq.Lps) limvl = 6
                     do i = 1,5
                        do j = 1,limvl
                           Ybody(i,j,k4,k5) = Ybody(i,j,k3,k5)
                           end do
                        end do
                  endif
                  if(ngo.lt.3) call
     .                 YPRTCD(Body(1,k1,k5),Ybody(1,1,k3,k5),k2)
                  call PRTERP(Ybody(1,1,1,k5),Dbody)
               endif
               if(Nptrp.lt.9) then
c
c perform venus through pluto interpolations
                  Nptrp = Nptrp + 1
               else
c
c set up mercury interpolation
                  if(Kiss(1).lt.0) return
                  ntab  = nn/2
                  n     = nn - 2*ntab
                  P(1)  = n
                  P(1)  = (P(1) + fract1)*0.5_10
                  ntab  = ntab + 1
                  Nptrp = 1
                  goto 700
               endif
               end do
         else
c
c rest of setup for moon, nutation, libration interpolation
            Nptrp = 10
            if(ntab.ne.ntmoon) then
               if(mtab.eq.ntmoon) then
                  k1 = ntab + 36
                  k3 = 1
               else if(ntab.eq.ntmoon+1) then
                  k1 = mtab + 36
                  k3 = 2
               else
                  k1 = ntab + 36
                  k2 = 2
                  k3 = 1
                  goto 640
               endif
               k2 = 1
               k4 = 3 - k3
c
c shift y-vectors for moon, nutation libration interpolation
               if(Kiss(10).ge.0) then
                  limvl = Limvel
                  if(Ler.lt.0 .and. Nptrp.eq.Lps) limvl = 6
                  do i = 1,5
                     do j = 1,limvl
                        Ymoon(i,j,k4) = Ymoon(i,j,k3)
                        end do
                     end do
               endif
               if(Kiss(11).gt.0) then
                  do i = 1,3
                     do j = 1,Limnut
                        Yeps(i,j,k4) = Yeps(i,j,k3)
                        Ypsi(i,j,k4) = Ypsi(i,j,k3)
                        end do
                     end do
               endif
               if(Kiss(12).gt.0) then
                  do i = 1,3
                     do j = 1,Limlib
                        do k = 1,3
                           Ylib(i,j,k4,k) = Ylib(i,j,k3,k)
                           end do
                        end do
                     end do
               endif
c
c calculate y-vectors for moon nutation libration interp.
  640          ntmoon = ntab
               if(Kiss(10).ge.0)
     .             call YPRTCD(Mon(1,k1),Ymoon(1,1,k3),k2)
               if(Kiss(11).gt.0) then
                  call YNUTLB(Psid(k1),Ypsi(1,1,k3),k2,Kiss(11))
                  call YNUTLB(Epsd(k1),Yeps(1,1,k3),k2,Kiss(11))
               endif
               if(Kiss(12).gt.0) then
                  do i = 1,3
                     call YNUTLB(Librt(k1,i),Ylib(1,1,k3,i),k2,
     .                           Kiss(12))
                  end do
               endif
            endif
         endif
         goto 800
  700 end do
c
c perform moon, nutation, libration interpolations
  800 if(Kiss(10).ge.0) call PRTERP(Ymoon,Dmoon)
      if(Kiss(11).gt.0) then
         call NUTERP(Ypsi,Dpsi,Dmoon,Kiss(11))
         call NUTERP(Yeps,Deps,Dmoon,Kiss(11))
      endif
      if(Kiss(12).gt.0) then
         do i = 1,3
            call NUTERP(Ylib(1,1,1,i),Librat(1,i),Dmoon,Kiss(12))
            end do
      endif
c
c set up venus through pluto interpolation
  900 ntab  = nn/4
      n     = nn - 4*ntab
      P(1)  = n
      P(1)  = (P(1) + fract1)*0.25_10
      ntab  = ntab + 1
      Nptrp = 2
      goto 600
 
      end
