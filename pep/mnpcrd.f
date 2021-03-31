      subroutine MNPCRD(jd,fract)
 
      implicit none
c
c j.f.chandler   2015 april   subroutine mnpcrd
c determination of moon coordinates for n-body integration that doesn't
c include the moon as an integrated body
c based on prtcrd
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
      include 'aprtbf.inc'
      include 'bddtaint.inc'
      integer*2 Iparm
      equivalence (Iparm,Nplbd)
      include 'bdydta.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'prtcod.inc'
      include 'prtpin.inc'
      include 'scdtaint.inc'

c local
      integer   i,ivl,j,jdsim,k,k1,k2,k3,k4,l,l1,l2,
     . limvl,m1,mm,mtab,mvl,n,nn,ntab,ntmoon
      character*56 gmess
      character*10 errcod
c
c see if jd is on perturbing planet tape
      if(jd.gt.Jdbd1. and. jd.lt.Jdbd2) goto 300
 
c not on tape
  100 errcod=' NOT FOUND'
      goto 220
c
c error message for moon tape
  200 errcod='READ ERROR'
  220 write(gmess,230) jd,errcod,Ipert2
  230 format('JD=',i7,a10,' ON MOON DATA SET',i3,', STOP IN MNPCRD')
      call SUICID(gmess,14)
c
c get correct three records of perturbing planet tape into
c storage
  300 nn = Ibdsgn*(jd - Jdbd(2))
      if(nn.lt.0) then
c
c correct records are behind on the tape
         nn = 4 - (nn + 1)/4
         do i = 1,nn
            backspace Ipert2
         end do
      else if(nn.eq.0) then
         goto 500
      else
c
c correct records are ahead on tape
         n = nn/4
         if(n.eq.0) goto 500
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on the tape
c storage must be shifted
            l2 = (3-n)*8
            m1 = n*8
            nn = 4-n
            do i=1,3-n
               Jdbd(i) = Jdbd(i+n)
               Frct(i) = Frct(i+n)
               Mvel(i) = Mvel(i+n)
            end do
            mvl = Mvel(1)
            do j = 1,l2
               if(j.eq.9) mvl=Mvel(2)
               mm = j + m1
               do i = 1,mvl
                  Mon(i,j) = Mon(i,mm)
               end do
               Psid(j) = Psid(mm)
               Epsd(j) = Epsd(mm)
               do i = 1,3
                  Librt(j,i) = Librt(mm,i)
               end do
            end do
            goto 400
         else if(mm.ne.0) then
            do i = 1,mm
               read(Ipert2,end=100,err=310)
               goto 320
  310          read(Ipert2)
c second read of error record might not be needed if system changes
  320       end do
         endif
         nn = 1
         goto 400
      endif
c
c read moon tape into nbody data storage
      entry MNREDP(jd)
      nn = 1
  400 do l  = nn,3
         l2 = l*8
         l1 = l2 - 7
 
c read a record from tape
         read(Ipert2,end=100,err=200) Jdbd(l),Frct(l),mvl,
     .    (((Mon(i,MIN(j,2)+k-1),i=1,mvl),j=1,Iparm),k=l1,l2),
     .           (Psid(i),Epsd(i),i=l1,l2),
     .           ((Librt(j,i),i=1,3),j=l1,l2)
         Mvel(l) = mvl
      end do
      ntmoon = -99
      if(jd.eq.0) return
c
c
c some setup for moon, nutation, libration interpolation
      nn     = Ibdsgn*(jd - Jdbd(2))
  500 ntab = 2*nn + 1
      P(1) = fract - 0.5_10
      if(P(1).lt.0._10) then
         P(1) = fract
      else
         ntab = ntab + Ibdsgn
      endif
      P(1)  = P(1)/Dmoon(1)
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
c
c rest of setup for moon interpolation
      Nptrp = 10
      if(ntab.ne.ntmoon) then
         if(mtab.eq.ntmoon) then
            k1 = ntab + 4
            k3 = 1
         else if(ntab.eq.ntmoon+1) then
            k1 = mtab + 4
            k3 = 2
         else
            k1 = ntab + 4
            k2 = 2
            k3 = 1
            goto 640
         endif
         k2 = 1
         k4 = 3 - k3
c
c shift y-vectors for moon
         limvl = Limvel
         do i = 1,5
            do j = 1,limvl
               Ymoon(i,j,k4) = Ymoon(i,j,k3)
            end do
         end do
c
c calculate y-vectors for moon interp.
  640    ntmoon = ntab
         call YPRTCD(Mon(1,k1),Ymoon(1,1,k3),k2)
      endif
c
c perform moon interpolation
      call PRTERP(Ymoon,Dmoon)
      return
      end
