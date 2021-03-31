      subroutine MNPRD1(lice4)
 
      implicit none
c
c j.f.chandler    2015 april   subroutine mnprd1
c first five records of moon tape are read for n-body integration that
c doesn't integrate the moon -- alternative to reading n-body
c this routine therefore usurps the storage normally used by the
c n-body reading routine
c this routine is called if and only if ipert2 is positive
c
      integer*4 lice4
c lice4 =0 printout of data on first two records of moon tape
c lice4 =1 no such printout

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'aprtbf.inc'
      include 'bddtaint.inc'
      integer*2 Iparm
      equivalence (Iparm,Nplbd)
      integer*2 tempk(200),kimn(99),nkimn,km(u_nmprm)
      equivalence (Betabd,tempk,nkimn),(tempk(2),kimn),(tempk(101),km)
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'petina.inc'
      include 'prtpin.inc'
      include 'yvectrd1.inc'
      character*8 nm8(3)
      equivalence (Nammon,nm8)

c local
      integer*4 i,intmx,rec
      integer*2 np10/10/,klmn/-2/,lice
      real*10 frm1(2)

      Ntab1 = 1
      Mtab1 = 2
c
c some moon constants already read from disk
c restarted moon integrations not supported
      Jdxx9=0
c
c read first two records of moon peripheral data set
      lice=lice4
      call XXRD1(lice,np10,Ipert2,klmn,Jdbd1,Jdbd2,
     . Iparm,i_mxplprt+1,intmx,Ibdsgn,km,Dmoon(1),frm1,nkimn,kimn)

      if(Kkxx(70).ne.Jct(13)) call SUICID('MOON IPERT2 HAS WRONG REFEREN
     .CE FRAME, STOP IN MNPRD1   ',14)

      nm8(1)=Tpname
      nm8(2)=Title(1)
      nm8(3)=Title(10)
      Nbdy=10
      do i=1,6
         Betabd(i,Nbdy)=Mcond(i)
      end do
      rec=0
      do i=1,3
         Recbeg(i)=rec+1
         rec=rec+8
         Recend(i)=rec
      end do
      Dir = Ibdsgn
      Dmoon(1) = Dmoon(1)*Dir
      Dmoon(2) = Dmoon(1)**2
      do i = 1,9
         Kiss(i) = -1
      end do
      Kiss(10) = 1
      Kiss(11) = 0
      Kiss(12) = 0
      Ler    = 1
      Lps    = 0
      Limvel = 6
      Limnut = 0
      Limlib = 0
c
c read first three data records of moon peripheral data set
      call MNREDP(0)
      return
      end
