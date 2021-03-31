      subroutine PUNCH
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, is, j, k, klam, n, ncard, nd, num
 
c*** end of declarations inserted by spag
 
 
c
c ash / amuchastegui  - october 1969 - subroutine punch
c punch adjusted parameters at end of program run
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'empcnd.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      character*1 datel(8)
      equivalence (Date,datel)
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'plnhar.inc'
      include 'scoef4.inc'
      common/WRKCOM/ Bnum((u_mxtes*(u_mxtes+1))/2-1),Pun(100),Ipun(100)
      character*8 Bnum
      real*10 Pun
      integer*4 Ipun
      include 'zeroes.inc'
c
c quantities internal to this subroutine
      character*1 tem(9)/'9','9','(','9','9',')',3*' '/
      character*4 blank/'    '/
c
c punch title card
      ncard = 0
      write(Ipunch,100) Heding,datel(1),datel(2),datel(4),datel(5)
     . ,datel(7),datel(8),Jct(13)
  100 format(18A4,6A1/' &NMLST1  JCT(13)=',I2,',  JCT(27)=1,')
      ncard = ncard + 2
c
c punch adjusted parameters in &nmlst1
c
c search for positive lprm
      if(Lprm(1).gt.0) then
         num = 1
         do i = 2, 100
            if(Lprm(i).le.0) goto 150
            num = num + 1
         end do
c
c punch positive lprm
  150    write(Ipunch,200) (Lprm(i),i = 1,num)
  200    format(' LPRM=', (3x,15(i3,',')))
         ncard = ncard + (num - 1)/15 + 1
c
c punch adjusted parameters
         do i = 1, num
            j = Lprm(i)
            Pun(i)  = prmter(j)
            Ipun(i) = j
         end do
         write(Ipunch,250) (blank,Ipun(i),Pun(i),i = 1,num)
  250    format(2(a1,'PRMTER(',i3,')=',1pd26.19,','))
         ncard = ncard + (num - 1)/2 + 1
      endif
c
c punch adjusted dt table
      call DTPUN(ncard)
c
c set up harmonic labels
      nd = 1
      k  = 1
      call MVC(')  ',1,3, tem,6)
      do n = 2, 25
         is = 3 - nd
         call EBCDI(n,tem,2)
         do i = 1, n
            call EBCDI(i,tem(4),2)
            call MVC(tem,is,8,Bnum,k)
            k = k + 8
         end do
         if(n.eq.9) nd = 2
      end do
c start 'expediency mode' for higher than 25
c use CH(nnnn) instead of Ci(j),
c  where nnnn = i*(i-1)/2 - 1 + j
      tem(2) = 'H'   
      tem(8) = ')'
      is = k/8
      do n = 26,u_mxtes
         do i = 1, n
            is=is+1
            call EBCDI(is,tem(4),4)
            call MVC(tem,2,8,Bnum,k)
            k = k + 8
         end do
      end do
 
c punch &nmlst2 for embary
      call BODPUN(Nplnt(-3),Npcent(-3),Aplnt(-3),Iem,Jdem0,Lem,Econd,
     .            Nezone, Lezhar, Ezhar, Netess, Lechar, Echar, Leshar,
     .            Eshar, 1, 1, Izr2, 0, 0., ncard)
c
c punch &nmlst2 for moon
      call BODPUN(Nplnt(-2),Npcent(-2),Aplnt(-2),Imn,Jdmn0,Lmn,Mcond,
     .            Nmzone, Lmzhar, Mzhar, Nmtess, Lmchar, Mchar, Lmshar,
     .            Mshar, 1, 1, Izr2, 0, 0., ncard)
c
c punch &nmlst2 for earth rotation
      call BODPUN(Nplnt(-1),Npcent(-1),Aplnt(-1),Inut,Jder0,Ler,Ercond,
     .            Izr2, Izr2, 0._10, Izr2, Izr2, 0._10, Izr2,
     .            0._10, 1, 1, Izr2, 0, 0., ncard)
c
c punch &nmlst2 for moon rotation
      call BODPUN(Nplnt(0),Npcent(0),Aplnt(0),Ilib,Jdmr0,Lmr,Mrcond,
     .            Izr2, Izr2, 0._10, Izr2, Izr2, 0._10, Izr2,
     .            0._10, 1, 1, Izr2, 0, 0., ncard)
c
c punch &nmlst2 for planets, asteroids,satellites,
c probes and planet rotation
      if(Numpln.gt.0) then
         do j = 1, Numpln
            if(Nplnt(j).eq.0) goto 400
            do i = 1, 4
               if(Nplhar(i).eq.Nplnt(j)) then
                  klam = i
c
c with harmonics
                  call BODPUN(Nplnt(j),Npcent(j),Aplnt(j),Iplnt(j),
     .                        Jdpl0(j),Lpl(1,j),Pcond(1,j),
     .                        Npzone(klam),Lpzhar,Pzhar,Nptess(klam),
     .                        Lpchar, Pchar, Lpshar, Pshar, 4, klam,
     .                        Icnd(j),Nshape(klam),Scontl,ncard)
                  goto 300
               endif
            end do
c
c no harmonics
            call BODPUN(Nplnt(j),Npcent(j),Aplnt(j),Iplnt(j),
     .            Jdpl0(j),Lpl(1,j),Pcond(1,j),Izr2,Izr2,0._10,Izr2,
     .            Izr2, 0._10, Izr2, 0._10, 1, 1, Icnd(j),0,0.,ncard)
  300    end do
      endif
c
c
c punch adjusted pulsar parameters
  400 call PSRPUN(ncard)
c
c punch adjusted site coordinates
      call SITPUN(ncard)
c
c punch adjusted spot coordinates
      call SPTPUN(ncard)
c
c punch adjusted radar observation biases
      call RBSPUN(ncard)
c
c punch adjusted equinox-equator-latitude coorections
      call EQNPUN(ncard)
c
c punch adjusted planetary phase corrections
      call PHSPUN(ncard)
c
c punch adjusted star catalog error coefficients
      call SKYPUN(ncard)
c
c completion message
      write(Ipunch,500)
  500 format('**END'/'/*')
      ncard = ncard + 1
      write(Iout,600) ncard
  600 format(/i5,' INPUT DATA CARDS PUNCHED')
 
      return
      end
