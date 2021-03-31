      subroutine PLRD1(lice)
 
      implicit none

c
c m.e.ash   march 1967    subroutine plrd1
c reference to s-body tape added 1977 jul - j.f.chandler
c first five records of planet tape are read
c also set up planet shape quantities
c
c parameters 
      integer*2 lice
c lice =0 printout of data on first two records of planet tape
c lice =1 no such printout
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'lcntrl.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'param.inc'
      include 'pemctl.inc'
      include 'plnhar.inc'
      include 'rotcom.inc'
      include 'scoef4.inc'
      include 'shpcom.inc'
      include 'shphar.inc'
      real*4    tlatin, tlonin
      equivalence (Scntrl(5),tlatin),(Scntrl(6),tlonin)
      include 'tapdta.inc'
      include 'trpcom.inc'
      include 'yvectrd1.inc'
 
c local variables
      real*10 frp1(2)
      integer*4 i,ibrnch,ict66,ictl,ishp,j,k,lpli,nt1,nz1

c
c read body constants from disk
c test if coordinates available from n-body tape and,
c if so, set up controls as if from individual tape
      call XXRDBD(Nplnt(Klan),Jplnt,Klan,Jdp,Intpx,Idirpl,Kp,Pintx,
     . Nkipl,Kipl,4,Pcom)
      if(Kp(88).ne.-8) then
c
c read first two records of planet peripheral data set
         if(Jtest.eq.0) call XXRD1(lice,Nplnt(Klan),Jplnt,Klan,
     .    Jdp1,Jdp2,Iparp,i_mxplprt+1,Intpx,Idirpl,Kp,Pintx,frp1,
     .    Nkipl,Kipl)
      endif
 
      Intp5 = 5*Intpx
      Pint  = Pintx
      if(Idirpl.lt.0) Pint = -Pintx
      T0svpl = Jdxx0
 
c setup planet shape and rotation quantities
      Pradls = Pcond(7,Klan)/Ltvel
      Venm   = Pcond(11,Klan)*Convd
      Psim   = Pcond(12,Klan)*Convd
      Omegm  = 0.0_10
      if(Pcond(13,Klan).ne.0.0_10) Omegm = Twopi/Pcond(13,Klan)
      Alphc(1) = (90.0_10 + Pcond(15,Klan))*Convd
      Alphc(2) = 0.0_10
      Salphc   = SIN(Alphc(1))
      Calphc   = COS(Alphc(1))
      Deltc(1) = (90.0_10 - Pcond(14,Klan))*Convd
      Deltc(2) = 0.0_10
      Sdeltc   = SIN(Deltc(1))
      Cdeltc   = COS(Deltc(1))
      Psic0    = Pcond(12,Klan)*Convd
c
c mars rotation setups
      ict66 = Ict(66)
      if(Nplnt(Klan).eq.4 .and. mod(ict66,2).ne.0) then
c
c        decide if psi0,i0 or alpha0,delta0 needs to be calculated.
c        this was done in bodred, but needs to be done again if
c        iterating and changing parameter values.  warning -- if i0
c        and psi0 are being adjusted, new r.a. and dec. will be
c        calculated.  this means that you cannot iterate if adjusting
c        i0,psi0 simultaneously with r.a. and dec.
c
         ictl = 0
         do i = 1, 30
            lpli = Lpl(i,Klan)
            if(lpli.eq.10 .and. lpli.eq.11) then
               ictl = 1
               goto 50
            endif
         end do
 
   50    call ROCHNG(Pcond(16,Klan),Pcond(17,Klan),Pcond(14,Klan),
     .               Pcond(15,Klan),Pcom(7),Pcom(8),ictl)
         call ROTSET(Trig,Pcond(15,Klan),Pcond(14,Klan),Pcom(7),
     .               Pcom(8),Pcond(16,Klan),Pcond(17,Klan),
     .               Pcond(12,Klan),Pcond(13,Klan),Pcond(18,Klan),
     .               Pcond(19,Klan))
         call QSET
      endif
c
c
c setup planet shape
      Npshp1 = 0
 
c this loop also zeroes out lz,lc,ls
      do i = 1, 500
         Zone(i) = 0._10
      end do
      do i = 1, 1000
         Lz(i) = 0
      end do

      Nshp = 0

      if(Nmphar.gt.0) then
         do i = 1, Nmphar
            if(Nplnt(Klan) + Nplhar(i).eq.0) then
               Npshp1 = i
               Nz     = Npzone(Npshp1)
               nz1    = Nz - 1
               Nt     = Nptess(Npshp1)
               nt1    = Nt*(Nt + 1)/2 - 1
               Nshp   = Nshape(Npshp1)
               Shpzer = Szero(Npshp1)
               do j = 1, 9
                  Scntrl(j) = Scontl(Npshp1,j)
               end do
               call FLATS(Pcond(8,Klan))
 
c branch for different shape models
               ibrnch = Nshp + 1
               if(ibrnch.eq.1) then
               else if(ibrnch.eq.2) then
c fourier series shape model (nshp=1)
c at present 122 coefs. used
                  nz1 = 122
                  nt1 = 0
               else if(ibrnch.eq.3) then
 
c altitude grid local shape model (nshp=2)
                  ishp  = 0
                  Grdf1 = 1._10/tlonin
                  Grdf1 = Grdf1/tlatin
                  Grdf2 = -2._10*Grdf1/Ltvel
                  do j = 1,1000/u_stdsz
 
c move two (or more) coefficients as one
                     Zone(j) = Pzhar(Npshp1,j)
                     do k = 1,u_stdsz
                        ishp     = ishp + 1
                        Lz(ishp) = Lpzhar(Npshp1,ishp)
                        if(ishp.ge.ngdpts) goto 100
                     end do
                  end do
                  goto 100
               else if(ibrnch.eq.4) then
c external shape model
                  goto 100
               else
                  call SUICID('NSHP OUT OF RANGE. STOP IN PLRD1', 8)
               endif
 
c spherical harmonic shape model (nshp=0)
               if(nz1.gt.0) then
                  do j = 1, nz1
                     Zone(j) = Pzhar(Npshp1,j)
                     Lz(j)   = Lpzhar(Npshp1,j)
                  end do
               endif
               if(nt1.gt.0) then
                  do j = 1, nt1
                     Ctess(j) = Pchar(Npshp1,j)
                     Stess(j) = Pshar(Npshp1,j)
                     Lc(j)    = Lpchar(Npshp1,j)
                     Ls(j)    = Lpshar(Npshp1,j)
                  end do
               endif
               goto 100
            endif
         end do
      endif
c
c count partials and reset jd0 for next iteration
  100 call XXRDCK(Lparp,Kipl,Klan,1)
      if(Kp(88).ne.-8) then
 
c see if any i.c. partials needed
         if(Jdpl0(Klan).eq.0) goto 200
c
c initial condition partials to be gotten from elliptic
c orbit partials
         if(Kipl(1).ne.0) goto 200
      endif
      Kipl(1) = -2
      i = Nplnt(Klan)
      call IMITL(Gauss,Mass(i),Pcond(1,Klan),1,3)
c
c read first three data records of planet data set
  200 if(Jtest.eq.0 .and. Kp(88).ne.-8) call PLRED1(0)
      return
      end
