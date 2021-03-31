      subroutine RTRD1(lice,nplntr)
 
      implicit none

c
c m.e. ash, r.w. king september 1972    subroutine rtrd1
c first five records of earth,moon,or planet rotation tape are read
c
c parameters
      integer*4 nplntr
      integer*2 lice
c lice =0 printout of data on first two records of planet tape
c lice =1 no such printout

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'rotdta.inc'
      include 'trpcom.inc'
      include 'yvectrd1.inc'
 
      integer   i,inter,intpr
      integer*2 klpr, nppr, kler/-1/, nper/-3/
c
c
c determine rotating body tape to be read
      if(nplntr.eq.10) then
         klpr = 0
         nppr = -10
      else if(nplntr.eq.3) then
c
c read first record of earth rotation data set
c read earth rotation constants from disk
         call XXRDBD(nper,Jer,kler,Jder,inter,Idirer,Ker,Erintx,
     .               Nkier,Kier,0,Ercom)
 
         call XXRD1(lice,nper,Jer,kler,Jder1,Jder2,Iparer,i_mxplprt+1,
     .              inter,Idirer,Ker,Erintx,Frer1,Nkier,Kier)
 
         Erint = Erintx
         if(Idirer.lt.0) Erint = -Erintx
 
c count partials and reset jd0 for next iteration
         call XXRDCK(Lparer,Kier,kler,1)
         goto 100
      else
         klpr = Klanr
         nppr = Nplnt(Klanr)
      endif
 
c read rotation constants from disk
      call XXRDBD(nppr,Jpr,klpr,Jdpr,intpr,Idirpr,Kpr,Printx,
     .            Nkipr,Kipr,0,Mrcom)
 
      call XXRD1(lice,nppr,Jpr,klpr,Jdpr1,Jdpr2,Iparpr,i_mxplprt+1,
     .           intpr,Idirpr,Kpr,Printx,Frpr1,Nkipr,Kipr)
 
      Print = Printx
      if(Idirpr.lt.0) Print = -Printx
      if(Kpr(100).eq.1) then
         do i = 1, 12
            Mrcom(i) = Cn1x(i + 24)
         end do
      endif
c
c count partials and reset jd0 for next iteration
      call XXRDCK(Lparpr,Kipr,klpr,1)
c
c read first three records of earth, moon, or planet rotation
c data set
  100 call RTRED1(0,nplntr)
      return
      end
