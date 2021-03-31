      subroutine MOON
 
      implicit none
c
c m.e.ash   april 1967    subroutine moon
c main program for moon numerical integration
c modified august 1977 by r.king and r.cappallo to perform
c simultaneous integration of orbit and rotation
c
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'aacoff.inc'
      include 'bdctrl.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'metuna.inc'
      include 'plndta.inc'
c
c quantities internal to this routine
      real*4    epsmr(6)
      real*10 rd8
      integer*4 i, iperts, ngo
      integer*2 id2
      character*20 typmes/'  MOON  INTEGRATION '/
c
c printout moon title page
      call PAGSET(typmes, 5)
      call OPRMSG(typmes, 5)
      call NEWPG
      write(Iout, 100)
      write(Iout, 200)
  100 format('-'/ '-'/ '-'/ '-'//35x,
     . '**            **    *********     *********    **           **'
     . /35x,
     . '***          ***   ***********   ***********   ***          **'
     . /35x,
     . '****        ****   **       **   **       **   ****         **'
     . /35x,
     . '** **      ** **   **       **   **       **   ** **        **'
     . /35x,
     . '**  **    **  **   **       **   **       **   **  **       **'
     . /35x,
     . '**   **  **   **   **       **   **       **   **   **      **'
     . /35x,
     . '**    ****    **   **       **   **       **   **    **     **'
     . )
  200 format(35x,
     . '**     **     **   **       **   **       **   **     **    **'
     . /35x,
     . '**            **   **       **   **       **   **      **   **'
     . /35x,
     . '**            **   **       **   **       **   **       **  **'
     . /35x,
     . '**            **   **       **   **       **   **        ** **'
     . /35x,
     . '**            **   **       **   **       **   **         ****'
     . /35x,
     . '**            **   ***********   ***********   **          ***'
     . /35x,
     . '**            **    *********     *********    **           **'
     . )
      call PLINK
c
c determine control and data constants
      read(Iplcon) Econ1
      read(Iplcon) Mcon1,Epsm,Km,Jdm1,Jdm2,Intmn,i,i,id2,id2,rd8,
     . Kkm,id2,Tconm,Numkim,(Kim(i),i=1,Numkim)
      read(Iplcon) Ercon1
      read(Iplcon) Mrcon1,epsmr,Kmr,Jdmr1,Jdmr2,Intmr,i,i,id2,id2,rd8,
     . Kkmr,id2,Tconmr,Numkir,(Kir(i),i=1,Numkir)
      rewind Iplcon
c
c set planet numbers
      Mplnt  = 10
      Mplntr = -10
      Mcentr = 3
 
      iperts = Ipert
      if(Nbody.gt.0 .and. Ict(9).gt.0) Ipert = Ibody

c reference epoch for coordinate systems
      Kepoch= 1
      if(Jct(13).gt.0) Kepoch=2
 
c new code calls prtrd1 in morset
      call MORSET(0)
      ngo = 8
c
c perform numerical integration
      if(Km(88).lt.2) then
         call NUMINT(ngo)
      else if(Km(88).eq.2) then
         call ADAM(ngo)
      else
         call RROAD(ngo)
      endif
 
      rewind Ipert
      Itrwnd(Ipert) = 0
      Ipert = iperts
      call EMPRWD
      if(Imn.gt.0) rewind Imn
      if(Ilib.gt.0) rewind Ilib
      return
      end
