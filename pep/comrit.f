      subroutine COMRIT(marta)
 
      implicit none

c m.e.ash    oct 1969     subroutine comrit
c observed minus theory and partial derivative tape is written

c arguments
      integer*4 marta
c marta=0 comrit called in midst of observing series
c marta=1 comrit called at end of observing series

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'comdat.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      integer*2 ildt(6)
      equivalence (ildt(1), Ildt1)
      include 'nutprc.inc'
      real*10 nutprc(18)
      equivalence (nutprc, Nutpr)
      include 'obscrd.inc'
      real*10 ctrecf
      equivalence (ctrecf, Dstf)
      include 'prpgat.inc'
      include 'sitcrd.inc'
      include 'zeroes.inc'
c
c local
      integer   i, icntob, ict20, j, jdu, jdu1, lprnt
      integer*2 numprs
      character*1 blank/' '/
      integer*4 nnutsav
c
c
c saved quantities to go on output observation library tape
      integer*2 nmp2/1/
      real*10 dvx(2, 2)/4*0._10/
      data icntob/0/
c
c single precision quantities for nout output
      real*4    hr, xnout(3, 2)
c
c           jct(45)= 0   only td&dop regular print in compar
c           jct(45)< 0   also print out of az&el, ra&decl for receive
c                        times of td&dop   this done in angot
c           jct(45)> 0   also arecibo bcd site tape   this done in angot
      if(Jct(45).ne.0) call ANGOT(marta)
c
c indicate no longer first point of series
      if(Nk1.le.0) Nk1 = 1
c
c set constants for end of observation series
      if(marta.gt.0) then
         if(Ict(20).lt.0) icntob = 0
         Ncode  = 0
         Numobs = 1
         Num1   = 1
         Num2   = 1
         numprs = Numpar
         Numpar = 1
         Ncal   = 1
         if(Iabs2.le.0 .or. (Idumob.eq.1 .and. Ict(3).le.0)) then
c
c restore numpar if end of observation series
            Numpar = numprs
            return
         else
            Numsav = 1
            goto 100
         endif
c
c write observed minus theory and partial derivative tape
      else if(Iabs2.gt.0) then
         if(Idumob.ne.1 .or. Ict(3).gt.0) then
            do i = 1,8
               Save(i) = Dstf(i+1)
            end do
 
c note: nutprc not restored in obsred
            nnutsav=19
            if(Ncodf.ge.4 .and. Ncodf.le.9) nnutsav=9
            do i=1,nnutsav
               Save(i+8) = nutprc(i)
            end do
            if(nnutsav.ge.19) Save(27) = Pc(1)
            if(Numsav.lt.nnutsav+8) Numsav = nnutsav+8
            if((Ncodf.ge.4 .and. Ncodf.le.6) .or. Ncodf.gt.9) then
               Save(18) = Salph
               Save(19) = Calph
               Save(20) = Tdelt
               Save(21) = Estf(10)
               if(Numsav.lt.21) Numsav=21
            endif
            goto 100
         endif
      endif
 
      goto 200
 
100   continue
      write(Iabs2) Ncode, Ihr, Imin, Sec,
     .             (Result(j), Error(j), j = 1, 2), Atuts, Ututs,
     .             Clamp, Limb, Observ, Imonth, Iday, Iyear, Jds,
     .             Jd, Ctat, ctrecf, Numsav,
     .             (Save(i), i = 1, Numsav), Numobs, Numpar,
     .             ((Deriv(i,j),i=1,Numpar), j = Num1, Num2), ildt,
     .             nmp2, ((dvx(i,j),i=1,nmp2), j = 1, Numobs),
     .             Ncal, (Cal(i), Scal(i), Ical(i), i = 1, Ncal),
     .             (Sumcor(i), i = 1, 2), Lnshob, Lixshp,
     .             (Lshobs(i), i = 1, Lnshob), izero
 
200   continue
 
c
c write extra listing on nout
      if(marta.gt.0) then
         Numpar = numprs
      else
         if(Nout.gt.0) then
            jdu = Jds - 2440000
            if(Nobs(1).le.1 .and. Nobs(2).le.1) jdu1 = jdu
            if((jdu-jdu1) .gt. 1) jdu1 = jdu
            hr = (Ihr*3600 + Imin*60 + Sec)/3600. + (jdu - jdu1)*24
            do j = Num1, Num2
               xnout(1, j) = Result(j)
               xnout(2, j) = Deriv(1, j)
               xnout(3, j) = Deriv(2, j)
            end do
            write(Nout, 20) jdu, Ihr, Imin, Sec, hr,
     .                      ((xnout(i,j),i=1,3), j = Num1, Num2)
   20       format(i4, 2I3, f5.1, f9.5, 1pd10.2, d8.1, 2D10.2, d8.1,
     .             d10.2)
         endif
c
c printout partials if ict(20) so indicates
         if(Ict(1).gt.0) then
            ict20 = Ict(20)
            ict20 = iabs(ict20)
            if(ict20.gt.icntob) then
               icntob = icntob + 1
               lprnt  = (Numpar + 8)/6
               do j = Num1, Num2
                  if(Line + lprnt.gt.58) call OBSPAG
                  write(Iout, 30) ildt, Result(j), j,
     .                            (blank, Deriv(i,j), i = 1, Numpar)
   30             format(' ILDT=', 6I4, 1pd22.14, 5x, 'DERIV(', i1,
     .                   ')=', 3(a1,1pd21.14), a1/ (1x,6(d21.14,a1)))
                  Line = Line + lprnt
               end do
            endif
         endif
      endif
      return
      end
