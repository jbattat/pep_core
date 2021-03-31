      subroutine ATMDEN(rho,drhodr,atmexp)
 
      implicit none
c
c n.beebe/r.reasenberg  july 1970  subroutine atmden
c
 
c computation of atmospheric density of planet by formula or table
c lookup, with or without interpolation between taple points.

c arguments
      real*10 rho,drhodr,atmexp

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'lothrf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'sbthng.inc'
c
c     kp(82)= +-1 -- compute from formula
c     kp(82)= +-2 -- table lookup with four-point interpolation
c     kp(82)= +-3 -- table lookup without interpolation
c
      integer*2 kap
      real*10 emin/-60._10/,eterm,rsbkm

      rsbkm = Rsb*Aultsc*Ltvel
      kap   = Kp(82)
      if(kap.lt.0) kap = -kap
      if(kap.gt.1) then
c
c determine rho from table
         call SUICID('STOP IN ATMDEN, TABLE LOOKUP NOT CODED  ', 10)
         if(kap.eq.3) then
c use four-point interpolation
c
c no interpolation
         endif
      else
c
c compute rho from formula
         eterm  = (Rhzkm - rsbkm)/Sh
         rho    = 0.0_10
         drhodr = 0._10
         atmexp = 0._10
         if(eterm.ge.emin) then
            atmexp = EXP(eterm)
            rho    = Rhoz*atmexp
            drhodr = (-1._10*rho/Sh)*Aultsc*Ltvel
            if(eterm.gt.Drgeps) then
               call PAGCHK(60,1,0)
               write(Iout,10) eterm
   10          format(' IN ATMDEN CALCULATION, VALUE OF EXPONENT IS ',
     .                d22.15)
            endif
         endif
      endif
 
      return
      end
