      subroutine YHERMT(tab, yh, ntype, mtype, hc2, nb)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   mb, mtype, nb, ntype, ntype1, ntype2
 
c*** end of declarations inserted by spag
 
 
c
c m.ash,smith,friedman - subroutine yhermt-hermite coeff. -sept.1969
c evaluate y-vectors for hermite interpolation
c
      real*10 tab(6, 8), yh(6)
      real*10 hc2(6)
      real*10 c1, c2, c3, hh
c
c          ntype=1 position interp.    ntype=2 velocity interp.
c          mtype=0 coordinate itself   mtype=1 partial of corrdinate
c
c   hemite 2 point, 2nd derivitives, fifth deg. polynomial interpolation
c               same method used for all four types
c
      mb     = nb + 1
      ntype1 = ntype + 1
      ntype2 = ntype + 2
c
c first three interpolation coefficients
      yh(1) = tab(nb, ntype)
      yh(2) = tab(nb, ntype1)
      yh(3) = tab(nb, ntype2)*0.5_10
c
c temporary quantities
      c3 = tab(mb, ntype2) - tab(nb, ntype2)
      c2 = tab(mb, ntype1) - tab(nb, ntype1) - hc2(nb)*tab(nb, ntype2)
      c1 = tab(mb, ntype) - tab(nb, ntype) - hc2(nb)
     .     *(tab(nb,ntype1) + hc2(nb)*tab(nb,ntype2)*0.5_10)
c
c last three interpolation coeffiecients
      hh    = 1.0_10/hc2(nb)**3
      yh(6) = (0.5_10*c3 - (3._10*c2-6._10*c1/hc2(nb))/hc2(nb))*hh
      yh(5) = (-c2 + hc2(nb)*0.5_10*c3)*0.5_10*hh - 2.5_10*hc2(nb)*yh(6)
      yh(4) = c1*hh - hc2(nb)*(yh(5) + hc2(nb)*yh(6))
 
      return
      end
