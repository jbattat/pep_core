      subroutine SIDTIM(jd,ctut1,sidtm0,fract,gmstera)
 
      implicit none

c
c m.e.ash and r.a.goldstein   july 1964  subroutine sidtim

c parameters
      real*10 ctut1,sidtm0,fract,gmstera
      integer jd
c          jd    =input julian day number
c          ctut1 =input offset of CT from UT1 (s)
c          sidtm0=output mean sidereal time at 0 hr uninversal time UT1
c                =mean sidereal time at julian date jd-.5 (in radians)
c          fract =output ratio between mean sidereal time and universal
c                 time multiplied by twopi divided by 86,400
c          gmstera=difference between Greenwich mean sidereal time and
c                 Earth rotation angle at epoch (in radians)
c                 (only if 2000 or later IAU model selected)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'

c local variables
      real*10 t,tt,tu,era0
      integer*4 j
      logical*4 init/.false./

      if(Jct(21).lt.2) then
         t = jd - 2415020
         gmstera=0._10
      else
         t = jd - 2451545
      endif
      t = t - 0.5_10
c
c     logic added mar 1980 to select old or new iau definition of
c     sidereal time.  the new definition was recommended by iau
c     commissions 4, 19, and 31 in aug 1980 so that there will be
c     no change in the rate of ut1 when the iau(1976) precession
c     constant if used, or a change in the value of ut1 when the fk5
c     equinox is used.
c
c     selection logic now superseded by Jct(21)
c
      if(Jct(21).eq.0) then
c
c use IAU 1980 expression for sidtm
c changed by r.king dec 1980 to final iau 1982 values
c Corresponding expressions in seconds of time and julian centuries
c relative to J2000 are:
c excess fract= 8640184.81284857 +0.186209355 t -1.864121e-05 t^2
c sidtm0=
c  24110.5483974 +8640184.81285379 t +0.09310441 t^2 -6.197924e-06 t^3

         fract  = 7.2921158548774E-5_10 + t*(1.1750609E-19_10 -
     .    t*3.22E-28_10)
         sidtm0 = 1.739935359516_10 + t*(1.72027914345E-2_10 +
     .    t*(5.076246E-15_10-t*9.25E-24_10))

      else if(Jct(21).lt.0) then
         if(Ercond(28).ne.0._10 .and. Ercond(28).ne.5026.75_10) then
c if a non-standard value of the precession constant is used,
c print out a warning message and use the old definitin of sidtm
            if(.not. (init)) then
               init = .true.
               call SUICID(
     .'***WARNING:  NON-STANDARD VALUE OF PRECESSION CONSTANT USED, SIDT
     .M WILL CAUSE A CHANGE IN UT1***', -24)
            endif
         endif
c
c use old expression for sidtm
c corresponding expressions in seconds of time and julian centuries
c relative to J2000 are:
c excess fract= 8640184.72800493 + 0.185874259 t
c sidtm0= 24110.4708965 + 8640184.72779644 t + 0.09290001 t^2
         fract  = 7.2921158546827E-5_10 + t*1.1727115E-19_10
         sidtm0 = 1.73993589372_10 + t*(1.7202791266E-2_10 +
     .    t*5.06409E-15_10)
      else
c use IAU 2006 model
c This procedure adopts the new formula for true sidereal time (GST) as
c Earth rotation angle (ERA) minus Equation of origins (EO) instead of
c the old mean sidereal time (GMST) plus projected nutation in longitude.
c Instead of returning the rate of GMST, it now returns the rate of ERA.
c However, to preserve compatibility, the operation still calculates and
c returns GMST at 0h UT and also now the offset of the latter from ERA.
c This offset is to be saved by the caller, and the EO is added to it
c later to get the offset of GST from the approximation to GMST computed
c by updating GMST at 0h UT with the ERA rate.  Thus:
c ERA-EO = GMST(0) + UT*ERA' - (GMST(0)-ERA(0)) - EO
c Note: the formulation here comes from Capitaine, Wallace, & Chapront
c 2003, but it follows the expression in terms of angles, rather than
c time, the two expressions being inconsistent.
c In terms of time seconds and Julian centuries, this formulation is:
c excess fract= 8640184.79448082 + 0.185544227 tt -8.833E-8 tt^2
c -7.98827E-6 tt^3 -1.227E-8 tt^4  (not actually used)
c sidtm0 = 24110.5493771 + 8639877.31737601 tu + 307.4771023 tt
c +0.092772113 tt^2 -2.933E-8 tt^3 -1.99707E-6 tt^4 -2.453E-9 tt^5
         tu = t/36525._10
         tt = (t+ctut1/86400._10)/36525._10
         gmstera =(0.014506_10 + tt*(4612.156534_10 + tt*(1.3915817_10
     .    + tt*(-4.4e-7_10 + tt*(-2.9956e-5_10 + tt*(-3.68e-8_10))))))
     .    *Convds
         era0 = (0.2790572732640_10 + t*.00273781191135448_10)*Twopi
         sidtm0 = era0 + gmstera
c         fract = ((129602771.91721231_10 + tt*(2.7831634_10 + tt*
c     .    (-1.3249e-6_10 + tt*(-1.19824e-4_10 + tt*(-1.84e-7_10)))))
c     .    /(36525._10*1296000._10) + 1._10)*Convhs
         fract=1.00273781191135448_10*Convhs
      endif
 
      j = sidtm0/Twopi
      if(sidtm0.lt.0) j = j - 1
      t = j
      sidtm0 = sidtm0 - t*Twopi
      return
      end
