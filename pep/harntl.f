      subroutine HARNTL(nplnt, zhar, char, shar, nzone, ntess)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash    aug 1971    subroutine harntl
c initialize gravitational potential or shape harmonics
c
c
c        only harmonics corresponding to nzone=20,ntess=10
c        are brought in here for initialization. although
c        model has been increased to a maximum
c        capacity of nzone=25,ntess=20
c        if a priori values for the higher order harmonics
c        are to be coded, the following specification &
c        equiv. statements must be enlarged to accomodate
c        them. (see bodred or bodntl to compare)
c        also, do loops 993 & 995 would have to be increased
c        s.brody 5/76   r.b. goldstein feb78
      real*10 zhar(19), char(54), shar(54)
      integer*2 nzone, ntess, nplnt
      real*10 zhr(19), chr(54), shr(54)
 
c 3 c11(11),c12(12),c13(13),c14(14),c15(15),
      real*10 j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14,
     .          j15, j16, j17, j18, j19, j20, c2(2), c3(3), c4(4), 
     .          c5(5), c6(6), c7(7), c8(8), c9(9), c10(10), s2(2), 
     .          s3(3), s4(4), s5(5), s6(6), s7(7), s8(8), s9(9), s10(10)
 
c 5,s11(11),s12(12),s13(13),s14(14),s15(15)
      equivalence(zhr, j2), (zhr(2), j3), (zhr(3), j4), (zhr(4), j5)
     .            , (zhr(5), j6), (zhr(6), j7), (zhr(7), j8),
     .            (zhr(8), j9), (zhr(9), j10), (zhr(10), j11),
     .            (zhr(11), j12), (zhr(12), j13), (zhr(13), j14),
     .            (zhr(14), j15), (zhr(15), j16), (zhr(16), j17),
     .            (zhr(17), j18), (zhr(18), j19), (zhr(19), j20)
      equivalence(chr, c2), (chr(3), c3), (chr(6), c4),
     .            (chr(10), c5), (chr(15), c6), (chr(21), c7),
     .            (chr(28), c8), (chr(36), c9), (chr(45), c10)
      equivalence(shr, s2), (shr(3), s3), (shr(6), s4),
     .            (shr(10), s5), (shr(15), s6), (shr(21), s7),
     .            (shr(28), s8), (shr(36), s9), (shr(45), s10)
c
c initialize all harmonics to zero
      do i = 1, 19
         zhr(i) = 0.0_10
      end do
      do i = 1, 54
         chr(i) = 0.0_10
         shr(i) = 0.0_10
      end do
c
c     zonal harmonics are unnormalized, but tesserals are normalized to
c     four pi (for gravitational harmonics).
c     shape harmonics are all normalized to one, except that shape
c     harmonics (nplnt negative) might be fourier coefficients depending
c     on shape subroutines that are used in compar link.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           initialize earth gravitational potential harmonics
      if( nplnt .ne. 3 ) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c initialize moon gravitational potential harmonics
         if( nplnt .ne. 10 ) return
c
c j.p.l. tech report 32-1306
c constants and related information
c for astrodynamic calculations, 1968 - pag. 24
         j2    = 2.054E-4_10
         c2(2) = 0.231E-04_10
 
c unnormalized tesseral for moon link (to be changed eventually)
         nzone = 2
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c no standard values yet for other body harmonics
         ntess = 2
      else
c
c smithsonian astrophysical observatory special report no.353
c 1973 smithsonian standard earth (iii), edited by e.m.
c gaposchkin, p.276 (zonal harmonics 1973 i)
         j2  = 1082.637E-6_10
         j3  = -2.541E-6_10
         j4  = -1.618E-6_10
         j5  = -0.228E-6_10
         j6  = 0.552E-6_10
         j7  = -0.352E-6_10
         j8  = -0.205E-6_10
         j9  = -0.154E-6_10
         j10 = -0.237E-6_10
         j11 = 0.312E-6_10
         j12 = -0.192E-6_10
         j13 = -0.339E-6_10
         j14 = 0.105E-6_10
         j15 = 0.105E-6_10
         j16 = 0.034E-6_10
         j17 = -0.220E-6_10
         j18 = -0.102E-6_10
         j19 = 0.099E-6_10
         j20 = -0.119E-6_10
c     j21 =  -0.083E-6_10
c     j22 =   0.092E-6_10
c     j23 =   0.145E-6_10
c     j35 =  -0.134E-6_10
c     j36 =   0.199E-6_10
         nzone = 11
c
c          earth's gravitational field's tesseral
c          harmonics derived from a combination of
c          satellite data with gravity anomalies-
c          smithsonian astrophysical observatory special report no.353
c          1973 smithsonian standard earth (iii), edited by e.m.
c          gaposchkin, pp.290-293
         c2(2)   = 2.3799E-06_10
         s2(2)   = -1.3656E-06_10
         c3(1)   = 1.9977E-06_10
         s3(1)   = 2.2337E-07_10
         c3(2)   = 7.7830E-07_10
         s3(2)   = -7.5519E-07_10
         c3(3)   = 4.9011E-07_10
         s3(3)   = 1.5283E-06_10
         c4(1)   = -5.1748E-07_10
         s4(1)   = -4.8140E-07_10
         c4(2)   = 3.4296E-07_10
         s4(2)   = 6.7174E-07_10
         c4(3)   = 1.0390E-06_10
         s4(3)   = -1.1923E-07_10
         c4(4)   = -1.0512E-07_10
         s4(4)   = 3.5661E-07_10
         c5(1)   = -5.3667E-08_10
         s5(1)   = -7.9973E-08_10
         c5(2)   = 5.9869E-07_10
         s5(2)   = -3.9910E-07_10
         c5(3)   = -5.8429E-07_10
         s5(3)   = -1.6338E-07_10
         c5(4)   = -1.1583E-07_10
         s5(4)   = -4.5393E-08_10
         c5(5)   = 1.3956E-07_10
         s5(5)   = -8.6841E-07_10
         c6(1)   = -7.2166E-08_10
         s6(1)   = 1.7756E-08_10
         c6(2)   = 2.4670E-08_10
         s6(2)   = -4.0654E-07_10
         c6(3)   = 4.4139E-09_10
         s6(3)   = 2.9055E-08_10
         c6(4)   = -1.0003E-07_10
         s6(4)   = -3.0297E-07_10
         c6(5)   = -1.3504E-07_10
         s6(5)   = -6.0964E-07_10
         c6(6)   = -2.9136E-08_10
         s6(6)   = -2.6327E-07_10
         c7(1)   = 2.3532E-07_10
         s7(1)   = 5.5634E-08_10
         c7(2)   = 2.0425E-07_10
         s7(2)   = 1.7321E-07_10
         c7(3)   = 2.1994E-07_10
         s7(3)   = -3.4644E-07_10
         c7(4)   = -2.8617E-07_10
         s7(4)   = -2.7738E-07_10
         c7(5)   = 3.4727E-08_10
         s7(5)   = 8.7014E-08_10
         c7(6)   = -2.7496E-07_10
         s7(6)   = 8.5865E-08_10
         c7(7)   = -2.4856E-08_10
         s7(7)   = -8.8968E-09_10
         c8(1)   = 1.0946E-08_10
         s8(1)   = 4.8429E-08_10
         c8(2)   = 1.1084E-07_10
         s8(2)   = 1.0359E-07_10
         c8(3)   = -8.8578E-08_10
         s8(3)   = -5.0715E-08_10
         c8(4)   = -2.2315E-07_10
         s8(4)   = 2.6511E-07_10
         c8(5)   = 1.5318E-07_10
         s8(5)   = 8.1158E-08_10
         c8(6)   = -9.7542E-08_10
         s8(6)   = 2.8082E-07_10
         c8(7)   = 2.0498E-07_10
         s8(7)   = 2.4592E-07_10
         c8(8)   = 1.6967E-07_10
         s8(8)   = 9.3261E-08_10
         c9(1)   = 1.8099E-07_10
         s9(1)   = 4.1091E-08_10
         c9(2)   = -2.2013E-08_10
         s9(2)   = 2.4215E-08_10
         c9(3)   = -9.9252E-08_10
         s9(3)   = -2.3085E-08_10
         c9(4)   = -4.0867E-08_10
         s9(4)   = -3.8525E-08_10
         c9(5)   = -5.8957E-08_10
         s9(5)   = 3.6834E-09_10
         c9(6)   = 4.8812E-08_10
         s9(6)   = 1.1115E-07_10
         c9(7)   = -1.9880E-07_10
         s9(7)   = -1.4978E-07_10
         c9(8)   = 2.3523E-07_10
         s9(8)   = 9.6355E-09_10
         c9(9)   = -3.4533E-08_10
         s9(9)   = 5.9502E-08_10
         c10(1)  = 8.9008E-08_10
         s10(1)  = -6.0157E-08_10
         c10(2)  = -3.7256E-08_10
         s10(2)  = -6.3676E-08_10
         c10(3)  = -1.3307E-07_10
         s10(3)  = -7.2728E-08_10
         c10(4)  = -2.1887E-08_10
         s10(4)  = -7.8408E-08_10
         c10(5)  = -6.1509E-09_10
         s10(5)  = -1.1904E-07_10
         c10(6)  = -9.4142E-08_10
         s10(6)  = -1.1728E-08_10
         c10(7)  = 1.8525E-07_10
         s10(7)  = 2.1656E-08_10
         c10(8)  = 1.0887E-09_10
         s10(8)  = 7.0781E-09_10
         c10(9)  = 7.8473E-08_10
         s10(9)  = 5.6381E-09_10
         c10(10) = 1.3321E-07_10
         s10(10) = 9.8839E-08_10
c     c11( 1) = -1.2194E-08_10
c     s11( 1) =  7.5463E-08_10
c     c11( 2) = -2.0255E-08_10
c     s11( 2) = -6.2998E-08_10
c     c11( 3) = -1.0988E-09_10
c     s11( 3) = -3.8098E-08_10
c     c11( 4) =  1.5676E-08_10
c     s11( 4) = -1.9551E-07_10
c     c11( 5) = -1.8591E-09_10
c     s11( 5) =  6.1113E-08_10
c     c11( 6) =  6.3601E-08_10
c     s11( 6) = -2.6457E-08_10
c     c11( 7) = -3.3761E-08_10
c     s11( 7) = -1.2825E-07_10
c     c11( 8) = -1.3634E-08_10
c     s11( 8) =  4.5229E-08_10
c     c11( 9) =  2.1256E-08_10
c     s11( 9) =  6.6721E-08_10
c     c11(10) =  5.2555E-08_10
c     s11(10) = -7.7401E-08_10
c     c11(11) =  8.6996E-08_10
c     s11(11) = -2.5691E-08_10
c     c12( 1) = -5.6935E-08_10
c     s12( 1) = -6.6159E-08_10
c     c12( 2) = -9.7424E-08_10
c     s12( 2) =  4.6341E-08_10
c     c12( 3) =  1.1555E-07_10
c     s12( 3) = -4.8666E-08_10
c     c12( 4) = -5.0379E-08_10
c     s12( 4) =  5.3568E-08_10
c     c12( 5) =  8.1834E-08_10
c     s12( 5) =  2.7932E-08_10
c     c12( 6) = -2.1177E-08_10
c     s12( 6) =  3.5034E-08_10
c     c12( 7) =  2.9751E-08_10
c     s12( 7) =  3.1783E-08_10
c     c12( 8) =  4.0190E-08_10
c     s12( 8) =  5.6877E-08_10
c     c12( 9) = -1.1503E-07_10
c     s12( 9) =  1.4508E-08_10
c     c12(10) = -4.5921E-08_10
c     s12(10) = -4.3264E-08_10
c     c12(11) = -7.8443E-09_10
c     s12(11) = -4.7858E-08_10
c     c12(12) = -2.7617E-08_10
c     s12(12) = -1.6808E-08_10
c     c13( 1) =  8.6136E-09_10
c     s13( 1) = -3.2401E-08_10
c     c13( 2) = -1.0679E-08_10
c     s13( 2) = -9.0670E-08_10
c     c13( 3) = -3.2361E-08_10
c     s13( 3) =  4.9286E-08_10
c     c13( 4) =  3.9852E-08_10
c     s13( 4) = -1.0608E-07_10
c     c13( 5) =  4.0047E-08_10
c     s13( 5) =  3.8114E-08_10
c     c13( 6) = -2.1906E-08_10
c     s13( 6) = -1.1321E-08_10
c     c13( 7) = -7.6933E-08_10
c     s13( 7) =  1.1140E-08_10
c     c13( 8) = -2.7448E-09_10
c     s13( 8) =  1.4309E-08_10
c     c13( 9) = -1.1588E-08_10
c     s13( 9) =  7.2989E-08_10
c     c13(10) =  4.1979E-09_10
c     s13(10) =  7.6769E-09_10
c     c13(11) = -5.4381E-08_10
c     s13(11) =  1.3450E-08_10
c     c13(12) = -4.6633E-08_10
c     s13(12) =  7.9963E-08_10
c     c13(13) = -6.8944E-08_10
c     s13(13) =  7.1891E-08_10
c     c14( 1) = -1.4359E-08_10
c     s14( 1) =  5.2390E-08_10
c     c14( 2) = -1.5908E-08_10
c     s14( 2) =  2.7374E-08_10
c     c14( 3) =  9.6915E-08_10
c     s14( 3) = -2.5631E-08_10
c     c14( 4) = -2.9864E-08_10
c     s14( 4) = -3.8189E-09_10
c     c14( 5) = -1.3828E-09_10
c     s14( 5) = -5.8680E-08_10
c     c14( 6) = -1.3872E-08_10
c     s14( 6) = -2.7976E-08_10
c     c14( 7) =  7.1056E-08_10
c     s14( 7) =  2.4043E-09_10
c     c14( 8) = -1.8779E-08_10
c     s14( 8) = -5.8750E-08_10
c     c14( 9) = -2.4322E-08_10
c     s14( 9) =  6.0461E-08_10
c     c14(10) =  2.8985E-08_10
c     s14(10) = -3.4224E-08_10
c     c14(11) =  8.2611E-08_10
c     s14(11) = -1.9627E-09_10
c     c14(12) =  1.1751E-09_10
c     s14(12) = -3.0967E-08_10
c     c14(13) =  3.0793E-08_10
c     s14(13) =  4.7620E-08_10
c     c14(14) = -6.5969E-08_10
c     s14(14) =  3.3030E-09_10
c     c15( 1) =  2.9358E-08_10
c     s15( 1) = -1.6691E-08_10
c     c15( 2) = -1.2291E-08_10
c     s15( 2) = -6.8963E-08_10
c     c15( 3) = -5.8921E-08_10
c     s15( 3) =  4.4772E-08_10
c     c15( 4) =  1.4876E-08_10
c     s15( 4) =  7.0359E-09_10
c     c15( 5) =  3.6806E-08_10
c     s15( 5) = -8.4051E-09_10
c     c15( 6) =  1.0081E-08_10
c     s15( 6) = -3.0473E-08_10
c     c15( 7) =  3.0439E-08_10
c     s15( 7) =  1.5775E-08_10
c     c15( 8) = -6.8884E-08_10
c     s15( 8) =  6.0808E-08_10
c     c15( 9) = -4.5169E-08_10
c     s15( 9) =  5.5556E-08_10
c     c15(10) =  6.2126E-08_10
c     s15(10) = -7.1799E-09_10
c     c15(11) = -4.4724E-08_10
c     s15(11) = -3.4391E-09_10
c     c15(12) = -4.2025E-08_10
c     s15(12) =  5.9072E-09_10
c     c15(13) = -4.1654E-08_10
c     s15(13) = -5.5892E-09_10
c     c15(14) =  9.5654E-09_10
c     s15(14) = -2.7145E-08_10
c     c15(15) = -5.6358E-08_10
c     s15(15) =  3.4895E-08_10
c     c16( 1) = -9.9588E-09_10
c     s16( 1) =  5.4160E-08_10
c     c16( 2) =  5.5086E-09_10
c     s16( 2) =  4.9455E-08_10
c     c16( 3) =  5.4189E-08_10
c     s16( 3) =  5.4887E-09_10
c     c16( 4) =  4.6176E-08_10
c     s16( 4) =  3.6270E-08_10
c     c16( 5) = -2.4432E-08_10
c     s16( 5) =  2.9671E-08_10
c     c16( 6) = -3.7203E-09_10
c     s16( 6) = -2.0786E-08_10
c     c16( 7) = -2.2794E-09_10
c     s16( 7) =  3.0609E-09_10
c     c16( 8) = -1.0459E-07_10
c     s16( 8) = -4.4731E-08_10
c     c16( 9) =  2.4845E-08_10
c     s16( 9) = -8.6262E-08_10
c     c16(10) = -3.9928E-08_10
c     s16(10) = -4.5058E-09_10
c     c16(11) = -2.0848E-08_10
c     s16(11) =  2.9738E-08_10
c     c16(12) =  1.5930E-08_10
c     s16(12) = -1.2703E-08_10
c     c16(13) =  2.5280E-08_10
c     s16(13) =  6.6240E-09_10
c     c16(14) = -1.4852E-08_10
c     s16(14) = -8.1713E-09_10
c     c16(15) = -7.7425E-08_10
c     s16(15) = -2.6491E-08_10
c     c16(16) = -1.8538E-08_10
c     s16(16) = -2.2310E-08_10
c     c17( 1) =  8.6593E-09_10
c     s17( 1) = -4.1093E-08_10
c     c17( 2) = -9.0769E-09_10
c     s17( 2) = -2.7205E-08_10
c     c17( 3) = -7.7864E-09_10
c     s17( 3) = -1.7913E-08_10
c     c17( 4) = -4.3231E-08_10
c     s17( 4) =  6.8203E-08_10
c     c17( 5) =  4.1513E-08_10
c     s17( 5) = -2.5453E-08_10
c     c17( 6) = -4.5453E-08_10
c     s17( 6) = -1.7273E-08_10
c     c17( 7) =  1.6938E-08_10
c     s17( 7) = -3.3752E-08_10
c     c17( 8) =  4.1231E-08_10
c     s17( 8) =  5.8792E-09_10
c     c17( 9) = -4.3119E-08_10
c     s17( 9) = -1.5974E-08_10
c     c17(10) = -1.0844E-08_10
c     s17(10) =  5.5628E-08_10
c     c17(11) = -4.4136E-08_10
c     s17(11) = -4.3123E-09_10
c     c17(12) =  3.1661E-08_10
c     s17(12) =  6.2982E-09_10
c     c17(13) =  2.5147E-08_10
c     s17(13) =  9.7728E-09_10
c     c17(14) = -5.5945E-09_10
c     s17(14) =  7.2604E-09_10
c     c17(15) =  4.9113E-08_10
c     s17(15) =  3.1958E-08_10
c     c17(16) = -2.3540E-08_10
c     s17(16) = -1.5882E-08_10
c     c17(17) = -9.0191E-08_10
c     s17(17) = -9.4775E-09_10
c     c18( 1) = -2.3557E-08_10
c     s18( 1) = -7.4536E-08_10
c     c18( 2) = -9.4249E-09_10
c     s18( 2) =  3.0353E-08_10
c     c18( 3) = -3.5003E-08_10
c     s18( 3) = -2.0464E-08_10
c     c18( 4) =  2.9433E-08_10
c     s18( 4) = -4.4672E-08_10
c     c18( 5) =  1.7511E-09_10
c     s18( 5) = -6.0367E-09_10
c     c18( 6) =  2.3931E-08_10
c     s18( 6) = -4.4966E-09_10
c     c18( 7) = -7.8040E-10_10
c     s18( 7) = -8.2010E-09_10
c     c18( 8) =  5.3819E-08_10
c     s18( 8) = -2.2106E-08_10
c     c18( 9) = -3.6120E-10_10
c     s18( 9) = -5.0562E-09_10
c     c18(10) =  4.2146E-08_10
c     s18(10) =  7.8924E-09_10
c     c18(11) =  2.4981E-08_10
c     s18(11) =  2.3183E-08_10
c     c18(12) = -6.2242E-09_10
c     s18(12) =  6.6025E-09_10
c     c18(13) = -2.6685E-08_10
c     s18(13) = -4.2500E-08_10
c     c18(14) =  9.1191E-09_10
c     s18(14) = -3.3129E-08_10
c     c18(15) = -4.1521E-08_10
c     s18(15) = -1.7610E-08_10
c     c18(16) =  2.4850E-08_10
c     s18(16) = -4.8182E-09_10
c     c18(17) =  3.5357E-08_10
c     s18(17) = -4.7166E-08_10
c     c18(18) = -3.4701E-10_10
c     s18(18) =  5.0554E-08_10
c     c19(12) =  3.6058E-08_10
c     s19(12) = -3.4421E-09_10
c     c19(13) =  9.6876E-09_10
c     s19(13) = -6.6095E-08_10
c     c19(14) =  7.6389E-09_10
c     s19(14) = -2.7649E-08_10
c     c20(13) =  2.7630E-08_10
c     s20(13) =  3.2389E-08_10
c     c20(14) =  3.3687E-08_10
c     s20(14) = -6.5741E-08_10
c     c21(13) = -1.9799E-08_10
c     s21(13) = -3.0711E-08_10
c     c21(14) =  1.6623E-08_10
c     s21(14) =  8.7215E-09_10
c     c22(13) = -7.9435E-09_10
c     s22(13) =  4.1452E-09_10
c     c22(14) =  2.8516E-09_10
c     s22(14) = -4.2148E-08_10
c     c23(13) = -1.3236E-08_10
c     s23(13) = -4.8892E-09_10
c     c23(14) = -2.1148E-08_10
c     s23(14) =  2.2010E-08_10
c     c24(14) =  3.4668E-09_10
c     s24(14) =  2.2983E-08_10
         ntess = 5
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c move harmonics that were set up
      do i = 1, 19
         zhar(i) = zhr(i)
      end do
      do i = 1, 54
         char(i) = chr(i)
         shar(i) = shr(i)
      end do
 
      return
      end
