      real*10 function A1WWV(jd, fract)
c
c Ash/Amuchastegui - June 69  double precision function A1WWV
c  7/17/98 updated for 1999 Jan 1 quantum jump, plus a guess at the next
c  1/10/01 re-revised guess of the next leap second
c  1/15/03 removed guess of the next leap second
c  7/18/05 updated for 2006 Jan 1 leap second (at last)
c  7/04/08 updated for 2009 Jan 1 leap second
c  1/05/12 updated for 2012 Jul 1 leap second
c  1/05/15 updated for 2015 Jul 1 leap second
c  7/06/16 updated for 2017 Jan 1 leap second
c
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ilp, int, int1, j, jd, jd0, m1, mjd, nchang
 
c*** end of declarations inserted by spag
 
 
c
c     Calculate A.1-WWV time from 1956.
c     Linear interpolation from final bulletins at ten day intervals
c     from January 1956 to January 1966.
c     Linear formula used thereafter.
c$$$$ routine to be revised as history progresses
      real*10 t, fract, a1, a2, result, fract0
      data fract0/0.79166666666666666666_10/, nchang/302/
c
c table of A.1-WWV at 10 day intervals
      integer*2 a(375)
c
c frequency offset of WWV from A.1 in units of 10**-10
c 1956 Jan 07 to 1957 Aug 19 frequency offset varied between -36 to -98
c 1957 Aug 29 to 1958 Dec 22 frequency offset varied between -83 to -127
c 1959 Jan 01 to 1959 Dec 31 frequency offset varied between -92 to -112
c 1960 to 1961  frequency offset was kept near -150
c 1962 to 1963  frequency offset was kept near -130 (change in slope
c               1962 Jan 1 ignored between tabular points)
c 1964 to 1965  frequency offset was kept near -150 (change in slope
c               1964 Jan 1 ignored between tabular points)
c 1966 to 1971  frequency offset was kept at   -300
c 1972 to ----  frequency offset was kept at   0
c               linear formulae used from from 1966 Jan 1 to
c                  1972 Jan 1
c               step function used after 1972 Jan 1
c
c        tabular points from 1956 jan 7 to 1958 dec 22
c        julian day number 2435480 to 2436560
      integer*2 a5658(109)/ -8188, -8050, -7918, -7796, -7669, -7539,
     1 -7425, -7521, -7410, -7514, -7433, -7361, -7291, -7223, -7155,
     2 -7087, -7011, -6944, -6896, -6855, -6821, -6582, -6552, -6325,
     3 -6309, -6279, -6034, -5984, -5938, -5907, -5672, -5638, -5393,
     4 -5337, -5281, -5209, -5132, -5052, -4968, -4687, -4605, -4531,
     5 -4468, -4408, -4145, -4080, -4022, -3959, -3886, -3609, -3545,
     6 -3490, -3238, -2983, -2928, -2672, -2419, -2368, -2304, -2031,
     7 -1939, -1857, -1764, -1661, -1556, -1255, -1155,  -875,  -803,
     8  -708,  -598,  -297,  -200,  -103,   186,   276,   362,   644,
     9   926,  1010,  1098,  1182,  1272,  1561,  1649,  1736,  1817,
     1  1890,  1979,  2270,  2361,  2653,  2744,  3034,  3122,  3211,
     2  3302,  3389,  3477,  3565,  3652,  3740,  4024,  4109,  4192,
     3  4275,  4556,  4637,  4724/
c        tabular points from 1959 jan 1 to 1962 mar 6
c        julian day number 2436570 to 2437730
      integer*2 a5862(117)/  5011,  5102,  5193,  5486,  5579,  5671,
     1  5968,  6060,  6154,  6247,  6342,  6435,  6527,  6618,  6708,
     2  6798,  6892,  6983,  7073,  7162,  7250,  7335,  7622,  7709,
     3  7995,  8082,  8168,  8257,  8541,  8630,  8720,  9011,  9100,
     4  9389,  9476,  9766,  9855,  9965, 10095, 10223, 10351, 10480,
     5 10608, 10735, 10862, 10987, 11113, 11242, 11370, 11497, 11625,
     6 11751, 11878, 12004, 12130, 12258, 12385, 12514, 12641, 12767,
     7 12892, 13018, 13144, 13271, 13397, 13523, 13650, 13777, 13904,
     8 14031, 14161, 14290, 14419, 14549, 14728, 14858, 14987, 15116,
     9 15245, 15373, 15502, 15631, 15760, 15889, 16017, 16146, 16275,
     1 16404, 16532, 16661, 16789, 16917, 17045, 17173, 17302, 16931,
     2 17060, 17190, 17321, 17451, 17580, 17709, 17838, 17967, 18097,
     3 18226, 18356, 18486, 18616, 18746, 18869, 18981, 19093, 19205,
     4 19317, 19429, 19541/
c        tabular points from 1962 mar 16 to 1965 jun 28
c        julian day number 2437740 to 2438940
      integer*2 a6265(121)/ 19654, 19767, 19879, 19992, 20104, 20216,
     1 20328, 20440, 20552, 20664, 20775, 20887, 20999, 21111, 21223,
     2 21335, 21447, 21559, 21671, 21783, 21894, 22006, 22119, 22231,
     3 22343, 22455, 22567, 22679, 22791, 22903, 23015, 23126, 23238,
     4 23350, 23462, 23574, 23687, 23799, 23912, 24024, 24136, 24248,
     5 24361, 24473, 24586, 24698, 24810, 24923, 25035, 25147, 25260,
     6 25372, 25484, 25595, 25707, 25819, 25932, 26043, 26155, 26267,
     7 27379, 27491, 27603, 27715, 27827, 27939, 28058, 28188, 28318,
     8 28447, 28577, 28707, 28837, 28966, 29096,   225,   355,   484,
     9   613,   743,   872,  1002,  1131,  1260,  1390,  1519,  1649,
     1  1778,  1907,  2036,  3166,  3295,  3425,  3565,  3694,  3824,
     2  3953,  4083,  4212,  4342,  4471,  4601,  4730,  5859,  5989,
     3  6118,  6248,  6377,  6507,  7637,  7766,  7896,  8026,  8155,
     4  8285,  8415,  8544,  8674,  8803,  8933,  9062/
c        tabular points from 1965 jul 8 to 1966 jan 4
c        julian day number 2438950 to 2439130
      integer*2 a6566( 19)/ 10192, 10322, 10451, 10581, 10710, 10840,
     1 11969, 12099, 12229, 12358, 12488, 12617, 12747, 12876, 13006,
     2 13136, 13265, 13395, 13525/
c$$$$ last entry is what A.1-WWV would have been on 1966 Jan 4
c     if slope had not changed on 1966 Jan 1
 
      equivalence (a, a5658), (a(110), a5862),
     .            (a(227), a6265), (a(348), a6566)
c$$$$ last part of equivalence statement could change if routine revised
 
c points at which discontinuities occurred
c                    amount   between   Julian
c  instant of jump   of jump  tabular   day
c year month day hr  (sec)    points    number
c 1956 Jan    4  19  0.061              2435477   not in range of table
c      March  7  19 -0.020     7   8    2435540
c      March 28  19 -0.020     9  10    2435561
c      July  26  19  0.020    21  22    2435681
c      Aug   22  19  0.020    23  24    2435708
c      Sept  19  19  0.020    26  27    2435736
c      Oct   31  19  0.020    30  31    2435778
c      Nov   14  19  0.020    32  33    2435792
c 1957 Jan   23  19  0.020    39  40    2435862
c      Mar   13  19  0.020    44  45    2435911
c      May   01  19  0.020    49  50    2435960
c      June  05  19  0.020    52  53    2435995
c      June  19  19  0.020    53  54    2436009
c      July  03  19  0.020    55  56    2436023
c      July  17  19  0.020    56  57    2436037
c      Aug   14  19  0.020    59  60    2436065
c      Oct   16  19  0.020    65  66    2436128
c      Nov   06  19  0.020    67  68    2436149
c      Dec   11  19  0.020    71  72    2436184
c 1958 Jan   15  19  0.020    74  75    2436219
c      Feb   05  19  0.020    77  78    2436240
c      Feb   19  19  0.020    78  79    2436254
c      April 09  19  0.020    83  84    2436303
c      June  11  19  0.020    89  90    2436366
c      July  02  19  0.020    91  92    2436387
c      July  16  19  0.020    93  94    2436401
c      Oct   22  19  0.020   102 103    2436499
c      Nov   26  19  0.020   106 107    2436534
c      Dec   24  19  0.020   109 110    2436562
c 1959 Jan   28  19  0.020   112 113    2436597
c      Feb   25  19  0.020   115 116    2436625
c      Aug   05  19  0.020   131 132    2436786
c      Aug   26  19  0.020   133 134    2436807
c      Sept  30  19  0.020   137 138    2436842
c      Nov   04  19  0.020   140 141    2436877
c      Nov   18  19  0.020   142 143    2436891
c      Dec   16  19  0.020   144 145    2436919
c 1961 Jan   01  00  0.005   183 184    2437301
c      Aug   01  00 -0.050   204 205    2437513
c 1963 Nov   01  00  0.100   286 287    2438335
c 1964 April 01  00  0.100   301 302    2438487
c      Sept  01  00  0.100   316 317    2438640
c      Oct   01  00  0.001   319 320    2438670
c 1965 Jan   01  00  0.100   329 330    2438762
c      March 01  00  0.100   335 336    2438821
c      July  01  00  0.100   347 348    2438943
c      Sept  01  00  0.100   353 354    2439005
c 1966 Dec   01  00 -0.0007             2439461   linear formulae
c 1968 Feb   01  00 -0.100              2439888
c 1972 Jan   01  00  0.1076             2441318   step function
c      Jul   01  00  1.0                2441500
c 1973 Jan   01  00  1.0                2441684
c 1974 Jan   01  00  1.0                2442049
c 1975 Jan   01  00  1.0                2442414
c 1976 Jan   01  00  1.0                2442779
c 1977 Jan   01  00  1.0                2443145
c 1978 Jan   01  00  1.0                2443510
c 1979 Jan   01  00  1.0                2443875
c 1980 Jan   01  00  1.0                2444240
c 1981 Jul   01  00  1.0                2444787
c 1982 Jul   01  00  1.0                2445152
c 1983 Jul   01  00  1.0                2445517
c 1985 Jul   01  00  1.0                2446248
c 1988 Jan   01  00  1.0                2447162
c 1990 Jan   01  00  1.0                2447893
c 1991 Jan   01  00  1.0                2448258
c 1992 Jul   01  00  1.0                2448805
c 1993 Jul   01  00  1.0                2449170
c 1994 Jul   01  00  1.0                2449535
c 1996 Jan   01  00  1.0                2450084
c 1997 Jul   01  00  1.0                2450631
c 1999 Jan   01  00  1.0                2451180
c 2006 Jan   01  00  1.0                2453737
c 2009 Jan   01  00  1.0                2454833
c 2012 Jul   01  00  1.0                2456110
c 2015 Jul   01  00  1.0                2457205
c 2017 Jan   01  00  1.0                2457755
c$$$$ above comment table to be extended as routine is revised
c
c      Julian dates of jumps in time           (jd-2430000)
      integer*2 jdtabl(46)/ 5540, 5561, 5681, 5708, 5736, 5778, 5792,
     .          5862, 5911, 5960, 5995, 6009, 6023, 6037, 6065, 6128,
     .          6149, 6184, 6219, 6240, 6254, 6303, 6366, 6387, 6401,
     .          6499, 6534, 6562, 6597, 6625, 6786, 6807, 6842, 6877,
     .          6891, 6919, 7301, 7513, 8335, 8487, 8640, 8670, 8762,
     .          8821, 8943, 9005/
c
c left hand tabular index of jump
      integer*2 jump(46)/7, 9, 21, 23, 26, 30, 32, 39, 44, 49, 52, 53,
     .          55, 56, 59, 65, 67, 71, 74, 77, 78, 83, 89, 91, 93, 102,
     .          106, 109, 112, 115, 131, 133, 137, 140, 142, 144, 183,
     .          204, 286, 301, 316, 319, 329, 335, 347, 353/
c
c amount of jump
      integer*2 dif(46)/2*-200, 34*200, 50, -500, 3*1000, 10, 4*1000/
c
c dates of leap seconds
      integer*4 nleaps/27/, jdleap(27)/    2441500, 2441684, 2442049,
     .          2442414, 2442779, 2443145, 2443510, 2443875, 2444240,
     .          2444787, 2445152, 2445517, 2446248, 2447162, 2447893,
     .          2448258, 2448805, 2449170, 2449535, 2450084, 2450631,
     .          2451180, 2453737, 2454833, 2456110, 2457205, 2457755/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if(jd.lt.2435480) call SUICID(
     .  'CANNOT CALCULATE A1WWV BEFORE JAN 7, 1956. STOP IN A1WWV', 14)
 
      if(jd.ge.2439127) then
c
c linear formula after 1966 Jan 1 (frequency offset-300.E-10)
c and before 1972 Jan 1 (Julian day number 2441318)
         if(jd.lt.2441318) then
            result = 6.2398697_10+0.002592_10*((jd+fract)-2439857._10)
c (jd+fract) = Julian date + 0.5
c 2439856.5 = epoch of linear formula
            if(jd.lt.2439461) then
 
c correction before 1966 Dec 1
               result = result+0.0007_10
            else if(jd.ge.2439888) then
 
c correction after 1968 Feb 1
               result = result - 0.1_10
            endif
         else
c
c constant value after 1972 Jan 1 (Jul.day num. 2441318)
c with one second quantum jumps at the following dates
            result = 10.0343817_10
            ilp    = 0
            do i = 1, nleaps
               if(jd.lt.jdleap(i)) goto 20
               ilp = i
            end do
   20       result = result + ilp
         endif
         goto 300
      else
c
c interpolation formula before 1966 Jan 1
c calculate int
         t   = jd - 2435480
         t   = (t + fract)/10.0_10
         int = t
         t   = t - int
         int = int + 1
c
c check if int is irregular point
         int1 = int + 1
         mjd  = jd - 2430000
         jd0  = 5470 + 10*int
         m1   = 1
         if(jd.ge.2437000) m1 = 37
         do j = m1, 46
            if(jd0.eq.jdtabl(j)) then
               if(j.ge.37) goto 100
               if(mjd.ne.jd0 .or. fract.ge.fract0) goto 40
            else if(int.lt.jump(j)) then
               if(int1.lt.jump(j)) goto 100
               if(int1.ne.jump(j)) goto 50
               if(j.lt.37) goto 100
            else if(int.eq.jump(j)) then
               if(mjd.lt.jdtabl(j)) then
               else if(mjd.eq.jdtabl(j)) then
                  if(fract.gt.fract0 .or. mjd.gt.7000) goto 40
               else
                  goto 40
               endif
            else
               goto 50
            endif
            a1 = a(int)
            a2 = a(int1) - dif(j)
            goto 200
   40       a1 = a(int) + dif(j)
            a2 = a(int1)
            goto 200
   50    end do
  100    a1 = a(int)
         a2 = a(int1)
      endif
  200 if(int.ge.nchang) a1  = a1 + 3E4_10
      if(int1.ge.nchang) a2 = a2 + 3E4_10
c
c interpolation
      result = (a1 + (a2-a1)*t)*1E-4_10
 
  300 A1WWV = result
      return
      end
