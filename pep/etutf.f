      real*10 function ETUTF(jd, fract)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   int, intmax, jd
 
c*** end of declarations inserted by spag
 
 
c  M.E.Ash   April 1967    real*10 function ETUTF
c  Determination of ephemeris time minus universal time in seconds
c  by table interpolation or extrapolation.
c
c     before 1799.5  constant value
c     1799.5-1954.5 Clemence table
c     1955.5-2014.5 32.15+A1-UT2
c     2014.5-       extrapolation
c
      real*10 fract, t
c jd   = input julian day number
c fract= input fraction of a day from midnight (universal time)
      real*4    etutc(216), a1800(100), a1900(100), a2000(16)
      data intmax/216/
      equivalence (a1800,etutc), (a1900,etutc(101)), (a2000,etutc(201))
c The first item in a1800 is for 1799.5...
      data a1800/
     .  5.90, 5.80,  5.70, 5.60,  5.50, 5.40,  5.30, 5.20,  5.10, 5.00,
     1  4.90, 4.80,  4.70, 4.76,  4.83, 4.90,  4.97, 5.04,  5.11, 5.18,
     2  5.25, 5.32,  4.80, 4.32,  3.88, 3.49,  3.15, 2.85,  2.59, 2.37,
     3  2.20, 2.01,  1.77, 1.49,  1.16, 0.79,  0.36,-0.04, -0.28,-0.37,
     4 -0.29,-0.06,  0.25, 0.56,  0.85, 1.13,  1.39, 1.65,  1.89, 2.12,
     5  2.33, 2.53,  2.72, 2.90,  3.06, 3.22,  3.35, 3.46,  3.50, 3.51,
     6  3.43, 3.32,  3.14, 2.90,  2.61, 2.26,  1.85, 1.39,  0.87, 0.20,
     7 -0.68,-1.78, -3.09,-4.48, -5.65,-6.57, -7.24,-7.67, -7.87,-8.00,
     8 -8.09,-8.14, -8.17,-8.17, -8.14,-8.07, -7.97,-7.84, -7.67,-7.58,
     9 -7.58,-7.67, -7.85,-8.04, -8.07,-7.93, -7.63,-7.19, -6.57,-5.80/
      data a1900/
     . -4.87,-3.79, -2.54,-1.13, +0.35, 1.80,  3.26, 4.69,  6.11, 7.51,
     1  8.90,10.28, 11.64,12.95, 14.18,15.31, 16.39,17.37, 18.27,19.08,
     2 19.83,20.48, 21.06,21.56, 21.97,22.29, 22.55,22.72, 22.82,22.92,
     3 23.05,23.18, 23.34,23.50, 23.60,23.64, 23.63,23.58, 23.63,23.76,
     4 23.99,24.30, 24.71,25.15, 25.61,26.08, 26.57,27.08, 27.61,28.15,
     5 28.94,29.42, 29.66,30.29, 30.96,31.09, 31.23,31.47, 31.88,32.40,
     6 32.88,33.34, 33.79,34.23, 34.74,35.39, 36.15,36.99, 37.87,38.75,
     7 39.69,40.69, 41.67,42.81, 43.94,44.98, 45.96,46.98, 48.02,49.08,
     8 50.09,50.96, 51.79,52.55, 53.42,54.07, 54.61,55.08, 55.56,56.08,
     9 56.55,57.20, 57.94,58.72, 59.57,60.38, 61.23,61.98, 62.64,63.29/
      data a2000/
     . 63.64,63.96, 64.19,64.39, 64.53,64.63, 64.77,64.97, 65.32,65.61,
     1 65.93,66.22, 66.46,66.75, 67.11,67.47/
c
c 1799.5 = 2378314
      if(jd.ge.2378314) then
c
c calculate et-ut within table (linear interpolation)
         t   = jd - 2378314
         t   = (t + fract)/365.25_10
         int = t
         t   = t - int
         int = int + 1
         if(int.ge.intmax) then
c linear extrapolation beyond table
            t   = t + (int - intmax + 1)
            int = intmax - 1
         endif
         ETUTF = etutc(int) + (etutc(int+1) - etutc(int))*t
      else
         ETUTF = 6.0_10
      endif
      return
      end
