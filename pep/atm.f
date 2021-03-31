      subroutine ATM(sitf, month, day, wetz, dryz)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      arg, c, delt1, delt2, delt3, delt4, delt5, delt6, delt7,
     .          delt8, delt9, delta, e, frac, fric, gam, gam1,
     .          gam2, gam3
      real      gam4, gam5, gamma, p0, p1, p2, p3, p4, p5, pres, rh,
     .          rh1, rh2, rh3, rh4, rh5, rhum, scale
      real      t0, t1, t2, t3, t4, t5, temp
      integer   i, ii, ik, im, jj, k, l
 
c*** end of declarations inserted by spag
 
 
c
c g.slater and s.synnott   aug 1972   subroutine atm
c program to compute standard zenith correction for neutral
c atmosphere delay given month day,and site names
c ref: jpl iom 391.3-352 25may 1971
c
c program assumes monthly values at 15th of month and linearly
c interpolates  to obtain value for day.
c
      character*4 sitf(2,2),sitn(2),sitnam
      integer*2 month, day
      real*4 wetz(2), dryz(2)
      dimension p1(12), p2(12), p3(12), p4(12), p5(12)
      dimension t1(12), t2(12), t3(12), t4(12), t5(12)
      dimension rh1(12), rh2(12), rh3(12), rh4(12), rh5(12)
      dimension gam1(12), gam2(12), gam3(12), gam4(12), gam5(12)
      dimension p0(12, 5), gam(12, 5), rh(12, 5), t0(12, 5), delta(4, 9)
      dimension delt1(4), delt2(4), delt3(4), delt4(4), delt5(4),
     .          delt6(4), delt7(4), delt8(4), delt9(4)
      equivalence(delt1, delta(1,1)), (delt2, delta(1,2)),
     .            (delt3, delta(1,3)), (delt4, delta(1,4)),
     .            (delt5, delta(1,5)), (delt6, delta(1,6)),
     .            (delt7, delta(1,7)), (delt8, delta(1,8)),
     .            (delt9, delta(1,9))
      equivalence(p1, p0(1,1)), (p2, p0(1,2)), (p3, p0(1,3)),
     .            (p4, p0(1,4)), (p5, p0(1,5)), (t1, t0(1,1)),
     .            (t2, t0(1,2)), (t3, t0(1,3)), (t4, t0(1,4)),
     .            (t5, t0(1,5)), (rh1, rh(1,1)),
     .            (rh2, rh(1,2)), (rh3, rh(1,3)),
     .            (rh4, rh(1,4)), (rh5, rh(1,5)),
     .            (gam1, gam(1,1)), (gam2, gam(1,2)),
     .            (gam3, gam(1,3)), (gam4, gam(1,4)),
     .            (gam5, gam(1,5))
      character*4 ds12/'12DS'/, ds13/'13DS'/, ds14/'14DS'/,
     .          ds41/'41DS'/,
     .          ds51/'51DS'/, ds61/'61DS'/, ds62/'62DS'/, ds42/'42DS'/,
     .          ds11/'11DS'/, ds43/'43DS'/, ds63/'63DS'/
      data delt1/ - 32.0, 0., -2.3, 0./
      data delt2/ - 30.2, 0., -2.2, 0./
      data delt3/ - 38.9, 0., -2.8, 0./
      data delt4/ - 30.2, 0., -2.2, 0./
      data delt5/1.4, 0., 0.1, 0./
      data delt6/ - 6.7, 0., -0.5, 0./
      data delt7/ - 18.0, 0., -1.2, 0./
      data delt8/ - 18.5, 0., -1.2, 0./
      data delt9/ - 43.5, 0., -2.9, 0./
c     p0(i,k) avg. static pressure month i, station k
c          gam,rh,t0 are temp. gradient,rel. hum. and temp.
c          k=1 edwards afb
c          k=2 woomera
c          k=3 pretoria
c          k=4 madrid
c          k=5 wagga
c
c        delta's are corrections to station values for each of dsn sites
c                  j=1 ds11
c                  j=2 ds12
c                  j=3 ds13
c                  j=4 ds14
c                  j=5 ds41
c                  j=6 ds51
c                  j=7 ds61
c                  j=8 ds62
c                  j=9 ds42
c
c     edwards afb data
      data p1/936.92, 936.86, 932.90, 930.23, 928.88, 925.66, 931.11,
     .     930.74, 930.04, 934.58, 933.82, 932.78/
      data gam1/6.201, 6.745, 6.526, 6.076, 7.008, 6.401, 7.110, 7.533,
     .     6.253, 6.230, 6.011, 6.178/
      data t1/292.28, 290.32, 288.38, 288.36, 292.98, 295.16, 302.70,
     .     303.80, 295.61, 298.16, 291.48, 284.74/
      data rh1/.2195, .1969, .3241, .2316, .2891, .2544, .3005, .2986,
     .     .3920, .2088, .3908, .3126/
c
c woomera data
      data p2/991.06, 993.32, 1002.45, 1000.12, 1002.50, 1003.65,
     .     999.40, 995.21, 1000.52, 999.30, 996.20, 995.41/
      data gam2/7.2571, 5.794, 5.591, 7.103, 5.856, 6.2762, 5.4954,
     .     5.2133, 5.6135, 5.7513, 5.3609, 6.1352/
      data t2/304.84, 298.58, 291.44, 297.56, 288.82, 290.03, 283.61,
     .     284.00, 288.68, 291.37, 289.70, 294.97/
      data rh2/.3058, .3029, .4555, .3292, .4709, .3734, .4378, .3467,
     .     .3190, .3529, .3257, .3335/
c
c pretoria data
      data p3/862.24, 866.00, 864.05, 862.71, 867.30, 870.81, 869.80,
     .     870.52, 869.10, 865.33, 865.50, 865.73/
      data gam3/4.823, 6.069, 5.627, 6.299, 5.843, 5.554, 6.562, 6.8898,
     .     6.102, 6.250, 6.237, 4.724/
      data t3/292.61, 293.22, 289.21, 292.03, 288.52, 285.95, 289.75,
     .     291.20, 289.25, 292.11, 292.80, 288.00/
      data rh3/.4559, .4786, .5294, .4119, .3833, .3239, .2484, .3325,
     .     .4294, .5285, .6214, .4810/
c
c madrid data
      data p4/949.87, 946.76, 945.23, 941.62, 942.93, 946.08, 945.05,
     .     944.43, 944.68, 946.24, 944.17, 947.74/
      data gam4/6.663, 6.598, 6.528, 6.722, 5.839, 6.250, 6.552, 6.230,
     .     6.011, 5.945, 6.040, 6.089/
      data t4/285.34, 283.87, 287.29, 285.17, 286.29, 291.28, 298.33,
     .     295.7, 292.29, 290.27, 285.12, 282.62/
      data rh4/.4723, .5164, .4064, .5508, .4994, .4287, .3477, .3542,
     .     .3866, .5222, .7448, .8157/
c
c wagga data
      data p5/990.00, 994.30, 992.05, 987.00, 993.15, 994.65, 995.87,
     .     988.01, 994.87, 992.70, 988.62, 989.54/
      data t5/295.14, 296.34, 293.89, 288.19, 285.11, 287.12, 283.32,
     .     282.42, 286.12, 289.00, 290.72, 293.10/
      data rh5/.3867, .3412, .3318, .4397, .5101, .4413, .4880, .6227,
     .     .4835, .4774, .3559, .3875/
      data gam5/5.991, 5.692, 5.709, 5.246, 5.5118, 6.135, 6.004, 6.303,
     .     6.289, 6.007, 5.561, 6.3419/
      data c/2.997925E+8/
c
c stdn data
      real*4    ndds11
 
c nd is surface refractivity nd(imonth,1site)
      character*4 stdn(12)/'BDA3','CYI3','CRO3','HAW3','TEX3','MAD8',
     .          'GWM3', 'HSK8', 'GDS8', 'MIL3', 'ACN3', 'ETC3'/
      real*4    nd(12, 12)/335., 333., 335., 343., 353., 372., 377.,
     .          378., 371., 359., 345., 338., 333., 336., 336., 335.,
     .          338., 346., 343., 352., 354., 350., 343., 336., 346.,
     .          348., 344., 334., 327., 329., 325., 322., 325., 326.,
     .          329., 339., 299., 298., 297., 299., 301., 305., 308.,
     .          310., 309., 309., 304., 302., 328., 333., 340., 350.,
     .          362., 377., 395., 375., 372., 353., 336., 334., 295.,
     .          293., 295., 294., 298., 302., 299., 300., 306., 304.,
     .          299., 296., 372., 377., 378., 383., 381., 377., 373.,
     .          379., 381., 380., 375., 379., 306., 312., 309., 293.,
     .          300., 297., 296., 295., 298., 301., 300., 301., 279.,
     .          281., 278., 279., 279., 277., 279., 279., 278., 277.,
     .          276., 279., 338., 339., 341., 348., 362., 378., 385.,
     .          386., 381., 361., 349., 339., 361., 364., 366., 363.,
     .          361., 356., 353., 353., 354., 356., 355., 357., 265.,
     .          261., 259., 255., 257., 266., 283., 286., 273., 271.,
     .          265., 265./
c
c i=1 receive site
c i=2 send site
      do i = 1, 2
         sitn(i) = sitf(1, i)
      end do
      if( day .ge. 15. ) then
         im = month
         ik = month + 1
         if( ik .eq. 13 ) ik = 1
         frac = (day - 15.)/30.
      else
 
c day less than 15
         im = month - 1
         ik = month
         if( im .eq. 0 ) im = 12
         frac = (day + 15.)/30.
      endif
      fric = 1. - frac
      do i = 1, 2
         sitnam = sitn(i)
         scale  = 1.0
         k  = 1
         jj = 1
         if( sitnam .ne. ds11 ) then
            jj = 2
            if( sitnam .ne. ds12 ) then
               jj = 3
               if( sitnam .ne. ds13 ) then
                  jj = 4
                  if( sitnam .ne. ds14 ) then
                     k  = 2
                     jj = 5
                     if( sitnam .ne. ds41 ) then
                        k  = 3
                        jj = 6
                        if( sitnam .ne. ds51 ) then
                           k  = 4
                           jj = 7
                           if( sitnam .ne. ds61 ) then
                              if( sitnam .ne. ds63 ) then
                                 jj = 8
                                 if( sitnam .ne. ds62 ) then
                                    k  = 5
                                    jj = 9
                                    if( sitnam .ne. ds42 ) then
                                       if( sitnam .ne. ds43 ) then
c
c set up goldstone values for stdn or non-standard site
                                         e = 7.4475*(t0(im,1) - 273.16)
     .                                      /(t0(im,1) - 35.46)
                                         e = 6.1*rh(im, k)*10**e
                                         ndds11 = 77.6*p0(im, 1)
     .                                      /t0(im, 1)
     .                                      + 77.6*4810.*e/t0(im, 1)**2
                                         k  = 1
                                         jj = 1
 
c see if stdn site
                                         do ii = 1, 12
                                         if( sitnam .eq. stdn(ii) )
     .                                     then
                                         l     = ii
                                         scale = nd(im, l)/ndds11
                                         go to 50
                                         endif
c
c non-standard site , use goldstone values (scale=1.0)
                                         end do
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
   50    pres    = fric*p0(im, k) + frac*p0(ik, k) + delta(1, jj)
         gamma   = fric*gam(im, k) + frac*gam(ik, k) + delta(2, jj)
         temp    = fric*t0(im, k) + frac*t0(ik, k) + delta(3, jj)
         rhum    = fric*rh(im, k) + frac*rh(ik, k) + delta(4, jj)
         dryz(i) = 2.2757E-3*pres/c
         arg     = (17.1486*temp - 4684.1331)/(temp - 38.45)
         wetz(i) = 0.5667*(rhum/gamma)*(1. - 38.45/temp)**2*EXP(arg)/c
         dryz(i) = dryz(i)*scale
         wetz(i) = wetz(i)*scale
         if( sitn(1) .eq. sitn(2) ) then
            dryz(2) = dryz(1)
            wetz(2) = wetz(1)
            return
         endif
      end do
      return
      end
