      real*10 function EXANM(v,e)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 de, dm, e, ecose, esine, v
 
c*** end of declarations inserted by spag
 
 
c
c w.taylor   may 1974   real*10 function exanm
c
 
      EXANM = v
      do while( .true. )
         esine = e*SIN(EXANM)
         ecose = e*COS(EXANM)
 
         dm = v - EXANM + esine
         de = dm/(1._10 - ecose + ((.5_10*dm)/(1._10-ecose))*esine)
 
         if(ABS(de).gt.1.) de = ABS(de)/de
         EXANM = EXANM + de
 
         if(ABS(de).le.2.E-11_10) return
      end do
      end
