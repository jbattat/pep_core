      subroutine KOUTXPR(t,acc,neqtot,dim,xprctl)
      implicit none
c subroutine KOUTXPR - j.f.chandler - 2015 Jul 1
c  print extra information on kout as indicated by control integer

c arguments
c t      - time argument in days
c acc    - array of time derivatives of integrated quantities
c neqtot - 6 times number of integrated vectors
c dim    - shape of acc array (3 or 6 elements per vector)
c xprctl - control integer
      integer neqtot,dim,xprctl
      real*10 t,acc(dim,*)

c commons
      include 'inodta.inc'
      include 'xprcom.inc'

c local
      integer npr,i,k,king,l

      if(MOD(xprctl/32,2).eq.1 .and. MOD(t-1._10,100._10).ne.0._10)
     . return
      if(MOD(xprctl,2).eq.1) then
         npr=neqtot/6
         write(Kout,270) t-.5_10,((acc(i,k),i=dim-2,dim),k=1,npr)
  270    format(f14.5/(1p3d26.19))
      else if(MOD(xprctl,32).gt.0) then
         write(Kout,280) t-.5_10
  280    format(f14.5,' EXTRA PRINT')
      endif
      if(MOD(xprctl/2,2).eq.1)
     . write(Kout,310) Relacc,Relpar
  310 format(' RELATIVISTIC ACCELERATION AND FIRST PARTIAL'/
     . (t13,1p3d22.15))
      if(MOD(xprctl/4,2).eq.1)
     . write(Kout,320) Xprnam,Pcorsav,
     . Xprnam,(((Dadxsav(k,l,1,king,1),k=1,3),l=1,3),
     . king=1,3),Xprnam,((((Dadxsav(k,l,i,king,2),k=1,3),l=1,3),
     . i=1,2),king=1,3)
  320 format(' COORDINATES (POS+VEL) FOR ',3A8/
     . 2(4X,1p3d22.15/),2(8X,1p3d22.15/),2(12X,1p3d22.15/),
     . ' NEWTONIAN GRAVITY GRADIENT W.R.T. ',3A8,' POSITIONS'/
     . 3(4X,1p3d22.15/),3(8X,1p3d22.15/),3(12X,1p3d22.15/),
     . ' RELATIVISTIC GRAVITY GRADIENTS FOR ',3A8,' POS+VEL'/
     . 6(4X,1p3d22.15/),6(8X,1p3d22.15/),(12X,1p3d22.15))
      if(MOD(xprctl/8,2).eq.1)
     . write(Kout,330) Xprnam,Parsav
  330 format(' PARTIALS OF ',3A8,' POS+VEL W.R.T. 1ST PARAMETER'/
     . 2(4X,1p3d22.15/),2(8X,1p3d22.15/),(12X,1p3d22.15))
      if(MOD(xprctl/16,2).eq.1)
     . write(Kout,340) Xprnam,Indpsav
  340 format(' NEWTONIAN INDIRECT PARTIALS DUE TO POS+VEL OF ',3A8/
     . 2(4X,1p3d22.15/),2(8X,1p3d22.15/),2(12X,1p3d22.15/),
     . ' RELATIVISTIC DITTO'/
     . 2(4X,1p3d22.15/),2(8X,1p3d22.15/),(12X,1p3d22.15))
      return
      end

