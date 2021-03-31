      block data BLKDATA9
 
      implicit none
 
 
      include 'rotcom.inc'
c        r. goldstein   july 1977
c        mars rotation model constants
c        ref:  reasenberg and king, the rotation of mars
c           to be submitted to jgr
c
c        coeficients for  mars rotation have been calculated for
c        eccentricity=0.093384664 obtained from esae(1961)
c        pg.113 for a data jan 1 1978 (jd2443510)
c
c
      data Mmot/.009146100372_10/
      data Qqdot/2.474515959_10, 6.5412E-7_10/
      data Ecof/1.0132250881_10, 0.2829314855_10, 0.03951222901_10,
     . 0.005417067576_10, 7.330557876E-4_10, 9.837E-5_10, 1.313E-5_10,
     . 0._10, -0.04664152593_10, 0.9782599720_10, 0.3206129166_10,
     . .07267658075_10, 0.01403569802_10, 2.476119213E-3_10/
      end
