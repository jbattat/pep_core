      block data BLKDATA8
 
      implicit none
 
c
c r.reasenberg/d.white   jan 1974   block data for analiz link
c
c Set up B and Bz matrices and dimensions thereof
c Maximum allowable parameters = 820
c 820 * (820+1) / 2 = 336610
      include 'numnum.inc'
      data Nrmsiz/336610/
      include 'matzro.inc'
      real*10 Bz(336610)
      equivalence (Bz,Coeff)
      common /NRMMAT/ B
      real*10 B(336610)
      end
