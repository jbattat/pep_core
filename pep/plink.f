      subroutine PLINK
 
      implicit none
 
 
      include 'inodta.inc'
      include 'loadid.inc'
 
      write(Iout, 100)
  100 format('-'/'-'/42x,
     .       '**           **    **           **   **       **'/42x,
     .       '**           **    ***          **   **      ** '/42x,
     .       '**           **    ****         **   **     **  '/42x,
     .       '**           **    ** **        **   **    **   '/42x,
     .       '**           **    **  **       **   **   **    '/42x,
     .       '**           **    **   **      **   *** **     '/42x,
     .       '**           **    **    **     **   *****      '/42x,
     .       '**           **    **     **    **   ** **      '/42x,
     .       '**           **    **      **   **   **  **     '/42x,
     .       '**           **    **       **  **   **   **    '/42x,
     .       '**           **    **        ** **   **    **   '/42x,
     .       '**           **    **         ****   **     **  '/42x,
     .       '***********  **    **          ***   **      ** '/42x,
     .       '***********  **    **           **   **       **')
      write(Iout, 200) Lnkdat, Lnkdsn
  200 format('- PEP VERSION= ', A8, 1x, A44)
      return
      end
