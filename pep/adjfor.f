      subroutine ADJFOR(ind,ll,aplnt,npl)
 
      implicit none
c
c        r.b. goldstein  sept 1978
c        called by adjust for fourier shape model parameters
c
c arguments
      character*8 aplnt
      integer   ind,ll
      integer*2 npl
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'inodta.inc'
      include 'scoef4.inc'
      real*10 cfour(4,20,6)
      integer*2 lfour(4,20,6)
      equivalence (Pzhar,cfour),(Lpzhar,lfour)
 
c local
      integer   i,k
      character*1 astrik(3)/'*','&',' '/
      character*2 wrds(6)/'A ','B ','C ','D ','AP','CP'/
      character*2 pp/'PP'/
      real*10 fact/1._10/
      character*2 stw(4)/'  ',' C','OE','FS'/
 
      do k = 1,6
         call LINCHK
         do i = 1,20
            if(lfour(ind,i,k).gt.0) then
               call ADJAST(cfour(ind,i,k),fact)
               write(Iout,10) astrik(Ntype),N,wrds(k),i,wrds(k),
     .          i,aplnt,npl,cfour(ind,i,k),Adj,Nwv,Sig,Fract
   10          format(1x,a1,i4,'. LF',a2,'(',i2,')   = 1 F',a2,
     .                '(',i3,') OF ',a8,i3,8x,1pd22.15,1pd16.8,
     .                1pd22.15,1pd10.3,0pf8.3)
               if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,
     .          wrds(k),i,wrds(k),i,aplnt,npl,cfour(ind,i,k),
     .          Adj,Nwv,Sig,Fract
               if(Keepit) cfour(ind,i,k) = Nwv
            endif
         end do
         stw(1) = wrds(k)
         call FRSTAT(1,2,stw)
      end do
c
c
      do i = 121,122
         if(Lpzhar(ind,i).gt.0) then
            call LINCHK
            call ADJAST(Pzhar(ind,i),fact)
            write(Iout,20) astrik(Ntype),N,wrds(i - 120),
     .                      wrds(i - 120),aplnt,npl,Pzhar(ind,i),
     .                      Adj,Nwv,Sig,Fract
   20       format(1x,a1,i4,'. LF',a1,'PP      = 1 F',a1,
     .             'PP     OF ',a8,i3,8x,1pd22.15,1pd16.8,
     .             1pd22.15,1pd10.3,0pf8.3)
            if(Jout.gt.0) write(Jout,20) astrik(Ntype),N,
     .                              wrds(i - 120),wrds(i - 120),
     .                              aplnt,npl,Pzhar(ind,i),Adj,
     .                              Nwv,Sig,Fract
            if(Keepit) Pzhar(ind,i) = Nwv
         endif
      end do
      stw(1) = pp
      call FRSTAT(1,2,stw)
      call FRSTAT(2,3,'ALL COEF')
 
      return
      end
