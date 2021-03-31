      subroutine ADJHAR(lhar,har,nn,kk,ll,jj,chars,chart,name,
     .                  npl,fcthar)
 
      implicit none
c
c m.e.ash   june 1969   subroutine adjhar
c adjust gravitational potential harmonic coefficients
c
c arguments
      integer*4 nn,kk,ll,jj
      integer*2 lhar(nn,100),npl
      real*10 har(nn,100),fcthar(100)
      character*8 name
      character*(*) chars
      character*4 chart
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adjstf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'

c local
      integer   j,j1
      character*1 astrik(3)/'*','&',' '/
 
      do j = 1,ll
         j1 = j + jj
         if(lhar(kk,j).gt.0) then
            if(Jct(78).gt.0) then
               call ADJAST(har(kk,j),1._10)
            else
               call ADJAST(har(kk,j),fcthar(j))
            endif
            if(jj.gt.0) then
               write(Iout,10) astrik(Ntype),N,chars,
     .                         j1,chart,lhar(kk,j),j1,name,npl,
     .                         har(kk,j),Adj,Nwv,Sig,Fract
   10          format(1x,a1,i4,'. L',a2,i2,a4,'  =',i2,' J',
     .                i2,5x,' OF ',1A8,i3,8x,1pd22.15,1pd16.8,
     .                1pd22.15,1pd10.3,0pf8.3)
               if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,
     .            chars,j1,chart,lhar(kk,j),j1,
     .            name,npl,har(kk,j),Adj,Nwv,Sig,Fract
            else
               write(Iout,20) astrik(Ntype),N,chars,j1,chart,
     .                         lhar(kk,j),chars(2:4),j,
     .                         name,npl,har(kk,j),Adj,Nwv,Sig,
     .                         Fract
   20          format(1x,a1,i4,'. L',a5,i2,a3,'=',i2,
     .                1x,a3,'(',i2,')  OF ',1A8,i3,8x,
     .                1pd22.15,1pd16.8,1pd22.15,1pd10.3,0pf8.3)
               if(Jout.gt.0) write(Jout,20) astrik(Ntype),N,
     .            chars,j1,chart,lhar(kk,j),chars(2:4),
     .            j,name,npl,har(kk,j),Adj,Nwv,Sig,Fract
            endif
            if(Keepit) har(kk,j) = Nwv
         endif
      end do
 
      return
      end
