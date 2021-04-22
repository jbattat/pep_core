      subroutine PSRPUN(ncard)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, nc, nl
 
c*** end of declarations inserted by spag
 
 
c psrpun - j.f.chandler - 1984 jul
c punch adjusted pulsar parameters

c arguments
      integer*4 ncard

c array dimensions
      include 'globdefs.inc'

c
c common
      include 'inodta.inc'
      include 'psrstf.inc'
 
      if(Numpsr.le.0) return
      do j = 1, Numpsr
         if(Lpsrcn(1,j).gt.0) then
            nl = u_nmpsr
            do k = 2, u_nmpsr
               if(Lpsrcn(k,j).le.0) then
                  nl = k - 1
                  goto 20
               endif
            end do
   20       nc = u_nmpsr
            do i = 1, u_nmpsr
               if(Psrcn(nc,j).ne.0._10) goto 40
               nc = nc - 1
            end do
   40       write(Ipunch,60) Sptpsr(j),Jdpsr0(j),Plspr(j),Ntypsr(j)
     .                        , (Lpsrcn(i,j),i = 1,nl)
   60       format('*OBJECT'/' NAME=''', a4, ''', JD0=', i8,
     .             ', CON1(2)=', 1pe24.17, ', NTYPE=',
     .             i2/' PULSAR=T,'/' L=', 20(i2,','))
            ncard = ncard + 4
            if(nc.gt.0) then
               write(Ipunch,70) (Psrcn(i,j),i = 1,nc)
   70          format(' CON=', (t6,3(1pe24.17,',')))
               ncard = ncard + 1 + (nc - 1)/3
            endif
         endif
      end do
 
      return
      end
