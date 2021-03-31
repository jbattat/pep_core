      subroutine ADJPHS(nphase,nplnt,name)
 
      implicit none
c
c ash/forni october 1967  subroutine adjphs
c adjust planetary optical observation phase corrections
c
c arguments
      integer   nphase
      character*8 name
      integer*2 nplnt
c           nphase= phase counter
c           nplnt = planet number
c           name = eight character planet name
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'inodta.inc'
      include 'phase.inc'
c
c internal to subroutine adjphs
      real*4    snwv
      integer   i,n2
      real*10 old
      real*10 fctphs/1._10/
      character*1 astrik(3)/'*','&',' '/
 
      do while( nphase.lt.Numphs )
 
         if(nplnt.ge.0) then
 
            if(nplnt.ne.Nplphs(nphase+1)) return
            nphase = nphase + 1
         else
            nphase = nphase + 1
            if(Nplphs(nphase).ge.0) goto 100
         endif
 
         call LINCHK
         n2 = Ncphs(nphase)
         if(n2.gt.0) then
            do i = 1,n2
               if(Lphs(i,nphase).gt.0) then
                  old = Aphase(i,nphase)
                  call ADJAST(old,fctphs)
                  snwv = Nwv
                  write(Iout,10) astrik(Ntype),N,i,nphase,
     .                            Lphs(i,nphase),i,name,
     .                            Phsit(nphase),Phser(nphase),
     .                            Aphase(i,nphase),Adj,snwv,Sig,
     .                            Fract
   10             format(1x,a1,i4,'. LPHS(',i1,',',i3,')= ',i1,
     .                   ' APHASE(',i1,') ',2x,a8,1x,a4,1x,a4,
     .                   6x,1pe12.5,5x,1pd16.8,5x,1pe12.5,5x,
     .                   1pd10.3,0pf8.3)
                  if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,
     .               i,nphase,Lphs(i,nphase),i,name,
     .               Phsit(nphase),Phser(nphase),Aphase(i,nphase),
     .               Adj,snwv,Sig,Fract
                  if(Keepit) Aphase(i,nphase) = snwv
               endif
            end do
 
            call FRSTAT(1,2,'OP PHASE')
         endif
  100 end do
 
      return
      end
