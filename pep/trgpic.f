      subroutine TRGPIC(kick)
 
      implicit none
c
c r.j. cappallo     october 1972     sr trgpic
c
c
c     trgpic calculates partial derivatives of observations of
c     planets wrt target planet initial conditions

c arguments
      integer*4 kick

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdat.inc'
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      include 'mtrapx.inc'
      include 'number.inc'
      include 'partcm.inc'
      include 'pemctl.inc'
      include 'tapdta.inc'

c local
      integer*4 i,iflag,j,l,ll,lt,lte,ltg0,ltype,mm,mt
      integer*2 lplext
 
      ll     = 1
      mm     = 1
      lt     = 8
      lte    = 8
      Kembry = Lparem
      Kplnt  = Lparp
      do while( .true. )
 
         if(ll.gt.Numtar) goto 900
 
         ltg0 = Ntrg(ll)*100
c
c see if target body is on obs.lib tape
         mt = 0
         if(Iabs1.gt.0 .and. mm.le.Mumtar) then
            if(Mtrg(mm).eq.Ntrg(ll)) mt = 1
         endif
c
c start of loop for target body ic partials
         l = 0
         call PCOPS(7,'TBOD',mt)
         do while( .true. )
            if(l.lt.6) then
               iflag = 0
               call PCOPY(l,6,iflag,1,Ltbod(1,ll),Mtbod(1,mm))
               if(iflag.le.0) then
                  ltype  = l
                  lplext = ltg0 + ltype
c
c interpolate for target body ic partial from embary tape if its
c there. note that if target body affects earth motion the partial
c must be integrated for correct results.
                  iflag = -1
                  call PBDPRM(Nkiem,Kiem,lte,Kembry,lplext,iflag)
                  if(iflag.gt.0) then
c
c correct partial not on earth tape, assume to be zero
                     do i = 1,Mouse
                        do j = 1,Index
                           derem(j,i) = 0._10
                        end do
                     end do
                  else
                     call CPARTL(1,1,kick)
                  endif
c
c interpolate for target body ic partial from planet tape, if any
                  if(Klan.gt.0) then
                     iflag = -1
                     call PBDPRM(Nkipl,Kipl,lt,Kplnt,lplext,iflag)
                     if(iflag.gt.0) call SUICID(
     .' TARGET BODY PARTIAL NOT ON PLANET TAPE, STOP IN TRGPIC ',14)
                     call CPARTL(3,1,kick)
                     if(kick.ne.3) then
                        do i = 1,Mouse
                           do j = 1,Index
                              derem(j,i) = derem(j,i) - Derpl(j)
                           end do
                        end do
                        Ivze(1) = 1
                     endif
                  endif
                  call CPARTC(kick)
               endif
            endif
 
            if(l.ge.6) then
               if(mt.gt.0) mm = mm + 1
               ll = ll + 1
               goto 100
            endif
         end do
  100 end do
 
  900 return
      end
