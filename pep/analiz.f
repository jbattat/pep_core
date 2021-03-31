      subroutine ANALIZ
 
      implicit none
 
c     r.reasenberg/d.white   jan 1974   subroutine analiz
c     main program for maximum liklihood estimator (least squares
c     analysis) or kalman-bucy filter
c         main program for forming, solving, and using normal equations
c
      include 'aprtbf.inc'
      include 'ciptrx.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'nrmmat.inc'
      include 'obsdta.inc'
      include 'scail.inc'
      integer*2 iptr(1000)
      equivalence (Scale,iptr)
c
c local
      integer*2 icv(12)
c
c initialize
c hedskp = .true. causes nrmfrm (via frmhed) to skip m-vector
c records
      Hedskp = .false.
      Pprdon = .false.
c
c write title page for least squares analysis
      call PAGSET('LEAST SQUARES ANALYSIS  ', 6)
      call NEWPG
      call ANAZNM
c
c shall analiz link computations be performed
      if(Ict(15).lt.0) then
         write(Iout,100) Ict(15)
  100    format(//' ANALIZ LINK CANCELLED BY ICT(15)=', i3)
         Ict(1) = Iterat - 1
         return
      endif
c
c translate multi-parameter & constraint controls
      if(Iterat.eq.1) then
         call PRMTRN
      else if(Ibuf5.gt.0 .and. Itrwnd(Ibuf5).gt.0) then
         rewind Ibuf5
         Itrwnd(Ibuf5)=-2
      endif
c
c set up parameter counts and pointers
      call ANSET
c
c clear normal equations
      call NRMSET(0)
      if(Iobcon.gt.0) rewind Iobcon
c
c prepare normal equations
      if(Iterat.gt.1 .or. Ict(5).le.0) call NRMICT
c
c filter normal equations
      if(Ict(42).gt.0) then
         call FILCTL(B)
         return
      endif
c
c add other info in saved norm eqn form and solve
      call SOLCTL(B,.false.)
      if(Pprdon) return
      if(Ict(15).le.0) then
c
c adjust parameters and prepare other applications of solut.
         call SOLUSE(B,Jct(51).le.0)
c
c find new rms for each set of normal equations
         if(Jct(51).lt.0) then
            call ZFILL(icv,2*12)
            icv(2) = -1
            icv(6) = 8
            icv(7) = 2
            if(Ict(5).eq.0 .or. (Ict(5).gt.0 .and. Iterat.le.1)) then
               if(Iobcon.gt.0) icv(2) = 0
               call NRMFRM(icv,iptr,B)
            else if(Imat1.gt.0) then
               icv(1) = 1
               call NRMFRM(icv,iptr,B)
            endif
c
c need to restore arrays if going to do prediction
            if(Ict(10).ge.-1) then
               if(Ict(14).le.0) Ict(14)=1
               call SOLUS2(B,.true.)
            endif
         endif
      endif
c
c initiate multiple parameter set runs
      if(Jct(51).gt.0) call ANALIX
      return
      end
