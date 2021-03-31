      subroutine PERCHK(nstop)
      implicit none

c arguments
      integer*4 nstop
c nstop is a cumulative count of errors

c array dimensions
      include 'globdefs.inc'

c common
      include 'aprtbf.inc'
      integer*4 ijkprt(3,3),ibufa(5),ibufb(5)
      equivalence (ijkprt,Ipert0),(ibufa,Ibuf1),(ibufb,Ibuf6)
      include 'bdctrl.inc'
      include 'fcntrl.inc'
      include 'filtda.inc'
      include 'filtds.inc'
      include 'france.inc'
      include 'inodta.inc'
      integer*4 ijkp(3),ijke(3)
      equivalence (ijkp,Ipert),(ijke,Ieng)
      include 'obsdta.inc'
      integer*4 iob12(10,2)
      equivalence (iob12,Iobs1)
      include 'plndta.inc'
      common/WRKCOM/ Idta(200),Namdta(200)
      character*9 Namdta
      integer*4 Idta

c local variables
      integer*4 i,j,n,ncomp,ntrash,nwrite
      character*1 ijetc(6)/'I','J','K','L','M','N'/,
     . digs(10)/'0','1','2','3','4','5','6','7','8','9'/

      n=0
c datasets that are written and not to be re-read by other parts of PEP
      Idta(n+1)  = Iout
      Namdta(n+1)='IOUT'
      n=1
      Idta(n+1)  = Jout
      Namdta(n+1)='JOUT'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Kout
      Namdta(n+1)='KOUT'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Lout
      Namdta(n+1)='LOUT'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Mout
      Namdta(n+1)='MOUT'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Nout
      Namdta(n+1)='NOUT'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Ipunch
      Namdta(n+1)='IPUNCH'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Ict(77)
      Namdta(n+1)='ICT(77)'
      if(Idta(n+1).gt.1 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Jpunch
      Namdta(n+1)='JPUNCH'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
      Idta(n+1)  = Kpunch
      Namdta(n+1)='KPUNCH'
      if(Idta(n+1).gt.0 .and. Idta(n+1).ne.Iout) n=n+1
c datasets that are written and perhaps read back in immediately, but
c are not to be re-read by other parts or later invocations of PEP
      Idta(n+1)  = Iplcon
      Namdta(n+1)='IPLCON'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Iobcon
      Namdta(n+1)='IOBCON'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Intern
      Namdta(n+1)='INTERN'
      if(Idta(n+1).gt.0) n=n+1
c      Idta(n+1)  = 3
c      Namdta(n+1)='ITHING'
c      if(Idta(n+1).gt.0) n=n+1
c      Idta(n+1)  = in0
c      Namdta(n+1)='IN0'
c      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = MOD(Ibuf,100)
      Namdta(n+1)='IBUF'
      if(Idta(n+1).gt.0 .and. (Ibuf.lt.100 .or. Kkp(13,4).le.0)) n=n+1
      if(Ibuf.gt.99) then
         Idta(n+1)  = Ibuf/100
         Namdta(n+1)='IBUF'''
         if(Idta(n+1).gt.0) n=n+1
      endif
      Idta(n+1)  = Kkp(12,4)
      Namdta(n+1)='KKMR(12)'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Kkp(13,4)
      Namdta(n+1)='KKMR(13)'
      if(Idta(n+1).gt.0) n=n+1
      do i=1,5
         Idta(n+1)  = ibufa(i)
         Namdta(n+1)='IBUF'//digs(i+1)
         if(Idta(n+1).gt.0) n=n+1
      end do
      do i=6,10
         Idta(n+1)  = ibufb(i-5)
         Namdta(n+1)='IBUF'//digs(i+1)
         if(i.ge.10) call EBCDI(i,Namdta(n+1)(5:6),2)
         if(Idta(n+1).gt.0) n=n+1
      end do
      if(Ict(42).gt.0) then
         Idta(n+1)  = Kfile
         Namdta(n+1)='KFILE'
         if(Idta(n+1).gt.0) n=n+1
         Idta(n+1)  = Lfile
         Namdta(n+1)='LFILE'
         if(Idta(n+1).gt.0) n=n+1
         Idta(n+1)  = Mfile
         Namdta(n+1)='MFILE'
         if(Idta(n+1).gt.0) n=n+1
         Idta(n+1)  = Filter
         Namdta(n+1)='FILTER'
         if(Idta(n+1).gt.0) n=n+1
         Idta(n+1)  = Insne
         Namdta(n+1)='INSNE'
         if(Idta(n+1).gt.0) n=n+1
         Idta(n+1)  = Outsne
         Namdta(n+1)='OUTSNE'
         if(Idta(n+1).gt.0) n=n+1
      endif
      ntrash=n

c datasets that are written and may be read back by other parts of PEP
      Idta(n+1)  = Ibody
      Namdta(n+1)='IBODY'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Iem
      Namdta(n+1)='IEM'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Imn
      Namdta(n+1)='IMN'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Inut
      Namdta(n+1)='INUT'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Ilib
      Namdta(n+1)='ILIB'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Kkp(11,4)
      Namdta(n+1)='KKMR(11)'
      if(Idta(n+1).gt.0) n=n+1
      do i=1,u_mxpl
         Idta(n+1)  = Iplnt(i)
         Namdta(n+1)='IPLNT(..)'
         call EBCDI(i,Namdta(n+1)(7:8),2)
         if(Idta(n+1).gt.0) n=n+1
      end do
      Idta(n+1)  = Imat
      Namdta(n+1)='IMAT'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Imat1
      Namdta(n+1)='IMAT1'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Imat2
      Namdta(n+1)='IMAT2'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Imat3
      Namdta(n+1)='IMAT3'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Imat4
      Namdta(n+1)='IMAT4'
      if(Idta(n+1).gt.0) n=n+1
      do i=1,Numobt
         do j=1,2
            Idta(n+1)  = iob12(i,j)
            Namdta(n+1)='IOBS'//digs(j+1)//'(..)'
            call EBCDI(i,Namdta(n+1)(7:8),2)
            if(Idta(n+1).gt.0) n=n+1
         end do
      end do
      if(Ict(42).gt.0) then
         Idta(n+1)  = Iconof
         Namdta(n+1)='ICONOF'
         if(Idta(n+1).gt.0) n=n+1
         Idta(n+1)  = Wzero
         Namdta(n+1)='WZERO'
         if(Idta(n+1).gt.0) n=n+1
      endif
      nwrite=n

c read-only datasets
      do i=1,3
         Idta(n+1)  = ijke(i)
         Namdta(n+1)=ijetc(i)//'ENG'
         if(Idta(n+1).gt.0) n=n+1
      end do
      Idta(n+1)  = Ictat
      Namdta(n+1)='ICTAT'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Libhoc
      Namdta(n+1)='LIBHOC'
      if(Idta(n+1).gt.0) n=n+1
      do i=1,3
         Idta(n+1)  = ijkp(i)
         Namdta(n+1)=ijetc(i)//'PERT'
         if(Idta(n+1).gt.0) n=n+1
         do j=1,3
            Idta(n+1)  = ijkprt(j,i)
            Namdta(n+1)=ijetc(i)//'PERT'//digs(j)
            if(Idta(n+1).gt.0) n=n+1
         end do
      end do
      do i=1,Numobt
         Idta(n+1)  = Iobs0(i)
         Namdta(n+1)='IOBS0(..)'
         call EBCDI(i,Namdta(n+1)(7:8),2)
         if(Idta(n+1).gt.0) n=n+1
      end do
      do i=1,Nummt0
         Idta(n+1)  = Imat0(i)
         Namdta(n+1)='IMAT0(..)'
         call EBCDI(i,Namdta(n+1)(7:8),2)
         if(Idta(n+1).gt.0) n=n+1
      end do
      Idta(n+1)  = Iobs
      Namdta(n+1)='IOBS'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Jmat
      Namdta(n+1)='JMAT'
      if(Idta(n+1).gt.0) n=n+1
      Idta(n+1)  = Ict(42)
      Namdta(n+1)='ICT(42)'
      if(Idta(n+1).gt.1) n=n+1
      Idta(n+1)  = Ict(44)
      Namdta(n+1)='ICT(44)'
      if(Idta(n+1).gt.1) n=n+1
      Idta(n+1)  = Jct(33)
      Namdta(n+1)='JCT(33)'
      Idta(n+2)  = Jct(33)+1
      Namdta(n+2)='JCT(33)+1'
      if(Idta(n+1).gt.0) n=n+2
      if(Ict(42).gt.0) then
         Idta(n+1)  = Smat
         Namdta(n+1)='SMAT'
         if(Idta(n+1).gt.0) n=n+1
      endif

      ncomp=1
      do i=2,n
         do j=1,ncomp
            if(Idta(i).eq.Idta(j)) then
               write(Iout,100) Namdta(i),Namdta(j),Idta(i)
  100          format(' **** DATASET CONFLICT: ',A,' AND ',A,
     .          ' ARE BOTH ',i4)
               if(Mout.gt.0) write(Mout,100) Namdta(i),Namdta(j),Idta(i)
               nstop=nstop+1
            endif
         end do
         if(i.lt.nwrite) then
            ncomp=i
         else
            ncomp=ntrash
         endif
      end do
      return
      end
