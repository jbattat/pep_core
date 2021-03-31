      subroutine FILCTL(b)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, lepoch, m, mode, n, ncol, np, npnpx2
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine filctl
c paul macneil january, 1978 modified for iteration
c
 
c parameter is b matrix - dim set in block data
      real*10 b(1)
c
c common
      include 'aprtbf.inc'
      include 'fcntrl.inc'
      include 'filptr.inc'
      include 'filtda.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'rtside.inc'
c
c local
      real*10 jd1, jd2, jda, jdb
      real*10 buff(1000)
      equivalence (buff, Sigma)
      logical*4 keepit
      logical*4 pvoflg
c
c form filter
      npnpx2 = 2*Npnp
c during pvo extended starting procedure, iconof now contains
c forward differences only.  backward differences will be
c needed by nomtrs, so extend iconof
      if(Jct(56).eq.1 .and. (Filflg(1).or.Fict(8).gt.0)) then
         pvoflg = .true.
         call WUPDAT(Iconof, Mfile, buff, npnpx2, Nepoch, Itrwnd, Side,
     .               pvoflg)
         call FLICPM(Iconof, Mfile, buff, npnpx2, Nepoch, Itrwnd)
      endif
      rewind Insne
      rewind Filter
      rewind Smat
 
c extended precision, double amount zeroed
      Nparam = Nparam*2
      call NRMSET(1)
      Nparam = Nparam/2
      call FORM(b)
c
c propagate residual, write saved norm eqn out
      rewind Insne
      rewind Filter
      rewind Outsne
 
c extended precision, double amount zeroed
      Nparam = Nparam*2
      call NRMSET(1)
      Nparam = Nparam/2
      call USE(b)
      rewind Outsne
c
c find epoch for which solutions are to be kept
      do i = 1, Nepoch
         if(Lfilt(i).eq.1) go to 100
         end do
  100 lepoch = i
c
c read heading for filtered normal eqn data set
      call FILHED(np, Npnp, Nepoch)
c
c loop through epochs
c order of solutions on outsne is from last to first epoch
c only one solution for forward filtering
      keepit = .false.
      n = Nepoch
      if(Fict(2).eq.0) n = 1
      do i = 1, n
         Ithsep = n - i + 1
c
c zero norm eqn
         call NRMSET(1)
c
c restore filtered norm eqn for this epoch
         jd1 = Fep(Ithsep)
         jd2 = Fep(Ithsep + 1)
         m   = Nepoch - i + 1
         k   = 3
         if(m.eq.1) k = 2
         if(m.eq.Nepoch) k   = 1
         if(Lfilt(m).eq.0) k = -k
         call FILSTR(b, Side, buff, np, m, jd1, jd2, k)
         if(Lfilt(m).ne.0) then
c
c get solution
            call SOLCTL(b, .false.)
c
c store data for solution summary printout
            if(Jct(52).gt.0) write(Ibuf6) Ithsep,
     .         (Side(j), Sigma(j), j = 1, Nparam)
c
c
c do adjustment
            if(m.eq.lepoch) keepit = .true.
            call SOLUSE(b, keepit)
c
c predict
            if(Iterat.eq.Ict(1)) then
               if(Ict(10).ge.-1) then
c
c write out side,jd1,jd2, of direct access if complete saved
c solution will be done in prdict
                  if(Fict(5).eq.1) write(Kfile, rec = m) jd1, jd2,
     .               (Side(j), j = 1, Nparam)
 
                  if(Lprdct(m).ne.0) then
                     jda = jd1
                     jdb = jd2
 
c mode=1  means do a time span prdict
                     mode = 1
                     call PRDICT(jda, jdb, mode)
                  endif
               endif
            endif
         endif
 
         end do
c
c print solution summary
      if(Jct(52).gt.0) then
         rewind Ibuf6
         ncol = 2*Nepoch + 2
         call PRNSUM(Nepoch, b, Nparam, ncol)
      endif
c
c
c mode=2 means do a complete saved solution predict over all
c epochs
      if(Iterat.eq.Ict(1)) then
         mode = 2
         if(Fict(5).eq.1) call PRDICT(0._10, 3E7_10, mode)
      endif
c
c
c transfer sum os w's from previous iterations
c to mfile
      pvoflg = .false.
      if((Iterat.gt.1) .or. (Jct(56).ge.2) .or. (Fict(8).gt.0)
     .   .or. (Filflg(1)) ) call WUPDAT(Iconof, Mfile, buff, npnpx2,
     .   Nepoch, Itrwnd, Side, pvoflg)
c
c transfer ic's and offsets to permanent file
      if(Iconof.gt.0) call FLICPM(Iconof, Mfile, buff, npnpx2,
     .                                Nepoch, Itrwnd)
c
c
c
c rewind for next iteration
      rewind Insne
      rewind Filter
      rewind Outsne
      if(Iconof.gt.0) rewind Iconof
      return
      end
