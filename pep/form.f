      subroutine FORM(b)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, i1, i2, j, k, LEG, m, n, nd, ne, np, npnp1
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine form
c paul macneil january, 1978 modified for iteration
c
 
c parameter is b matrix - dim set in block data
      real*10   b(1)
c
c common
      include 'fcntrl.inc'
      include 'filptr.inc'
      include 'filstf.inc'
      real*10   wtrans(300, 7)
      equivalence (G(1,1), wtrans(1,1))
      include 'filtda.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'rtside.inc'
 
c temporary work arrays
      common /WRKCOM/ Bthts(300, 7), Fstuff(1000), H(20, 20), S(20, 20),
     .        A6(20, 20), Sh(20, 20)
      real*10   Bthts, H, Sh, has(20, 20)
      real*10 Fstuff, S, A6
      equivalence (Sh(1,1), has(1,1))
c
c local
      character*16 nm(1000)
      equivalence (nm, G(1,1))
      integer*2 iptr(1000)
      equivalence (iptr, Sigma)
      real*10 wr(20), wi(20)
      character*80 title
      logical*4 shoinf
      character*8 label/'I   /  '/
c
c write page heading
      call PAGSET('FORM FILTER ', 3)
      call NEWPG
      i = Fict(1)
      shoinf = mod(i, 2) .eq. 1
c
c read smat into direct access
      npnp1 = Npnp
      call RSMAT(S, S, npnp1)
c
c transfer insne title to filter
      read(Insne) title, np, ne
      write(Filter) title, np, ne, Npnp
      write(Iout, 100) title, np, ne
  100 format('-INPUT SNE:             TITLE= ', A80/' NPARAM=', i4,
     .       '  NEPOCH=', i4)
c
c insne residual summary
      read(Insne)
c
c names
      read(Insne) (nm(i), i = 1, np)
      write(Filter) (nm(i), i = 1, np)
c
c         nominals
c
c         calculate pointers to process noise params (nptr) and
c         deterministic params (iptr).  i indexes all params, j noise,
c         k determs.
      i = 0
      j = 1
      k = 1
      do while(.true.)
         i = i + 1
         if(i.gt.np) then
            nd    = np - Npnp
            Lnptr = Npnp
            if((j-1).ne.Npnp )
     .          call SUICID('PARM NAMES NOT MATCH IN FORM', 7)
c
c write out process noise pointer
            write(Filter) (Nptr(i), i = 1, Npnp)
c
c a priori for forwards form (delta b(1))
            call RINSNE(b, Side, Fstuff, np, 1, 1)
            call EBCDIX(1, label, 3, 2)
            call EBCDIX(1, label, 6, 2)
            if(shoinf) call XPMPO(np, b, 1, 1, label, 8)
c
c loop on number of epochs
            n = ne - 1
            do i = 1, n
c
c read in appropriate s
               read(Lfile, rec = i) ((S(i1,i2),i1=1,Npnp), i2 = 1, Npnp)
c
c         smear covariances
c        pages 3-6 of memo
c        p** -1 is replaced by m** -1 (ie i-(k-1)/(k-1) is replaced by
c        i(k)/(k-1) addsmr used methods outlined on pgs 3-6 to execute
c        eq 5 of memo
               call ADDSMR(b, S, H, Sh, G, np, Npnp, Nptr, iptr, Smat,
     .                     Bthts, nd, A6, wr, wi)
               call EBCDIX(i + 1, label, 3, 2)
               if(shoinf) call XPMPO(np, b, 1, 0, label, 8)
c
c form wtrans
               if(Iconof.gt.0) call FRMWTR(b, wtrans, Npnp, np, Nptr)
c
c         save lhs in direct access for use
c        mfile now has i(k/k-1) for forward only rhs is also written but
c        overwritten in use write so that record calc done correctly
c        (i(k/k)-delta b)
               if(i.lt.n) call WFILDE(b, Side, i, np, Fstuff, 1, 1)
c
c add delta b(i+1)
c to form ik/k as in eq 6 of memo
               call RINSNE(b, Side, Fstuff, np, 1, i + 1)
               call EBCDIX(i + 1, label, 6, 2)
               if(shoinf) call XPMPO(np, b, 1, 0, label, 8)
c
c write b, bthts, and wtrans for this epoch to filter data set
c filter data set now has i(k/k) and i(k/k-1)*s
               call WFILTR(b, Fstuff, Bthts, np, 3, i + 1, Npnp, wtrans,
     .                     has, nd)
               end do
c
c
c
c smoothing backwards
            if(Fict(2).eq.0) return
c
c check insne for backward data
            read(Insne) i
            if(i.ne.0) call SUICID('NO BACKWARDS DATA IN FORM   ',7)
            write(Filter) i
c
c reset normal equations
            Nparam = Nparam*2
            call NRMSET(1)
            Nparam = Nparam/2
c
c a priori for backwards form (delta b(n))
            call RINSNE(b, Side, Fstuff, np, 1, ne)
            call EBCDIX(ne, label, 3, 2)
            call EBCDIX(ne, label, 6, 2)
            if(shoinf) call XPMPO(np, b, 1, 1, label, 8)
c
c loop again
            do i = 1, n
               m = n - i + 1
c
c read in appropriate s
               read(Lfile, rec = m) ((S(i1,i2),i1=1,Npnp), i2 = 1, Npnp)
c
c smear covariances
               call ADDSMR(b, S, H, Sh, G, np, Npnp, Nptr, iptr, Smat,
     .                     Bthts, nd, A6, wr, wi)
               call EBCDIX(m, label, 3, 2)
               if(shoinf) call XPMPO(np, b, 1, 0, label, 8)
c
c form wtrans
               if(Iconof.gt.0) call FRMWTR(b, wtrans, Npnp, np, Nptr)
c
c         make residual gain matrix
c        call to makeg deleted
c
c         add delta b(m)
               call RINSNE(b, Side, Fstuff, np, 1, m)
               call EBCDIX(m, label, 6, 2)
               if(shoinf) call XPMPO(np, b, 1, 0, label, 8)
c
c write b, bthts, and wtrans for this epoch to filter data set
               call WFILTR(b, Fstuff, Bthts, np, 3, i + 1, Npnp, wtrans,
     .                     has, nd)
               end do
            return
         else
            if(j.le.Npnp) then
               if(LEG(16,1,nm(i),1,Pnames(1,j)).eq.0) then
                  Nptr(j) = i
                  j = j + 1
                  go to 200
               endif
            endif
            iptr(k) = i
            k = k + 1
         endif
  200    end do
      end
