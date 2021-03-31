      subroutine CNSTRN(icntrl, nrec, neqn, b, side, iptr)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 bdmst
      integer   i, icntrl, ii, im, imst, is, islv, jboth, jdmst, jdslv,
     .          jmst, jslv, km, ks, neqn, nrec
 
c*** end of declarations inserted by spag
 
 
      real*10 b(1), side(neqn)
      integer*2 ktype, npm, iptr(neqn)
c       subroutine cnstrn - j.f.chandler - 1979 nov 9
c       apply parameter constraints to a set of normal equations
c
c  icntrl - unit number of data set with constraints
c  nrec -   number of records (on return)
c  neqn -   total number of equations
c  b -      triangular matrix form of equation coefficients
c  side -   right-hand side vector
c  iptr -   work space used in applying constraints
c
c       entry cnstri - added 1980 april
c       redistribute solutions, scales, & covariances
c
      real*10 soln(neqn), scale(neqn)
c
c soln -   solution vector
c scale -  vector of scale factors
c
      include 'fcntrl.inc'
 
      nrec = 0
      do while( .true. )
         read(icntrl) ktype, npm, (iptr(i), i = 1, npm)
         Itrwnd(icntrl) = 2
         if( ktype .ne. 3 ) then
            if( ktype .ne. 2 ) call SUICID('BAD CONSTRAINT RECORD   ',
     .          6)
            if( npm .eq. 0 ) return
 
c count usable records
            nrec = nrec + 1
 
c pointers for 'master' parameter
            imst  = iptr(1)
            jmst  = (imst*(imst-1))/2
            jdmst = jmst + imst
c
c
c loop over 'slave' parameters
c note: islv.ne.imst guaranteed by prmtrn
            do i = 2, npm
               islv  = iptr(i)
               jslv  = (islv*(islv-1))/2
               jdslv = jslv + islv
c add together rows and columns
c symmetric matrix, so row = column, one operation
c (only special case: diagonal elts. & their corners)
               jboth = jmst + islv
               if( islv .gt. imst ) jboth = jslv + imst
 
c take care of b(m,m) b(m,s) b(s,m) b(s,s)
               b(jdmst) = b(jdmst) + 2._10*b(jboth) + b(jdslv)
               b(jboth) = 0._10
               b(jdslv) = 0._10
 
c add right-hand sides
               side(imst) = side(imst) + side(islv)
               side(islv) = 0._10
 
c initialize pointers for rows
               km = jmst
               ks = jslv
               im = 1
               is = 1
               do ii = 1, neqn
                  km    = km + im
                  ks    = ks + is
                  b(km) = b(km) + b(ks)
                  b(ks) = 0._10
                  if( ii .ge. imst ) im = ii
                  if( ii .ge. islv ) is = ii
               end do
            end do
         endif
      end do
 
      entry CNSTRI(icntrl, neqn, b, soln, scale, iptr)
 
      if( icntrl .le. 0 .or. Itrwnd(icntrl) .ne. 2 ) return
      read(icntrl) ktype, npm, (iptr(i), i = 1, npm)
      Itrwnd(icntrl) = 2
      if( ktype .ne. 4 ) call SUICID('BAD CONSTRAINT RECORD IN CNSTRI '
     .                               , 8)
      if( npm .eq. 0 ) return
      do i = 1, npm
         if( soln(i) .eq. 0._10 .and. scale(i) .eq. 1._10 ) then
 
c pointers for 'master' and 'slave' parameter
            imst = iptr(i)
            islv = i
            jmst = (imst*(imst-1))/2
            jslv = (islv*(islv-1))/2
 
c distribute solution and scale factor
            soln(islv)  = soln(imst)
            scale(islv) = scale(imst)
 
c save variance of master
            bdmst = b(jmst + imst)
 
c initialize pointers for rows
            km = jmst
            ks = jslv
            im = 1
            is = 1
c
c copy master row & column to slave
            do ii = 1, neqn
               km    = km + im
               ks    = ks + is
               b(ks) = b(km)
               if( ii .ge. imst ) im = ii
               if( ii .ge. islv ) is = ii
            end do
 
c insert variance for slave
            b(jslv + islv) = bdmst
         endif
      end do
      return
      end
