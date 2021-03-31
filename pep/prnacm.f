      subroutine PRNACM(nstop)
 
      implicit none

c subroutine prnacm - j.f.chandler - 1979 april
c print out a priori matrix, if any
c
c parameter
      integer*4 nstop
c nstop - cumulative count of input errors, possibly incremented here

c common
      include 'aprtbf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
 
c shared external work space
      common/WRKCOM/ Buff(1000),Aptitl,Bufl,Iptr(1000)
      real*10 Buff
      character*80 Aptitl
      character*128 Bufl
      integer*2 Iptr
c local variables
      real*10 sumsq
      integer*4 i,j,n,nsave
      character*8 phdr(8)/'      PA','RAMETER ','NAME    ',' NOMINAL',
     .          ' (SCALED', ')', ' SCALE F', 'ACTOR'/,
     .          rhdr(3)/'ROW   CO', 'EFFICIEN', 'TS'/
 
      if(Lout.le.0) return
      if(Ict(44).eq.0 .and. Ibuf1.le.0) return
      rewind Intern
c
c write out saved summary of adjusted parameters
      if(Line.gt.50) call NEWPG
      write(Iout,100) phdr
  100 format('-SUMMARY OF ADJUSTABLE PARAMETERS'/'0', 10A8)
      Line = Line + 5
      call PAGSET(phdr,-16)
      do while( .true. )
         read(Intern,150,end=200) Bufl
  150    format(a128)
         call PAGCHK(60,1,1)
         write(Iout,150) Bufl
         end do
 
c finished copying parameter summary
  200 rewind Intern
 
      if(Ibuf2.le.0) return
      if(Ict(44).eq.0) goto 400
 
c read back and print a priori saved normal equations
      read(Ibuf2,end=400,err=400) Aptitl
      if(Line.gt.52) call NEWPG
      write(Iout,250) Aptitl,Ibuf2
  250 format('-TITLE= ', a80, 5x, 'ON A PRIORI DATA SET', i3)
c
c read error statistics
      read(Ibuf2) nsave,sumsq
      write(Iout,260) sumsq
  260 format(' EFFECTIVE A PRIORI SUM-SQUARED O-C/ERROR',1pd13.5)
      Line = Line + 4
c
c read nsave rows plus pointer vector
      read(Ibuf2) nsave,(Iptr(i),i=1,nsave)
      n = (nsave - 1)/30 + 3
      call PAGCHK(60,n,0)
      write(Iout,300) nsave,(Iptr(i),i=1,nsave)
  300 format('0NSAVE    POINTERS'/ i6, (t9,30I4))
c
c read right hand side contribution
      read(Ibuf2) nsave,(Buff(i),i=1,nsave)
      n = (nsave - 1)/10 + 3
      call PAGCHK(60,n,0)
      write(Iout,350) (Buff(i),i=1,nsave)
  350 format('0RHS'/ (7x,1p,10D12.5))
 
      n = (nsave - 1)/10 + 1
 
c set flag to require header for 1st line
      call PAGSET(rhdr,-5)
      if(Line + n.gt.56) call NEWPG
      call PAGHED(0)
      do while( .true. )
c
c read rows of compressed coeff array
         read(Ibuf2) i,(Buff(j),j=1,nsave)
         call PAGCHK(60,n,1)
         write(Iout,360) i,(Buff(j),j=1,nsave)
  360    format(i4,(t8,1p,10D12.5))
         if(i.ge.nsave) goto 365
      end do
c
c read estimate vector
  365 read(Ibuf2) nsave,(Buff(i),i=1,max0(1,nsave))
      if(nsave.eq.-1) then
         call PAGCHK(60,4,0)
         write(Iout,370)
  370    format('0', 76('*')/
     .          ' * WARNING, A PRIORI NOMINALS CANNOT BE ',
     .          'COMPUTED WITH OFFSET A PRIORI INPUT *'/
     .          1x,76('*'))
      else
         n = (nsave - 1)/10 + 3
         call PAGCHK(60,n,0)
         write(Iout,380) (Buff(i),i=1,nsave)
  380    format('0ESTIMATES'/ (7x,1p,10D12.5))
      endif
      goto 600
 
c probably no non-zero sigmas, bad ibuf2
  400 call PAGCHK(60,2,0)
      write(Iout,500) Ibuf2
  500 format('0 NO SAVED A PRIORI MATRIX ON DATA SET', i3)
  600 rewind Ibuf2
      Itrwnd(Ibuf2) = 0
      return
      end
