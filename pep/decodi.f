      subroutine DECODI(text, pos, end, num1)
 
      implicit none
 
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  J.F.Chandler  --  1978 July 8  (Assembler version)
c     Converted to Fortran 77 1992 Jun 16
c
c  Fortran-callable routine to decode text strings containing integers.
c  The text string may be padded with any number of blanks between
c  the integers, but there may be no imbedded blanks in the integers
c  themselves.   Any character other than a digit, sign or blank will
c  be interpreted as a separator for an omitted number.
c
c  Calling sequence --  call DECODI(text,pos,end,num1)
c        TEXT    Array of text input to DECODI.
c        POS     Starting character position in string for decoding.
c        END     Last character position in string for decoding.
c        NUM1    Integer to be read from string.
c  NOTE: If the text string runs out before NUM1 is filled, NUM1 will
c  be zeroed.  Originally, multiple integers could be specified, but
c  the Fortran version can handle only one.
c
c  On return from DECODI, 'POS' will contain the position of the next
c  character remaining undecoded.  Any non-numeric character other than
c  a sign used to terminate a number will be considered "decoded" when
c  encountered.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c---routines
c   called:      ichar
c
 
      character*(*) text
      integer      pos, end, num1, ICHAR, digit
      logical      negate
 
c set default value for num1 and check for legal value of pos
 
      num1 = 0
      negate = .false.
 
c find something non-blank
      do while (pos .le. end)
        if( text(pos:pos) .ne. ' ' ) goto 100
        pos = pos + 1
        end do
 
      return
 
c note sign, if any
  100 if( text(pos:pos) .eq. '+' ) then
         pos = pos + 1
      else if( text(pos:pos) .eq. '-' ) then
         pos = pos + 1
         negate = .true.
         end if
 
c decode digits
      do while (pos .le. end)
         digit = ICHAR(text(pos:pos)) - ICHAR('0')
         if( digit .lt. 0 .or. digit .gt. 9) then
           if( text(pos:pos) .ne. '-' .and. text(pos:pos) .ne. '+')
     .         pos = pos + 1
           goto 200
           endif
         num1 = 10*num1 + digit
         pos = pos + 1
         enddo
 
c apply sign, if any
  200 if( negate ) num1 = -num1
      return
      end
