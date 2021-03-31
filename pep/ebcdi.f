c-----------------------------------------------------------------------
c
c  EBCDI
c
c  j.f.chandler, 1977 may
c
c  modified 11/25/91  mam (sao)
c
c    num       i*4 value to be encoded
c    store     any output location or array
c    iwid      i*4 width of field desired
c
c  output field is padded on left with blanks and
c  right justified.  if the width 'iwid' is too narrow
c  to hold the output, the field is flagged with asterisks as in
c  Fortran formatted writes.
c
c-----------------------------------------------------------------------
      subroutine EBCDI(num, store, iwid)
 
      implicit none
 
      integer      num, iwid
      character*(*) store
      character*10 fmt
 
      if(iwid.lt.10) then
         write(fmt, '( ''(I'',i1,'')'')' ) iwid
      else if(iwid.lt.100) then
         write(fmt, '( ''(I'',i2,'')'')' ) iwid
      else
         write(fmt, '( ''(I'',i3,'')'')' ) iwid
      endif
 
      write(store(1:iwid), fmt) num
 
      return
      end
 
 
 
c-----------------------------------------------------------------------
c
c  EBCDIX
c
c     same parameters as ebcdi except
c          ibeg - start formatting at ibeg'th character of storage
c
c-----------------------------------------------------------------------
      subroutine EBCDIX(num, store, ibeg, iwid)
 
      implicit none
 
      integer      num, ibeg, iwid
      character*(*) store
      character*10 fmt
 
      if(iwid.lt.10) then
         write(fmt, '( ''(I'',i1,'')'')' ) iwid
      else if(iwid.lt.100) then
         write(fmt, '( ''(I'',i2,'')'')' ) iwid
      else
         write(fmt, '( ''(I'',i3,'')'')' ) iwid
      endif
 
      write(store(ibeg:ibeg+iwid-1), fmt) num
 
      return
      end
