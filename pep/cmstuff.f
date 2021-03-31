      subroutine cmstuff
      implicit none

c Fortran replacement for PEP assembly-language routines for reading
c disk files into the input stream.

c        Based on PDSTUFF routines - J.F. Chandler - 1978 April 
c        Modified for CMS      1980 December
c        Converted to Fortran  1992 July

c        Calling sequences:

c call CMSOPN(fileid,nid,icode)
c    Open file with name 'fileid' (which must be a complete filespec);
c    name length is 'nid' bytes.
c    Returns indication of success in 'icode'
c call CMSGET(array,recno,icode)
c    Read a record from last-opened 'fileid' into 'array', success
c    returned in 'icode'.  If the latest CMSOPN failed, the same 
c    error code is returned again in 'icode'.
c    Note: the assumed record length is 80 bytes.  Longer records 
c    are truncated at 80, and shorter ones are padded to 80 with 
c    blanks.
c    Note: the value of icode passed to CMSGET has meaning if the 
c    data set has an uncorrectable i/o error -- then, if icode=6,
c    the bad record is accepted and passed back to the caller with
c    a warning -- any other value of icode causes the bad record 
c    to be rejected.
c call CMSCLS
c    Close latest file (if any) and release buffers (if any).
c    This is not necessary for program operation  -- the old file
c    is automatically closed whenever CMSOPN is called with a 
c    new one.  CMSCLS may, however, be useful when no further 
c    reading is to be done.
c call PDSOPN(fileid,member,icode)
c    Open file named 'fileid' and find member 'member'.  If that file
c    is already open, start the scan from the current point.
c    Return indicator of success in 'icode'.
c    Partitioned files must have separators between members in the
c    format: "./ ADD NAME=member".  The end should be marked by a
c    record containing "./ ENDUP".
c call PDSGET(array,icode)
c    Read a record from last-opened 'member' into 'array', success
c    returned in 'icode'.  If the latest PDSOPN failed, the same 
c    error code is returned again in 'icode'.
c    Note: the assumed record length is 80 bytes.  Longer records 
c    are truncated at 80, and shorter ones are padded to 80 with 
c    blanks.
c    Note: the value of icode passed to PDSGET has meaning if the 
c    data set has an uncorrectable i/o error -- then, if icode=6,
c    the bad record is accepted and passed back to the caller with
c    a warning -- any other value of icode causes the bad record 
c    to be rejected.
c call PDSCLS 
c    Close latest file (if any) and release buffers (if any).
c    This is not necessary for program operation  -- the old file
c    is automatically closed whenever PDSOPN is called with a 
c    new one.  PDSCLS may, however, be useful when no further 
c    reading is to be done.

c Composite description of calling arguments

c    fileid - filespec for the desired file (in a form usable
c             by the Fortran OPEN statement).
c    nid    - length in bytes of 'fileid' (only for CMSOPN: length is
c             always 8 for PDSOPN).
c    member - 8-character name of pds member to be read 
c    recno  - record number to read. If 0, then next sequential.
c             This argument is not yet implemented and doesn't exist
c             for PDSGET at all.
c    array  - 80-byte array for returned record from file.
c    icode  - integer*4 variable for return code 
c         -1  reached end of file (no info in array)
c          0  action completed ok 
c          1  file (or member) not found on CMSOPN or PDSOPN
c          2  file could not be opened 
c          3  improper record format 
c          4  CMSOPN or PDSOPN never called (or not since close)
c          5  uncorrectable i/o error (no info in array)
c          6  uncorrectable i/o error (but array has something)
c          7  invalid fileid 

      integer*4 fidlen, nid, i, icode, recno, itype, currec, recsav,
     .     ISCAN, oldtyp
      character*80 fileid
      character*80 array, cntrl
      character*8  member, oldfil

      integer*4 iunit/3/, istat/4/
      logical*4 exist, d/.false./

      entry CMSOPN(fileid,nid,icode)

      fidlen = nid
      itype  = 1

 50   if(istat.le.0) then
         close(iunit)
         istat = 4
      endif
      if(d) write(6,fmt='('' >>Inquire '',a)') fileid(1:fidlen)
      inquire(file=fileid(1:fidlen),exist=exist)
      if(.not.exist) then
c file doesn't exist
         istat = 1
         icode = 1
         return
      endif
      if(d) write(6,fmt='('' >>Open '',a)') fileid(1:fidlen)
      open(iunit,file=fileid(1:fidlen),form='FORMATTED',status='OLD',
     .     err=100)
      istat = 0
      icode = 0
      currec= 0
      oldfil= fileid(1:8)
      oldtyp= itype
      if(itype.eq.2) goto 120
      return

c OPEN failure
 100  istat = 2
      icode = 2
      return
c-----------------------------------------------------------------------
      entry PDSOPN(fileid,member,icode)

      itype = 2
      fidlen= 8
      if((istat.le.0.and.fileid(1:8).ne.oldfil) .or.
     .    istat.gt.0 .or. itype.ne.oldtyp) goto 50
 120  recsav = currec
      if(d) write(6,fmt='('' >> Starting scan at recno'',i6)') currec
c the current "./" card might be the one we need
      if(currec.gt.0) goto 135
 130  read(iunit,200,err=300,end=140) cntrl
      currec = currec+1
c check for PDS format if starting
      if(currec.eq.1 .and. cntrl(1:3).ne.'./ ') then
         istat = -1
         icode = 3
         if(d) write(6,fmt='('' >> Bad start:'',a20)') cntrl(1:20)
         return
      endif
c check for complete wrap
      if(currec.eq.recsav) then
         istat = -1
         icode = 1
         if(d) write(6,fmt='('' >> Wrapped at rec'',i5)') currec
         return
      endif
 135  if(cntrl(1:3).ne.'./ ') goto 130
      if(d) write(6,fmt='('' >> Control='',a40)') cntrl(1:40)
      if(cntrl(4:8).eq.'ENDUP') then
         goto 140
      else if(cntrl(4:12).eq.'ADD NAME=') then
         i=ISCAN(cntrl(13:20), 8, ',')
         if(i.gt.0) cntrl(12+i:18+i) = '       '
         if(d) write(6,fmt='('' >> ?? '',a8,''='',a8)') cntrl(13:20),
     .        member
         if(cntrl(13:20).ne.member) goto 130
c found member
         istat = 0
         icode = 0
         return
      endif

c anomalous control card
      istat = -1
      icode = 5
      return

c reached end during scan for member
 140  rewind iunit
      currec = 0
      if(d) write(6,fmt='('' >> Hit end'')')
      goto 130
c-----------------------------------------------------------------------
      entry CMSGET(array,recno,icode)

      itype = 1
      goto 150

      entry PDSGET(array,icode)
      itype = 2

c check if ok to read a record
 150  if(istat.gt.0) then
         icode = 4
         return
      else if(istat.lt.0) then
         icode = -1
         return
      endif

c make sure same type of request (CMS vs PDS)
      if(itype.ne.oldtyp) then
         istat = -1
         icode = 4
         return
      endif

c grab a record
      read(iunit,200,err=300,end=400) array
 200  format(a80)
      currec = currec+1
      cntrl = array
      icode = 0

c check for logical EOF if PDS
      if(itype.eq.2 .and. array(1:3).eq.'./ ') goto 400
      return

c I/O error (pretend end of file for future calls)
 300  istat = -1
      icode = 5
      return

c EOF
 400  istat = -1
      icode = -1
      return

c-----------------------------------------------------------------------
      entry CMSCLS
      entry PDSCLS

      if(istat.le.0) then
         close(iunit)
         istat = 4
      endif
      return
      end
