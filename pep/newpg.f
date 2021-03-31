      subroutine NEWPG

      implicit none


c*** start of declarations inserted by spag
      integer   i, ipage, itf, lfmth, linc, mpage, nln, nst
c*** end of declarations inserted by spag


c  j.f.chandler - 1976 sep
c           extra headers added for integrations, 1978 sep
c           merged with nrmpg, 1980 oct
c
c           write new page top when necessary
c
c        parameters
      character*8 pname
      character*4 lhed(1),fmth(34)
      integer*4 ioutt, itfi, jd1, jd2, jdz, linci, lmax, n, nn
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'

      character*4 hed1(10)/10*'****'/, blank/'    '/, to/' TO '/
      character*8 hedn
      equivalence (hed1, hedn)
      character*4 hed2(33),hedc(43)
      equivalence (hedc, hed1), (hedc(11), hed2)
      character*4 iter(3)/' (IT','ERAT','...)'/,
     .            soln(3)/' (SO','LUTN','...)'/
      integer*4 izero/0/, ineg1/-1/
      character*4 czero, cneg1
      equivalence (izero,czero),(ineg1,cneg1)

      character*4 fmtloc(34)
      integer*4 iout2/0/, nextra
c
c unconditional new page, no secondary header.
c
      itf  = 0
      linc = 0
      goto 100

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry PAGCHK(lmax, linci, itfi)
c           test current line counter and start new page if needed.
c           increment counter before returning.
c           optionally print secondary header if new page.
c  lmax - maximum lines allowed on this page
c  linc - number of lines about to be printed
c  itf  - non-zero to select optional secondary header
c
c           copy input values to local storage
      linc = linci
      itf  = itfi

      if(Line+linc.le.lmax) goto 500
      goto 100

      entry PAGE(linci, itfi)
c unconditional new page
c optionally print secondary header as well
c
c copy input values to local storage
      linc = linci
      itf  = itfi

  100 write(Iout, 200) hed1, Heding, Date, Npage
  200 format('1', 10A4, 1x, 18A4, 1x, 2A4, ' PAGE', i5)
      if(iout2.gt.0) write(iout2, 200) hed1, Heding, Date, Npage
      Npage = Npage + 1
      Line  = 1
      if(itf.eq.0) goto 500
      goto 300

      entry PAGHED(linci)
c enter here for just the secondary header.  increment line
c
c copy input value to local storage
      linc = linci

      itf = 1
      if(Line+max0(2,linc).gt.58) goto 100
  300 write(Iout, 400) hed2
  400 format('0', 33A4)
      if(iout2.gt.0) write(iout2, 400) hed2
      Line = Line + 2

  500 Line = Line + linc
c
c Normalize page and line numbers
      nextra=(Line-1)/60
      if(nextra.gt.0) then
         Npage=Npage+nextra
         Line=Line-60*nextra
      endif

      return

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry PAGSET(lhed, n)
c           save header string for subsequent calls
c           set up extra output unit (or lack thereof)
c  lhed - text string for header, or lhed(1)=-1 for setting up iout2
c  n    - iabs(n)=length of lhed in words, + if top line, - if 2nd line
c         if lhed(1)=-1, then n is second output unit
c
c           clear second output unit
      iout2 = 0
      if(lhed(1).eq.cneg1) then

c set extra output unit from 'n'
         iout2 = n
         return
      endif
      if(n.le.0) then

c copy into hed2 (extra header line)
         nn  = -n
         nst = 10
         nln = 33
      else

c copy into hed1 (top line left-hand text)
         nn  = n
         nst = 0
         nln = 10
      endif

c limit copy to available array size
      if(nn.gt.nln) nn = nln
      do i = 1, nn
         hedc(i + nst) = lhed(i)
         end do
      if(n.le.0 .or. nn+3.gt.nln) goto 700

c copy iteration number into title
  600 if(Nsoltn.le.0) then
         call EBCDI(Iterat, iter(3), 3)
         do i = 1, 3
            hedc(nn + i) = iter(i)
            end do
      else
         call EBCDI(Nsoltn, soln(3), 3)
         do i = 1, 3
            hedc(nn + i) = soln(i)
            end do
      endif
      nn = nn + 3
  700 if(nn.ge.nln) return
      nn = nn + 1
      do i = nn, nln
         hedc(i + nst) = blank
         end do
      return

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry NEWPGS(pname, jd1, jd2, fmth)
c set up page top
c pname  - name of integrated object (r*8)
c jd1,2  - range of integration
c fmth   - array containing format(if any) for optional 2nd line
c          (MUST be either a binary 0 or a 136-character string)
      hedn = pname
      call EBCDI(jd1, hed1(3), 8)
      hed1(5) = to
      call EBCDI(jd2, hed1(6), 8)
      nn  = 7
      nln = 10
      nst = 0
c save either 4 bytes (if =0) or 136 bytes (if format string) of fmth
      lfmth = 4
      if(fmth(1).ne.czero) lfmth = 136
      call MVC(fmth,1,lfmth, fmtloc,1)
      goto 600

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry NEWPGT(ioutt, mpage, jdz)
c print new page top + maybe second line
c ioutt  - output print unit number
c mpage  - next available page number on ioutt
c jdz    - optional extra line control, if .gt.0 print extra line
c          (format saved on latest call to NEWPGS)
      ipage = mpage

c caution: argument copying conventions require  npage assignment
      if(ioutt.eq.Iout) ipage = Npage
      mpage = ipage + 1
      if(ioutt.eq.Iout) Npage = ipage + 1
      write(ioutt, 200) hed1, Heding, Date, ipage
      if(jdz.gt.0.and.fmtloc(1).ne.czero) write(ioutt,fmtloc)
      return
      end
