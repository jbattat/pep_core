      subroutine NAMRED(in0,nm,nappar,n)
 
      implicit none

c     d. white  april 1973  subroutine namred
c         reads one or two part parameter names in free field format
c         using blanks as delimiters, one name per record, until eof or
c         nappar have been read.
c        r.reasenberg jan 77   a name must not start with '
c        an internal ' is ok if there are no internal blanks
c        for names with internal blanks, use ' as the delimiters
c
c         parameters
      integer*4 in0,nappar,n
      character*8  nm(2,1)
c in0   - unit number to read for names
c nm    - array of names to fill
c nappar- expected or maximum number of names to read(size of nm)
c n     - returned number of names found
c
c array dimensions
      include 'globdefs.inc'

c common
c note: numpln,nplnt,aplnt should really be passed as arguments
      include 'namtim.inc'
      include 'prmnms.inc'

c external functions
      integer*4 ISCAN,LEG,NSCAN
c
c locals
      character*8 buff(9)
      character*8 blanks(2)/2*'        '/, mssbdy(2)/'MASS    ','BODY'/
      character*1 blank /' '/
      character*1 apos/''''/,dlm,cbuff(72)
      equivalence (buff,cbuff)
      integer   i,j,k,np
c
c zero n,the number of pairs found
      n = 0
  100 do while( .true. )
c
c read a card
         read(in0,150,end=400) buff
  150    format(9A8)
 
c find 1st non-blank character
         i = NSCAN(buff,72,blank) + 1
         if(i.gt.0) then
            n = n + 1
            if(n.gt.nappar) goto 400
c
c have another param, so first blank names
            call MVC(blanks,1,16,nm(1,n),1)
            dlm = blank
            if(LEG(1,i,buff,1,apos).eq.0) then
               dlm = apos
               i   = i + 1
            endif
            j = ISCAN(cbuff(i),73 - i,dlm) + i
            if(j.le.0) j = 72
            k = min0(j - i,8)
c
c move all non-blank char of first name
            call MVC(buff,i,k,nm(1,n),1)
            if(j.lt.72) then
               i = NSCAN(cbuff(j+1),72 - j,blank) + j + 1
               if(i.gt.0) then
                  dlm = blank
                  if(LEG(1,i,buff,1,apos).eq.0) then
                     dlm = apos
                     i   = i + 1
                  endif
                  j = ISCAN(cbuff(i),73 - i,dlm) + i
                  if(j.le.0) j = 72
                  k = min0(j - i,8)
 
c move all non-blank char of second name
                  call MVC(buff,i,k,nm(2,n),1)
                  goto 160
               endif
            endif
c
c check for special single names
            if(LEG(6,1,nm(1,n),1,Qprmtr).eq.0) then
               i = 7
               call DECODI(nm(1,n),i,8,np)
               if(i.ne.9) goto 200
               if(np.le.30) goto 300
               do i = 1, Prmsmx
                  if(np.eq.Iprms(i)) then
 
c found in list
                     nm(1,n) = Prms(i)
                     goto 200
                  endif
               end do
               goto 200
            endif
c
c check for masses
  160       if(LEG(12,1,nm(1,n),1,mssbdy).eq.0) then
               i = 5
               call DECODI(nm(2,n),i,8,np)
               if(i.eq.8) then
                  if(LEG(2,7,nm(2,n),1,blanks).eq.0) goto 300
               endif
            endif
         endif
  200 end do
  300 if(np.gt.10 .and. np.le.30) then
         if(Numpln.gt.0) then
            do i = 1, Numpln
               if(np.eq.Nplnt(i)) then
 
c found body
                  if(Aplnt(i).ne.blanks(1)) then
                     nm(1,n) = mssbdy(1)
                     nm(2,n) = Aplnt(i)
                  endif
                  goto 100
               endif
            end do
         endif
      endif
      goto 100
c
c end of input
  400 rewind in0
      return
      end
