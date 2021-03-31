      subroutine PLPNDX(nctl,ictrl,nlist,klist,nplp,nq,kq,kr,kp,
     .                  numki,ki)
 
      implicit none
 
c     routine plpndx - j.f.chandler - 1982 feb
c     cross-reference partials requested for integration with list of
c     partials available on central-body or target-body tape
c     quantities of interest ('possible partials' + motion) are (in
c     the order of the output arrays):
c             1:   motion
c             2-7: partials w.r.t. initial conditions
c             8... any other partials
c
c arguments
      integer*4 nctl,ictrl(nctl),nlist,klist(nlist),nq
      integer*2 nplp,kq(8),kr(8),kp(nctl),numki,ki(99)
c  nctl = number of partials to be integrated
c  ictrl= list of descriptors for integrated partials
c  nlist= number of possible partials besides central/target-body i.c.'s
c  klist= list of descriptors to be compared with 'ictrl'
c  nplp=  planet number of central or target body
c  nq =   (on return) number of partials (+ motion) to be interpolated
c  kq =   array of 'nq' ptrs to list of possible partials (maps the
c         interpolation y-vectors to the output arrays)
c  kr =   array of 'nq' ptrs into partials on integration tape
c  kp =   array of 'nctl' ptrs from integrated partials to output arrays
c  numki= number of integration partials controls on integration tape
c  ki =   array of integration controls
c
c array dimensions
      include 'globdefs.inc'
c        common
      include 'inodta.inc'
c
c local
      character*4 prmnm(2)/'I.C.', 'PARM'/
      character*14 source(2)/'  CALCULATED  ','READ FROM TAPE'/
      integer*2 ktmp(i_mxplp)
      integer i,ierr,ip,isrc,it,j,jcp,jp,k,king,kong,nic,nn,nlp
c
c form total number of possible partials (+ motion)
      nn = 7 + nlist
      if(nn.gt.i_mxplp)
     . call SUICID('TOO MANY INDIRECT PARTIALS, STOP IN PLPNDX  ',11)
 
c initialize ptr arrays
      ierr   = 0
      nq     = 1
      kq(1)  = 1
      ktmp(1)= 1
      do i = 2, nn
         kq(i)  = 0
         kr(i)  = 0
         ktmp(i)= 0
      end do
      if(nctl.gt.0) then
 
c search for i.c. partials in integration controls
         do i = 1, nctl
            king = (ictrl(i)-1)/100
 
c check for parameter of central or target body
            if(king.eq.nplp) then
               kong = ictrl(i) - king*100
               if(kong.gt.6) goto 50
 
c found an initial condition partial
               nq     = nq + 1
               kq(nq) = kong + 1
            endif
         end do
 
c search for any 'extra' possible partials
   50    do k = 1, nlist
            do i = 1, nctl
               if(ictrl(i).eq.klist(k)) then

c found one, index it
                  nq     = nq + 1
                  kq(nq) = k + 7
                  goto 60
               endif
            end do
   60    end do
      endif
c
c check for partials on tape
      if(nplp.le.0) return
      if(numki.le.0) goto 150
      if(nq.le.1) goto 100
      ip = 1
      if(ki(1).ne.0) then
         do i = 2, 7
            if(ki(i).ge.0) then
               ip     = ip + 1
               ktmp(i)= ip
            endif
         end do
      endif

c check for 'extra' partials on tape
      if(nlist.gt.0) then
         i=8
         do while (i.le.numki)
            ip = ip + 1
            do it=1,nlist
               if(ki(i).eq.klist(it)) ktmp(it+7) = ip
            end do
            king=(ki(i)-1)/100
            if(king.gt.0) then
               kong=ki(i)-100*king
               if(kong.gt.30) then
                  if(kong.ge.40) then
                     ip=ip+(ki(i+3)*(ki(i+3)+1))/2+ki(i+4)
     .                -(ki(i+1)*(ki(i+1)+1))/2-ki(i+2)
                     i=i+4
                  else
                     ip=ip+ki(i+2)-ki(i+1)
                     i=i+2
                  endif
               endif
            endif
            i=i+1
         end do
      endif
 
c combine into output reference array
  100 do i = 1, nq
         ip    = kq(i)
         kr(i) = ktmp(ip)
         if(kr(i).le.0) then
 
c partial not available on tape, remove corresponding ptr to output
            do j=1,nctl
               if(kp(j).eq.ip) kp(j)=0
            end do
            if(ip.gt.7) then
 
c body parameter partial missing
               jp  = klist(ip - 7)
               jcp = 2
            else
 
c missing i.c. partial
               jp  = ip - 1
               jcp = 1
            endif
            write(Iout,120) prmnm(jcp),jp,nplp
  120       format('0PARTIAL W.R.T. ', a4, i4,
     .             ' MISSING FROM TAPE FOR BODY',i2)
            ierr = 1
         endif
      end do
c print index
  150 isrc=2
      if(numki.le.0) isrc=1
      call PAGCHK(60,4+(nq-1)/15*(isrc+1),0)
      write(Iout,200) nq,source(isrc),nplp,(kq(i),i = 1,nq)
  200 format('0', i5,' QUANTITIES TO BE ',a14,' FOR BODY',
     . i2,':',(t51,15i4))
      nic=nq-nlist
      nlp=min0(nlist,nq-1)
      write(Iout,250) 0,(nplp*100+kq(i)-1,i=2,nic),(klist(i),i=1,nlp)
  250 format(30x,'CODES FOR THE ABOVE:', (t51,15i4))
      if(isrc.eq.2) write(Iout,300) (kr(i),i = 1,nq)
  300 format(24x,'RELATIVE POSITION ON TAPE:', (t51,15i4))
      if(ierr.gt.0) then
         call PAGCHK(60,2+(nn-1)/20+(nctl-1)/20,0)
         write(Iout,350) (ktmp(i),i=1,nn)
  350    format(' INDEX OF STORAGE ARRAY:',(1x,20i4))
         write(Iout,370) (ictrl(i),i=1,nctl)
  370    format(' INTEGRATION CONTROLS:'(1x,20i4))
         if(nplp.ne.10)
     .    call SUICID('MISSING PARTIALS, STOP IN PLPNDX',8)
      endif
      return
      end

