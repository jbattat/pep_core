      subroutine PCHECK(nstop)
 
      implicit none
c
c m.e.ash  june 1969    subroutine pcheck
c consistency check for embary,moon,erotat,mrotat,planets,n-body
c
c arguments
      integer*4 nstop

c array dimensions
      include 'globdefs.inc'
c common
      include 'bdctrl.inc'
      include 'fcntrl.inc'
      include 'france.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      integer*2 lll(u_nmbod,u_mxpl+4)
      equivalence (Lem(1),lll(1,1))
      include 'namtim.inc'
      integer*4 jd0(u_mxpl+4)
      equivalence (Jdem0,jd0)
 
c local
      real*10 frp,t0,t1,t2
      integer*4 i,ifrp,j,jj,k,klip,kold,ll1,ll2,ncp,nki,nn1,nn2,npl
      character*8 name
 
c
c check embary,moon,erotat,mrotat,planet quantities
      do i = 1,u_mxpl+4
         nki = Ndumki(i)
         if(Nplnt(i-4).eq.0) goto 600
         name = Aplnt(i-4)
         npl  = Nplnt(i-4)
         ncp  = Npcent(i-4)
         if(npl.eq.3 .and. Ict(40).gt.0) ncp = -1
c
c see if target planets have appeared in the input
         if(npl.gt.0) then
            do j = 1,i_mxtrg
               if(Kkk(j,i).gt.0 .and. Kkk(j,i).le.30) then
                  if(Kkk(j,i).eq.npl) then
                     write(Iout,10) name,j,Kkk(j,i),'IS ITSELF'
                     if(Mout.gt.0) write(Mout,10) name,j,Kkk(j,i),
     .                'IS ITSELF'
                     nstop=nstop+1
                  else if(Kkk(Kkk(j,i)+30,i).lt.0) then
                     write(Iout,10) name,j,Kkk(j,i),'IS IGNORED'
                     if(Mout.gt.0) write(Mout,10) name,j,Kkk(j,i),
     .                'IS IGNORED'
                     nstop=nstop+1
                  endif
               endif
               if(Kkk(j,i).gt.0 .and. Kkk(j,i).ne.3 .and.
     .         Kkk(j,i).ne.10) then
                  do jj = 1,u_mxpl
                     if(Nplnt(jj).eq.Kkk(j,i)) goto 20
                  end do
                  write(Iout,10) name,j,Kkk(j,i),'DOES NOT EXIST'
                  if(Mout.gt.0) write(Mout,10) name,j,Kkk(j,i),
     .             'DOES NOT EXIST'
   10             format(' ERROR FOR ',a8,
     .             ' BECAUSE TARGET PLANET K(',i2,')=',i3,1x,a)
                  nstop = nstop + 1
               endif
   20       end do
         endif
c
c check consistency of beginning,initial,ending times
         if(jd0(i).eq.0) goto 250
 
c check if there is a possibility of integration with jd0.le.0
         if(jd0(i).le.-3000000) goto 250
 
c check integration epoch not midnight
         frp = 0._10
         if(Intp1(i).lt.0) then
         else if(Intp1(i).eq.0) then
            goto 100
         else
            ifrp = Intp2(i)
            if(ifrp.lt.0 .and. ifrp.ge.-30) then
               frp = Intp1(i)
               frp = frp*2._10**ifrp
               if(frp.lt.1._10) goto 100
            endif
         endif
         write(Iout,50) name,Intp1(i),Intp2(i)
         if(Mout.gt.0) write(Mout,50) name,Intp1(i),Intp2(i)
   50    format(' ERROR FOR ',a8,' BECAUSE INT1=',i11,',  INT2=',i5)
         nstop = nstop + 1
         goto 300
  100    if(jd0(i).lt.0) frp = Dumeps(6,i)
         t0 = iabs(jd0(i)) + frp
         t1 = Jd1(i)
         t1 = t1 + Dumeps(1,i)
         t2 = Jd2(i)
         t2 = t2 + Dumeps(2,i)
         if(jd0(i).eq.-1) t0 = t1
         if(jd0(i).lt.-1) t1 = t0
         call TCHECK(name,t1,t0,t2,Int(i),nstop)
         goto 300
 
c error message for jd1,jd2 inconsistency
  150    write(Iout,200) name
         if(Mout.gt.0) write(Mout,200) name
  200    format(' ERROR FOR ',a8,
     .' BECAUSE VALUES OF JD1,JD2 NOT CONSISTENT FOR INTEGRATION WHICH W
     .ILL OCCUR ON SECOND ITERATION')
         nstop = nstop + 1
         goto 300
  250    if(Ict(1).gt.1) then
            do k = 1,Nbody
               if(npl.eq.Nplbdy(k)) goto 300
            end do
            if(Kkk(88,i).le.-8) goto 300
            do j = 1,6
               if(lll(j,i).gt.0) then
                  if(Jd1(i).le.0 .or. Jd2(i).le.0) goto 150
                  if(Jd2(i).ne.Jd1(i)) goto 300
                  goto 150
               endif
            end do
         endif
c
c check consistency of controls for integration of partial derivatives
  300    kold = 0
         j    = 7
  350    j    = j + 1
         if(Kdumi(j,i).lt.0) then
            if(Kdumi(j,i).lt.kold) goto 450
         else if(Kdumi(j,i).eq.0) then
 
c check for zeros
            do k = j,nki
               if(Kdumi(k,i).ne.0) then
                  write(Iout,360) name,k,Kdumi(k,i),j
                  if(Mout.gt.0) write(Mout,360) name,k,Kdumi(k,i),j
  360             format(' ERROR FOR ',a8,' BECAUSE KI(',i2,')=',
     .                   i3,' WHEREAS KI(',i2,')=0')
                  nstop = nstop + 1
               endif
            end do
            goto 550
         else
 
c see if effect is included for which there is partial

            if(Kdumi(j,i).le.100) then
               jj = Kdumi(j,i)
               if(jj.gt.50) then
                  jj = jj + 50
               else
                  jj = jj + 30
               endif
               if(Kkk(jj,i).le.0) then
                  write(Iout,370) name,j,Kdumi(j,i),jj,Kkk(jj,i)
                  if(Mout.gt.0) write(Mout,370) name,j,Kdumi(j,i),jj,
     .             Kkk(jj,i)
  370             format(' ERROR FOR ',a8,' BECAUSE KI(',i2,')=',i3,
     .                   ' NOT CONSISTENT WITH K(',i3,')=',i3)
                  nstop = nstop + 1
               endif
            endif
 
c check for order
            if(Kdumi(j,i).gt.kold) goto 450
         endif
         write(Iout,400) name,j,Kdumi(j,i)
         if(Mout.gt.0) write(Mout,400) name,j,Kdumi(j,i)
  400    format(' ERROR FOR ',a8,' BECAUSE KI(',i2,')=',i3,
     .          ' IS NOT IN ORDER')
         nstop = nstop + 1
 
c define new old value
  450    kold = Kdumi(j,i)
         klip = kold/100
         if(klip.gt.0) then
            klip = kold - klip*100
 
c increment counter for zonal or tesseral harmonics partials
            if(klip.ge.31) then
               if(klip.eq.31) then
                  nn1 = Kdumi(j + 1,i) - 1
                  nn2 = Kdumi(j + 2,i) - 1
                  ll1 = j
                  j   = j + 2
                  ll2 = j
                  if((nn1.gt.0) .and. (nn2.ge.nn1)) goto 500
               else
                  nn1 = (Kdumi(j+1,i)*(Kdumi(j+1,i)-1))
     .                  /2 + Kdumi(j + 2,i) - 1
                  nn2 = (Kdumi(j+3,i)*(Kdumi(j+3,i)-1))
     .                  /2 + Kdumi(j + 4,i) - 1
                  ll1 = j
                  j   = j + 4
                  ll2 = j
                  if((nn1.gt.0) .and. (nn2.ge.nn1)) goto 500
               endif
               write(Iout,460) name,ll1,ll2
               if(Mout.gt.0) write(Mout,460) name,ll1,ll2
  460          format(' ERROR FOR ',a8,' BECAUSE KI(',i2,
     .                ') TO KI(',i2,') NOT CONSISTENT')
               nstop = nstop + 1
            endif
         endif
  500    if(j.lt.nki) goto 350
         if(j.ne.nki) then
            write(Iout,520) name,j
            if(Mout.gt.0) write(Mout,520) name,j
  520       format(' ERROR FOR ',a8,' BECAUSE KI VECTOR COUNTER IS',i3)
            nstop = nstop + 1
         endif
c
c check consistency for variable output nordsieck integration
  550    if(Kkk(88,i).eq.0) then
            if((Kdumi(1,i).lt.0) .or. (Kkk(100,i).ge.0)) then
               write(Iout,560) name,Kdumi(1,i),Kkk(100,i),Kkk(88,i)
               if(Mout.gt.0) write(Mout,560) name,Kdumi(1,i),
     .          Kkk(100,i),Kkk(88,i)
  560          format(' ERROR FOR ',a8,' BECAUSE KI(1)=',i3,
     .                ' AND/OR K(100)=',i3,
     .                ' NOT CONSISTENT WITH K(88)=',i3)
               nstop = nstop + 1
            endif
         endif
 
      end do
c
c check n-body beginning,initial,ending times
  600 if(Nbody.gt.0) then
         name = ' N-BODY '
         if(Jdbdy0.lt.0) then
            if(Jdbdy0.le.-3000000) goto 700
         else if(Jdbdy0.eq.0) then
            goto 700
         endif
         t0 = iabs(Jdbdy0)
         t1 = Jdbdy1
         t2 = Jdbdy2
         if(Jdbdy0.eq.-1) t0 = t1
         if(Jdbdy0.lt.-1) t1 = t0
         call TCHECK(name,t1,t0,t2,Intbdy,nstop)
      endif
      return
  700 if(Ict(1).le.1) then
      endif
c check consistency of jdbdy1,jdbdy2 for integration on second
c iteration could be inserted here
c
      return
      end
