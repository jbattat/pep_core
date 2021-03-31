      subroutine TRGLIC(nkbb,kbb,onlyic)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, kbbll, klantg, klnh, klnhr, lgg, ll, lm, lpp, ltg,
     .          ltg1, ltg2, ltt, mm, mtg1, n1, n2, n9, ngo, ntg1
 
c*** end of declarations inserted by spag
 
 
      integer*2 nkbb, kbb(30)
      logical*4 onlyic
c           subroutine trglic - j.f.chandler - 1980 oct
c           form partial derivative controls for target bodies
c           initial conditions plus (if onlyic is false) parameters
c           and harmonic coefficients, from integrated partials
c           indicated by kbb plus saved obslib partials (if any)
c           (based on  m.e.ash   subr. trgcnt)
c
c array dimensions
      include 'globdefs.inc'
c common
      include 'comdat.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'ltrapx.inc'
      include 'mtrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'plnhar.inc'
      include 'scoef4.inc'
c
c local
      integer*2 mmtz, mmtt
      integer*2 izr2/0/
c           local pointers:
c ll:  next unfinished in kbb
c mm:  next uncopied old target body (1-mumtar)
c lm:  next uncopied in mtbod (+mt.har)
c ltt: last output to ltbod
c lpp: next unchecked in lpl
c ngo: indicates latest item from kbb: 1-nothing, 2-zonal,
c           3-tesseral cosine, 4-tesseral sine
c
c*  start=100
      mm = 1
      ll = 8
  100 do while( .true. )
c
c see if there are target body partials on integration tape
         ltg1 = 0
         if(ll.gt.nkbb) goto 200
         kbbll = kbb(ll)
         if(kbbll.eq.0) goto 200
         ltg1 = kbbll/100
         if(ltg1.gt.0) then
            if(ltg1.ne.Ncp0) then
               if(ltg1.ne.3 .and. ltg1.ne.10) then
c*  start=200
c find harmonics, if any, for body 'ltg1'
                  klnhr = 0
                  if(.not. (Nmphar.le.0 .or. onlyic)) then
                     do i = 1, Nmphar
                        if(Nplhar(i).eq.ltg1) then
                           klnhr = i
                           goto 200
                        endif
                     end do
                  endif
                  goto 200
               endif
            endif
            if(kbbll.gt.(ltg1*100+30)) ll = ll + 2
            if(kbbll.gt.(ltg1*100+40)) ll = ll + 2
         endif
         ll = ll + 1
      end do
c
c*  start=300
c see if there are target body partials on input observation
c library tape
  200 mtg1 = 0
      if(Iabs1.gt.0 .and. mm.le.Mumtar) then
         mtg1 = Mtrg(mm)
         mmtz = Mtzone(mm)
         mmtt = Mttess(mm)
         if(onlyic .and. (mmtz.gt.0 .or. mmtt.gt.0)) goto 900
      endif
      if(ltg1.le.0 .and. mtg1.le.0) return
      Numtar = Numtar + 1
      if(Numtar.gt.i_mxtrg) goto 700
      ntg1 = ltg1
      klnh = klnhr
      if(mtg1.gt.0) ntg1 = mtg1
      if(ltg1.gt.0 .and. ltg1.lt.ntg1) ntg1 = ltg1
 
c find corresponding input planet
      klantg = 0
      do while( .true. )
         klantg = klantg + 1
         if(klantg.gt.Numpln) then
c
c*  start=9000
c found no klan belonging to target body
            write(Iout,220) ntg1
  220       format('0* * * TARGET BODY', i3, ' NOT INPUT, ERROR * * *')
            call SUICID('TARGET BODY NOT INPUT, STOP TRGLIC  ', 9)
            goto 700
         else if(ntg1.eq.Nplnt(klantg)) then
            Klant(Numtar) = klantg
            Ntrg(Numtar)  = ntg1
            if(ltg1.eq.ntg1) then
c
c*  start=1000
c target body is from integration (and maybe obslib)
               lm  = 7
               ltt = 6
               lpp = 7
               if(mtg1.ne.ntg1) then
c
c*  start=2000
c target body is from integration tape only
c indicate no copying from mtbod, etc.
                  lm   = 31
                  mmtz = 0
                  mmtt = 0
               else
c target body is from both integration and obslib tapes
c copy i.c. controls from old obslib tape
                  do i = 1, 6
                     Ltbod(i,Numtar) = Mtbod(i,mm)
                  end do
               endif
               goto 300
            else
c
c target body is from obslib tape only
               klnh = 0
               do i = 1, 30
                  Ltbod(i,Numtar) = Mtbod(i,mm)
               end do
 
c copy harmonic controls, if any
               lm  = 31
               ngo = 1
               goto 600
            endif
         endif
      end do
  300 do while( .true. )
 
c start of ll loop
         ngo = 1
         if(ll.gt.nkbb) goto 500
         ltg2 = kbb(ll)/100
         if(ltg2.ne.ltg1) goto 500
         ltg = kbb(ll) - ltg2*100
         if(ltg.gt.6) then
 
c parameters
            if(onlyic) goto 900
            if(ltg.ge.31) then
c*  start=1200
c harmonics
               ngo = (ltg - 1)/10 - 1
               if(ltg.gt.31) then
c*  start=1300
c tesseral harmonics
                  n1 = (kbb(ll+1)*(kbb(ll+1)-1))/2 + kbb(ll + 2) - 1
                  n2 = (kbb(ll+3)*(kbb(ll+3)-1))/2 + kbb(ll + 4) - 1
                  ll = ll + 4
               endif
               goto 500
            else
               lgg = ltg - 6
               do while( lm.le.30 )
                  if(Mtbod(lm,mm).le.0) then
                     lm = 31
                     goto 310
                  else if(Mtbod(lm,mm).lt.lgg) then
                     ltt = ltt + 1
                     Ltbod(ltt,Numtar) = Mtbod(lm,mm)
                     lm = lm + 1
                  else if(Mtbod(lm,mm).eq.lgg) then
                     lm = lm + 1
                     goto 320
                  else
                     goto 310
                  endif
               end do
  310          do while( lpp.le.30 )
                  if(Lpl(lpp,klantg).le.0) then
                     goto 400
                  else if(Lpl(lpp,klantg).lt.lgg) then
                     lpp = lpp + 1
                  else if(Lpl(lpp,klantg).eq.lgg) then
                     goto 320
                  else
                     goto 400
                  endif
               end do
 
c end of ll loop
               goto 400
            endif
  320       ltt = ltt + 1
            Ltbod(ltt,Numtar) = lgg
         else
 
c initial conditions
            if(Lpl(ltg,klantg).gt.0) Ltbod(ltg,Numtar) = 1
         endif
400      ll = ll + 1
      end do
c
c*  start=1700
c no more integration, take rest from obslib tape
c or, just starting new type of partials
  500 if(.not. (onlyic)) then
         n9 = ltt + 1
         if(n9.le.30) then
 
c copy remaining old obslib partial controls
            do i = n9, 30
               if(lm.gt.30) goto 600
               Ltbod(i,Numtar) = Mtbod(lm,mm)
               lm = lm + 1
            end do
         endif
      endif
  600 if(ngo.eq.2) then
 
c zonal harmonics
         n1 = kbb(ll + 1) - 1
         n2 = kbb(ll + 2) - 1
         ll = ll + 2
         call LVTHAR(Ltzhar(1,Numtar),Mtzhar(1,mm),Lpzhar,4,klnh,
     .               n1, n2, Ntzone(Numtar),mmtz,4)
         lm = 41
         ll = ll + 1
         goto 300
      else
         if(lm.lt.41) then
            call LVTHAR(Ltzhar(1,Numtar),Mtzhar(1,mm),izr2,1,0,
     .                  0, 0, Ntzone(Numtar),mmtz,4)
            lm = 41
         endif
         if(ngo.eq.3) then
 
c tesseral cosine harmonics
            call LVTHAR(Ltchar(1,Numtar),Mtchar(1,mm),Lpchar,4,klnh,
     .                  n1, n2, Nttess(Numtar),mmtt,5)
            lm = 51
            ll = ll + 1
            goto 300
         else
            if(lm.lt.51) then
               call LVTHAR(Ltchar(1,Numtar),Mtchar(1,mm),izr2,1,0,
     .                     0, 0, Nttess(Numtar),mmtt,5)
               lm = 51
            endif
            if(ngo.eq.4) then
 
c tesseral sine harmonics
               call LVTHAR(Ltshar(1,Numtar),Mtshar(1,mm),Lpshar,4,
     .                     klnh, n1, n2, Nttess(Numtar),mmtt,5)
               lm = 61
               ll = ll + 1
               goto 300
            else
               if(lm.lt.61) then
                  call LVTHAR(Ltshar(1,Numtar),Mtshar(1,mm),izr2,1,
     .                        0, 0, 0, Nttess(Numtar),mmtt,5)
                  lm = 61
               endif
               if(mtg1.eq.ntg1) mm = mm + 1
               goto 100
            endif
         endif
      endif
 
  700 write(Iout,800) i_mxtrg
  800 format('0* * * NO MORE THAN', i3, ' TARGETS ALLOWED * * *')
      call SUICID('TOO MANY TARGET BODIES, STOP IN TRGLIC  ', 10)
 
  900 call SUICID(
     .'ONLY TARGET BODY INITIAL CONDITION PARTIALS ALLOWED, STOP IN TRGL
     .IC ', 17)
c
c*  start=9990
      return
      end
