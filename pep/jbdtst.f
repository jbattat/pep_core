      integer function JBDTST(nql)
 
      implicit none

c m.e.ash   sept 1968    integer function jbdtst

c           jbdtst=0 if planet nql is not on n-body tape
c           jbdtst.gt.0 if planet nql is on n-body tape, in which case
c                    we have nql=npl(jbdtst)
c           jbdtst.lt.0 if planet nql is on s-body tape, in which case
c                    nql=np2(-jbdtst)

c arguments
      integer*2 nql
c ict(39)=-2 take nothing from n-body data set (n-body tape exists
c            only for center of mass of solar system in compar link)
c ict(39)=-1 take only moon from n-body tape if there is one
c ict(39)= 0 take everything from n-body tape that is there (compar l.)
c ict(39)= 1 use individual body tape whenever it is available
c
c        jct(54) = -2  take nothing from s-body data set (compar link)
c        jct(54) =  0  take anything from s-body data set
c        jct(54) =  1  take from s-body data set unless the
c                      individual body tape is available

c array dimensions
      include 'globdefs.inc'

c commons
      include 'bdctrl.inc'
      include 'bddta.inc'
      include 'b2dta.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      integer*4 i2bod
      equivalence (i2bod,Jpert)
      include 'namtim.inc'
      include 'plndta.inc'
c local variables
      integer   i, j, k, num4
 
      JBDTST = 0
      if(nql.gt.10) then
 
c check s-body tape
         if(i2bod.gt.0 .and. Jct(54).gt.-2) then
            do i = 1, Nast
               if(Np2(i).eq.nql) then
                  JBDTST = -i
                  if(Jct(54).eq.1) goto 200
                  goto 100
               endif
            end do
         endif
 
c satellites, probes, etc. not allowed on n-body tape
      else if(Nbody.gt.0) then
         if(Ict(39).gt.-2) then
            if((nql.eq.10) .or. (Ict(39).ge.0)) then
 
c look for nql on n-body tape
               do i = 1, Nbdy
                  if(Npl(i).eq.nql) then
                     JBDTST = i
                     if(Ict(39).eq.1) goto 200
                     goto 100
                  endif
               end do
            endif
         endif
      endif
  100 return
 
c body was found, now see if individual tape exists
  200 num4 = Numpln + 4
      do j = 1, num4
 
c get index into nplnt, iplnt arrays
         k = j - 4
         if(nql.eq.Nplnt(k)) then
c found nql in list of planets
            if(Iplnt(k).gt.0) JBDTST = 0
            goto 100
         endif
      end do
      goto 100
      end
