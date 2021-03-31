      subroutine RNGMOD
 
      implicit none
c
c        r.b. goldstein  march 1978
c        this routine processes the range observation that is in
c        modulo form. it does the following:
c             a. demods the range observation
c             b. "fixes" the observed range based on the
c                assumption that it was in error by an integral
c                number of "subcodes"
c        the routine can do the above tasks a and b independently
c
c        some variables use and their meaning:
c             cdlmax: equivalenced to relult(2). the maximum
c                     codelength. the obs. range is given modulo this
c                     value.
c             ict(25):alim=2**(-ict(25))
c             alim:    the smallest value of cdlmax for which demodding
c                      is done (seconds).
c             ict(67):blim=2**(-abs(ict(67)))
c                     if(ict(67).lt.0, use the predicted
c                     residual instead of the current resid.
c             blim  : maximum allowed range residual (seconds).
c
c        print out:
c             a. if no alteration to obs, no p/o
c             b. if alim>cdlen, no p/o
c             c. demod info: nrng, cdlmax,raw range
c             d. fix info: subcode, binary and decimal m, prefix resid.

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      include 'obscrd.inc'
      real*10 obs,comp,cdlmax
      equivalence (Result,obs),(comp,Deriv(2,1)),
     .            (cdlmax,Result(2))

c local variables
      real*10 alim,alt,blim,blim2,cmo,dblim,
     .          rawrng,resid,resid2,s,sc
      integer   i,ia67,ibit,ict25,ict67,ixor,m,mmn,mmo,mnew,
     .          mold,ndif,nrng
      logical*4 prdrs1/.false./, init/.false./, prdres/.false./
c
c
      if(Ict(25).ne.0 .or. Ict(67).ne.0) then
         if(.not. (init)) then
c
c initialize some numbers
c
            init  = .true.
            ict25 = Ict(25)
            ict67 = Ict(67)
            ia67  = iabs(ict67)
            alim  = 2._10**(-ict25)
            blim  = 2._10**(-ia67)
            blim2 = blim/2._10
            dblim = 2._10*blim
            if(Ict(67).lt.0) prdres = .true.
         endif
 
         prdrs1 = prdres
         if(Numsav.lt.47 .or. Save(47).eq.0._10) prdrs1 = .false.
c
c
c get the residual to use in the computations
         resid = obs - comp
         if(prdrs1) resid = Save(47)
         cmo = -resid
c
c now demod if requested
c
         if(ict25.ne.0) then
            if(alim.lt.cdlmax) then
               nrng   = cmo/cdlmax + SIGN(0.5_10,cmo)
               rawrng = obs
               alt    = nrng*cdlmax
               obs    = obs + alt
               if(prdrs1) Save(47) = Save(47) + alt
c
c print out altered observable
               if(nrng.ne.0) then
                  write(Iout,10) nrng,cdlmax,rawrng
   10             format(100x,i10,6p,2F10.3)
                  Line = Line + 1
               endif
            endif
         endif
c
c
c
c
c        now fix the observed if necessary
c        this must use a de-modded residual
c
         if(ict67.ne.0) then
            if(blim.lt.cdlmax) then
               resid2 = obs - comp
               if(prdrs1) resid2 = Save(47)
               if(ABS(resid2).ge.blim2) then
c
c        a fix is necessary. first get subcode
c        subcode is defined as the smallest integral power of
c        two divisor of cdlmax that is larger than blim
c
                  sc = cdlmax
                  do i = 1,20
                     if(sc.lt.dblim) goto 20
                     sc = sc/2._10
                  end do
                  call SUICID('RNGMOD ERROR',3)
c
c now get the multiple of sc to alter obs by
   20             m = ABS(resid2)/sc + 0.5_10
c
c now alter obs and new predicted residual
                  s    = SIGN(1._10,resid2)
                  alt  = s*m*sc
                  mold = obs/sc + .5_10
                  mnew = (obs - alt)/sc + .5_10
c
c form exclusive 'or' of old and new
c also count bits that are different
                  mmo  = mold
                  mmn  = mnew
                  ibit = 1
                  ixor = 0
                  ndif = 0
                  do i = 1,28
                     if(mod(mmo-mmn,2).ne.0) then
                        ndif = ndif + 1
                        ixor = ixor + ibit
                     endif
                     mmo  = mmo/2
                     mmn  = mmn/2
                     ibit = ibit*2
                  end do
 
c reject the change if too many bits are different
                  if(ndif.le.4) then
 
                     obs = obs - alt
                     if(prdrs1) Save(47) = Save(47) - alt
                  endif
c
c now printout the results
                  write(Iout,30) m,ndif,ixor,sc,resid2
   30             format(83x,2I6,1x,z7,2x,6p,2F10.3)
                  Line = Line + 1
               endif
            endif
         endif
      endif
c
c
c
      return
      end
