      subroutine EMTRP(idumm,jdt,frt,nvel,lcntl,ip)
 
      implicit none
c
c        emtrp - dec. 1975 - r. goldstein
c        velocities from numerical differentiation added  (nvel= -1)
c        by j.f.chandler, 1977 may
c        this is a control routine for interpolating the earth-moon
c        barycenter from its integration tape
c        input parameters:
      integer*4 idumm,jdt,nvel,lcntl,ip
      real*10 frt
c          idumm     ignored
c          jdt       integer value of the coord. time of the desired
c                    point
c          frt       fractional part of the coordinate time of the point
c          nvel  = 0 positions determined
c                  1 velocities determined
c                 -1 velocities determined from position y-vectors
c          lcntl = 0 use results of last call if possible
c                  1 do not use results (must recalculate y vectors)
c                 -1 use old results only if jdt.frt is in same tab pt.
c          ip    = 1 use recieve arrays
c          ip    = 2 use send arrays
c        output is put in the xem(1-6,ip) vector
c
c        entry point 'etrp' is identical, except that it also calls
c        mntrp for the same epoch and stores the earth position in
c        'xemlsc(1-6,ip)' - control ixctl governs whether this
c        output is scaled:
c        ixctl=1 leave units as au, au/day
c             =2 convert to lt.sec., lt.sec./sec.
c             =3 compute position of site 'ip' w.r.t. sun
      integer*4 ixctli,ixctl
c
c        at present, only everett interpolation is allowed
c        hermite method would be indicated by kem(88)=0

c array dimensions
      include 'globdefs.inc'
c
c        commons
      include 'comdateq.inc'
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'pqind.inc'
      include 'sitcrd.inc'
      include 'tapdta.inc'
      include 'tapdte.inc'
      include 'tapdtm.inc'
      include 'yvect.inc'

c local
      real*10 xx
      integer*4 i,ipsv,j,jm,
     .          jo,jy,nvela
      integer*4 icall/1/
      integer*2 jtest(4)/4,8,16,32/
      integer*2 iwrn/0/
 
c argument array for acceleration
      real*10 xema(3)
 
      ixctl = 0
      goto 100
 
      entry ETRP(idumm,jdt,frt,nvel,lcntl,ip,ixctli)
      ixctl = ixctli
 
c everett interpolation
  100 call EVTRP(jdt,frt,nvel,lcntl,icall,Yempos,Xem(1,ip),
     .           Embary,i_mxplprt+1,Jdem,Fem,Pem(1,ip))
      if(jdt.gt.0) then
         Jdsvem(ip)=jdt
         Frsvem(ip)=frt
         ipsv = ip
         if(ixctl.gt.0) then
c
c compute earth position by correcting for offset from embary
            call MNTRP(idumm,jdt,frt,nvel,lcntl,ip)
            if(jdt.gt.0) then
               jo = 3*iabs(nvel)
               do jy = 1,3
                  j  = jy + jo
                  xx = Xem(j,ip) - Xm(j,ip)*Mnfct
                  if(ixctl.gt.1) then
 
c convert from au to light-sec units
                     if(jo.eq.0) xx = xx*Aultsc
                     if(jo.eq.3) xx = xx*Aultvl
 
c add on offset of observing site
                     if(ixctl.ge.3) xx = xx + Xsite(j,ip)
                  endif
                  Xemlsc(j,ip) = xx
               end do
 
c extra printout
               nvela = iabs(nvel)
               jm    = Jct(6)/jtest(nvela + 1)
               if(mod(jm,2).ne.0) then
                  jo = 3*nvela
                  if(Line.gt.56) call OBSPAG
                  write(Iout,110) ip,ixctl,
     .                             (Xemlsc(i+jo,ip),i = 1,3)
  110             format(' ETRP:    IP=',i1,' IXCTL=',i1,t52,'X=',
     .                   1p,3D23.15)
                  Line = Line + 1
               endif
            endif
         endif
      endif
      return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry EMTRPA(xema)
c enter here to get acceleration for same (jdt,frt) as last
c call. for everett method, uses whichever y vectors are
c available and corresponding differentiating interpolation.
      if(ipsv.ne.1 .and. iwrn.ne.1) then
         iwrn = 1
         write(Iout,150)
  150    format('0',8('* '),'EMTRPA MAY BE WRONG - SEE J.F.CHANDLER')
      endif
c
c everett acceleration
      call EVTRP(0,0._10,-2,0,icall,Yempos,xema,Embary,i_mxplprt+1,
     .           Jdem,Fem,Pem(1,ipsv))
      return
      end
