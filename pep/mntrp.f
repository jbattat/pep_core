      subroutine MNTRP(idumm,jdt,frt,nvel,lcntl,ipi)
 
      implicit none

c
c        mntrp - dec. 1975 - r. goldstein
c        velocities from numerical differentiation added  (nvel= -1)
c        by j.f.chandler, 1977 may
c        this is a control routine for interpolating the moon from its
c        integration tape (or n-body tape)

c input parameters:
      real*10 frt
      integer*4 idumm, jdt, nvel, lcntl, ipi
 
c parameter for acceleration
      real*10 xmna(3)

c          idumm     ignored
c          jdt       integer value of the coord. time of the desired pt
c          frt       fractional part of the coordinate time of the point
c          nvel  = 0 positions determined
c                  1 velocities determined
c                 -1 velocities determined from position y-vectors
c          lcntl = 0 use results of last call if possible
c                  1 do not use results (must recalculate y vectors)
c                 -1 use old results only if jd.fract is in same tab pt.
c          ip    = 1 use recieve arrays
c          ip    = 2 use send arrays
c        output is put in the xm(1-6,ip) vector
c
c        at present, only everett interpolation is allowed
c        hermite method would be indicated by kmn(88)=0
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'empcnd.inc'
      include 'param.inc'
      include 'pqind.inc'
      include 'tapdta.inc'
      include 'tapdtm.inc'
      include 'yvect.inc'

c local 
      real*10 ang, cang, omcang, sang
      integer   i, i0, icall, ip, j
      data icall/2/
      real*10 xmt(3),scwt(3,3)
 
      ip = ipi
c
c everett interpolation
      call EVTRP(jdt,frt,nvel,lcntl,icall,Ymnpos,Xm(1,ip),Moon,
     .           i_mxplprt+1, Jdm, Fm, Pm(1,ip))
c
c ad hoc precession of node
      Tgdp(ip) = jdt + frt - 0.5_10 - prmter(97)
      if(Dgdprc(3,1).ne.Soblq0 .or. Dgdprc(2,1).ne.Coblq0) then
         scwt(2,2) = Soblq0**2
         scwt(2,3) = -Soblq0*Coblq0
         scwt(3,2) = scwt(2,3)
         scwt(3,3) = Coblq0**2
 
c initialize matrices for rate=0
         call ZFILL(Gdprc,16*27)
         do j = 1, 2
            do i = 1, 3
               Gdprc(i,i,j) = 1._10
            end do
         end do
         Dgdprc(1,2) = -Coblq0
         Dgdprc(1,3) = -Soblq0
         Dgdprc(2,1) = Coblq0
         Dgdprc(3,1) = Soblq0
      endif
      if(Mcond(12).ne.0._10 .and. jdt.gt.0) then
         ang  = Mcond(12)*Tgdp(ip)
         cang = COS(ang)
         sang = SIN(ang)
         Gdprc(1,1,ip) = cang
         omcang = 1._10 - cang
         do i = 2, 3
            do j = 2, 3
               Gdprc(i,j,ip) = omcang*scwt(i,j)
            end do
            Gdprc(i,i,ip) = Gdprc(i,i,ip) + cang
            Gdprc(1,i,ip) = sang*Dgdprc(1,i)
            Gdprc(i,1,ip) = -Gdprc(1,i,ip)
         end do
         i0 = 3*iabs(nvel)
         do i = 1, 3
            xmt(i) = Xm(i + i0,ip)
         end do
c           ignore extra term in velocity (has order of rate)
c           anyway, not guaranteed to have position here
c     if(nvel.eq.0) go to 990
c     call prodct(dgdprc,xm(1,ip),xmt,3,3,1)
c     do 130 i=1,3
c 130 xm(i+3,ip)=xm(i+3,ip)+mcond(12)*xmt(i)
         call PRODCT(Gdprc(1,1,ip),xmt,Xm(i0+1,ip),3,3,1)
      endif
      return
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry MNTRPA(xmna)
c           enter here to get acceleration for same (jd,fract) as last
c           call. for everett method, uses whichever y vectors are
c           available and corresponding differentiating interpolation.
c
c           everett acceleration
      call EVTRP(0,0._10,-2,0,icall,Ymnpos,xmna,Moon,i_mxplprt+1,Jdm,
     .           Fm, Pm(1,ip))
 
      return
      end
