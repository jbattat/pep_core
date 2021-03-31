      subroutine SJ2SET
      implicit none
c sj2set - j.f.chandler - 2017/6/9
c set up for solar j2 effect on planets or probes or moon
c code copied from bodset and sbset

c array dimensions
      include 'globdefs.inc'
c common
      include 'funcon.inc'
      include 'intstf.inc'
      include 'param.inc'

c local
      real*10 cebar,comegs,cosis,ebar,omegs,sebar,si,sinis,somegs

      omegs  = Seqasc*Convd
      ebar   = Ecinc*Convd
      si     = Seqinc*Convd
      comegs = COS(omegs)
      somegs = SIN(omegs)
      cebar  = COS(ebar)
      sebar  = SIN(ebar)
      cosis  = COS(si)
      sinis  = SIN(si)
      C3(1)  = somegs*sinis
      C3(2)  = -comegs*sinis*cebar - cosis*sebar
      C3(3)  = -comegs*sinis*sebar + cosis*cebar
      Sunrd2 = (Sunrad/Ltvel/Aultsc)**2
      Shar2  = Sunhar*Sunrd2
      return
      end

