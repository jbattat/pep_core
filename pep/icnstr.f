      subroutine ICNSTR(solut)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, iep, ind, j, k, lnptr2
 
c*** end of declarations inserted by spag
 
 
c
c paul macneil january, 1978
c
 
c store adjusted process noise ic's on direct access
      real*10 solut(1)
 
      include 'filptr.inc'
      include 'filtda.inc'
 
      real*10 wbuff(40)/40*0.0_10/
 
      lnptr2 = 2*Lnptr
      do i = 1, Lnptr
         k   = Nptr(i)
         ind = i
 
c update parameter value
         wbuff(ind) = solut(k)
      end do
 
      iep = 2*Ithsep - 1
      write(Mfile, rec = iep) (wbuff(j), j = 1, lnptr2)
 
      return
      end
