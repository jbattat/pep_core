      subroutine REDCPY(nprmt, nplntt, sert, sitt, prmt, lprmt,
     .                  ndim, ndims, ndimp, l, npl,
     .                  nprm, nplnt, ser, site, prm, lprm)
 
      implicit none
 
c
c J.F. Chandler - 1992 June - subroutine REDCPY
c
c---purpose:     to copy input series/site-related parameters from
c                temporary storage to final storage, sorted by planet
c
c---input:
c       nprmt= number of relevant parameters in temp storage
c                (integer*2)
c       sert = name of the observation series in temp storage
c                (character*4)
c       sitt = name(s) of the observing site(s) in temporary storage
c                (character*4)
c       prmt = array of parameters in temporary storage
c                (real*4)
c       lprmt= flags designating whether to adjust the parameters
c                (integer*2)
c       ndim = maximum number of parameters in prm/prmt
c                (integer*4)
c       ndims= number of site names to copy
c                (integer*4)
c       ndimp= number of sets of parameters in storage
c                (integer*4)
c       l    = index for copying to final storage
c                (integer*4)
c       npl  = current planet number for copying
c                (integer*4)
c
c---output:
c       nprm = number of relevant parameters in final storage
c                (integer*2)
c       ser  = name of the observation series in final storage
c                (character*4)
c       site = name(s) of the observing site(s) in final storage
c                (character*4)
c       prm  = array of parameters in final storage
c                (real*4)
c       lprm = flags designating whether to adjust the parameters
c                (integer*2)
c
      integer*4 i, l, n, ndim, ndims, ndimp, npl
      real*4    prm(ndim,ndimp), prmt(ndim,ndimp)
      integer*2 nprm(ndimp), nprmt(ndimp), nplnt(ndimp), nplntt(ndimp),
     .          lprm(ndim,ndimp), lprmt(ndim,ndimp)
      character*4 ser(ndimp), sert(ndimp), site(ndims,ndimp),
     .            sitt(ndims,ndimp)
 
      do n=1,ndimp
        if(nplntt(n) .eq. npl) then
 
c found a match, copy this set
          l=l+1
          if(ndim.gt.2) nprm(l) = nprmt(n)
          nplnt(l) = nplntt(n)
          nplntt(n) = -1
          ser(l) = sert(n)
          do i = 1, ndims
             site(i,l) = sitt(i,n)
             end do
          do i = 1, ndim
             prm(i,l) = prmt(i,n)
             lprm(i,l) = lprmt(i,n)
             end do
          end if
        end do
 
      return
      end
