      subroutine USE(b)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, kode, m, n, nd, nef, nei, nf, ni, npf, npi
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine use
c paul macneil january, 1978 modified for iteration
c
 
c parameter is b matrix - dim set in block data
      real*10   b(1)
c
c common
      include 'fcntrl.inc'
      include 'filstf.inc'
      real*10   bthts(300,7)
      equivalence (G(1,1),bthts(1,1))
      include 'filtda.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'revptr.inc'
      include 'rtside.inc'
      include 'scail.inc'
 
c temporary work arrays
      common/WRKCOM/ Wtrans(300,7),Has(20,20)
      real*10   Wtrans, Has
c
c local
      character*80 tin,tfi
      equivalence (tin,G(1,1)),(tfi,G(11,1))
      character*16 nmi(1000),nmf(1000)
      equivalence (nmi,G(1,1))
      equivalence (nmf,G(1,3))
      real*10 buff(1000)
      equivalence (buff,Scale)
      integer*2 nptr(20)
      logical*4 shorhs
      character*8 label/'U   /  '/
c
c page heading
      call PAGSET('USE FILTER MATRICES ', 5)
      call NEWPG
      i = Fict(1)
      shorhs = mod(i/2,2).eq.1
c
c read headers from insne and filter
      read(Insne) tin,npi,nei
      read(Filter) tfi,npf,nef,Npnp
      write(Iout,100) tin,npi,nei
  100 format('-INPUT SNE:             TITLE= ', a80/' NPARAM=', i4,
     .       '  NEPOCH=', i4)
      write(Iout,200) tfi,npf,nef,Npnp
  200 format('-INPUT FILTER DATA SET: TITLE= ', a80/' NPARAM=', i4,
     .       '  NEPOCH=', i4, '  NPRNOISE=', i3)
c
c assuming for now that npi = npf
      if(npi.ne.npf) call SUICID('NPI.NE.NPF, STOP IN USE ', 6)
c
c start outsne
      write(Outsne) npi,nei,Npnp
      write(Outsne) tin
      write(Outsne) tfi
c
c transfer residual stats from insne
      read(Insne) i,(buff(j),j = 1,3)
      write(Outsne) i,(buff(j),j = 1,3)
c
c names
      read(Insne) (nmi(i),i = 1,npi)
      read(Filter) (nmf(i),i = 1,npf)
c
c compare names
      ni = 1
      do nf = 1, npf
 
cdw   iptr(ni) = nf
         if(nmi(ni).eq.nmf(nf)) ni = ni + 1
         end do
      if(ni-1.ne.npi) call SUICID(' PARAMETER MIXUP IN USE ', 6)
      write(Outsne) (nmi(i),i = 1,npi)
c
c nominals
      if(Iconof.gt.0) Itrwnd(Iconof) = 1
 
c skip epochs for nomtrs
      if(Iconof.gt.0) read(Iconof)
c
c process pointer
      read(Filter) (nptr(i),i = 1,Npnp)
      write(Outsne) (nptr(i),i = 1,Npnp)
c
c         list of output epochs
c
c         construct reverse pointer for noise parameters
c         rptr(i) = position in noise matrix of parameter i
      do i = 1, npi
         Rptr(i) = 0
      end do
      do i = 1, Npnp
         Rptr(nptr(i)) = i
      end do
      nd = npi - Npnp
c
c a priori for forward data (delta u(1))
      call RINSNE(b,Side,buff,npi,0,1)
      call EBCDIX(1,label,3,2)
      call EBCDIX(1,label,6,2)
      if(shorhs) call GPMPO(npi,Side,0,1,label,8)
c
c loop over epochs from insne
      n = nei - 1
      do i = 1, n
c
c b=i(k/k)
c g=i(k/k-1)*s
c read wtrans, has, bthts, and b from filter data set
         call RFILTR(b,buff,bthts,npf,Npnp,i+1,Wtrans,Has,nd)
c
c propagate residual
c perform part of equation 20 calc
c (e-i(k/k-1)*s) u(k-1)
         call PROPGU(bthts,Side,buff,nptr,Rptr,Npnp,npi,Has,nd)
         call EBCDIX(i+1,label,3,2)
         if(shorhs) call GPMPO(npi,Side,0,0,label,8)
c
c subtract i(k/k-1)*w
         if(Iconof.gt.0) call NOMTRS(Side,buff,nptr,npi,Wtrans)
c
c write rhs to direct access for later use
c m file now has above quantity
c (u(k)- delta u)
         j = i*(npi+1) + 1
         if(i.lt.n) write(Mfile,rec = j) (Side(m),m = 1,npi)
c
c add delta u(i+1)
         call RINSNE(b,Side,buff,npi,0,i+1)
         call EBCDIX(i+1,label,6,2)
         if(shorhs) call GPMPO(npi,Side,0,0,label,8)
         end do
c
c write out nth epoch
c i(ne/ne), u(ne)
c nth epoch has only forward data
      call WOTSNE(b,Side,buff,npi,1,nei)
c
c
c
c doing smoothed multi-epoch solutions?
      if(Fict(2).eq.0) return
c
c check data sets for backward data
      read(Insne) i
      read(Filter) j
      if(i.ne.0 .or. j.ne.0) call SUICID('NO BACKWARD DATA IN USE ', 6)
c
c reset normal equations
      Nparam = Nparam*2
      call NRMSET(1)
      Nparam = Nparam/2
c
c a priori for backward data (delta u(n))
      call RINSNE(b,Side,buff,npi,0,nei)
      call EBCDIX(nei,label,3,2)
      call EBCDIX(nei,label,6,2)
      if(shorhs) call GPMPO(npi,Side,0,1,label,8)
c
c loop backward
      kode = 3
      do i = 1, n
         if(i.eq.n) kode = 2
         m = n - i + 1
c
c b=i(k/k)
c g=i(k/k-1)*s
c read wtrans, has, bthts, and b from filter data set
         call RFILTR(b,buff,bthts,npf,Npnp,i + 1,Wtrans,Has,nd)
c
c propagate residual
c perform part of eq 20
c (e-g) u(k-1)
         call PROPGU(bthts,Side,buff,nptr,Rptr,Npnp,npi,Has,nd)
         call EBCDIX(m,label,3,2)
         if(shorhs) call GPMPO(npi,Side,0,0,label,8)
c
c subtract i(k/k-1)*w
         if(Iconof.gt.0) call NOMTRS(Side,buff,nptr,npi,Wtrans)
c
c add delta u(m)
c and thereby complete eq 20 to form u(k)
         call RINSNE(b,Side,buff,npi,0,m)
         call EBCDIX(m,label,6,2)
         if(shorhs) call GPMPO(npi,Side,0,0,label,8)
c
c write out mth epoch
c write out i(k/k) (backwards) on outsne, reads in (u(k/k)-delta)
c (forward) from mfile and write it out outsne
         call WOTSNE(b,Side,buff,npi,kode,m)
      end do
      if(Iconof.gt.0) rewind Iconof
      if(Iconof.gt.0) Itrwnd(Iconof) = 0
      return
      end
