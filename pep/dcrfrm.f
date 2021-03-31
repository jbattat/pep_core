      subroutine DCRFRM(imats, b, lenb, iptr, jptrc, jptrd)
 
      implicit none

c       j.f.chandler  aug 1986   subroutine dcrfrm
c     routine for loading d-inv*f matrix from ppr nrmeqs
c           notation:
c 'mparam': saved parameter set, composed of 'ncparm' interesting
c           and 'ndparm' uninteresting parameters
c 'nparam': current parameter set; some overlap with 'mparam', 'ncparm',
c           and 'ndparm'.  the latter intersections are 'ncr' and 'ndr'
c
c        parameters
      real*10 b(1)
      integer*4 imats,lenb
      integer*2 iptr(1),jptrc(1),jptrd(1)
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'iptrst.inc'
      include 'mcnfrm.inc'
      include 'restor.inc'
      include 'rtsidesl.inc'
c
c local
      real*10 ai
      integer   i,ia,ip,ippr,j,jp,ncr,ndr,nparm
      equivalence (ai,ia)
c
c read header from imats
      call FRMHED(imats,'IMATS','    ',2,ippr,1)
      if(ippr .ne. 1) call SUICID(
     . 'NORMAL EQUATIONS NOT PREREDUCED, STOP IN DCRFRM ',12)
c
c skip series header
      read(imats)
 
c read pointer group for pre-reduced sne
      read(imats) nparm,Ncparm,Ndparm,Znsqpp,Nserpp,Nauxpp,
     .            (Icrest(i),i = 1,Ncparm),
     .            (Idrest(i),i = 1,Ndparm)
c
c get ptrs from 'mparam' to 'nparam'
      call ZFILL(Sigma,16*Nparam)
      call ZFILL(iptr,2*Mparam)
      do i = 1,Mparam
         ia     = i
         Sav(i) = ai
         end do
      call FRMMVE(Sigma,Nparam)
      do i = 1,Nparam
         ai = Sigma(i)
         if(ia .gt. 0) iptr(ia) = i
         end do
c
c     at this point:
c           icrest(1-ncparm) --> mparam
c           idrest(1-ndparm) --> mparam
c           iptr(1-mparam) --> nparam
c
c     now get ptrs from 'ncparm', 'ndparm', 'nparam' to 'ncr','ndr'
c           icrest(1-ncparm) --> ncr
c           idrest(1-ndparm) --> ndr
c           jptrc(1-nparam) --> ncr
c           jptrd(1-nparam) --> ndr
      ncr = 0
      ndr = 0
      call ZFILL(jptrc,2*Nparam)
      call ZFILL(jptrd,2*Nparam)
      do i = 1,Ncparm
         jp = iptr(Icrest(i))
         Icrest(i) = 0
         if(jp .gt. 0) then
            ncr = ncr + 1
            Icrest(i) = ncr
            jptrc(jp) = ncr
         endif
         end do
      do i = 1,Ndparm
         jp = iptr(Idrest(i))
         Idrest(i) = 0
         if(jp .gt. 0) then
            ndr = ndr + 1
            Idrest(i) = ndr
            jptrd(jp) = ndr
         endif
         end do
      if(ncr*ndr .gt. lenb)
     .     call SUICID('FBAR ARRAY TOO LARGE, STOP IN DCRFRM',9)
c
c skip v-bar and c-bar
      call BSKIP(imats,Ncparm)
c
c restore partial solution (z-bar)
      call ZFILL(Solut,16*ndr)
      call QREAD(imats,i,Sav,Ndparm)
      do i = 1,Ndparm
         ip = Idrest(i)
         if(ip .gt. 0) Solut(ip) = Sav(i)
         end do
c
c restore transformation matrix (f-bar adjoint)
      call ZFILL(b,16*ncr*ndr)
      do while( .true. )
 
c read and transpose a row
         call QREAD(imats,i,Sav,Ncparm)
         ip = Idrest(i)
         if(ip .gt. 0) then
            do j = 1,Ncparm
               jp = Icrest(j)
               if(jp .gt. 0) b(ndr*(jp-1) + ip) = Sav(j)
               end do
         endif
 
c check if done reading
         if(i .ge. Ndparm) then
            rewind imats
            Itrwnd(imats) = 0
            Ncparm = ncr
            Ndparm = ndr
 
            return
         endif
         end do
      end
