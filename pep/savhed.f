      subroutine SAVHED(imats, namhed)
 
      implicit none
 
c
c ash - amuchastegui / december 1969 / subroutine savhed
c
c arguments
      integer*4 imats
      character*8 namhed

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'eqenox.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'lcntrl.inc'
      include 'loadid.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'phase.inc'
      include 'plnhar.inc'
      integer*4 ngdpts(4)
      equivalence (ngdpts, Scontl(1,9))
      include 'psrstf.inc'
      include 'rdbias.inc'
      include 'scoef4.inc'
      integer*2 numphr
      equivalence (numphr, Nmphar)
      include 'skystf.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'
      include 'zeroes.inc'
c
c local
      integer*4 i,ippr,j,n,ngd,ngd2,nsh
      integer*2 nzone1, ntess1, numdt1
      integer*2 npln3/3/, ncntr3/0/
      integer*2 npln10/10/, ncnt10/3/
      integer*2 mpln3/-3/, mcntr3/3/
      integer*2 mpln10/-10/, mcnt10/10/
      integer*2 nsitcr/6/
      character*8 nppr/'PPRNRMEQ'/,
     1    embary/' EMBARY '/,erotat/' EROTAT '/,mrotat/' MROTAT '/
 
c ::::::::: temporary :::::::::::
      integer*2 nphr1/6/, nshp2x/-1/, nshap(6)/6*0/
 
c :::::::::
      do i = 1, 4
         nshap(i+2) = Nshape(i)
      end do
c ::::::::: end :::::::::::::::::
c
      write(imats) namhed, Heding, Date, Lnklvl
      numdt1 = Numdt
      if(Numdt.le.0) numdt1 = 1
      if(Jddt0.le.0 .and. Numdt.gt.0) numdt1 = Numdt + 400
      ippr = 0
      if(namhed.eq.nppr) ippr = 1
      write(imats) u_nmprm,u_nmbod,prmter,Lprm,Ict,Nparam,Iterat,Npage,
     .             Numpln,numphr,Numsit,Numspt,Numstr,Numpsr,Numrbs,
     .             Numeqn,Numphs,Numdt,numdt1,
     .             (Dt(i),Jddt(i),Ldt(i),i = 1,numdt1),ippr,
     .             nphr1,nshp2x,nshap,Jddt0,Lnklvl,nsitcr,Izero
c
c earth-moon barycenter parameters
      write(imats) npln3,ncntr3,Jdem0,Econd,Lem,embary
c
c moon i.c. and parameters
      write(imats) npln10,ncnt10,Jdmn0,Mcond,Lmn,Aplnt(17)
c
c earth rotation i.c. and parameters
      write(imats) mpln3,mcntr3,Jder0,Ercond,Ler,erotat
c
c moon rotation initial conditions and parameters
      write(imats) mpln10,mcnt10,Jdmr0,Mrcond,Lmr,mrotat
c
c planet initial conditions and parameters
      if(Numpln.gt.0) then
         do i = 1,Numpln
            write(imats) Nplnt(i),Npcent(i),Jdpl0(i),
     .       (Pcond(j,i),j=1,u_nmbod),(Lpl(j,i),j=1,u_nmbod),
     .       Aplnt(i)
         end do
      endif
c
c write earth harmonics
      nzone1 = Nezone - 1
      if(nzone1.le.0) nzone1 = 1
      ntess1 = (Netess*(Netess+1))/2 - 1
      if(ntess1.le.0) ntess1 = 1
      write(imats) npln3, Nezone, nzone1, (Ezhar(i),i=1,nzone1),
     .             (Lezhar(i),i=1,nzone1), Netess, ntess1,
     .             (Echar(i),i=1,ntess1), (Lechar(i),i=1,ntess1),
     .             (Eshar(i),i=1,ntess1), (Leshar(i),i=1,ntess1)
c
c write moon  harmonics
      nzone1 = Nmzone - 1
      if(nzone1.le.0) nzone1 = 1
      ntess1 = (Nmtess*(Nmtess+1))/2 - 1
      if(ntess1.le.0) ntess1 = 1
      write(imats) npln10, Nmzone, nzone1, (Mzhar(i),i=1,nzone1),
     .             (Lmzhar(i),i=1,nzone1), Nmtess, ntess1,
     .             (Mchar(i),i=1,ntess1), (Lmchar(i),i=1,ntess1),
     .             (Mshar(i),i=1,ntess1), (Lmshar(i),i=1,ntess1)
c
c write planet harmonics
      if(numphr.gt.0) then
         do j = 1, numphr
            nsh = Nshape(j) + 1
            if(nsh.eq.2) then
 
c fourier
               write(imats) Nplhar(j), (Pzhar(j,i),i=1,122),
     .                      (Lpzhar(j,i),i=1,122)
            else if(nsh.eq.3) then
 
c grid
               ngd  = ngdpts(j)
               ngd2 = (ngd+u_stdsz-1)/u_stdsz
               write(imats) Nplhar(j), ngd, ngd2, (Scontl(j,i),i=1,9),
     .                     (Pzhar(j,i),i=1,ngd2), (Lpzhar(j,i),i=1,ngd)
            else
c
c spherical harmonics
               nzone1 = Npzone(j) - 1
               if(nzone1.le.0) nzone1 = 1
               ntess1 = (Nptess(j)*(Nptess(j)+1))/2 - 1
               if(ntess1.le.0) ntess1 = 1
               write(imats) Nplhar(j), Npzone(j), nzone1,
     .               (Pzhar(j,i),i=1,nzone1), (Lpzhar(j,i),i=1,nzone1),
     .               Nptess(j), ntess1, (Pchar(j,i),i=1,ntess1),
     .               (Lpchar(j,i),i=1,ntess1), (Pshar(j,i),i=1,ntess1),
     .               (Lpshar(j,i),i=1,ntess1)
            endif
 
         end do
      endif
c
c sites
      if(Numsit.gt.0) then
         write(imats) ((Site(j,i),j=1,2), Kscrd(i), T0site(i),
     .    (Scord(j,i),j=1,nsitcr), (Lscrd(j,i),j=1,nsitcr),
     .    i=1,Numsit)
      endif
c
c spots
      if(Numspt.gt.0) write(imats) (Spot(i), Nsplnt(i),
     .             (Spcord(j,i),j=1,3), (Lspcrd(j,i),j=1,3), i=1,Numspt)
c
c stars
      if(Numstr.gt.0) then
         n = 1
         do i = 1, Numstr
            if(Nskycf(i).gt.n) n = Nskycf(i)
         end do
         write(imats) (Ctlgnm(i), n, Nskycf(i), (Skycf(j,i),j=1,n),
     .                (Lskycf(j,i),j=1,n), i = 1, Numstr)
      endif
c
c pulsar parameters
      if(Numpsr.gt.0) write(imats)
     .                          (Sptpsr(i), Jdpsr0(i), Plspr(i),
     .                          Ntypsr(i), (Psrcn(j,i),j=1,16),
     .                          (Lpsrcn(j,i),j=1,16), i = 1, Numpsr)
c
c radar biases
      if(Numrbs.gt.0) write(imats)
     .                          ((Rdbsit(j,i),j=1,2), Rdbser(i),
     .                          (Rbias(j,i),j=1,2), (Lrbs(j,i),j=1,2),
     .                          Nplrbs(i), i = 1, Numrbs)
c
c equinox equator
      if(Numeqn.gt.0) write(imats)
     .                          (Eqnsit(i), Eqnser(i), Denox(i),
     .                          Dequat(i), Dlat(i), (Leqn(j,i),j=1,3),
     .                          i = 1, Numeqn)
c
c phase corrections
      if(Numphs.gt.0) write(imats)
     .                          (Phsit(i), Phser(i), (Aphase(j,i),j=1,
     .                          9), Ncphs(i), Nplphs(i),
     .                          (Lphs(j,i),j=1,9), i = 1, Numphs)
 
      return
      end
