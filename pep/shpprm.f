      subroutine SHPPRM(model,cycle,prefix,k,names,nshp,ngrid,xnom)
 
      implicit none
c
c subroutine shpprm       6/76   s.brody
c generates names for shape models with nshp gt 0
c (spherical harmonics names (nshp=0) done in hrmprm)
      real*10 xnom(1)
      character*8 names(2,1),prefix
      integer*4 model,k,nshp,ngrid,cycle
c
c array dimensions
      include 'globdefs.inc'
 
c shape coeffs and l-vectors
      common/SCOEF4/ Shape(2000),Lshape(4000)
      real*10 Shape
      integer*2 Lshape
 
c equivalencing to use real*4 shape coeffs
      real*4    shape4(u_stdsz,4000/u_stdsz)
      equivalence (Shape,shape4)
 
c locals
      integer   i,ifs,index,inum,j,l,nend,nmax,
     .          nstop,num
      character*8    fours(8)/'PA (  )', 'PB (  )', 'PC (  )', 
     .         'PD (  )', 'PAP(  )', 'PCP(  )', 'PAPP', 'PBPP'/
      character*8    grid/'GRD(   )'/, temp
      character*1 tem(8)
      equivalence (temp,tem)
      character*1 nums(10)/'1', '2', '3', '4', '5', '6', '7', '8',
     .                     '9', '0'/
c
      if(nshp.eq.1) then
c        fourier series model, nshp=1
c        122 coefficients are used, stored in shape
c        as shape(i),i=1,122*cycle,cycle
c
         nmax = 20
         inum = 0
         ifs  = 1
         nend = 122*cycle
         do i = model, nend, cycle
            inum = inum + 1
            if(inum.gt.nmax) then
               inum = 1
               ifs  = ifs + 1
               if(ifs.eq.7) nmax = 1
            endif
            if(Lshape(i).gt.0) then
               k = k + 1
               names(1,k) = prefix
               xnom(k)     = Shape(i)
               temp = fours(ifs)
               if(ifs.lt.7) then
                  if(inum.ge.10) then
                     j = mod(inum,10)
                     if(j.eq.0) j = 10
                     tem(6) = nums(j)
                     tem(5) = nums(inum/10)
                  else
                     tem(6) = nums(inum)
                  endif
               endif
               names(2,k) = temp
            endif
         end do
      else if(nshp.eq.2) then
c        altitude grid-local model, nshp=2
c        there are ngrid coefficients stored as
c        real*4 variables in the real*8 shape array
c        a pair of r*4 coefs for a given planet
c        can be extracted at every applicable element
c        of shape (each pair separated by 'cycle' elements)
         inum  = model
         nend  = ((ngrid+u_stdsz-1)/u_stdsz)*cycle
         nstop = (ngrid - 1)*cycle + model
         do i = model, nend, cycle
            do j = 1,u_stdsz
               if(inum.gt.nstop) return
               if(Lshape(inum).gt.0) then
                  k = k + 1
                  names(1,k) = prefix
                  xnom(k)    = shape4(j,i)
                  temp  = grid
                  num   = (inum-1)/cycle + 1
                  index = 7
                  do while( num.ge.10 )
                     l = mod(num,10)
                     if(l.eq.0) l = 10
                     tem(index) = nums(l)
                     index = index - 1
                     num   = num/10
                  end do
                  tem(index)  = nums(num)
                  names(2,k) = temp
               endif
               inum = inum + cycle
            end do
 
         end do
      endif
c
c names for additional shape models not coded
      return
      end
