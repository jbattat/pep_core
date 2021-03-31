      subroutine RG(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
C***BEGIN PROLOGUE  RG
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A2
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes eigenvalues and, optionally, eigenvectors of a
C            real general matrix.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     To find the eigenvalues and eigenvectors (if desired)
C     of a REAL GENERAL matrix.
C
C     On Input
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A.
C
C        A  contains the real general matrix.
C
C        MATZ  is an integer variable set equal to zero if
C        only eigenvalues are desired.  Otherwise it is set to
C        any non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        WR  and  WI  contain the real and imaginary parts,
C        respectively, of the eigenvalues.  Complex conjugate
C        pairs of eigenvalues appear consecutively with the
C        eigenvalue having the positive imaginary part first.
C
C        Z  contains the real and imaginary parts of the eigenvectors
C        if MATZ is not zero.  If the J-th eigenvalue is real, the
C        J-th column of  Z  contains its eigenvector.  If the J-TH
C        eigenvalue is complex with positive imaginary part, the
C        J-th and (J+1)-th columns of  Z  contain the real and
C        imaginary parts of its eigenvector.  The conjugate of this
C        vector is the eigenvector for the conjugate eigenvalue.
C
C        IERR  is an integer output variable set equal to an
C        error completion code described in section 2B of the
C        documentation.  The normal completion code is zero.
C
C        IV1  and  FV1  are temporary storage arrays.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  BALANC,BALBAK,ELMHES,ELTRAN,HQR,HQR2
C***END PROLOGUE  RG
C
      implicit none
      integer n,nm,is1,is2,ierr,matz
      real*10 a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
      integer iv1(n)
C
C***FIRST EXECUTABLE STATEMENT  RG
      if(n .le. nm) goto 10
      ierr = 10 * n
      goto 50
c
   10 call BALANC(nm,n,a,is1,is2,fv1)
      call ELMHES(nm,n,is1,is2,a,iv1)
      if(matz .ne. 0) goto 20
c     .......... find eigenvalues only ..........
      call HQR(nm,n,is1,is2,a,wr,wi,ierr)
      goto 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call ELTRAN(nm,n,is1,is2,a,iv1,z)
      call HQR2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if(ierr .ne. 0) goto 50
      call BALBAK(nm,n,is1,is2,fv1,n,z)
   50 return
      end
      subroutine BALANC(nm,n,a,low,igh,scale)
C***BEGIN PROLOGUE  BALANC
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1A
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Balances a general REAL matrix and isolates eigenvalue
C            whenever possible.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BALANCE,
C     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
C     HANDBOOk FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
C
C     This subroutine balances a REAL matrix and isolates
C     eigenvalues whenever possible.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        A contains the input matrix to be balanced.
C
C     On OUTPUT
C
C        A contains the balanced matrix.
C
C        LOW and IGH are two integers such that A(I,J)
C          is equal to zero if
C           (1) I is greater than J and
C           (2) J=1,...,LOW-1 or I=IGH+1,...,N.
C
C        SCALE contains information determining the
C           permutations and scaling factors used.
C
C     Suppose that the principal submatrix in rows LOW through IGH
C     has been balanced, that P(J) denotes the index interchanged
C     with J during the permutation step, and that the elements
C     of the diagonal matrix used are denoted by D(I,J).  Then
C        SCALE(J) = P(J),    for J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     The order in which the interchanges are made is N to IGH+1,
C     then 1 TO LOW-1.
C
C     Note that 1 is returned for IGH if IGH is zero formally.
C
C     The ALGOL procedure EXC contained in BALANCE appears in
C     BALANC  in line.  (Note that the ALGOL roles of identifiers
C     K,L have been reversed.)
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BALANC
C
      implicit none
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      real*10 a(nm,n),scale(n)
      real*10 c,f,g,r,s,b2,radix
      logical noconv
C
C***FIRST EXECUTABLE STATEMENT  BALANC
      radix = 16
c
      b2 = radix * radix
      k = 1
      l = n
      goto 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if(j .eq. m) goto 50
c
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
c
   50 goto (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if(l .eq. 1) goto 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if(i .eq. j) goto 110
            if(a(j,i) .ne. 0.0_10) goto 120
  110    continue
c
         m = l
         iexc = 1
         goto 20
  120 continue
c
      goto 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if(i .eq. j) goto 150
            if(a(i,j) .ne. 0.0_10) goto 170
  150    continue
c
         m = k
         iexc = 2
         goto 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0_10
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0_10
         r = 0.0_10
c
         do 200 j = k, l
            if(j .eq. i) goto 200
            c = c + abs(a(j,i))
            r = r + abs(a(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if(c .eq. 0.0_10 .or. r .eq. 0.0_10) goto 270
         g = r / radix
         f = 1.0_10
         s = c + r
  210    if(c .ge. g) goto 220
         f = f * radix
         c = c * b2
         goto 210
  220    g = r * radix
  230    if(c .lt. g) goto 240
         f = f / radix
         c = c / b2
         goto 230
c     .......... now balance ..........
  240    if((c + r) / f .ge. 0.95_10 * s) goto 270
         g = 1.0_10 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
c
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
c
  270 continue
c
      if(noconv) goto 190
c
  280 low = k
      igh = l
      return
      end
      subroutine BALBAK(nm,n,low,igh,scale,m,z)
C***BEGIN PROLOGUE  BALBAK
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of REAL general matrix from
C            eigenvectors of matrix output from BALANC.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BALBAK,
C     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
C     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
C
C     This subroutine forms the eigenvectors of a REAL GENERAL
C     matrix by back transforming those of the corresponding
C     balanced matrix determined by  BALANC.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH are integers determined by  BALANC.
C
C        SCALE contains information determining the permutations
C          and scaling factors used by  BALANC.
C
C        M is the number of columns of Z to be back transformed.
C
C        Z contains the real and imaginary parts of the eigen-
C          vectors to be back transformed in its first M columns.
C
C     On OUTPUT
C
C        Z contains the real and imaginary parts of the
C          transformed eigenvectors in its first M columns.
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BALBAK
C
      implicit none
      integer i,j,k,m,n,ii,nm,igh,low
      real*10 scale(n),z(nm,m)
      real*10 s
C
C***FIRST EXECUTABLE STATEMENT  BALBAK
      if(m .eq. 0) goto 200
      if(igh .eq. low) goto 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
c
  110 continue
c     ......... for i=low-1 step -1 until 1,
c               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if(i .ge. low .and. i .le. igh) goto 140
         if(i .lt. low) i = low - ii
         k = scale(i)
         if(k .eq. i) goto 140
c
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine ELMHES(nm,n,low,igh,a,int)
C***BEGIN PROLOGUE  ELMHES
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B2
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduces REAL general matrix to upper Hessenberg form
C            stabilized elementary similarity transformations.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMHES,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     Given a REAL GENERAL matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     LOW through IGH to upper Hessenberg form by
C     stabilized elementary similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  BALANC.  If  BALANC  has not been used,
C          set LOW=1, IGH=N.
C
C        A contains the input matrix.
C
C     On OUTPUT
C
C        A contains the Hessenberg matrix.  The multipliers
C          which were used in the reduction are stored in the
C          remaining triangle under the Hessenberg matrix.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction.
C          Only elements LOW through IGH are used.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ELMHES
C
      implicit none
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real*10 a(nm,n)
      real*10 x,y
      integer int(igh)
C
C***FIRST EXECUTABLE STATEMENT  ELMHES
      la = igh - 1
      kp1 = low + 1
      if(la .lt. kp1) goto 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0_10
         i = m
c
         do 100 j = m, igh
            if(ABS(a(j,mm1)) .le. ABS(x)) goto 100
            x = a(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if(i .eq. m) goto 130
c    .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
c
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c    .......... end interchange ..........
  130    if(x .eq. 0.0_10) goto 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            y = a(i,mm1)
            if(y .eq. 0.0_10) goto 160
            y = y / x
            a(i,mm1) = y
c
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
c
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
c
  160    continue
c
  180 continue
c
  200 return
      end
      subroutine ELTRAN(nm,n,low,igh,a,int,z)
C***BEGIN PROLOGUE  ELTRAN
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Accumulates the stabilized elementary similarity
C            transformations used in the reduction of a REAL general
C            matrix to upper Hessenberg form by ELMHES.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine accumulates the stabilized elementary
C     similarity transformations used in the reduction of a
C     REAL GENERAL matrix to upper Hessenberg form by  ELMHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  BALANC.  If  BALANC  has not been used,
C          set LOW=1, IGH=N.
C
C        A contains the multipliers which were used in the
C          reduction by  ELMHES  in its lower triangle
C          below the subdiagonal.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction by  ELMHES.
C          only elements LOW through IGH are used.
C
C     On OUTPUT
C
C        Z contains the transformation matrix produced in the
C          reduction by  ELMHES.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ELTRAN
C
      implicit none
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real*10 a(nm,igh),z(nm,n)
      integer int(igh)
C
C***FIRST EXECUTABLE STATEMENT  ELTRAN
      do 80 i = 1, n
c
         do 60 j = 1, n
   60    z(i,j) = 0.0_10
c
         z(i,i) = 1.0_10
   80 continue
c
      kl = igh - low - 1
      if(kl .lt. 1) goto 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
c
         i = int(mp)
         if(i .eq. mp) goto 140
c
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0_10
  130    continue
c
         z(i,mp) = 1.0_10
  140 continue
c
  200 return
      end
      subroutine HQR(nm,n,low,igh,h,wr,wi,ierr)
C***BEGIN PROLOGUE  HQR
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C2B
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes eigenvalues of a REAL upper Hessenberg matrix
C            using the QR method.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure HQR,
C     NUM. MATH. 14, 219-231(1970) by Martin, Peters, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
C     This subroutine finds the eigenvalues of a REAL
C     UPPER Hessenberg matrix by the QR method.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  BALANC.  If  BALANC  has not been used,
C          set LOW=1, IGH=N.
C
C        H contains the upper Hessenberg matrix.  Information about
C          the transformations used in the reduction to Hessenberg
C          form by  ELMHES  or  ORTHES, if performed, is stored
C          in the remaining triangle under the Hessenberg matrix.
C
C     On OUTPUT
C
C        H has been destroyed.  Therefore, it must be saved
C          before calling  HQR  if subsequent calculation and
C          back transformation of eigenvectors is to be performed.
C
C        WR and WI contain the real and imaginary parts,
C          respectively, of the eigenvalues.  The eigenvalues
C          are unordered except that complex conjugate pairs
C          of values appear consecutively with the eigenvalue
C          having the positive imaginary part first.  If an
C          error exit is made, the eigenvalues should be correct
C          for indices IERR+1,...,N.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  HQR
C
      implicit none
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      real*10 h(nm,n),wr(n),wi(n)
      real*10 p,q,r,s,t,w,x,y,zz,norm,s1,s2
      logical notlas
C
C***FIRST EXECUTABLE STATEMENT  HQR
      ierr = 0
      norm = 0.0_10
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + ABS(h(i,j))
c
         k = i
         if(i .ge. low .and. i .le. igh) goto 50
         wr(i) = h(i,i)
         wi(i) = 0.0_10
   50 continue
c
      en = igh
      t = 0.0_10
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if(en .lt. low) goto 1001
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if(l .eq. low) goto 100
         s = ABS(h(l-1,l-1)) + ABS(h(l,l))
         if(s .eq. 0.0_10) s = norm
         s2 = s + ABS(h(l,l-1))
         if(s2 .eq. s) goto 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if(l .eq. en) goto 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if(l .eq. na) goto 280
      if(itn .eq. 0) goto 1000
      if(its .ne. 10 .and. its .ne. 20) goto 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = ABS(h(en,na)) + ABS(h(na,enm2))
      x = 0.75_10 * s
      y = x
      w = -0.4375_10 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = ABS(p) + ABS(q) + ABS(r)
         p = p / s
         q = q / s
         r = r / s
         if(m .eq. l) goto 150
         s1 = ABS(p) * (ABS(h(m-1,m-1)) + ABS(zz) + ABS(h(m+1,m+1)))
         s2 = s1 + ABS(h(m,m-1)) * (ABS(q) + ABS(r))
         if(s2 .eq. s1) goto 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0_10
         if(i .eq. mp2) goto 160
         h(i,i-3) = 0.0_10
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if(k .eq. m) goto 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0_10
         if(notlas) r = h(k+2,k-1)
         x = ABS(p) + ABS(q) + ABS(r)
         if(x .eq. 0.0_10) goto 260
         p = p / x
         q = q / x
         r = r / x
  170    s = SIGN(SQRT(p*p+q*q+r*r),p)
         if(k .eq. m) goto 180
         h(k,k-1) = -s * x
         goto 190
  180    if(l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
c     .......... row modification ..........
         do 210 j = k, en
            p = h(k,j) + q * h(k+1,j)
            if(.not. notlas) goto 200
            p = p + r * h(k+2,j)
            h(k+2,j) = h(k+2,j) - p * zz
  200       h(k+1,j) = h(k+1,j) - p * y
            h(k,j) = h(k,j) - p * x
  210    continue
c
         j = MIN(en,k+3)
c     .......... column modification ..........
         do 230 i = l, j
            p = x * h(i,k) + y * h(i,k+1)
            if(.not. notlas) goto 220
            p = p + zz * h(i,k+2)
            h(i,k+2) = h(i,k+2) - p * r
  220       h(i,k+1) = h(i,k+1) - p * q
            h(i,k) = h(i,k) - p
  230    continue
c
  260 continue
c
      goto 70
c     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.0_10
      en = na
      goto 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0_10
      q = p * p + w
      zz = SQRT(ABS(q))
      x = x + t
      if(q .lt. 0.0_10) goto 320
c     .......... real pair ..........
      zz = p + SIGN(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if(zz .ne. 0.0_10) wr(en) = x - w / zz
      wi(na) = 0.0_10
      wi(en) = 0.0_10
      goto 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      goto 60
c     .......... set error -- no convergence to an
c                eigenvalue after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine HQR2(nm,n,low,igh,h,wr,wi,z,ierr)
C***BEGIN PROLOGUE  HQR2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C2B
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes eigenvalues and eigenvectors of real upper
C            Hessenberg matrix using QR method.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure HQR2,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a REAL UPPER Hessenberg matrix by the QR method.  The
C     eigenvectors of a REAL GENERAL matrix can also be found
C     if  ELMHES  and  ELTRAN  or  ORTHES  and  ORTRAN  have
C     been used to reduce this general matrix to Hessenberg form
C     and to accumulate the similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH are integers determined by the balancing
C          subroutine  BALANC.  If  BALANC  has not been used,
C          set LOW=1, IGH=N.
C
C        H contains the upper Hessenberg matrix.
C
C        Z contains the transformation matrix produced by  ELTRAN
C          after the reduction by  ELMHES, or by  ORTRAN  after the
C          reduction by  ORTHES, if performed.  If the eigenvectors
C          of the Hessenberg matrix are desired, Z must contain the
C          identity matrix.
C
C     On OUTPUT
C
C        H has been destroyed.
C
C        WR and WI contain the real and imaginary parts,
C          respectively, of the eigenvalues.  The eigenvalues
C          are unordered except that complex conjugate pairs
C          of values appear consecutively with the eigenvalue
C          having the positive imaginary part first.  If an
C          error exit is made, the eigenvalues should be correct
C          for indices IERR+1,...,N.
C
C        Z contains the real and imaginary parts of the eigenvectors.
C          If the I-th eigenvalue is real, the I-th column of Z
C          contains its eigenvector.  If the I-th eigenvalue is complex
C          with positive imaginary part, the I-th and (I+1)-th
C          columns of Z contain the real and imaginary parts of its
C          eigenvector.  The eigenvectors are unnormalized.  If an
C          error exit is made, none of the eigenvectors has been found.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C
C     Calls CDIV for complex division.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  CDIV
C***END PROLOGUE  HQR2
C
      implicit none
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn
      integer igh,itn,its,low,mp2,enm2,ierr
      real*10 h(nm,n),wr(n),wi(n),z(nm,n)
      real*10 p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,s1,s2
      logical notlas
C
C***FIRST EXECUTABLE STATEMENT  HQR2
      ierr = 0
      norm = 0.0_10
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + ABS(h(i,j))
c
         k = i
         if(i .ge. low .and. i .le. igh) goto 50
         wr(i) = h(i,i)
         wi(i) = 0.0_10
   50 continue
c
      en = igh
      t = 0.0_10
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if(en .lt. low) goto 340
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if(l .eq. low) goto 100
         s = ABS(h(l-1,l-1)) + ABS(h(l,l))
         if(s .eq. 0.0_10) s = norm
         s2 = s + ABS(h(l,l-1))
         if(s2 .eq. s) goto 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if(l .eq. en) goto 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if(l .eq. na) goto 280
      if(itn .eq. 0) goto 1000
      if(its .ne. 10 .and. its .ne. 20) goto 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = ABS(h(en,na)) + ABS(h(na,enm2))
      x = 0.75_10 * s
      y = x
      w = -0.4375_10 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = ABS(p) + ABS(q) + ABS(r)
         p = p / s
         q = q / s
         r = r / s
         if(m .eq. l) goto 150
         s1 = ABS(p) * (ABS(h(m-1,m-1)) + ABS(zz) + ABS(h(m+1,m+1)))
         s2 = s1 + ABS(h(m,m-1)) * (ABS(q) + ABS(r))
         if(s2 .eq. s1) goto 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0_10
         if(i .eq. mp2) goto 160
         h(i,i-3) = 0.0_10
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if(k .eq. m) goto 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0_10
         if(notlas) r = h(k+2,k-1)
         x = ABS(p) + ABS(q) + ABS(r)
         if(x .eq. 0.0_10) goto 260
         p = p / x
         q = q / x
         r = r / x
  170    s = SIGN(SQRT(p*p+q*q+r*r),p)
         if(k .eq. m) goto 180
         h(k,k-1) = -s * x
         goto 190
  180    if(l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
c     .......... row modification ..........
         do 210 j = k, n
            p = h(k,j) + q * h(k+1,j)
            if(.not. notlas) goto 200
            p = p + r * h(k+2,j)
            h(k+2,j) = h(k+2,j) - p * zz
  200       h(k+1,j) = h(k+1,j) - p * y
            h(k,j) = h(k,j) - p * x
  210    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 230 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            if(.not. notlas) goto 220
            p = p + zz * h(i,k+2)
            h(i,k+2) = h(i,k+2) - p * r
  220       h(i,k+1) = h(i,k+1) - p * q
            h(i,k) = h(i,k) - p
  230    continue
c     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            if(.not. notlas) goto 240
            p = p + zz * z(i,k+2)
            z(i,k+2) = z(i,k+2) - p * r
  240       z(i,k+1) = z(i,k+1) - p * q
            z(i,k) = z(i,k) - p
  250    continue
c
  260 continue
c
      goto 70
c     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0_10
      en = na
      goto 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0_10
      q = p * p + w
      zz = SQRT(ABS(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if(q .lt. 0.0_10) goto 320
c     .......... real pair ..........
      zz = p + SIGN(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if(zz .ne. 0.0_10) wr(en) = x - w / zz
      wi(na) = 0.0_10
      wi(en) = 0.0_10
      x = h(en,na)
      s = ABS(x) + ABS(zz)
      p = x / s
      q = zz / s
      r = SQRT(p*p+q*q)
      p = p / r
      q = q / r
c     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
c     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
c     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
c
      goto 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      goto 60
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  340 if(norm .eq. 0.0_10) goto 1001
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if(q) 710, 600, 800
c     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0_10
         if(na .eq. 0) goto 800
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = h(i,en)
            if(m .gt. na) goto 620
c
            do 610 j = m, na
  610       r = r + h(i,j) * h(j,en)
c
  620       if(wi(i) .ge. 0.0_10) goto 630
            zz = w
            s = r
            goto 700
  630       m = i
            if(wi(i) .ne. 0.0_10) goto 640
            t = w
            if(t .ne. 0.0_10) goto 635
            t = norm
  632       t = 0.5_10*t
            if(norm + t .gt. norm) goto 632
            t = 2.0_10*t
  635       h(i,en) = -r / t
            goto 700
c     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if(ABS(x) .le. ABS(zz)) goto 650
            h(i+1,en) = (-r - w * t) / x
            goto 700
  650       h(i+1,en) = (-s - y * t) / zz
  700    continue
c     .......... end real vector ..........
         goto 800
c     .......... complex vector ..........
  710    m = na
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         if(ABS(h(en,na)) .le. ABS(h(na,en))) goto 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         goto 730
  720    call CDIV(0.0_10,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0_10
         h(en,en) = 1.0_10
         enm2 = na - 1
         if(enm2 .eq. 0) goto 800
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 790 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0_10
            sa = h(i,en)
c
            do 760 j = m, na
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
c
            if(wi(i) .ge. 0.0_10) goto 770
            zz = w
            r = ra
            s = sa
            goto 790
  770       m = i
            if(wi(i) .ne. 0.0_10) goto 780
            call CDIV(-ra,-sa,w,q,h(i,na),h(i,en))
            goto 790
c     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0_10 * q
            if(vr .ne. 0.0_10 .or. vi .ne. 0.0_10) goto 783
            s1 = norm * (ABS(w)+ABS(q)+ABS(x)+ABS(y)+ABS(zz))
            vr = s1
  782       vr = 0.5_10*vr
            if(s1 + vr .gt. s1) goto 782
            vr = 2.0_10*vr
  783       call CDIV(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,
     1                h(i,na),h(i,en))
            if(ABS(x) .le. ABS(zz) + ABS(q)) goto 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            goto 790
  785       call CDIV(-r-y*h(i,na),-s-y*h(i,en),zz,q,
     1                h(i+1,na),h(i+1,en))
  790    continue
c     .......... end complex vector ..........
  800 continue
c     .......... end back substitution.
c                vectors of isolated roots ..........
      do 840 i = 1, n
         if(i .ge. low .and. i .le. igh) goto 840
c
         do 820 j = i, n
  820    z(i,j) = h(i,j)
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zz = 0.0_10
c
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
c
            z(i,j) = zz
  880 continue
c
      goto 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine CDIV(ar,ai,br,bi,cr,ci)
C***BEGIN PROLOGUE  CDIV
C***REFER TO  EIS<DOC
C
C     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CDIV
      implicit none
      real*10 ar,ai,br,bi,cr,ci
c
      real*10 s,ars,ais,brs,bis
C***FIRST EXECUTABLE STATEMENT  CDIV
      s = ABS(br) + ABS(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
