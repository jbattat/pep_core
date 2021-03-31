      integer function ITYPOB(ncodf)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   iundef, ncodg
 
c*** end of declarations inserted by spag
 
 
c        j.f.chandler - 1979 march 1  - function itypob
c        return a code specifying the type of series for ncodf
c        1 - radar, 2 - optic, 3 - trnsit, 4 - fermtr, 5 - stmrdr
c        6 - undefined
c        this function should be used thruout pep where an observation
c        series must be classified.
      integer*2 ncodf
      data iundef/6/
      ncodg = ncodf
      if( ncodf .gt. 20 ) ncodg = ncodf - 20
      ITYPOB = (ncodg + 2)/3
      if( ncodg .le. 0 .or. ncodg .eq. 15 ) ITYPOB = iundef
      if( ncodg .eq. 20 ) ITYPOB = 2
      if( ncodg .eq. 19 ) ITYPOB = 1
      if( ncodg .eq. 18 ) ITYPOB = 1
      return
      end
