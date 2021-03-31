      subroutine FILHED(np, npnp, nepoch)
 
      implicit none
 
c
c d. white  subroutine filhed  september 1974
c
c read outsne heading
c n. b.  assuming np=nparam for now
c
c parameters
      integer*4 np, npnp, nepoch
c
c common
      include 'filtds.inc'
      include 'inodta.inc'
      include 'rtside.inc'
c
c local
      character*80 title
c
c page heading
      call PAGSET('READ FILTER SNE ', 4)
      call NEWPG
c
c first record
      read(Outsne) np, nepoch, npnp
c
c titles
      read(Outsne) title
      write(Iout, 100) title
  100 format('-INPUT SNE:             TITLE= ', A80)
      read(Outsne) title
      write(Iout, 200) title
  200 format('-INPUT FILTER DATA SET: TITLE= ', A80)
c
c error measurements of o-c
      read(Outsne) Measmt, Ermeas
c
c parameter names
      read(Outsne)
c nominals
c
c process noise pointer
      read(Outsne)
c
c list of output epochs
      return
      end
