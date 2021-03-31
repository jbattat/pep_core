      subroutine FILTPN(in0, nstop, init)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   in0, nprd, nstop
 
c*** end of declarations inserted by spag
 
 
      logical   init
c
c     z. goldberg march 1980 subroutine filtpn
c                            (from filtin 7/78)
c
c         subroutine to read in filter process noise parameters.
c         invoked by '*markov' in input stream.
c
c         the names of the process noise parameters follow '*markov'.
c         one parameter name goes on each card and the names must be in
c         the same order as the corresponding set of names input to the
c         smearing matrix program.  most parameters have two-part names;
c         these are described in subroutine namprm.  the names can be
c         placed on the name cards in free field format, at least one
c         blank separating the first part from the second.  the name
c         cards may be ended by a blank card or the next * card.
c
c         common
      include 'filnit.inc'
      include 'filtim.inc'
      include 'inodta.inc'
 
      if(  .not. init ) then
c
c spool process noise parameter names until 1 blank card
         call PEPTIC(In, Iout, in0, 8,
     .               'PROCESS NOISE PARAMETER NAMES   ', nstop, 1)
c
c error if &nmlst3 not yet read
         if( Nml3rd ) then
c
c read names ignoring blank records
            call NAMRED(in0, Pnames, Npnp, nprd)
            Pnpnrd = .true.
c
c check for mismatched number of parameters
            if( nprd .ne. Npnp ) then
               write(Iout, 10)
   10          format(82x, 'WRONG NUMBER OF PARAMETER NAMES ***')
               nstop = nstop + 1
            endif
 
            return
         else
            write(Iout, 20)
   20       format(
     .        ' *** PROCESS NOISE PARAMETER NAMES INPUT BEFORE &NMLST3'
     .        , ' AT LABEL=30 IN FILTPN ***')
            nstop = nstop + 1
            return
         endif
      else
c
c no process noise parm names input
         write(Iout, 50)
   50    format(
     .       ' *** NO PROCESS NOISE PARAMETER NAMES HAVE BEEN INPUT --'
     .       , ' ERROR AT LABEL=10 IN FILTPN ***')
         nstop = nstop + 1
         return
      endif
      end
