      subroutine MVC( string1, k1, n, string2, k2 )
 
c----------------------------------------------------------------------
c
c---purpose:     to move n characters beginning at k1 in character
c                variable string1 to string 2 beginning at position k2.
c
c---references:  none
c
c---input
c   arguments:   string1  = character array to move string from
c                           (character)
c                k1       = position of first character in string to be
c                           moved
c                           (integer)
c                n        = number of characters in string to be moved
c                           (integer)
c                k2       = position of first character in string2 where
c                           the string is to be inserted
c                           (integer)
c
c---output
c   arguments:   string2  = character array into which the string was
c                           inserted
c                           (character)
c
c---ver./date/
c   programmer:  v1.0/10-91/jlh (usno/omd)
c
c                11/25/91 mam (sao)  modified
c                06/15/92 jfc - changed to 256 limit
c
c----------------------------------------------------------------------
 
      implicit none
 
      character*(*) string1, string2
      integer  k1, k2, n, m
 
 
c don't move more than 256 bytes
      m = min( n, 256 )
 
c copy input string to output storage
      string2( k2 : k2+m-1) = string1( k1 : k1+m-1 )
 
      return
      end
