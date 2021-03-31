      double precision function random(iseed)
      implicit none
      integer iseed

c Return uniformly distributed random variable in the interval 0 to 1.

      iseed = iseed*117649 + 7
      random = iseed*2.328306437D-10
      if(random.lt.0.) random = random + 1.
      return
      end
