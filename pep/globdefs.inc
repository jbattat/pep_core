c PEP global parameter definitions

c universal parameters

c number of real*4 fields to fit in a standard variable
c (2 for r*8 or 4 for r*10 or r*16)
      integer u_stdsz
      parameter (u_stdsz=4)

c number of solar-system parameters
      integer*4 u_nmprm
      parameter (u_nmprm=100)

c number of body-specific parameters, including initial conditions
      integer*4 u_nmbod
      parameter (u_nmbod=30)

c maximum number of input bodies (excluding earth/moon, limited-scope
c asteroids, and pulsars)
      integer u_mxpl
      parameter (u_mxpl=16)

c maximum number of input limited-scope asteroids
      integer u_mxsml
      parameter (u_mxsml=200)

c maximum number of input pulsars
      integer u_mxpsr
      parameter (u_mxpsr=10)

c maximum order of zonal harmonics
      integer u_mxzon
      parameter (u_mxzon=50)

c maximum order of tesseral harmonics
      integer u_mxtes
      parameter (u_mxtes=50)

c maximum number of input observing sites
      integer u_mxsit
      parameter (u_mxsit=100)

c maximum number of input observed spots
      integer u_mxspt
      parameter (u_mxspt=100)

c maximum number of input radar biases
      integer u_mxrbs
      parameter (u_mxrbs=200)

c maximum number of input equator-equinox-latitude corrections
      integer u_mxeqx
      parameter (u_mxeqx=400)

c maximum number of input planetary phase corrections
      integer u_mxphs
      parameter (u_mxphs=200)

c maximum number of input star catalog error models
      integer u_mxsky
      parameter (u_mxsky=10)


c parameters for integrations

c maximum number of differential equations to be integrated simultaneously
c (motion and partials)
      integer i_mxeqn
      parameter (i_mxeqn=1000)

c maximum number of partial derivatives to be integrated for earth/moon/planet
      integer i_mxplprt
      parameter (i_mxplprt=98)

c maximum number of target bodies considered for a given integrated body
      integer i_mxtrg
      parameter (i_mxtrg=11)

c size of index for external indirect partials
c           see tapdtplp and yvectplp includes
      integer*4 i_mxplp
      parameter (i_mxplp=95)
