!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! parameters for mainfire.f90 - Tropical community                     !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module parafire

implicit none

! COLONIZATION AND MORTALITY RATES
! mortality=3/life_span (from 'boreal-species.csv')
! colonization is set so that the time to achieve the asumptotic state is
! the same as the time of max cover after burning in 'boreal-species.csv'

!for PFT1
!double precision,parameter :: c1  = 0.085d0     ! yr^-1
!double precision,parameter :: m1  = 0.035d0    ! yr^-1

!for PFT2
!double precision,parameter :: c1  = 0.13d0     ! yr^-1
!double precision,parameter :: m1  = 0.0015d0    ! yr^-1

!for PFT3
double precision,parameter :: c1  = 0.18d0     ! yr^-1
double precision,parameter :: m1  = 0.023d0    ! yr^-1

double precision,parameter :: c2 = 0.d0      ! yr^-1
double precision,parameter :: m2  = 0.d0     ! yr^-1
double precision,parameter :: c3  = 0.d0       ! yr^-1 
double precision,parameter :: m3  = 0.d0      ! yr^-1 

! TIME PARAMS
double precision,parameter   :: h   = 1.d0              ! time step of integration (day)
double precision,parameter   :: dtoyr = 1/365.d0        ! day to year converter
double precision,parameter   :: hy = h*dtoyr            ! time step in years (yr)
integer, parameter           :: maxT = 500            ! integration time (yr)
integer,parameter            :: NN = maxT*365           ! integration time steps (day)

end module parafire
