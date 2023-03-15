!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! parameters for mainfire.f90 - Tropical community                     !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module parafire

implicit none

! COLONIZATION AND MORTALITY RATES
double precision,parameter :: c1  = 0.085d0    ! yr^-1
double precision,parameter :: m1  = 0.035d0    ! yr^-1
double precision,parameter :: m2  = 0.015d0    ! yr^-1
double precision,parameter :: c3  = 0.17d0     ! yr^-1 
double precision,parameter :: m3  = 0.023d0    ! yr^-1 
double precision           :: c2               ! yr^-1

! TIME PARAMS
double precision,parameter   :: h   = 1.d0              ! time step of integration (day)
double precision,parameter   :: dtoyr = 1/365.d0        ! day to year converter
double precision,parameter   :: hy = h*dtoyr            ! time step in years (yr)
integer, parameter           :: maxT = 15000            ! integration time (yr)
integer,parameter            :: NN = maxT*365           ! integration time steps (day)
integer,parameter            :: statout = nint(0.8*NN)  !starting time of statistics (bav and tfav) calculation

!FIRE PARAMS
integer                     :: r1                               ! fire respoonse PFT1
integer,parameter           :: r2=55,r3=85                      ! fire response PFT2
double precision            :: R(3)                             ! fire respoonses
double precision,parameter  :: L(3) = (/0.004d0,0.0133d0,0.01d0/)  ! flammabilities (yr^-1)
double precision, parameter :: delt(3) = (/0.d0, 0.d0, 0.d0/)   ! minimum values of bi when fire
integer,parameter           :: minfirerettime = 2               ! minimum fire return time (yr)
double precision,parameter  :: eps = 0.0001d0                   ! maximum fire return time (yr)

end module parafire
