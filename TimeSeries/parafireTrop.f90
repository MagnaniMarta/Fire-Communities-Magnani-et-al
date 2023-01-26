!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! parameters for mainfire.f90 - Tropical community                     !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module parafire

implicit none

! COLONIZATION AND MORTALITY RATES
double precision,parameter :: c1  = 0.20d0     ! yr^-1
double precision,parameter :: m1  = 0.01d0    ! yr^-1
double precision,parameter :: c2 = 0.15d0      ! yr^-1
double precision,parameter :: m2  = 0.06d0     ! yr^-1
double precision,parameter :: c3  = 20d0       ! yr^-1 
double precision,parameter :: m3  = 3d0      ! yr^-1 

! TIME PARAMS
double precision,parameter   :: h   = 1.d0              ! time step of integration (day)
double precision,parameter   :: dtoyr = 1/365.d0        ! day to year converter
double precision,parameter   :: hy = h*dtoyr            ! time step in years (yr)
integer, parameter           :: maxT = 15000            ! integration time (yr)
integer,parameter            :: NN = maxT*365           ! integration time steps (day)
integer,parameter            :: statout = nint(0.8*NN)  !starting time of statistics (bav and tfav) calculation

!FIRE PARAMS
double precision            :: R(3) = (/0.10d0,0.70d0,0.85d0/)  ! fire respoonses
double precision            :: L(3) = (/0.001d0,0.2d0,0.5d0/)  ! flammabilities (yr^-1)
double precision, parameter :: delt(3) = (/0.d0, 0.d0, 0.d0/)   ! minimum values of bi when fire
integer,parameter           :: minfirerettime = 1               ! minimum fire return time (yr)
double precision            :: eps = 0.0001d0                   ! maximum fire return time (yr)

end module parafire
