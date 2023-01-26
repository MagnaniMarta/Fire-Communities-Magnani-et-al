!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! parameters for mainfire.f90  - Mediterranean community               !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module parafire

implicit none

! COLONIZATION AND MORTALITY RATES
double precision,parameter :: c1  = 0.047d0     ! yr^-1
double precision,parameter :: m1  = 0.0025d0    ! yr^-1
!double precision,parameter :: c2 = 0.053d0      ! yr^-1
double precision,parameter :: m2  = 0.008d0     ! yr^-1
double precision,parameter :: c3  = 0.3d0       ! yr^-1 
double precision,parameter :: m3  = 0.03d0      ! yr^-1 
double precision           :: c2                ! yr^-1

! TIME PARAMS
double precision,parameter   :: h   = 1.d0              ! time step of integration (day)
double precision,parameter   :: dtoyr = 1/365.d0        ! day to year converter
double precision,parameter   :: hy = h*dtoyr            ! time step in years (yr)
integer, parameter           :: maxT = 15000            ! integration time (yr)
integer,parameter            :: NN = maxT*365           ! integration time steps (day)
integer,parameter            :: statout = nint(0.8*NN)  !starting time of statistics (bav and tfav) calculation

!FIRE PARAMS
!double precision,parametr  :: R(3) = (/0.85d0,0.40d0,0.50d0/)  ! fire respoonses
integer                     :: r1                               ! fire respoonse PFT1
integer,parameter           :: r2=40, r3=50                     ! fire respoonse PFT2,3
double precision            :: R(3)                             ! fire respoonses
double precision,parameter  :: L(3) = (/0.002d0,0.05d0,0.1d0/)  ! flammabilities (yr^-1)
double precision, parameter :: delt(3) = (/0.d0, 0.d0, 0.d0/)   ! minimum values of bi when fire
integer,parameter           :: minfirerettime = 2               ! minimum fire return time (yr)
double precision,parameter  :: eps = 0.0001d0                   ! maximum fire return time (yr)

end module parafire
