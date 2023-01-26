!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                      PROGRAM TO MODEL FIRE ECOSYSTEMS                         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! THIS PROGRAMS MODELS THE DYNAMICS OF COMMUNITIES IN FIRE PRONE ENVIRONMENTS.  !
! IT INCLUDES 3 VEGETATION TYPES THAT HAVE DIFFERET FLAMMABILITIES (L) AND      !
! FIRE RESPONSES (R). DETERMINISTIC SUCCESSION (TILMAN,ECOLOGY 1994) PERTURBED  !
! BY STOCHASTIC FIRES, WITH AVERAGE RETURN TIME Tf. FIRE EVENTS REPRESENTED AS  !
! NONSTATIONARY POISSON PROCESS. FIRE DYNAMICS IMPLEMENTED USING BERNOULLI      !
! CONDITION.                                                                    !
!                                                                               !
! integration time step, dt : 1 day                                             !
! program input : 'parafire.f90', different for each biome                      !
! program output : statistics over last 20% simulation time ('stats.dat')       !
!                                                                               !
! authors: Magnani M., Baudena M.                                               !
! temporary reference: Fire responses shape plant communities in a minimal      !
!                      model for fire ecosystems worldwide. Magnani M.,         !
!                      Diaz-Sierra R., Sweneey L., Provenzale A., Baudena M.    !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program mainfire
use parafire

implicit none

double precision,dimension(3) ::  b                          ! plant covers
double precision,dimension(3) ::  f,fout                     ! variables for plant cover calculation
double precision              ::  tf,tfdum                   ! expected fire return time
double precision              ::  tfav                       ! algebric average fire return time
integer                       ::  lastF,numF                 ! time last fire and number of fires in the simulation 
integer                       ::  idummy                     ! variable for random number generator
real                          ::  dummy, ran3                ! output and function of random number generator 
integer                       ::  i                          ! time counter
double precision              ::  firevf                     ! fire free days counter
integer                       ::  iifire                     ! fire flag 
integer                       ::  ic                         ! iteration c2 parameter
integer                       ::  stfreq(8),k                ! frequency of states
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!open output file: iteration parameter values and frequency of each state before fire
open (22,file='stats.dat')
    
!initialization parameters
firevf=NN
idummy=12

!iteration over c2
do ic=235,949,17
!ic=235,949,17  ! for Med community
!ic=750,3000,56 ! for Trop community
!ic= 650,2650,50   ! for Bor community

 c2=0.0001d0*ic

 !iteration over r1
 do r1=1,90,2
 R=0.01d0*(/r1,r2,r3/)

 !initial covers: random covers between (0:1] 
 !renormalized to have total cover lower than 1 
 b(1)=ran3(idummy)
 b(2)=ran3(idummy)
 b(3)=ran3(idummy)
 b=b/sum(b)

 stfreq=0
 numF=0
 lastF=0
 iifire=0

  ! time loop
  do i=1,NN      
          
    ! integration deterministic succession
    f=b 
    call rk4(f,fout)
    f=fout
          
    ! stochastic fire dynamics:
    ! average retur time Tf(L_i,b_i)
    ! Bernoulli condition for Poisson events
    ! with occurrence probability P=dt/Tf

    dummy=ran3(idummy)             
    tf=1./(L(1)*b(1)+L(2)*b(2)+L(3)*b(3)+eps)
    tfdum=nint(dummy*365*tf)

    if (tfdum==nint(tf)*365-1 .and. int(firevf)>=minfirerettime) then
!    if (dummy<=hy/tf .and. int(firevf)>=minfirerettime) then  !alternative Bernoulli condition

      iifire=1 !fire
      firevf=0.d0

      if (i>=statout) then
        tfav=tfav+0.001*real(i-lastF)*dtoyr
        lastF=i
        numF=numF+1
        call decide(real(b),k)
        stfreq(k)=stfreq(k)+1
      end if
          
    else  
       iifire=0 !no fire
       firevf=firevf+hy
          
      end if

    !vegetation update
    call fireocc(f,iifire,fout)
    b=fout

  end do

  if(numF==0)then
    call decide(real(b),k)
    stfreq(k)=stfreq(k)+1
  end if

  write(22,'(11f15.4)')c2,R(1),real(stfreq),real(numF)

 end do
end do

close(22)


end program mainfire


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

    subroutine rk4(y,yout)
    
   ! THIS SUBROUTINE INTEGRATES THE DIFFERENTIAL EQUATIONS USING 
   ! RUNGE-KUPTA 4 INTEGRATION SCHEME
    
    use parafire

    implicit none

    double precision              :: h6
    double precision,dimension(3) :: y,yout
    double precision,dimension(3) :: k1,k2,k3,k4,y1,y2,y3

    h6=hy/6.d0

    call derivs(y,k1)
    y1=y+0.5*k1*hy
    call derivs(y1,k2)
    y2=y+0.5*k2*hy
    call derivs(y2,k3)
    y3=y+k3*hy
    call derivs(y3,k4)

    yout=y+h6*(k1+2*k2+2*k3+k4)

    return
    end


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

     subroutine derivs(v,dv)
     
     !SUBROUTINE CALCULATING THE R.H.S. OF MODEL DIFFERENTIAL EQUATIONS
     !v = vegetation cover, dv = derivatives
        
     use parafire

     implicit none
     double precision,dimension(3)  :: v,dv

        
     dv(1)=c1*v(1)*(1-v(1))-m1*v(1)
     dv(2)=c2*v(2)*(1-v(1)-v(2))-m2*v(2)-c1*v(1)*v(2)
     dv(3)=c3*v(3)*(1-v(1)-v(2)-v(3))-m3*v(3)-c1*v(1)*v(3)-c2*v(2)*v(3)

        return
     end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

     subroutine decide(y,i)

     !SUBROUTINE DECIDING THE TYPE OF STATE TO BE STORED IN stfreq()
             
     implicit none

     real,dimension(3):: y
     integer ::i
     real, parameter:: delta=0.03

     i=8

     if(y(1)>delta)then !combinations with b1

      if(y(2)>delta)then

       if(y(3)>delta)then
        i=7     !b1+b2+b3
       else
        i=4  !b1+b2
       end if

      else if(y(2)<=delta .and. y(3)>delta)then
       i=5 !b1+b3
      else
       i=1 !b1
      end if

     end if

     if(i==8 .and. y(2)>delta)then !combinations with b2

       if(y(3)>delta)then
        i=6 !b2+b3
       else
        i=2 !b2
       end if

     end if

     if(i==8 .and. y(3)>delta) then
      i=3 !b3
     end if

     return
     end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      subroutine fireocc(b,iifire,bout)
      !THIS SUBROUTINE UPDATE VEGETATION COVER AFTER THE FIRE CONDITION 
      ! IF iifire=1 THEN A FIRE OCCURRED AND b_i IS SET TO R_i*b_i 
      ! OTHERWISE THE COVER IS NOT CHANGED
      ! WHEN iifire=1 AND R_i*b_i IS LOWER THAN A MINIMUM COVER 'delt_i' (from parafire.f90)
      ! THE COVER IS SET TO delt_i

      use parafire
      
      implicit none
      
      double precision,dimension(3)       :: b,bout
      integer                             :: iifire


      bout=b
      bout(1)=b(1)*(1-iifire)+max(b(1)*R(1),delt(1))*iifire
      bout(2)=b(2)*(1-iifire)+max(b(2)*R(2),delt(2))*iifire
      bout(3)=b(3)*(1-iifire)+max(b(3)*R(3),delt(3))*iifire  



      end subroutine
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


      FUNCTION RAN3(IDUM)
      ! RANDOM NUMBER GENERATOR - FROM NUMERICAL RECIPES IN FORTRAN
      
      SAVE
      PARAMETER (MBIG=1000000000,MSEED=161843398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
