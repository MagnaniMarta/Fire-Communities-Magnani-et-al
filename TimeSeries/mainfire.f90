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
! program output : 1) time series of vegetation cover and expected fire return  !
!                  time ('tbTf.dat')                                            !
!                  2) statistics over last 20% simulation time ('stats.dati')   !
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
double precision,dimension(3) ::  bav                        ! average cover
double precision              ::  tf,tfdum                   ! expected fire return time
double precision              ::  tfav                       ! algebric average fire return time
integer                       ::  lastF,numF                 ! time last fire and number of fires in the simulation 
integer                       ::  idummy                     ! variable for random number generator
real                          ::  dummy, ran3                ! output and function of random number generator 
integer                       ::  i                          ! time counter
double precision              ::  firevf                     ! fire free days counter
integer                       ::  iifire                     ! fire flag 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!open output files
open (21,file='tbTf.dat')
open (22,file='stats.dat')
    
!initialization parameters
firevf=NN
iifire=0
idummy=12
lastF=0
numF=0

!initial covers: random covers between (0:1] 
!renormalized to have total cover lower than 1 
b(1)=ran3(idummy)
b(2)=ran3(idummy)
b(3)=ran3(idummy)
b=b/sum(b)
bav=0.d0*b
    
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
!  if (dummy<=hy/tf .and. int(firevf)>=minfirerettime) then  !alternative Bernoulli condition

    iifire=1 !fire
    firevf=0.d0

    if (i>=statout) then
      tfav=tfav+0.001*real(i-lastF)*dtoyr
      lastF=i
      numF=numF+1
      bav=bav+0.001d0*b
    end if
          
  else  
     iifire=0 !no fire
     firevf=firevf+hy
          
    end if

! vegetation update
  call fireocc(f,iifire,fout)
  b=fout

  if (mod(i,100)==0 .or. iifire>0) then
    write (21,'(1i15,4f15.4)') i,b,tf
  end if

end do

!compute statistics
if(numF>0)then !fires occurred
        bav=bav*1000/real(numF)
        tfav=tfav*1000/real(numF)
else !no fire occurred
        bav=b
        tfav=0.d0
end if

write(22,'(4f15.4)') bav,tfav

close(21)
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
