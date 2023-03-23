!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!            PROGRAM TO MODEL PLANT DYNAMICS, NO FIRE                  !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program mainconst
use parafire

implicit none

double precision,dimension(3) ::  b,f,fout
integer           ::  i
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! CODE STARTS

    open (21,file='tb.dat')

    !     variables initialization
    b=(/0.01,0.0,0.0/)

    !     temporal loop
    do i=1,NN

          ! integration deterministic succession
         f=b !b(t)
        call rk4(f,fout) !subroutine: fourth order Runge-Kutta scheme (see line 138)
        b=fout !b(t+dt)

          if (mod(i,100).eq.0) then
            write (21,'(1i15,3f15.4)') i,b
          end if
    end do


end program mainconst


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
