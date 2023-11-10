module rk_force2
  implicit none
  !--------------------------------------------------------------------
  !Declaration of private and public variables,functions and subroutines:
  private :: dp,r
  integer, parameter  :: NEQ=4   !This is a public variable that determines the number of equations that our ode system has
  public  :: NEQ,k1,k2,k3,k4     !These are public variables. Also the k1,k2,k3 and k4 are the constants for our 4 odes problem
  public  :: f,energy            !We have public also the functions f and energy
  public  :: rkstep              !We have public the subroutine rkstep
  integer, parameter  :: dp=kind(0.0_8) !We have double precision and also private parameter
  real(dp) :: k1,k2,k3,k4
  real(dp) :: r
  !---------------------------------------------------------------------
contains
  !---------------------------------------------------------------------
  !We have a function that takes as an input the time t and also the array X of size NEQ and gives us the result dXdt array
  function f(t,X)    result(dXdt)
    !-------------------------------------------------------------------
    !Declaration of local variables:
    real(dp)  :: t,X(:)             !We have as an input the time t and the array X, which we know it's size so we don't need the dimension statement
    real(dp)  :: dXdt(sizeof(X))    !We have that the result is going to have the same size with the input array X
    real(dp)  :: x1,x2,v1,v2,a1,a2  !Some useful local variables
    !-------------------------------------------------------------------
    !Calculations:
    x1=X(1)
    x2=X(2)
    v1=X(3)
    v2=X(4) !We take the useful values of our array
    r=(x1*x1 + x2*x2)**(-1.5_dp)
    a1= k1*x1*r
    a2= k1*x2*r
    dXdt(1)=v1
    dXdt(2)=v2
    dXdt(3)=a1
    dXdt(4)=a2  !We give the solutions of the ode system to the result array
  end function f
  !---------------------------------------------------------------------
  real(dp) function  energy(t,X)
    !-------------------------------------------------------------------
    !Declaration of local variables:
    real(dp)  :: t,X(:)
    real(dp)  :: x1,x2,v1,v2
    !-------------------------------------------------------------------
    !Calculations:
    x1=X(1)
    x2=X(2)
    v1=X(3)
    v2=X(4) !We take the useful values of our array
    r=sqrt(x1*x1 + x2*x2)
    energy=(0.5_dp)*(v1*v1 + v2*v2) + k1/r
  end function energy
  !---------------------------------------------------------------------
  subroutine rkstep(t,X,dt)
    !-------------------------------------------------------------------
    !Declaration of local variables:
    real(dp),dimension(NEQ)   :: X    !We declare here the size because we have a subroutine
    real(dp)                  :: t,dt
    real(dp),dimension(NEQ)   :: k1,k2,k3,k4
    real(dp)                  :: h,h2,h6
    !-------------------------------------------------------------------
    !Calculations:
    h=dt
    h2=h/(2.0_dp)
    h6=h/(6.0_dp)
    k1=f(t,X)
    k2=f(t+h2,X+(h2*k1))
    k3=f(t+h2,X+(h2*k2))
    k4=f(t+h,X+(h*k3))
    t=t+h
    X=X + h6*(k1+ (2.0_dp * k2)+ (2.0_dp *k3) +k4)
  end subroutine rkstep
  !---------------------------------------------------------------------
end module rk_force2
!=======================================================================
