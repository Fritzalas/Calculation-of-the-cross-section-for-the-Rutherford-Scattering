program planets
  use rk_force2   !We are going to use our module in this problem
  implicit none
  !-------------------------------------------------------------------
  !Declaration of variables:
  integer, parameter       :: dp=kind(0.0_8)  !We are going to have double precision to our problem
  real(dp),allocatable     :: t(:),x(:,:)     !We are going to have a matrix for x to store for every dimension the required data
  real(dp),dimension(NEQ)  :: x0              !We have initial value as many as the number of equations
  real(dp)                 :: ti,tf           !Initial and final time
  integer                  :: Nt              !The number of places between ti and tf
  integer                  :: i,j             !Useful integers for do loops
  real(dp)                 :: dum1,dum2       !Useful dummy indexes during our program for k1 and k2 (we have k1 and k2 for Runge-Kutta orders)
  real(dp)                 :: bmin,bmax,bstep !We have that b is going to take discrete values from bmin to bmax with step bstep
  real(dp)                 :: b               !The b every time
  real(dp)                 :: theta,thetaold  !The current angle theta and the previous angle
  real(dp)                 :: sigma           !The cross section that we calculate
  real(dp),allocatable     :: Q1(:),Q2(:,:)   !Dummy arrays for the first calculation
  real(dp)                 :: sigmat          !Theoretical value of s
  !-------------------------------------------------------------------
  !User Interface:
  print *, '#We are going to have 4th order Runge-Kutta method for number of ODEs:'
  print *, '#NEQ=',NEQ
  print *, '#Enter k1 and k2 for our problem constants:'
  read  *, k1,k2
  print *, '#k1,k2=',k1,k2
  print *, '#Enter the number of spaces between ti and tf:'
  read  *, Nt
  print *, '#Enter the initial time ti and the final time tf:'
  read  *, ti,tf
  print *, '#Enter the initial values of our problem:'
  read  *, x0  !It asks for values as many as the size of the array
  print *, '#Nt,ti,tf,x0=',Nt,ti,tf,x0
  print *, '#Give the values for the range of b, bmin and bmax:'
  read  *, bmin,bmax
  print *, '#Give also the step for b, bstep:'
  read  *, bstep
  print *, '#bmin,bmax,bstep=',bmin,bmax,bstep
  !------------------------------------------------------------------
  !Array Allocation:
  ALLOCATE(t(Nt))
  ALLOCATE(x(Nt,NEQ))
  ALLOCATE(Q1(Nt))
  ALLOCATE(Q2(Nt,NEQ))
  !------------------------------------------------------------------
  !Calculations:
  dum1=k1
  dum2=k2 !We keep the values of k1 and k2 to dummy indexes
  call rk  !We don't need to give input to our subroutine because we are going to have contains statement in a few lines
  Q1=t
  Q2=x
  !-------------------------------------------------------------------
  !Now we want to calculate the cross section of the Rutheford problem:
  open(unit=12,file='rutherford.dat')
  b=bmin  !We give to the distance b the initial value
  x0(2)=b !Despite the initial value we have the problem (x0,b) with initial position
  call rk !We call rk for our new calculations with the given velocity but with new initial position
  theta=atan(x(Nt,4),x(Nt,3)) !We calculate the theta angle
  thetaold=theta !We keep the previous theta
  b=b+bstep      !We go to the next value of b
  do while (b<=bmax)
     x0(2)=b
     call rk
     theta=atan(x(Nt,4),x(Nt,3))
     sigma =b/(sin(theta)*ABS((theta-thetaold)/bstep))
     sigmat=((k1**2)/(4.0_dp))*(x0(3)**(-4.0_dp))*(sin(0.5_dp * theta)**(-4.0_dp))
     write(12,*) b,theta,sigma,sigmat
     thetaold=theta
     b=b+bstep
  end do
  !-------------------------------------------------------------------
  !Output file:
  k1=dum1
  k2=dum2   !Now we have the correct k1 and k2 for the calculation of the energy
  open(unit=11,file='motion_energy.dat')
  do i=1,Nt
     write(11,*) Q1(i),Q2(i,:),energy(Q1(i),Q2(i,:))  !We have the ith line and all the collumns over there
  end do
  !-------------------------------------------------------------------
  close(11)
  close(12)
  !-------------------------------------------------------------------
contains
  !-------------------------------------------------------------------
  subroutine rk
    !-----------------------------------------------------------------
    !Declaration of local values:
    real(dp)   :: dt,ts,xs(NEQ) !We having the dt about space, ts that holds the values of t and also the vector of size NEQ to hold the values of the quantities because of dimensions
    integer    :: i,j           !Useful integers for do loops
    !-----------------------------------------------------------------
    !Calculations:
    dt=(tf-ti)/(Nt-1)
    t(1)=ti
    x(1,:)=x0 !We have our initial values for our problem
    ts=ti
    xs=x0     !We keep these values to make the next calculation with Runge-Kutta
    do i=2,Nt
       call rkstep(ts,xs,dt)
       t(i)=ts
       x(i,:)=xs
    end do
    !----------------------------------------------------------------
  end subroutine rk
  !-------------------------------------------------------------------
end program planets
!=====================================================================
