module lib_array2

  implicit none
  save

  !public :: locate
  !interface locate
  !   module procedure locate_dp
  !end interface locate

  !public :: interp1d_loglin
  !interface interp1d_loglin
  !   module procedure interp1d_loglin_dp
  !   module procedure interp1d_array_loglin_dp
  !end interface interp1d_loglin


contains

  
  real function interp1d_single_loglin_dp(x1,y1,x2,y2,xval) result(yval)
    real,intent(in) :: x1,y1,x2,y2,xval
    real :: frac
    frac = ( log10(xval) - log10(x1) ) / ( log10(x2) - log10(x1) )
    yval = y1 + frac * ( y2 - y1 )
  end function interp1d_single_loglin_dp


  real function interp1d_loglin_dp(x,y,xval,N,bounds_error,fill_value) result(yval)
    integer,intent(in) :: N
    real,dimension(N),intent(in) :: x,y
    real,intent(in) :: xval
    logical,intent(in),optional :: bounds_error
    real,intent(in),optional :: fill_value

    integer :: ipos
!    ! position of x value in x array

    logical :: bounds_error_tmp
    real :: fill_value_tmp

    if(present(bounds_error)) then
       bounds_error_tmp = bounds_error
    else
       bounds_error_tmp = .true.
    end if

    if(.not.bounds_error_tmp) then
       if(present(fill_value)) then
          fill_value_tmp = fill_value
       else
          fill_value_tmp = 0.
       end if
    end if


    ipos = locate_dp(x,xval,N)

!    ! --- First some error checking --- !

    if(ipos == -1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') xval,x(1),x(n)
          stop
       else
          yval = fill_value_tmp
          return
       end if
    end if

    if( ipos < n .and. ipos > 0) then
       yval = interp1d_single_loglin_dp(x(ipos), y(ipos), x(ipos+1), y(ipos+1), xval)
    else if(ipos == n) then
       yval = y(n)
    else if(ipos == 0) then
       yval = y(1)
    else
       write(0,'("ERROR: Unexpected value of ipos : ",I0)') ipos
       stop
    end if

  end function interp1d_loglin_dp

  function interp1d_array_loglin_dp(x,y,xval,N,bounds_error,fill_value)
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: x(N),y(N)
    real,intent(in) :: xval(N)
    logical,intent(in),optional :: bounds_error
    real,intent(in),optional :: fill_value
    real :: interp1d_array_loglin_dp(size(xval))
    integer :: i
    do i=1,size(xval)
       interp1d_array_loglin_dp(i) = interp1d_loglin_dp(x,y,xval(i),N,bounds_error,fill_value)
    end do
  end function interp1d_array_loglin_dp


    integer function locate_dp(xx,x,N)
!    ! Locate a value in a sorted array

    implicit none
    integer, intent(in) :: N
    real, dimension(N), intent(in) :: xx
    real, intent(in) :: x
    integer :: jl,jm,ju
    logical :: ascnd

    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate_dp = 1
    else if (x == xx(n)) then
       locate_dp = n-1
    else if(ascnd.and.x > xx(n)) then
       locate_dp = n
    else if (ascnd.and. x < xx(1)) then
       locate_dp = 0
!    else if(ascnd.and. (x > xx(n) .or. x < xx(1))) then
!       locate_dp = -1
    else if(.not.ascnd.and. (x < xx(n) .or. x > xx(1))) then
       locate_dp = -1
    else
       locate_dp = jl
    end if

  end function locate_dp

end module lib_array2
