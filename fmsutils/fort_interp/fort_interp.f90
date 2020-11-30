module fort_interp_mod

  use lib_array2, only : interp1d_array_loglin_dp
  
  implicit none

contains

  subroutine interp_data(nt, npx, npy, npz, data, p, newp,output)
    integer, intent(in)    :: nt, npx, npy, npz
    real,    intent(in) :: data(nt,npz,npy,npx)
    real,    intent(out) :: output(nt,npz,npy,npx)
    real,    intent(in)    ::    p(nt,npz,npy,npx)
    real,    intent(in)    :: newp(npz)

    integer :: t,x,y

    do t=1,nt
       do y=1,npy
          do x=1,npx
             output(t,:,y,x) = interp1d_array_loglin_dp(p(t,:,y,x), data(t,:,y,x), newp,npz)
          end do
       end do
    end do
    
  end subroutine interp_data
  
  subroutine height(temp,ph,nt,npz,npy,npx,R,g,h)
    integer, intent(in) :: nt,npz,npy,npx
    real,    intent(in) :: temp(nt,npz,npy,npx)
    real,    intent(in) :: ph(nt,npz+1,npy,npx)
    real,    intent(in) :: R,g
    real,    intent(out):: h(nt,npz+1,npy,npx)

    real :: delz
    integer :: t,z,y,x
    
    h = 0
    do t=1,nt
       do z=npz,1,-1
          do y=1,npy
             do x=1,npx
                delz = R/g * temp(t,z,y,x) * (log(ph(t,z+1,y,x)) - log(ph(t,z,y,x)))
                h(t,z,y,x) = h(t,z+1,y,x) + delz
             end do
          end do
       end do
    end do  
  end subroutine height
  
end module fort_interp_mod
