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
end module fort_interp_mod
