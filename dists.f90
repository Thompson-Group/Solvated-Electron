!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine to calculate distances to other atoms
!Written by Paul Burris and Ward Thompson
!January 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dists

use common_variables
implicit none

real(kind=dp) :: xtmp, ytmp, ztmp
! working variables
integer :: j, i

do j = 1, n_atoms-1
   do i = j+1, n_atoms
      xtmp = x(j) - anint( (x(j) - x(i))/Lx )*Lx
      ytmp = y(j) - anint( (y(j) - y(i))/Ly )*Ly
      ztmp = z(j) - anint( (z(j) - z(i))/Lz )*Lz

      rx(i,j) = x(i) - xtmp
      ry(i,j) = y(i) - ytmp
      rz(i,j) = z(i) - ztmp

      rx(j,i) = -rx(i,j)
      ry(j,i) = -ry(i,j)
      rz(j,i) = -rz(i,j)

   enddo
enddo

end subroutine dists
