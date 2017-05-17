subroutine dump(i)

use common_variables
use constants
implicit none

integer(kind=ip),intent(in) :: i
integer(kind=ip) :: j

	write(nxyz,*) n_atoms
	write(nxyz,*) 'Atoms. Timestep: ', i
	do j=1,n_atoms
		write(nxyz,'(I4.1,3F12.5)') a_type(j), x(j), y(j), z(j)
	end do

end subroutine dump
