subroutine dump(i)

use common_variables
use constants
implicit none

integer(kind=ip),intent(in) :: i
integer(kind=ip) :: j

	write(nxyz,*) n_atoms
	write(nxyz,*) 'Atoms. Timestep: ', i
	do j=1,n_atoms
		write(nxyz,'(I4.1,3F12.5)') elements(a_type(j)), x(j), y(j), z(j)
	end do
    if (run_style .eq. 'qm_nvt' .or 'qm_nve') then
        wrote(nxyz,'(I4.1,3F12.5)') 'I', r_e_avg(1), r_e_avg(2), r_e_avg(3)
    endif

end subroutine dump
