subroutine kinetic(nraw)
!--------------------------------------------------------------------------
!     Program to compute and store the one-dimensional kinetic energy
!       matrices using the -infinity to infinity form of Sinc-function DVR
!--------------------------------------------------------------------------

use common_variables
use quantum_variables
implicit none

real(kind=ip) :: idtmp, i , j
real(kind=dp) :: del,dtmp, pref

pref = 1d0/(2d0*me*del**2)

ke_x = 0d0
ke_y = 0d0
ke_z = 0d0
do i = 1, nraw
    ke_x(i,i) = pref*pi**2/3d0
    ke_y(i,i) = pref*pi**2/3d0
    ke_z(i,i) = pref*pi**2/3d0
enddo

do i = 1, nraw
    do j = 1, i - 1
        idtmp = i - j
        dtmp = dble(idtmp)
        ke_x(i,j) = (-1d0)**idtmp*pref*2d0/dtmp**2
        ke_x(j,i) = ke_x(i,j)
        
        ke_y(i,j) = (-1d0)**idtmp*pref*2d0/dtmp**2
        ke_y(j,i) = ke_y(i,j)

        ke_z(i,j) = (-1d0)**idtmp*pref*2d0/dtmp**2
        ke_z(j,i) = ke_z(i,j)
    enddo
enddo

end subroutine

