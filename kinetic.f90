subroutine kinetic
!--------------------------------------------------------------------------
!     Program to compute and store the one-dimensional kinetic energy
!       matrices using the -infinity to infinity form of Sinc-function DVR
!--------------------------------------------------------------------------

use common_variables
use quantum_variables
use constants
implicit none

integer(kind=ip) :: idtmp, i , j
real(kind=dp) :: dtmp, pref, del_tmp

del_tmp = del/angperau

pref = 1.0_dp/(2.0_dp*me*del_tmp**2)

kex = 0.0_dp
key = 0.0_dp
kez = 0.0_dp
do i = 1, nraw
    kex(i,i) = pref*pi**2/3.0_dp
    key(i,i) = pref*pi**2/3.0_dp
    kez(i,i) = pref*pi**2/3.0_dp
enddo

do i = 1, nraw
    do j = 1, i - 1
        idtmp = i - j
        dtmp = real(idtmp)
        kex(i,j) = (-1.0_dp)**idtmp*pref*2.0_dp/dtmp**2
        kex(j,i) = kex(i,j)
        
        key(i,j) = (-1.0_dp)**idtmp*pref*2.0_dp/dtmp**2
        key(j,i) = key(i,j)

        kez(i,j) = (-1.0_dp)**idtmp*pref*2.0_dp/dtmp**2
        kez(j,i) = kez(i,j)
    enddo
enddo

end subroutine

