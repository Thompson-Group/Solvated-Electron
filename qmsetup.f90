!***********************************************************************************
!     Subroutine to set up the quantum calculation.  
!
!  Written by A. Katiyar, P. Wimalasiri, and O. Mesele
!***********************************************************************************
subroutine qmsetup

use kinds
use common_variables
use quantum_variables
use ham_data
use rank3d
use constants
use grid_index

implicit none
real(kind=dp) , dimension(ng), intent(out) :: try
call grid(re)

allocate(Eigvec(niter,niter))
allocate(Eigval(niter))
allocate(krylov-vectors(0:niter,ng))
allocate(indxx(ng))

call kinetic
call rankyx
call rankzx
call looplims
call graphgrid
try = 1.0_dp/(sqrt(real(ng)))

call planczos(try)
!call direct_diag
end subroutine
