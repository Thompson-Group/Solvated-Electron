subroutine qm_dump(i)
   use common_variables
   use constants
   implicit none

   integer(kind=ip), intent(in) :: i
   integer(kind=ip) :: j

   write(nthermo, '(I4.1,3F12.5)') i, (Eigval(j),j=1,5)

end subroutine qm_dump
