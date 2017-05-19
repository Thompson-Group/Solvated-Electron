    subroutine qm_dump(i)

      use common_variables
      use quantum_variables
      use constants
      implicit none

   integer(kind=ip), intent(in) :: i
   integer(kind=ip) :: j

   write(neigs, '(21F12.5)') real(i)*dt, (Eigval(j)*evperau,j=1,10)
   write(nrg, '(2F12.5)') real(i)*dt, sqrt(r2_e_avg)*angperau
   write(nlanc, '(F12.5,I4)') real(i)*dt, actiter

end subroutine qm_dump
