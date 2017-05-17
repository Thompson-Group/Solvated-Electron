subroutine wrapper
!************************************************************************************************
!this subroutine calls other subroutine that computes the forces acting on atoms in the system
!Written by Mesele and Pubudu
!hackathon on Thursday Jan 12, 2017
!*************************************************************************************************
 use kinds
 use common_variables
  
 implicit none
  integer :: i 

    call dists
    call calculate_bond
    call calc_angle
    call calc_vdw
    call coulomb_dsf
     
     fx_tot = 0d0; fy_tot = 0d0; fz_tot = 0d0     
	
	do i = 1, n_atoms
	        fx_tot(i) = fx_a(i) + fx_v(i) + fx_c(i) + fx_b(i)
      		fy_tot(i) = fy_a(i) + fy_v(i) + fy_c(i) + fy_b(i)
       		fz_tot(i) = fz_a(i) + fz_v(i) + fz_c(i) + fz_b(i)
        enddo

    v_tot = v_a + v_b + v_v + v_c

end subroutine  wrapper
