subroutine forces
!************************************************************************************************
!this subroutine calls other subroutine that computes the forces acting on atoms in the system
!Written by Mesele and Pubudu
!hackathon on Thursday Jan 12, 2017
!*************************************************************************************************
 use kinds
 use common_variables
 use timings
  
 implicit none
  integer :: i 

  !Timing variables
  real(kind=dp) :: tinit,tfinal
  
  call cpu_time(tinit)

    call dists
    call calculate_bond
    call calc_angle
    call calc_vdw
    if(coul_flag.eq.0) then
       call coulomb
    elseif(coul_flag.eq.1) then
       call coulomb_sf
    elseif(coul_flag.eq.2) then
       call coulomb_dsf
    endif
     
    fx_tot = 0.0_dp; fy_tot = 0.0_dp; fz_tot = 0.0_dp     

    fx_tot = fx_a + fx_v + fx_c + fx_b
    fy_tot = fy_a + fy_v + fy_c + fy_b
    fz_tot = fz_a + fz_v + fz_c + fz_b
    
    !Add to total timing
    call cpu_time(tfinal)
    tforces = tforces + tfinal - tinit

end subroutine forces
