!  Subroutine to check the force calcualtion
!     An atom is chosen at random, moved by a small amount and                                          
!     numerical derivatives taken and compared to the forces         
!
  Subroutine force_check(fc_label)

    use kinds
    use common_variables
    implicit none
    
    character*3 :: fc_label

!Working variables
    integer :: iat
    integer, parameter :: nfc=29
    real(kind=dp) :: num
    real(kind=dp) :: fx0, fy0, fz0, fxnum, fynum, fznum
    real(kind=dp) :: fx0_a, fy0_a, fz0_a, fx0_b, fy0_b, fz0_b
    real(kind=dp) :: fx0_c, fy0_c, fz0_c, fx0_v, fy0_v, fz0_v
    real(kind=dp) :: fxnum_a, fynum_a, fznum_a, fxnum_b, fynum_b, fznum_b
    real(kind=dp) :: fxnum_c, fynum_c, fznum_c, fxnum_v, fynum_v, fznum_v
    real(kind=dp) :: v0, vp, vm
    real(kind=dp) :: v0_a, v0_b, v0_c, v0_v
    real(kind=dp) :: vm_a, vm_b, vm_c, vm_v, vp_a, vp_b, vp_c, vp_v


! Choose an atom at random

    call random_seed()
    call random_number(harvest = num)
    iat = nint( num*real(n_atoms-1) ) + 1

! Calculate the forces and potential

    call forces

    fx0 = fx_tot(iat); fy0 = fy_tot(iat); fz0 = fz_tot(iat)
    fx0_a = fx_a(iat); fy0_a = fy_a(iat); fz0_a = fz_a(iat)
    fx0_b = fx_b(iat); fy0_b = fy_b(iat); fz0_b = fz_b(iat)
    fx0_c = fx_c(iat); fy0_c = fy_c(iat); fz0_c = fz_c(iat)
    fx0_v = fx_v(iat); fy0_v = fy_v(iat); fz0_v = fz_v(iat)

    v0 = v_a + v_b + v_c + v_v
    v0_a = v_a; v0_b = v_b; v0_c = v_c; v0_v = v_v

! Calculate the numerical force in the x-direction
! shift the atom in the +x direction and recalculate the potential

    x(iat) = x(iat) + delta

    call forces

    vp = v_a + v_b + v_c + v_v
    vp_a = v_a; vp_b = v_b; vp_c = v_c; vp_v = v_v

! shift the atom in the -x direction and recalculate the potential

    x(iat) = x(iat) - 2.0_dp*delta

    call forces

    vm = v_a + v_b + v_c + v_v
    vm_a = v_a; vm_b = v_b; vm_c = v_c; vm_v = v_v

!   Calculate the forces numerically

    fxnum = -(vp - vm)/(2.0_dp*delta)
    fxnum_a = -(vp_a - vm_a)/(2.0_dp*delta)
    fxnum_b = -(vp_b - vm_b)/(2.0_dp*delta)
    fxnum_c = -(vp_c - vm_c)/(2.0_dp*delta)
    fxnum_v = -(vp_v - vm_v)/(2.0_dp*delta)
    
! Reset the x position

    x(iat) = x(iat) + delta

! Calculate the numerical force in the y-direction
! shift the atom in the +y direction and recalculate the potential

    y(iat) = y(iat) + delta

    call forces

    vp = v_a + v_b + v_c + v_v
    vp_a = v_a; vp_b = v_b; vp_c = v_c; vp_v = v_v

! shift the atom in the -y direction and recalculate the potential

    y(iat) = y(iat) - 2.0_dp*delta

    call forces

    vm = v_a + v_b + v_c + v_v
    vm_a = v_a; vm_b = v_b; vm_c = v_c; vm_v = v_v

!   Calculate the forces numerically

    fynum = -(vp - vm)/(2.0_dp*delta)
    fynum_a = -(vp_a - vm_a)/(2.0_dp*delta)
    fynum_b = -(vp_b - vm_b)/(2.0_dp*delta)
    fynum_c = -(vp_c - vm_c)/(2.0_dp*delta)
    fynum_v = -(vp_v - vm_v)/(2.0_dp*delta)
    
! Reset the y position

    y(iat) = y(iat) + delta

! Calculate the numerical force in the z-direction
! shift the atom in the +z direction and recalculate the potential

    z(iat) = z(iat) + delta

    call forces

    vp = v_a + v_b + v_c + v_v
    vp_a = v_a; vp_b = v_b; vp_c = v_c; vp_v = v_v

! shift the atom in the -z direction and recalculate the potential

    z(iat) = z(iat) - 2.0_dp*delta

    call forces

    vm = v_a + v_b + v_c + v_v
    vm_a = v_a; vm_b = v_b; vm_c = v_c; vm_v = v_v

!   Calculate the forces numerically

    fznum = -(vp - vm)/(2.0_dp*delta)
    fznum_a = -(vp_a - vm_a)/(2.0_dp*delta)
    fznum_b = -(vp_b - vm_b)/(2.0_dp*delta)
    fznum_c = -(vp_c - vm_c)/(2.0_dp*delta)
    fznum_v = -(vp_v - vm_v)/(2.0_dp*delta)
    
! Reset the z position

    z(iat) = z(iat) + delta

!   Write out results

    open(nfc,file='force_check_'//trim(fc_label))

    write(nfc,'(A,I4)') ' Atom chosen = ',iat
    write(nfc,*)
    write(nfc,*) ' ----- Total Forces -----'
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fx = ',fx0,' Numerical:  fx = ',fxnum
    write(nfc,'(A,F15.8)') ' % Error in fx = ',(fx0-fxnum)*100.0_dp/fxnum
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fy = ',fy0,' Numerical:  fy = ',fynum
    write(nfc,'(A,F15.8)') ' % Error in fy = ',(fy0-fynum)*100.0_dp/fynum
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fz = ',fz0,' Numerical:  fz = ',fznum
    write(nfc,'(A,F15.8)') ' % Error in fz = ',(fz0-fznum)*100.0_dp/fznum

    write(nfc,*)
    write(nfc,*) ' ----- Bond Stretching Forces -----'
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fx = ',fx0_b,' Numerical:  fx = ',fxnum_b
    write(nfc,'(A,F15.8)') ' % Error in fx = ',(fx0_b-fxnum_b)*100.0_dp/fxnum_b
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fy = ',fy0_b,' Numerical:  fy = ',fynum_b
    write(nfc,'(A,F15.8)') ' % Error in fy = ',(fy0_b-fynum_b)*100.0_dp/fynum_b
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fz = ',fz0_b,' Numerical:  fz = ',fznum_b
    write(nfc,'(A,F15.8)') ' % Error in fz = ',(fz0_b-fznum_b)*100.0_dp/fznum_b

    write(nfc,*)
    write(nfc,*) ' ----- Angle Bending Forces -----'
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fx = ',fx0_a,' Numerical:  fx = ',fxnum_a
    write(nfc,'(A,F15.8)') ' % Error in fx = ',(fx0_a-fxnum_a)*100.0_dp/fxnum_a
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fy = ',fy0_a,' Numerical:  fy = ',fynum_a
    write(nfc,'(A,F15.8)') ' % Error in fy = ',(fy0_a-fynum_a)*100.0_dp/fynum_a
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fz = ',fz0_a,' Numerical:  fz = ',fznum_a
    write(nfc,'(A,F15.8)') ' % Error in fz = ',(fz0_a-fznum_a)*100.0_dp/fznum_a

    write(nfc,*)
    write(nfc,*) ' ----- van der Waals Forces -----'
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fx = ',fx0_v,' Numerical:  fx = ',fxnum_v
    write(nfc,'(A,F15.8)') ' % Error in fx = ',(fx0_v-fxnum_v)*100.0_dp/fxnum_v
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fy = ',fy0_v,' Numerical:  fy = ',fynum_v
    write(nfc,'(A,F15.8)') ' % Error in fy = ',(fy0_v-fynum_v)*100.0_dp/fynum_v
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fz = ',fz0_v,' Numerical:  fz = ',fznum_v
    write(nfc,'(A,F15.8)') ' % Error in fz = ',(fz0_v-fznum_v)*100.0_dp/fznum_v

    write(nfc,*)
    write(nfc,*) ' ----- Electrostatic Forces -----'
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fx = ',fx0_c,' Numerical:  fx = ',fxnum_c
    write(nfc,'(A,F15.8)') ' % Error in fx = ',(fx0_c-fxnum_c)*100.0_dp/fxnum_c
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fy = ',fy0_c,' Numerical:  fy = ',fynum_c
    write(nfc,'(A,F15.8)') ' % Error in fy = ',(fy0_c-fynum_c)*100.0_dp/fynum_c
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fz = ',fz0_c,' Numerical:  fz = ',fznum_c
    write(nfc,'(A,F15.8)') ' % Error in fz = ',(fz0_c-fznum_c)*100.0_dp/fznum_c

    close(nfc)

    End subroutine force_check

    
