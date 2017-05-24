!  Subroutine to check the force calcualtion
!     An atom is chosen at random, moved by a small amount and                                          
!     numerical derivatives taken and compared to the forces         
!
  Subroutine pseudo_force_check

    use common_variables
    use quantum_variables
    use constants
    implicit none
    
!Working variables
    integer :: iat
    integer, parameter :: nfc=29
    real(kind=dp) :: num
    real(kind=dp) :: fx0, fy0, fz0, fxnum, fynum, fznum
    real(kind=dp) :: v0, vp, vm
    real(kind=dp) :: vtmp
    real(kind=dp), dimension(3) :: re
    real(kind=dp), dimension(n_atoms) :: fxtmp, fytmp, fztmp
    character(len=2) :: fc_ext

! Choose an atom at random

    call random_seed()
    call random_number(harvest = num)
    iat = nint( num*real(n_atoms-1) ) + 1
!    iat = 2

! Put the electron at the origin for simplicity

    re = 0.0_dp

! Calculate the forces and potential

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    fx0 = fxtmp(iat)*kcalperau/angperau; fy0 = fytmp(iat)*kcalperau/angperau
    fz0 = fztmp(iat)*kcalperau/angperau

    v0 = vtmp*kcalperau

! Calculate the numerical force in the x-direction
! shift the atom in the +x direction and recalculate the potential

    x(iat) = x(iat) + delta

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    vp = vtmp*kcalperau

! shift the atom in the -x direction and recalculate the potential

    x(iat) = x(iat) - 2.0_dp*delta

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    vm = vtmp*kcalperau

!   Calculate the forces numerically

    fxnum = -(vp - vm)/(2.0_dp*delta)
    
! Reset the x position

    x(iat) = x(iat) + delta

! Calculate the numerical force in the y-direction
! shift the atom in the +y direction and recalculate the potential

    y(iat) = y(iat) + delta

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    vp = vtmp*kcalperau

! shift the atom in the -y direction and recalculate the potential

    y(iat) = y(iat) - 2.0_dp*delta

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    vm = vtmp*kcalperau

!   Calculate the forces numerically

    fynum = -(vp - vm)/(2.0_dp*delta)
    
! Reset the y position

    y(iat) = y(iat) + delta

! Calculate the numerical force in the z-direction
! shift the atom in the +z direction and recalculate the potential

    z(iat) = z(iat) + delta

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    vp = vtmp*kcalperau

! shift the atom in the -z direction and recalculate the potential

    z(iat) = z(iat) - 2.0_dp*delta

    call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

    vm = vtmp*kcalperau

!   Calculate the forces numerically

    fznum = -(vp - vm)/(2.0_dp*delta)
    
! Reset the z position

    z(iat) = z(iat) + delta

!   Write out results

!    write(fc_ext,'(I0)') qfc_cnt
!    open(nfc,file='pseudo_force_check.'//trim(fc_ext))
    open(nfc,file='pseudo_force_check.dat')

    write(nfc,'(A,I4)')    ' Atom chosen = ',iat
    write(nfc,'(A,F12.5)') ' Increment   = ',delta
    write(nfc,*)
    write(nfc,*) ' ----- Total Forces -----'
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fx = ',fx0,' Numerical:  fx = ',fxnum
    write(nfc,'(A,F15.8)') ' % Error in fx = ',(fx0-fxnum)*100.0_dp/fx0
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fy = ',fy0,' Numerical:  fy = ',fynum
    write(nfc,'(A,F15.8)') ' % Error in fy = ',(fy0-fynum)*100.0_dp/fy0
    write(nfc,*)
    write(nfc,'(A,F15.8,A,F15.8)') ' Analytical: fz = ',fz0,' Numerical:  fz = ',fznum
    write(nfc,'(A,F15.8)') ' % Error in fz = ',(fz0-fznum)*100.0_dp/fz0

    close(nfc)

    End subroutine pseudo_force_check

    
