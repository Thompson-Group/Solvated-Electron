!*****************************************************************************
! subroutine to calculate the pseudo potential for the solvated electron
! without the corresponding forces on the particle
! Written by P.Wimalasiri and A.Katiyar
! Hackathon Wednesday May 17, 2017 
!******************************************************************************
  subroutine pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

  use common_variables
  use quantum_variables
  use constants
  use pseudo_constants
 implicit none

! scalars
  real(kind=dp)    :: dHe, dOe, e1, e2, e3, two_by_rtpi
  real(kind=dp)    :: Afac, B2fac, B3fac
  integer(kind=ip) :: i, k

! arrays
  real(kind=dp), dimension(n_atoms),intent(inout) :: fxtmp,fytmp,fztmp
  real(kind=dp),dimension(3)    :: L, rHe, rOe, r 
  real(kind=dp), dimension(3), intent(in) :: re
  real(kind=dp), intent(out) :: vtmp

  vtmp = 0.0_dp; fxtmp = 0.0_dp; fytmp = 0.0_dp; fztmp = 0.0_dp
  two_by_rtpi = 2.0_dp/sqrt(pi)

  do i = 1, n_atoms

        r(1) = x(i)
        r(2) = y(i)
        r(3) = z(i)

        L(1) = Lx
        L(2) = Ly
        L(3) = Lz

! wrap boundary conditions

        do k=1,3
           r(k) = r(k) - L(k)*anint((r(k) - re(k))/L(k))
        enddo

! calculate the distances

        dOe = 0.0_dp
        do k = 1, 3
           rOe(k) = re(k) - r(k)
           dOe = dOe + rOe(k)**2             
        enddo
        dOe = sqrt(dOe)
   
        if (dOe.le.r_cut) then

           if (n_a_type .eq.1) then    ! oxygen atom
      
              ! calculate the polarization part
              vtmp = vtmp - 0.5_dp*alpha/(dOe**2 + C1o**2)**2

              ! the force looks like
              fxtmp(i) = fxtmp(i) + 2.0_dp*alpha*rOe(1)/(dOe**2 + C1o**2)**3
              fytmp(i) = fytmp(i) + 2.0_dp*alpha*rOe(2)/(dOe**2 + C1o**2)**3
              fztmp(i) = fztmp(i) + 2.0_dp*alpha*rOe(3)/(dOe**2 + C1o**2)**3              

              ! calculate the electrostatics
  
              e1 = erf(A1o*dOe)  
              e2 = erf(B2o*dOe)  
              e3 = erf(B3o*dOe)  
  
              vtmp = vtmp - e1*qo/dOe + (e2 - e3)*B1o/dOe

              Afac = two_by_rtpi*A1o*exp(-(A1o*dOe)**2)
              B2fac = two_by_rtpi*B2o*exp(-(B2o*dOe)**2)
              B3fac = two_by_rtpi*B3o*exp(-(B3o*dOe)**2)

              fxtmp(i) = fxtmp(i) + rOe(1)*( qo*(e1/dOe**3 - Afac/dOe**2)      & 
                   - B1o*( (e2 - e3)/dOe**3 - (B2fac - B3fac)/dOe**2 ) )

              fytmp(i) = fytmp(i) + rOe(2)*( qo*(e1/dOe**3 - Afac/dOe**2)      &             
                   - B1o*( (e2 - e3)/dOe**3 - (B2fac - B3fac)/dOe**2) )

              fztmp(i) = fztmp(i) + rOe(3)*( qo*(e1/dOe**3 - Afac/dOe**2)      &             
                   - B1o*( (e2 - e3)/dOe**3 - (B2fac - B3fac)/dOe**2) )

           elseif (n_a_type .eq. 2) then   ! hydrogen atom

              dHe = dOe
              rHe = rOe
              e1 = erf(A1h*dHe)
              e2 = erf(B2h*dHe)
              e3 = erf(B3h*dHe)
              
              vtmp = vtmp - e1*qh/dHe + (e2 - e3)*B1h/dHe

              Afac = two_by_rtpi*A1h*exp(-(A1h*dHe)**2)
              B2fac = two_by_rtpi*B2h*exp(-(B2h*dHe)**2)
              B3fac = two_by_rtpi*B3h*exp(-(B3h*dHe)**2)

              fxtmp(i) = fxtmp(i) + rHe(1)*( qh*(e1/dHe**3 - Afac/dHe**2)      &             
                   - B1h*( (e2 - e3)/dHe**3 - (B2fac - B3fac)/dHe**2 ) )
              
              fytmp(i) = fytmp(i) + rHe(2)*( qh*(e1/dHe**3 - Afac/dHe**2)      &
                   - B1h*( (e2 - e3)/dHe**3 - (B2fac - B3fac)/dHe**2 ) )
              
              fztmp(i) = fztmp(i) + rHe(3)*( qh*(e1/dHe**3 - Afac/dHe**2)      &
                   - B1h*( (e2 - e3)/dHe**3 - (B3fac - B3fac)/dHe**2 ) )
           endif
        endif
  enddo

return
end subroutine pseudo_e_tb
