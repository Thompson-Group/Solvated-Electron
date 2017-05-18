!*****************************************************************************
! subroutine to calculate the pseudo potential for the solvated electron
! without the corresponding forces on the particle
! Written by P.Wimalasiri and A.Katiyar
! Hackathon Wednesday May 17, 2017 
!******************************************************************************
  subroutine pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)
  use kinds
  use common_variables
  use quantum_variables
  use constants
 implicit none

! scalars
  real(kind=dp)    :: qo,qh,A1o,A1h,B1o,B1h,B2o,B2h,B3o,B3h,C1o
  real(kind=dp)    :: dOe,e1,e2,e3
  real(kind=dp),parameter :: alpha_pol=9.7446_dp  !7.518E7
  integer(kind=ip) :: i,k
  real(kind=dp), external :: erf
! arrays
  real(kind=dp),allocatable,dimension(:),intent(inout) :: fxtmp,fytmp,fztmp
  real(kind=dp),dimension(3)    :: L,rOe,r 
  real(kind=dp), dimension(3), intent(in) :: re
  real(kind=dp), intent(out) :: vtmp

  allocate(fxtmp(n_atoms))
  allocate(fytmp(n_atoms))
  allocate(fztmp(n_atoms))

  qo  = -0.820_dp
  qh  = 0.410_dp
  A1o = 0.575_dp
  A1h = 0.750_dp
  B1o = 0.620_dp
  B1h = 0.150_dp
  B2o = 1.000_dp
  B2h = 0.500_dp
  B3o = 0.400_dp
  B3h = 0.350_dp
  C1o = 4.400_dp

!  A1o = 0.575/angperau
!  A1h = 0.750/angperau
!  B1o = 0.620/angperau
!  B1h = 0.150/angperau
!  B2o = 1.000/angperau
!  B2h = 0.500/angperau
!  B3o = 0.400/angperau
!  B3h = 0.350/angperau
!  C1o = 4.400*angperau

  vtmp = 0.0_dp
  fxtmp = 0.0_dp
  fytmp = 0.0_dp
  fztmp = 0.0_dp
 
   
  do i = 1, n_atoms

        r(1) = x(i)
        r(2) = y(i)
        r(3) = z(i)


        L(1) = Lx
        L(2) = Ly
        L(3) = Lz

! wrap boundary conditions

        do k=1,3
                r(k) = r(k) - L(k)*int((r(k) - re(k))/L(k))
        enddo

! calculate the distances
        dOe = 0_dp
  
        do k=1,3
                rOe(k) = re(k) - r(k)
                dOe = dOe + rOe(k)**2             
        enddo
        dOe = sqrt(dOe)
   
        if (dOe.le.r_cut) then
! calculate the polarization part
         if (n_a_type .eq.1) then    
      
                vtmp = vtmp - (0.5_dp*alpha/(dOe**2+C1o**2)**2)
! the force looks like
                fxtmp(i) = fxtmp(i) + (2_dp*alpha/(rOe(1)**2+C1o**2)**3)
                fytmp(i) = fytmp(i) + (2_dp*alpha/(rOe(2)**2+C1o**2)**3)
                fztmp(i) = fztmp(i) + (2_dp*alpha/(rOe(3)**2+C1o**2)**3)
                      
! the electrostatics
         
! oxygen atom
  
                e1 = erf(A1o*dOe)  
                e2 = erf(B2o*dOe)  
                e3 = erf(B3o*dOe)  
  
                vtmp = vtmp - (e1*qo/dOe + (e2 - e3)*B1o/dOe)
                fxtmp(i) = fxtmp(i) + qo*rOe(1)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1o*exp(-(A1o*dOe)**2)/dOe**2)      & 
                 + B1o*rOe(1)*(2_dp/sqrt(pi))*((B2o*exp(-(B2o*dOe)**2))&
                 - B3o*exp(-(B3o*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3 

                fytmp(i) = fytmp(i) + qo*rOe(2)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1o*exp(-(A1o*dOe)**2)/dOe**2)      &             
                 + B1o*rOe(1)*(2_dp/sqrt(pi))*((B2o*exp(-(B2o*dOe)**2))&
                 - B3o*exp(-(B3o*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3  

                fztmp(i) = fztmp(i) + qo*rOe(3)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1o*exp(-(A1o*dOe)**2)/dOe**2)      &             
                 + B1o*rOe(1)*(2_dp/sqrt(pi))*((B2o*exp(-(B2o*dOe)**2))&
                 - B3o*exp(-(B3o*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3  

          elseif (n_a_type .eq. 2) then
! H atom
                e1 = erf(A1h*dOe)
                e2 = erf(B2h*dOe)
                e3 = erf(B3h*dOe)

                vtmp = vtmp - (e1*qh/dOe + (e2 - e3)*B1h/dOe)
                fxtmp(i) = fxtmp(i) + qh*rOe(1)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*dOe)**2)/dOe**2)      &             
                 + B1h*rOe(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*dOe)**2))&
                 - B3h*exp(-(B3h*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3  

                fytmp(i) = fytmp(i) + qh*rOe(2)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*dOe)**2)/dOe**2)      &
                 + B1h*rOe(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*dOe)**2))&
                 - B3h*exp(-(B3h*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3

                fztmp(i) = fztmp(i) + qh*rOe(3)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*dOe)**2)/dOe**2)      &
                 + B1h*rOe(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*dOe)**2))&
                 - B3h*exp(-(B3h*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3
                endif 
        endif
  enddo

return
end subroutine pseudo_e_tb
