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
  use timings
 implicit none

! scalars
  real(kind=dp)    :: dHe, dOe, e1, e2, e3, two_by_rtpi
  real(kind=dp)    :: Afac, B2fac, B3fac, fmag, r_cut_au
  integer(kind=ip) :: i, k

! arrays
  real(kind=dp), dimension(n_atoms),intent(out) :: fxtmp,fytmp,fztmp
  real(kind=dp), dimension(3) :: L, rHe, rOe, r 
  real(kind=dp), dimension(3), intent(in) :: re
  real(kind=dp), intent(out) :: vtmp
  real(kind=dp) :: vtmp_h, vtmp_o, vtmp_pol

!External functions
  real(kind=dp), external :: derf

  !Timing variables
  real(kind=dp) :: tinit,tfinal

!  write(6,*) ' alpha_pol = ',alpha_pol
!  stop
  call cpu_time(tinit)

  vtmp = 0.0_dp; fxtmp = 0.0_dp; fytmp = 0.0_dp; fztmp = 0.0_dp
  two_by_rtpi = 2.0_dp/sqrt(pi)

  vtmp_h = 0.0_dp; vtmp_o = 0.0_dp; vtmp_pol = 0.0_dp

  L(1) = Lx/angperau
  L(2) = Ly/angperau
  L(3) = Lz/angperau
  r_cut_au = r_cut/angperau
!  write(6,*) Lx, Ly, Lz
!  write(6,*) ' L(1) = ',L(1)
!  write(6,*) ' L(2) = ',L(2)
!  write(6,*) ' L(3) = ',L(3)
!  write(6,*) ' rcut = ',r_cut,' Angs = ',r_cut_au,' au'
!  write(6,*) ' alpha = ',alpha_pol
!  write(6,*) ' A1o = ',A1o
!  write(6,*) ' B1o = ',B1o
!  write(6,*) ' B2o = ',B2o
!  write(6,*) ' B3o = ',B3o
!  write(6,*) ' qo = ',qo

!  write(6,*) ' qh = ',qh
!  write(6,*) ' A1h = ',A1h
!  write(6,*) ' B1h = ',B1h
!  write(6,*) ' B2h = ',B2h
!  write(6,*) ' B3h = ',B3h
 
  do i = 1, n_atoms

     r(1) = x(i)/angperau
     r(2) = y(i)/angperau
     r(3) = z(i)/angperau
     
! wrap boundary conditions

     do k = 1, 3
        r(k) = r(k) - L(k)*anint((r(k) - re(k))/L(k))
     enddo

! calculate the distances

     dOe = 0.0_dp
     do k = 1, 3
        rOe(k) = re(k) - r(k)
        dOe = dOe + rOe(k)**2             
     enddo
     dOe = sqrt(dOe)
!     write(35,*) i, dOe
     
     if (dOe.le.r_cut_au) then

        if (a_type(i) .eq.1) then    ! oxygen atom
           
           ! calculate the polarization part
           vtmp = vtmp - 0.5_dp*alpha_pol/(dOe**2 + C1o**2)**2
!           vtmp_pol = vtmp_pol - 0.5_dp*alpha_pol/(dOe**2 + C1o**2)**2
!           write(36,'(I4,4F12.5)') i, dOe, C1o, vtmp_pol, - 0.5_dp*alpha_pol/(dOe**2 + C1o**2)**2
           
           ! the force looks like
           fmag = 2.0_dp*alpha_pol/(dOe**2 + C1o**2)**3
           fxtmp(i) = fxtmp(i) + rOe(1)*fmag
           fytmp(i) = fytmp(i) + rOe(2)*fmag
           fztmp(i) = fztmp(i) + rOe(3)*fmag
           
           ! calculate the electrostatics
           
           e1 = derf(A1o*dOe)  
           e2 = derf(B2o*dOe)  
           e3 = derf(B3o*dOe)  
           
           vtmp = vtmp - e1*qo/dOe + (e2 - e3)*B1o/dOe
!           vtmp_o = vtmp_o - e1*qo/dOe + (e2 - e3)*B1o/dOe
           
           Afac = two_by_rtpi*A1o*exp(-(A1o*dOe)**2)
           B2fac = two_by_rtpi*B2o*exp(-(B2o*dOe)**2)
           B3fac = two_by_rtpi*B3o*exp(-(B3o*dOe)**2)
           fmag = qo*(e1/dOe**3 - Afac/dOe**2)      & 
                - B1o*( (e2 - e3)/dOe**3 - (B2fac - B3fac)/dOe**2 ) 
           
           fxtmp(i) = fxtmp(i) + rOe(1)*fmag
           fytmp(i) = fytmp(i) + rOe(2)*fmag
           fztmp(i) = fztmp(i) + rOe(3)*fmag
           
        elseif (a_type(i) .eq. 2) then   ! hydrogen atom

           dHe = dOe
           rHe = rOe
           e1 = derf(A1h*dHe)
           e2 = derf(B2h*dHe)
           e3 = derf(B3h*dHe)
!           write(35,*) i,dHe, e1, e2, e3
           
           vtmp = vtmp - e1*qh/dHe + (e2 - e3)*B1h/dHe
!           vtmp_h = vtmp_h - e1*qh/dHe + (e2 - e3)*B1h/dHe
           
           Afac = two_by_rtpi*A1h*exp(-(A1h*dHe)**2)
           B2fac = two_by_rtpi*B2h*exp(-(B2h*dHe)**2)
           B3fac = two_by_rtpi*B3h*exp(-(B3h*dHe)**2)
           fmag = qh*(e1/dHe**3 - Afac/dHe**2)      &
                - B1h*( (e2 - e3)/dHe**3 - (B2fac - B3fac)/dHe**2 )
           
           fxtmp(i) = fxtmp(i) + rHe(1)*fmag              
           fytmp(i) = fytmp(i) + rHe(2)*fmag
           fztmp(i) = fztmp(i) + rHe(3)*fmag

        endif ! differentiate O from H
     endif ! check if within cuttoff radius
  enddo

!  write(34,'(4F12.5)') vtmp, vtmp_pol, vtmp_o, vtmp_h
!  stop

  !Add to total timing
  call cpu_time(tfinal)
  tpseudo = tpseudo + tfinal - tinit

end subroutine pseudo_e_tb
