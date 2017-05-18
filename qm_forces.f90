
!************************************************************************************************
!this subroutine calls other subroutine that computes the forces acting on atoms in the system
!Written by Ward Thompson
!hackathon on Thursday May 18, 2017
!*************************************************************************************************

   subroutine qm_forces

     use common_variables
     use quantum_variables
     implicit none

     integer(kind=ip) :: i 
      real(kind=dp), dimension(ng) :: try
     real(kind=dp), dimension(n_atoms) :: fx_avg, fy_avg, fz_avg

     ! First calculate the classical forces
     call forces

     ! Next update the solvated electron potential on the grid
     call v_e_update
    
! Next calculate the eigenvalues and eigenfunctions
     
     ! First calculate the initial guess
     call initial_guess(try)

     ! Then diagonalize
               !     if(trim(diag_type).eq.'lanczos') then
     call planczos(try)
               !     elseif(trim(diag_type).eq.'direct') then
               !        call direct_diag(try)
               !     endif

! Then average the forces
     call qm_force_avg(fx_avg,fy_avg,fz_avg)

     ! Add the quantum forces to the classical ones (overwrite the latter)

     fx_tot = fx_tot + fx_avg
     fy_tot = fy_tot + fy_avg
     fz_tot = fz_tot + fz_avg

     ! Add the average potential of the electron to the classical electrostatic one (overwrite the latter)
     v_c = v_c + v_e_avg

   end subroutine qm_forces
   

!
!  This subroutine updates the solvated electron potential vector, v_e by
!  recalculating it on a pre-determined grid. It does NOT determine the grid
!

   subroutine v_e_update

     use common_variables
     use quantum_variables
     implicit none

     integer(kind=ip) :: i 
     real(kind=dp) :: vtmp
     real(kind=dp), dimension(3) :: re
     real(kind=dp), dimension(n_atoms) :: fxtmp, fytmp, fztmp

     do i = 1, ng
        re(:) = rg_e(i,:)
        call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)

        v_e(i) = vtmp
        fg_ex(i,:) = fxtmp(:)
        fg_ey(i,:) = fytmp(:)
        fg_ez(i,:) = fztmp(:)

     enddo
   end subroutine v_e_update



!
!  This subroutine calculates an initial guess to the solvated
!  electron wavefunction based on the previous solution
!

   subroutine initial_guess(try)

     use quantum_variables
     implicit none

     integer(kind=ip) :: i, j, k
     real(kind=dp), dimension(ng) :: try

      try = 0.0_dp
      do i = 1, ng
         do j = 1, n_eigtol

            do k = 0, actiter - 1 
               try(i) = try(i) + Eigvec(k+1,j)*Krylov_vectors(k,i)
            enddo

         enddo
      enddo
      try = try/real(n_eigtol)


   end subroutine initial_guess

!
!  This subroutine calculates the averages of the forces on the atoms due 
!   to the solvated electron by computing the Hellmann-Feynman forces
!   It also calculates the average electron position and the radius
!   of gyration.
!

   subroutine qm_force_avg(fx_avg,fy_avg,fz_avg)

     use common_variables
     use quantum_variables
     implicit none

     integer(kind=ip) :: i, j, k
     real(kind=dp) :: psi2tmp
     real(kind=dp), dimension(3) :: re
     real(kind=dp), dimension(ng) :: psi
     real(kind=dp), dimension(n_atoms) :: fx_avg, fy_avg, fz_avg

     ! First calculate the wavefunction of the Nst state

     psi = 0.0_dp; norm = 0.0_dp
     fx_avg = 0.0_dp; fy_avg = 0.0_dp; fz_avg = 0.0_dp
     v_e_avg = 0.0_dp; r_e_avg = 0.0_dp

     do i = 1, ng

        do k = 0, actiter - 1
           psi(i) = psi(i) + Eigvec(k+1,Nst)*Krylov_vectors(k,i)
        enddo
        psi2tmp = psi(i)**2

        ! Calculate averages
        norm = norm + psi2tmp
        v_e_avg = v_e_avg + v_e(i)*psi2tmp
        r_e_avg(:) = r_e_avg(:) + rg_e(i,:)*psi2tmp

        do j = 1, n_atoms
           fx_avg(j) = fx_avg(j) + fg_ex(i,j)*psi2tmp
           fy_avg(j) = fy_avg(j) + fg_ey(i,j)*psi2tmp
           fz_avg(j) = fz_avg(j) + fg_ez(i,j)*psi2tmp
        enddo

     enddo

     ! Calculate the radius of gyration (there must be a way to do this in the loop above)
     r2_e_avg = 0.0_dp
     do i = 1, ng
        r2_e_avg = r2_e_avg + psi(i)**2*( (rg_e(i,1) - r_e_avg(1))**2 &
             + (rg_e(i,2) - r_e_avg(2))**2 + (rg_e(i,3) - r_e_avg(3))**2 )
     enddo


   end subroutine qm_force_avg
