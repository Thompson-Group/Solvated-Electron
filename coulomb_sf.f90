! Subroutine to calculate the electrostatic interactions between atoms on different molecules
!   within a given cuttoff radius.  This is the *SHIFTED FORCE* Coulomb's law
!
!         V_c(r_ij) = C*q(i)*q(j)/r_ij
!
!         V_ij(r_ij) = V_c(r_ij) - V_c(r_cut) - (dV/dr)|_r_cut * (r_ij - r_cut)
!
! This routine assumes an array of the charge products qq(i,j) = q(i)*q(j) has been calculated
!
! The units are kcal/mol (energy), electron charge, e (charge), and Angstroms (distance)
!
! Jan. 12, 2017 - Paul Burris & Ward Thompson
!

  Subroutine coulomb_sf

  use common_variables
  use constants
  implicit none

!Working variables
  integer i, j
  real(kind=dp) :: rij, fxtmp, fytmp, fztmp

! Loop over pairs of atoms

  v_c = 0.0_dp; fx_c = 0.0_dp; fy_c = 0.0_dp; fz_c = 0.0_dp

  do j = 1, n_atoms - 1
     do i = j + 1, n_atoms

        rij = sqrt( rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2 )  ! Calculate the atom-atom distance

        if (rij.le.r_cut) then
           
           v_c = v_c + C_coul*qq(i,j)*(1.0_dp/rij - 1.0_dp/r_cut + (rij - r_cut)/r_cut**2)  ! Calculate the contribution to the potential

           fxtmp = C_coul*qq(i,j)*rx(i,j)*(1.0_dp/rij**2 - 1.0_dp/r_cut**2)/rij
           fytmp = C_coul*qq(i,j)*ry(i,j)*(1.0_dp/rij**2 - 1.0_dp/r_cut**2)/rij
           fztmp = C_coul*qq(i,j)*rz(i,j)*(1.0_dp/rij**2 - 1.0_dp/r_cut**2)/rij

           fx_c(i) = fx_c(i) + fxtmp 
           fy_c(i) = fy_c(i) + fytmp
           fz_c(i) = fz_c(i) + fztmp

           fx_c(j) = fx_c(j) - fxtmp 
           fy_c(j) = fy_c(j) - fytmp
           fz_c(j) = fz_c(j) - fztmp

        endif

     enddo
  enddo

 End Subroutine coulomb_sf

