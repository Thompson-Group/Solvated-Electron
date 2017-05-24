! Subroutine to calculate the electrostatic interactions between atoms on different molecules
!   within a given cuttoff radius.  This is the *DAMPED SHIFTED FORCE* Coulomb's law
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

  Subroutine coulomb_dsf

  use common_variables
  use constants
  use timings
  implicit none

!Working variables
  integer i, j
  real(kind=dp) :: rij, ftmp
  real(kind=dp) :: Cpref, pref, etmp, ecut, exptmp, expcut

!External functions
  real(kind=dp), external :: derfc

  !Timing variables
  real(kind=dp) :: tinit,tfinal
  
  call cpu_time(tinit)

! Loop over pairs of atoms

  v_c = 0.0_dp; fx_c = 0.0_dp; fy_c = 0.0_dp; fz_c = 0.0_dp

  pref = 2.0_dp*alpha/sqrt(pi)
  ecut = derfc(alpha*r_cut)/r_cut**2
  expcut = pref*exp(-(alpha*r_cut)**2)/r_cut

  do j = 1, n_atoms - 1
     do i = j + 1, n_atoms

        rij = sqrt( rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2 )  ! Calculate the atom-atom distance

        if (rij.le.r_cut) then
           
           !Define some quantities used multiple times
           Cpref = C_coul*qq(i,j)
           etmp = derfc(alpha*rij)/rij**2
           exptmp = pref*exp(-(alpha*rij)**2)/rij

           v_c = v_c + Cpref*( etmp*rij - ecut*r_cut + ( ecut + expcut)*(rij-r_cut) )  ! potential

           ftmp = ( Cpref/rij )*( (etmp + exptmp) - (ecut + expcut) )  

           fx_c(i) = fx_c(i) + ftmp*rx(i,j) 
           fy_c(i) = fy_c(i) + ftmp*ry(i,j)
           fz_c(i) = fz_c(i) + ftmp*rz(i,j) 

           fx_c(j) = fx_c(j) - ftmp*rx(i,j) 
           fy_c(j) = fy_c(j) - ftmp*ry(i,j)
           fz_c(j) = fz_c(j) - ftmp*rz(i,j)

        endif

     enddo
  enddo

  !Add to total timing
  call cpu_time(tfinal)
  tcoul = tcoul + tfinal - tinit

 End Subroutine coulomb_dsf

