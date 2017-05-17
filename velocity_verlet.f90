subroutine velocity_verlet
!************************************************************************
!
! Subroutine for updating positions and velocities using the velocity
! verlet algorithm
!
!************************************************************************

  use common_variables

  implicit none
  integer(kind=ip) :: i
  real(kind=dp), dimension(n_atoms) :: fx_o, fy_o, fz_o


  do i = 1, n_atoms
     ! First Step of the Velocity Verlet Algorithm
     x(i) = x(i) + vx(i)*dt + 0.5_dp*dt**2*fx_tot(i)/M(i)
     y(i) = y(i) + vy(i)*dt + 0.5_dp*dt**2*fy_tot(i)/M(i)
     z(i) = z(i) + vz(i)*dt + 0.5_dp*dt**2*fz_tot(i)/M(i)
     
     ! Wraps the coordinates into the box
     x(i) = x(i) - (Lx * anint(x(i)/Lx))
     y(i) = y(i) - (Ly * anint(y(i)/Ly))
     z(i) = z(i) - (Lz * anint(z(i)/Lz))
     
     ! Stores the Old Forces    
     fx_o(i)=fx_tot(i)
     fy_o(i)=fy_tot(i)
     fz_o(i)=fz_tot(i)
  enddo

  call forces

  do i = 1, n_atoms
     ! Updates the velocities
     vx(i)=vx(i) + 0.5_dp*dt*(fx_o(i)+fx_tot(i))/M(i)
     vy(i)=vy(i) + 0.5_dp*dt*(fy_o(i)+fy_tot(i))/M(i)
     vz(i)=vz(i) + 0.5_dp*dt*(fz_o(i)+fz_tot(i))/M(i)
  enddo

end subroutine
