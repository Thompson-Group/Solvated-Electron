subroutine integrator
!************************************************************************
!
! Subroutine for updating positions and velocities using the leapfrog
! verlet algorithm
!
!************************************************************************

      use common_variables

      implicit none
  
      real (kind = dp), allocatable, dimension(:) :: vtempx, vtempy, vtempz, vtemp2x, vtemp2y, vtemp2z
      integer (kind= ip) :: i

 allocate(vtempx(n_atoms))
 allocate (vtempy(n_atoms))
 allocate(vtempz(n_atoms))
 allocate(vtemp2x(n_atoms))
 allocate (vtemp2y(n_atoms))
 allocate(vtemp2z(n_atoms))
    
 do i = 1, n_atoms
    
    !updates velocity to t + 0.5dt
    vtempx(i) = vx(i) + dt/2_dp*fx_tot(i)/M(i)
    vtempy(i) = vy(i) + dt/2_dp*fy_tot(i)/M(i)
    vtempz(i) = vz(i) + dt/2_dp*fz_tot(i)/M(i)
    
    ! updates positions to t + dt
    x(i) = x(i) + vtempx(i)*dt
    y(i) = y(i) + vtempy(i)*dt
    z(i) = z(i) + vtempz(i)*dt
    
    ! wrap positions into box
    
    x(i) = x(i) - (Lx * anint(x(i)/Lx))
    y(i) = y(i) - (Ly * anint(y(i)/Ly))
    z(i) = z(i) - (Lz * anint(z(i)/Lz))

 enddo

 call forces
 
 do i = 1, n_atoms

    ! updates velocity to t + 3/2dt
    vtemp2x(i) = vtempx(i) + fx_tot(i)/M(i)*dt
    vtemp2y(i) = vtempy(i) + fy_tot(i)/M(i)*dt
    vtemp2z(i) = vtempz(i) + fz_tot(i)/M(i)*dt
    
    ! update current velocity
    vx(i) = 0.5*(vtemp2x(i) - vtempx(i)) 
    vy(i) = 0.5*(vtemp2y(i) - vtempy(i)) 
    vz(i) = 0.5*(vtemp2z(i) - vtempz(i)) 
               
 enddo


end subroutine
