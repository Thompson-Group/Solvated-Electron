subroutine integrator(Stage)
!************************************************************************
!
! Subroutine for updating positions and velocities using the velocity
! verlet algorithm
!
!************************************************************************

      use common_variables

      implicit none
  
      integer (kind = ip) :: Stage,i,j
    
! first stage of VV updates positions
   
      if(Stage .eq. 1) then

           do i = 1, n_atoms
     
! updates positions to t + dt

               x(i) = x_o(i) + vx_o(i)* dt +(fx_o(i)*dt**2)/2*M(i)
               y(i) = y_o(i) + vy_o(i)* dt +(fy_o(i)*dt**2)/2*M(i)
               z(i) = z_o(i) + vz_o(i)* dt +(fz_o(i)*dt**2)/2*M(i)
 
! wrap positions into box

               x(i) = x(i) - (Lx * anint(x(i)/Lx))
               y(i) = y(i) - (Ly * anint(y(i)/Ly))
               z(i) = z(i) - (Lz * anint(z(i)/Lz))
              
! save new positions to old positions

               x_o(i) = x(i)
               y_o(i) = y(i)
               z_o(i) = z(i)
          
           enddo

! second stage of VV updates velocities 

      else if(Stage .eq. 2) then 

            do j = 1, n_atoms

! calculate t + dt velocities
               
               vx(j) = vx_o(j) + ((fx_tot(j) +fx_o(j))/2)*dt*M(j)
               vy(j) = vy_o(j) + ((fy_tot(j) +fy_o(j))/2)*dt*M(j)
               vz(j) = vz_o(j) + ((fz_tot(j) +fz_o(j))/2)*dt*M(j)

! save new velocities to old velocities

               vx_o(j) = vx(j)
               vy_o(j) = vy(j)
               vz_o(j) = vz(j)

! save old forces to new forces
           
               fx_o(j) = fx_tot(j)
               fy_o(j) = fy_tot(j)
               fz_o(j) = fz_tot(j)
            
            enddo
     
       
      endif

end subroutine
