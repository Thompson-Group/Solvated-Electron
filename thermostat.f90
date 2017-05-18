   subroutine thermostat

     use common_variables
     use constants
     implicit none

     integer (kind = ip) :: i,j,step(8)
     integer :: n
     real (kind = dp) :: temp_inst,KE,vxb,vyb,vzb,numx,lambda,mass
     integer,dimension(:),allocatable :: seed
 
     call random_seed()

     
     if(trim(nvt_type) .eq. "rescale") then 

        ke = 0.0_dp
        do i = 1, n_atoms
           ke = ke + M(i)*(vx(i)**2 + vy(i)**2 + vz(i)**2)
        enddo
        ke = 0.5_dp*ke
        
        temp_inst = 2.0_dp*ke/(3.0_dp*real(n_atoms-2)*kb)
        
        if (temp_inst .gt. 0.0_dp) then
           lambda = sqrt(temp/temp_inst)
             
           vx = lambda*vx
           vy = lambda*vy
           vz = lambda*vz
        endif

     elseif(trim(nvt_type) .eq.  "andersen") then
           
        do j = 1, n_atoms
           
           call random_number(numx)
           if(numx .lt. nu*dt) then

              mass = M(j)
              call boltz_vel(vxb,vyb,vzb,mass)  
              vx(j)= vxb
              vy(j)= vyb
              vz(j)= vzb  
              
           endif
        enddo
     endif

end subroutine thermostat

       





 
