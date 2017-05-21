subroutine calc_angle
!******************************************************************************************
! this subroutine calculates the energetic cost of angle bending for the harmonic style
!Written by Mesele and Pubudu
!hackathon Thursday January 12, 2017
!*******************************************************************************************

 use common_variables
 use timings
 implicit none

!local variables
 integer(kind=ip) :: i,a,b,c,d
 real(kind=dp) :: dot_prod, pi, rad_equi,catmp,cfac,lenR1, lenR2, rad_ang

 !Timing variables
 real(kind=dp) :: tinit,tfinal

 call cpu_time(tinit)

!initialize forces and potential to zero
 fx_a = 0; fy_a = 0; fz_a= 0
 v_a = 0.0

!do loop to calculate the angle contribution to the potential energy and force experienced by each atom
	do i=1,n_angles

!initialize angle type and atoms involved in the angle 
		a = angle_table(i,2)
		b = angle_table(i,3)
		c = angle_table(i,4)
                d = angle_table(i,1)
!calculate parameters neccesary for bond angle evaluation

		dot_prod = rx(a,b)*rx(c,b) + ry(a,b)*ry(c,b) + rz(a,b)*rz(c,b)


		lenR1 = sqrt(rx(a,b)**2 + ry(a,b)**2 + rz(a,b)**2)
		lenR2 = sqrt(rx(b,c)**2 + ry(b,c)**2 + rz(b,c)**2)
! all angles used in the calculation of energy should be in radians as K is in kcall/mol/rad^2                
                catmp = dot_prod/(lenR1*lenR2)
		rad_ang = acos(catmp)
!                write(*,*) 'Angles'
!                write(*,"(4i2)") i , a , b , c
!                write(*,"(i2,f10.5)") i , rad_ang 
		pi = 4.0_dp*atan(1.0_dp)
		rad_equi = theta_eq(d)*pi/180_dp
   
!cfac is a common parameter for all the forces
                cfac = -k_ang(d)*(rad_ang - rad_equi)/sqrt(1.0_dp - catmp **2)
		v_a = v_a + 0.5*k_ang(d)*(rad_ang-rad_equi)**2 
!                write(*,*) v_a, theta_eq(d), rad_ang*180_dp/pi, lenR1, lenR2
                fx_a(a) = fx_a(a) - cfac*((rx(c,b))/(lenR1*lenR2) - catmp*rx(a,b)/lenR1**2)
                fy_a(a) = fy_a(a) - cfac*((ry(c,b))/(lenR1*lenR2) - catmp*ry(a,b)/lenR1**2)
                fz_a(a) = fz_a(a) - cfac*((rz(c,b))/(lenR1*lenR2) - catmp*rz(a,b)/lenR1**2)

		fx_a(b) = fx_a(b) + cfac*( (rx(a,b)+rx(c,b))/(lenR1*lenR2) - catmp*(rx(a,b)/lenR1**2 + rx(c,b)/lenR2**2) )
		fy_a(b) = fy_a(b) + cfac*( (ry(a,b)+ry(c,b))/(lenR1*lenR2) - catmp*(ry(a,b)/lenR1**2 + ry(c,b)/lenR2**2) ) 
		fz_a(b) = fz_a(b) + cfac*( (rz(a,b)+rz(c,b))/(lenR1*lenR2) - catmp*(rz(a,b)/lenR1**2 + rz(c,b)/lenR2**2) )
               
		fx_a(c) = fx_a(c) - cfac*(rx(a,b)/(lenR1*lenR2) - catmp*rx(c,b)/lenR2**2)
		fy_a(c) = fy_a(c) - cfac*(ry(a,b)/(lenR1*lenR2) - catmp*ry(c,b)/lenR2**2)
		fz_a(c) = fz_a(c) - cfac*(rz(a,b)/(lenR1*lenR2) - catmp*rz(c,b)/lenR2**2)

         enddo
!       write(1000,*) sum(fx_a(:)) + sum(fy_a(:)) + sum(fz_a(:))

         !Add to total timing
         call cpu_time(tfinal)
         tang = tang + tfinal - tinit
      
end subroutine calc_angle

         
