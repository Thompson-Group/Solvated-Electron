subroutine calc_vdw()
!************************************************************************
! 
! Subroutine for calculating LJ energy and force.
!
!************************************************************************

     use common_variables
     implicit none

     integer(kind = ip) :: i,j
     real(kind = dp) :: dist, ftmp

! initialize
    
     v_v = 0.0_dp; fx_v = 0.0_dp; fy_v = 0.0_dp; fz_v = 0.0_dp     

! loop over all atoms
     
     do i = 1, n_atoms-1

! loop over all later atoms

         do j = i + 1, n_atoms 

! distance
             
            dist = sqrt(rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2)

!            write(*,"(2i2,f10.5)") i , j , dist
     
            if (dist .lt. r_cut) then

!               write(*,*) "Distance less than cutoff!"

! energy
          
               v_v = v_v + 4.0_dp*ee(i,j)*((ss(i,j)/dist)**12-(ss(i,j)/dist)**6)

!               if ((i.eq.1).and.(j.eq.4)) then

!                    write(*,*) "oxygen-oxygen interaction"
!                    write(*,"(3f10.5)") ee(i,j) , ss(i,j) , dist

!               endif
!               write(*,*) "New VDW energy: " , v_v
! force (-du/dx)

               ftmp = 24.0_dp*ee(i,j)*(2.0_dp*(ss(i,j)/dist)**12 - (ss(i,j)/dist)**6)/dist**2
               fx_v(i) = fx_v(i) + ftmp*rx(i,j)
               fy_v(i) = fy_v(i) + ftmp*ry(i,j)
               fz_v(i) = fz_v(i) + ftmp*rz(i,j)

!               if ((i.eq.1).and.(j.eq.4)) then

!                    write(*,*) "oxygen-oxygen force"
!                    write(*,"(3f10.5)") fx_v(i) , fy_v(i) , fz_v(i)

!               endif

! equal and opposite force

               fx_v(j) = fx_v(j) - ftmp*rx(i,j)
               fy_v(j) = fy_v(j) - ftmp*ry(i,j)
               fz_v(j) = fz_v(j) - ftmp*rz(i,j)

            endif
         
          enddo

    enddo

end subroutine 
