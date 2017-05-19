subroutine initvel

       use common_variables
       use constants

       implicit none

       integer(kind = ip) :: clock(8),i
       integer,allocatable,dimension(:) :: seed
       real(kind = dp) :: RT,mtot,sx,sy,sz,num1,num2,r1,r2,ke
       real(kind = dp) :: rsq,fac,gset,gauss,vcm_x,vcm_y,vcm_z

       call random_seed()

       RT = kb * temp
       mtot = sum(M)

! loop over all atoms

       do i = 1, n_atoms

! get random number to assign sign of each velocity component

            call random_number(sx)
            call random_number(sy)
            call random_number(sz)

! call uniform numbers and transform to gaussian 
 
            transform:do

                 call random_number(num1)
                 call random_number(num2)

                 r1 = 2.0_dp * num1 - 1.0_dp
                 r2 = 2.0_dp * num2 - 1.0_dp
                 rsq = r1*r1 + r2*r2
 
                 if((rsq .lt. 1.0_dp) .and. (rsq .ne. 0.0_dp)) exit transform

            enddo transform

            fac = sqrt(-2.0_dp * log(rsq)/rsq)
            gset = r1*fac
            gauss = r2*fac
 
            vx(i) = gauss*sign(1.0_dp,2d0*sx-1.0_dp)*sqrt(RT/M(i))
            vy(i) = gauss*sign(1.0_dp,2d0*sy-1.0_dp)*sqrt(RT/M(i))
            vz(i) = gauss*sign(1.0_dp,2d0*sz-1.0_dp)*sqrt(RT/M(i))
 
            vcm_x = vcm_x + M(i)*vx(i)
            vcm_y = vcm_y + M(i)*vy(i)
            vcm_z = vcm_z + M(i)*vz(i)
          
       enddo

       vcm_x = vcm_x/mtot
       vcm_y = vcm_y/mtot
       vcm_z = vcm_z/mtot
!       vcm_x = vcm_x/(mtot*real(mol_id(n_atoms)))
!       vcm_y = vcm_y/(mtot*real(mol_id(n_atoms)))
!       vcm_z = vcm_z/(mtot*real(mol_id(n_atoms)))

!      Subtract the center of mass velocity

       vx = vx - vcm_x
       vy = vy - vcm_y
       vz = vz - vcm_z

       ke = 0.0_dp
       do i = 1, n_atoms
          ke = ke + M(i)*( vx(i)**2 + vy(i)**2 + vz(i)**2)
       enddo
       ke = ke*0.5_dp
       write(6,*) 'initvel: T_inst = ',ke*2.0_dp/(3.0_dp*real(n_atoms-2)*kb)
      
       return

end subroutine initvel
