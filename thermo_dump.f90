!  Subroutine to open and close the trajectory (xyz) and thermodynamic data files.
!
!    

  Subroutine output_setup(oflag)

    use common_variables
    use constants
    implicit none

!Working variables
    character*5 :: oflag

!  Check to see if we are initiating or 

    if(trim(oflag).eq.'start') then
       open(nxyz,file='traj.xyz')
       open(nthermo,file='eners.dat')
       open(neigs,file='eigs.dat')
       open(nrg,file='rg.dat')
       open(nlanc,file='iterations.dat')
    elseif(trim(oflag).eq.'end') then
       close(nxyz)
       close(nthermo)
       close(neigs)
       close(nrg)
       close(nlanc)
    endif

   End Subroutine output_setup


!  Subroutine to dump the thermodynamic properties 
!    input argument is the timestep, converted to time here!
!

  Subroutine thermo_dump(istep)

    use common_variables
    use constants
    implicit none

    integer(kind=ip) :: istep
    real(kind=dp) :: time

!Working variables
    integer :: i
    real(kind=dp) :: ke, e_tot, temp_inst

!  Calculate the time

    time = real(istep)*dt
    write(6,*) istep

!  Calculate the kinetic energy

    ke = 0.0_dp
    do i = 1, n_atoms
       ke = ke + 0.5_dp*M(i)*( vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i) )
    enddo
    temp_inst = ke*2.0_dp/( (3.0_dp*real(n_atoms) - 6.0_dp)*kb )

!  Calculate the potential energy

    v_tot = v_a + v_b + v_c + v_v

!  Calculate the total energy

    e_tot = ke + v_tot

!  Dump the output

   write(nthermo,'(2F12.5,7e15.6)') time, temp_inst, e_tot, ke, v_tot, v_b, v_a, v_v, v_c

   end subroutine thermo_dump

