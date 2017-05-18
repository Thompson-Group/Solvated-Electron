!   Subroutine to run an NVE trajectory for a given number of steps

  subroutine nve_run(nstep, df_thermo, df_xyz)

       use common_variables

       implicit none

       integer(kind = ip) :: df_xyz, df_thermo
       integer(kind = ip) :: i, nstep

       do i = 1 , nstep

          ! update positions in first stage of VV integrator
          call velocity_verlet
          
          ! calculate thermodynamic properties            
          if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)
          
          ! write traj                                                                                                                     
          if (mod(i,df_xyz) .eq. 0) call dump(i)

       enddo

     end subroutine nve_run


!   Subroutine to run an NVT trajectory for a given number of steps

  subroutine nvt_run(nstep, df_thermo, df_xyz)

       use common_variables

       implicit none

       integer(kind = ip) :: df_xyz, df_thermo
       integer(kind = ip) :: i, nstep

       do i = 1 , nstep

          ! update positions in first stage of VV integrator
          call velocity_verlet
          
          ! calculate thermodynamic properties            
          if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)

          ! apply thermostat
          if (mod(i,nvt_freq).eq.0) call thermostat
          
          ! write traj                                                                                                                     
          if (mod(i,df_xyz) .eq. 0) call dump(i)

       enddo

     end subroutine nvt_run
