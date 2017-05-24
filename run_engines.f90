!   Subroutine to run an NVE trajectory for a given number of steps

    subroutine nve_run

      use common_variables
      implicit none

      integer(kind = ip) :: i
      real(kind=dp) :: tinit,tfinal

      call cpu_time(tinit)

      do i = 1 , nstep
         
         ! update positions in first stage of VV integrator
         call velocity_verlet
         
         ! calculate thermodynamic properties            
         if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)
         
         ! write traj                                                                                                                     
         if (mod(i,df_xyz) .eq. 0) call dump(i)
         
         ! write restart         
         if (mod(i,df_rest) .eq. 0) call write_restart(i)
         
      enddo
      call cpu_time(tfinal)
      write(6,'(A,F10.1,A,F10.2,A)') ' Stage cpu time = ',tfinal-tinit,' s, = ',(tfinal-tinit)/60_dp,' min'

    end subroutine nve_run


!   Subroutine to run an NVT trajectory for a given number of steps

    subroutine nvt_run

      use common_variables
      implicit none
      
      integer(kind = ip) :: i
      real(kind=dp) :: tinit,tfinal

      call cpu_time(tinit)
      
      do i = 1 , nstep
         
         ! update positions in first stage of VV integrator
         call velocity_verlet
         
         ! calculate thermodynamic properties            
         if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)
         
         ! apply thermostat
         if (mod(i,nvt_freq).eq.0) call thermostat
         
         ! write traj                                                                                                                     
         if (mod(i,df_xyz) .eq. 0) call dump(i)
         
         ! write restart
         if (mod(i,df_rest) .eq. 0) call write_restart(i)
      enddo
      call cpu_time(tfinal)
      write(6,'(A,F10.1,A,F10.2,A)') ' Stage cpu time = ',tfinal-tinit,' s, = ',(tfinal-tinit)/60_dp,' min'      
    end subroutine nvt_run

! Subroutine to do carry out a QM trajectory of the solvated electron
!   with a thermostat on the waters

    subroutine qm_nvt_run

      use common_variables
      use quantum_variables
      implicit none

      integer(kind=ip) :: i
      real(kind=dp) :: tinit,tfinal

      call cpu_time(tinit)

      i = 0
      call qm_dump(i)

      do i = 1, nstep
         
         ! update positions in first stage of VV integrator
         call qm_velocity_verlet
         
         ! calculate thermodynamic properties            
         if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)
         if (mod(i,df_thermo) .eq. 0) call qm_dump(i)
         
         ! apply thermostat
         if (mod(i,nvt_freq).eq.0) call thermostat
         
         ! write traj                                                                                                                     
         if (mod(i,df_xyz) .eq. 0) call dump(i)
         
         ! write restart
         if (mod(i,df_rest) .eq. 0) call write_restart(i)
      enddo
      call cpu_time(tfinal)
      write(6,'(A,F10.1,A,F10.2,A)') ' Stage cpu time = ',tfinal-tinit,' s, = ',(tfinal-tinit)/60_dp,' min'

    end subroutine qm_nvt_run

! Subroutine to do carry out a QM trajectory of the solvated electron in NVE

    subroutine qm_nve_run

      use common_variables
      use quantum_variables
      implicit none

      integer(kind=ip) :: i
      real(kind=dp) :: tinit,tfinal

      call cpu_time(tinit)

      i = 0
      call qm_dump(i)

      do i = 1, nstep

         ! update positions in first stage of VV integrator
         call qm_velocity_verlet
         
         ! calculate thermodynamic properties            
         if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)
         if (mod(i,df_thermo) .eq. 0) call qm_dump(i)
        
         ! write traj                                                                                                                     
         if (mod(i,df_xyz) .eq. 0) call dump(i)

         ! write restart
         if (mod(i,df_rest) .eq. 0) call write_restart(i)
      enddo
      call cpu_time(tfinal)
      write(6,'(A,F10.1,A,F10.2,A)') ' Stage cpu time = ',tfinal-tinit,' s, = ',(tfinal-tinit)/60_dp,' min'

    end subroutine qm_nve_run

