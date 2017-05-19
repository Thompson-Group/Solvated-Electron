program classical_md
!************************************************************************
!
! Fortran program for running a classical molecular dynamics simulation.
! Use velocity verlet integrator.
!
! Built by the most awesomest group of scientists.
!
! Jan. 2017
!
!************************************************************************

       use common_variables
       use input_variables

       implicit none

       integer(kind = ip) :: i, istage, d, nargs

       character(len=50) :: input_filename, arg

       logical :: restart, change, fc_flag

! read inputs from command line

       nargs = command_argument_count()

       if (nargs .ne. 2) then

            write(*,*) "Usage: 'main -in input_file' where input_file is the input script for the molecular dynamics run"

            stop

       endif

       call getarg(2,arg) 
       read(arg,*) input_filename

! read simulation input parameters
! Note: this subroutine also calls read_data.f90 to read in initial
! configurations

       write(*,*) "Reading inputs"
       call read_input(input_filename,fc_flag)

! calculate force of initial configuration

       write(*,*) "Calculating initial forces"
       call forces

! do a force check

       write(*,*) "Doing a force check"
       if(fc_flag) call force_check

! open output files

       write(*,*) "Opening output files"
       call output_setup('start')

! check for quantum stages - if they exist, allocate arrays

       call qm_allocation

!*****************************START MD SIMULATION***********************!

       write(*,*) "Beginning Molecular Dynamics!"
       write(*,*) " n_stages = ",n_stages

       do istage = 1, n_stages
          dt = sdt(istage)
          nstep = snstep(istage)
          temp = stemp(istage)
          df_xyz = sdf_xyz(istage)
          df_thermo = sdf_thermo(istage)
          df_rest = sdf_rest(istage)
          nvt_freq = snvt_freq(istage)
          write(6,*) ' stage = ',istage

          ! set initial velocites if not a restart
          if (.not. restart) call initvel

          ! write out initial thermo properties       
          call thermo_dump(0)
          write(6,*) ' made it through thermo_dump '
          write(6,*) ' srun_style(istage) = ',srun_style(istage)
          

          if(trim(srun_style(istage)).eq.'nve') then

             call nve_run

          elseif(trim(srun_style(istage)).eq.'nvt') then

             nvt_type = snvt_type(istage)
             call nvt_run

          elseif(trim(srun_style(istage)).eq.'qm_nve') then

             call qmsetup
             call qm_nve_run

          elseif(trim(srun_style(istage)).eq.'qm_nvt') then

             call qmsetup
             nvt_type = snvt_type(istage)
             call qm_nvt_run
            
          endif

       enddo

       call qmsetup

! close output files

       call output_setup('end')

end program classical_md
