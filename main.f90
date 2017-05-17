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

       implicit none

       integer(kind = ip) :: df_xyz,df_thermo,df_rest,nstep,fc_flag
       integer(kind = ip) :: i,d,nargs

       character(len=50) :: input_filename,nvt_type,arg

       logical :: restart,change

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
       call read_input(input_filename,df_xyz,df_thermo,df_rest,nvt_type,fc_flag,&
                       nstep)

! set initial velocites if not a restart

       if (.not. restart) call initvel

! calculate force of initial configuration

       write(*,*) "Calculating initial forces"
       call forces

! do a force check

       write(*,*) "Doing a force check"
       call force_check

! open output files

       write(*,*) "Opening output files"
       call output_setup('start')

! write out initial thermo properties
       
       call thermo_dump(0)

!*****************************START MD SIMULATION***********************!

       write(*,*) "Beginning Molecular Dynamics!"
       do i = 1 , nstep

!            write(*,*) "Step: " , i 
! update positions in first stage of VV integrator

            call velocity_verlet

! calculate thermodynamic properties

            if (mod(i,df_thermo) .eq. 0) call thermo_dump(i)

! apply thermostat

            
            if ((trim(nvt_type) .eq. "rescale" .or. &
                trim(nvt_type) .eq. "anderson")) then

                 call thermostat(nvt_type)

            endif

! write traj

            if (mod(i,df_xyz) .eq. 0) call dump(i)

       enddo

! close output files

       call output_setup('end')

end program classical_md
