subroutine read_input(input_filename, fc_flag)
    use kinds
    use common_variables
    use quantum_variables
    integer(kind=ip) :: i, j, ign_int, done
    character(len=50) :: data_filename, word, input_filename, coul_tmp
    logical :: restart,fc_flag
    open(unit=ninput, file=input_filename)
     

    read(ninput,*) word 
    if(word .eq. 'read_data') then
       read(ninput,*) data_filename
       restart = .false.   
    else if(word .eq. 'read_restart') then
       read(ninput,*) data_filename
       restart = .true. 
    else
       write(*,*) "Error : Data file not found at top of input file"
    end if

    call read_data(data_filename,restart)

        read(ninput,*)
    loopread:do 
        read(ninput,*,iostat=done) word
        if(done .lt. 0) exit loopread

        if (word .eq. 'bond_style') then
            read(ninput,*) bond_style
        else if (word .eq. 'bond_coeffs') then
            do i = 1, n_b_type
            ! Provide the bond coeffs in bond type order            
                read(ninput,*)  ign_int , k_r(i), req(i)
            enddo

        else if (word .eq. 'angle_coeffs') then 
            do i = 1, n_angle_type
            ! Provide the angle coeffs in angle type order  
                read(ninput,*) ign_int, k_ang(i), theta_eq(i)
            enddo
        else if (word .eq. 'cutoff_distance') then
            read(ninput,*) r_cut
        !coul_flag = 0: off, =1: shifted force, =2: damped shifted force
        else if (word .eq. 'coulomb_style') then 
            read(ninput,*) alpha, coul_tmp
            if (coul_tmp .eq. 'none') then
                coul_flag = 0
            else if (coul_tmp .eq. 'sf') then
                coul_flag = 1
            else if (coul_tmp .eq. 'dsf') then
                coul_flag = 2
            else
                write(*,*) "Error (read_input.f90): Invalid electrostatics style"
                stop
            end if
        else if (word .eq. 'force_check') then
            read(ninput,*) fc_flag,del
        else if (word .eq. 'stages') then
            read(ninput,*) n_stages
            if (n_stages .gt. 5) then
                write(*,*) "Error (read_input.f90): You have selected way too
                many stages. Rethink your life choices, and choose a number that is five or
                less. Do not pass go. Do not collect 200 dollars."
                stop
            endif
            do i=1,n_stages
                read(ninput,*) run_style(i)
                if (run_style(i) .eq. 'nve') then
                    nvt_type(i)='none'
                    nvt_freq(i)=0
                else if (run_style(i) .eq. 'nvt') then
                    read(ninput,*) nvt_type(i), nvt_freq(i)
                    if (nvt_type(i) .eq. 'andersen') then
                        read(ninput,*) nu
                else if (run_style(i) .eq. 'qm_nve') then
                    nvt_type(i)='none'
                    nvt_freq(i)=0
                else if (run_style(i) .eq. 'qm_nvt') then
                    read(ninput,*) nvt_type(i), nvt_freq(i)
                    if (nvt_type(i) .eq. 'andersen') then
                        read(ninput,*) nu
                    endif
                else
                    write(*,*) "Error (read_input.f90): Couldn't find the
                    run_style in the input file."
                    stop
                endif
                if (nvt_type(i) .ne. 'rescale' .or. nvt_type(i) .ne. 'andersen'
                .or. nvt_type(i) .ne. 'none') then
                    write(*,*) "Error (read_input.f90): Incorrect nvt_type"
                    stop
                endif
            enddo
        else if (word .eq. 'timestep') then
            if (n_stages .eq. 0) then
                write(*,*) "Error (read_input.f90): Timestep defined before
                number of stages."
                stop
            endif
            do i=1,n_stages
                read(ninput,*) dt(i)
            enddo
            
        else if (word .eq. 'run') then 
            if (n_stages .eq. 0) then
                write(*,*) "Error (read_input.f90): Step defined before number
                of stages."
                stop
            endif
            do i=1,n_stages
                read(ninput,*) nstep(i)
            enddo

        else if (word .eq. 'dump_freq') then 
            if (n_stages .eq. 0) then
                write(*,*) "Error (read_input.f90): Dump Frequency defined
                before number of stages."
                stop
            endif
            do i=1,n_stages
                read(ninput,*) df_xyz(i), df_thermo(i), df_rest(i)
            enddo
        else if (word .eq. 'temperature') then 
            if (n_stages .eq. 0) then
                write(*,*) "Error (read_input.f90): Temperature defined before
                the number of stages"
                stop
            endif
            do i=1,n_stages
                read(ninput,*) temp(i)
            enddo
        else
            exit
        end if
        read(ninput,*)  
    end do loopread

       

end subroutine
