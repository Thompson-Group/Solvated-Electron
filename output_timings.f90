!   Subroutine to print out the summary of timings for the different pieces of the code
!
  subroutine output_timings

    use timings

    implicit none

    write(6,'(A,F10.1,A,F10.2,A)') ' Bond time     = ',tbonds,' s, = ',tbonds/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' Angle time    = ',tang,' s, = ',tang/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' vdWaals time  = ',tvdw,' s, = ',tvdw/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' Coulomb time  = ',tcoul,' s, = ',tcoul/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' Forces time   = ',tforces,' s, = ',tforces/60_dp,' min'
    write(6,*)
    write(6,'(A)') ' Quantum Components'
    write(6,'(A,F10.1,A,F10.2,A)') ' Grid time     = ',tgrid,' s, = ',tgrid/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' Pseudo time   = ',tpseudo,' s, = ',tpseudo/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' QM Force time = ',tqm_forces,' s, = ',tqm_forces/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' Lanczos time  = ',tlanc,' s, = ',tlanc/60_dp,' min'
    write(6,'(A,F10.1,A,F10.2,A)') ' Avgs time     = ',tavgs,' s, = ',tavgs/60_dp,' min'

  end subroutine output_timings
