subroutine write_restart()
    use kinds
    use common_variables
    implicit none
    integer(kind=ip) :: i
    real(kind=dp) :: m
    open(unit=nrest,file="traj.restart")
    

    !Dumps Lammps Header
    write(nrest,*) "Run Header"
    write(nrest,*)
    write(nrest,*) n_atoms, "atoms"
    write(nrest,*) n_bonds, "bonds"
    write(nrest,*) n_angles,"angles"
    write(nrest,*) n_dih, "dihedrals"
    write(nrest,*) n_imp, "impropers"
    write(nrest,*)
    write(nrest,*) n_a_type, "atom types"
    write(nrest,*) n_b_type, "bond types"
    write(nrest,*) n_angle_type, "angle types"
    if (n_dih_type .ne. 0) then
        write(nrest,*) n_dih_type, "dihedral types"
    end if
    if (n_imp_type .ne. 0) then
        write(nrest,*) "How did you get impropers into this?"
    end if
    write(nrest,*) 
    write(nrest,*) xlo,xhi "xlo xhi"
    write(nrest,*) ylo, yhi, "ylo yhi"
    write(nrest,*) zlo, zhi, "zlo zhi"
    write(nrest,*)
    ! Write out Masses
    write(nrest,*) "Masses"
    write(nrest,*)
    DO i=1,n_a_type
        m=M(i)*4.184*10**4
        write(nrest,*) i, m
    END DO
    write(nrest,*)
    ! Write out Pair Coeffs
    write(nrest,*) "Pair Coeffs"
    write(nrest,*)
    DO i=1,n_a_type
        write(nrest,*) i, ep(i), sig(i)
    END DO
    write(nrest,*)
    ! Writes the atom section
    write(nrest,*) "Atoms"
    write(nrest,*)
    DO i=1,n_atoms
        write(nrest,*) a_id(i), mol_id(i), a_type(i), q(i), x(i), y(i), z(i)
    END DO
    write(nrest,*)
    ! Writes the bonds section
    write(nrest,*) "Bonds"
    write(nrest,*)
    DO i=1, n_bonds
        write(nrest,*) i, bond_table(i,1), bond_table(i,2), bond_table(i,3)
    END DO
    write(nrest,*)
    ! Writes the angles section
    write(nrest,*) "Angles"
    write(nrest,*)
    DO i=1, n_angles
        write(nrest,*) i, angle_table(i,1), angle_table(i,2), angle_table(i,3), angle_table(i,4)
    END DO
    write(nrest,*)
    ! Writes the Dihedrals
    if (n_dih_type .ne. 0) then
        write(nrest,*) "Dihedrals"
        write(nrest,*)
        DO i=1, n_dih
            write(nrest,*) i, dih_table(i,1), dih_table(i,2), dih_table(i,3), dih_table(i,4), dih_table(i,5)
        END DO
        write(nrest,*)
    END IF
    ! Writes the impropers
    if (n_imp_type .ne. 0) then
        write(nrest,*) "Impropers"
        write(nrest,*)
        DO i=1, n_imp
            write(nrest,*) "huh? impropers? here?"
        END DO
    write(nrest,*)
    END IF
    ! Writes the Velocities
    write(nrest,*) "Velocities"
    write(nrest,*)
    DO i=1,n_atoms
        write(nrest,*) i, vx(i), vy(i), vz(i)
    END DO



end subroutine
