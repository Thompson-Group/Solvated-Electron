subroutine read_data(data_filename,restart) 
    use kinds
    use common_variables
    use constants
    implicit none
    character(len=50) :: ignore, data_filename
    integer(kind=ip) :: i,j
    real(kind=dp) :: Q_tot
    real(kind=dp), allocatable, dimension(:) :: M_bytype
    logical :: restart
    open(unit=ndata, file=data_filename)
    read(ndata,*)
    read(ndata,*)
    read(ndata,*) n_atoms, ignore
    read(ndata,*) n_bonds, ignore
    read(ndata,*) n_angles, ignore
    read(ndata,*) n_dih, ignore
    read(ndata,*) n_imp, ignore
    !Allocate Arrays that rely on n_atoms, n_bonds, n_angles,n_dih, and n_imp
    allocate(a_id(n_atoms))
    allocate(mol_id(n_atoms))
    allocate(a_type(n_atoms))
    allocate(q(n_atoms))
    allocate(x(n_atoms))
    allocate(y(n_atoms))
    allocate(z(n_atoms))
    allocate(rx(n_atoms,n_atoms))
    allocate(ry(n_atoms,n_atoms))
    allocate(rz(n_atoms,n_atoms))
    allocate(qq(n_atoms,n_atoms))
    allocate(ee(n_atoms,n_atoms))
    allocate(ss(n_atoms,n_atoms))
    allocate(vx(n_atoms))
    allocate(vy(n_atoms))
    allocate(vz(n_atoms))
    allocate(fx_tot(n_atoms))
    allocate(fy_tot(n_atoms))
    allocate(fz_tot(n_atoms))
    allocate(fx_b(n_atoms))
    allocate(fy_b(n_atoms))
    allocate(fz_b(n_atoms))
    allocate(fx_v(n_atoms))
    allocate(fy_v(n_atoms))
    allocate(fz_v(n_atoms))
    allocate(fx_a(n_atoms))
    allocate(fy_a(n_atoms))
    allocate(fz_a(n_atoms))
    allocate(fx_c(n_atoms))
    allocate(fy_c(n_atoms))
    allocate(fz_c(n_atoms))
    allocate(fx_q(n_atoms))
    allocate(fy_q(n_atoms))
    allocate(fz_q(n_atoms))
    allocate(M(n_atoms))
    allocate(bond_table(n_bonds,3))
    allocate(angle_table(n_angles,4))
    allocate(dih_table(n_dih, 5))

    read(ndata,*)
    read(ndata,*) n_a_type, ignore, ignore
    IF (n_bonds.ne.0) THEN
        read(ndata,*) n_b_type, ignore, ignore
    END IF
    IF (n_angles.ne.0) THEN
        read(ndata,*) n_angle_type, ignore, ignore
    END IF
    IF (n_dih.ne.0) THEN
        read(ndata,*) n_dih_type, ignore, ignore
    END IF
    IF (n_imp.ne.0) THEN
        read(ndata,*) n_imp_type, ignore, ignore
    END IF

    !Allocates basted on n_a_type, n_b_type, n_ang_type, n_dih_type, n_imp_type
    allocate(M_bytype(n_a_type))
    allocate(ep(n_a_type))
    allocate(sig(n_a_type))
    allocate(k_r(n_b_type))
    allocate(req(n_b_type))
    allocate(k_ang(n_angle_type))
    allocate(theta_eq(n_angle_type))
    !Reads in box size
    read(ndata,*)
    read(ndata,*) xlo, xhi, ignore, ignore
    read(ndata,*) ylo, yhi, ignore, ignore
    read(ndata,*) zlo, zhi, ignore, ignore
    read(ndata,*)
    Lx = 2.0_dp * xhi
    Ly = 2.0_dp * yhi
    Lz = 2.0_dp * zhi
    !Reads in the Masses
    read(ndata,*)
    read(ndata,*)
    DO i = 1, n_a_type
        read(ndata,*) ignore, M_bytype(i)
        M_bytype(i) = M_bytype(i)*mass_conv
    END DO

    !Reads in the Pair Coeffs
    read(ndata,*)
    read(ndata,*)
    read(ndata,*)
    DO i = 1, n_a_type
        read(ndata,*) ignore, ep(i), sig(i)
    END DO
    read(ndata,*)
    !Reads in the Atoms Section
    read(ndata,*)
    read(ndata,*)
    DO i = 1, n_atoms
        read(ndata,*) a_id(i), mol_id(i), a_type(i), q(i), x(i), y(i), z(i)
    END DO

    !Set the mass for each atom
    do i = 1, n_atoms
       M(i) = M_bytype(a_type(i))
    enddo
    deallocate(M_bytype)

    !Reads in the Bonds Section
    IF (n_bonds.ne.0) THEN
        read(ndata,*)
        read(ndata,*)
        read(ndata,*)
        DO i = 1, n_bonds
            read(ndata,*) ignore, bond_table(i,1), bond_table(i,2), bond_table(i,3)
        END DO
    END IF
    !Reads in the Angles Section
    IF (n_angles.ne.0) THEN
        read(ndata,*) 
        read(ndata,*) 
        read(ndata,*) 
        DO i = 1, n_angles
            read(ndata,*) ignore, angle_table(i,1), angle_table(i,2), angle_table(i,3), angle_table(i,4)
        END DO
    END IF
    !Reads in the Dihedrals Section
    IF (n_dih.ne.0) THEN
        read(ndata,*)
        read(ndata,*)
        read(ndata,*)
        DO i = 1, n_dih
            read(ndata,*) ignore, dih_table(i,1), dih_table(i,2), dih_table(i,3), dih_table(i,4), dih_table(i,5)
        END DO
    END IF
    IF (n_imp.ne.0) THEN
        write(*,*) "Error: Identified impropers, when code does not support them."
    END IF
    IF (restart) THEN
        DO i = 1, n_atoms
            read(ndata,*) ignore, vx(i), vy(i), vz(i)
        END DO

    END IF

    DO i=1,n_atoms
        DO j=1, n_atoms
            IF (mol_id(i).ne.mol_id(j)) THEN
                qq(i,j)=q(i)*q(j)
                ee(i,j)=sqrt(ep(a_type(i))*ep(a_type(j)))
                ss(i,j)=(sig(a_type(i))+sig(a_type(j)))*0.5
            ELSE
                qq(i,j)=0.0d0
                ee(i,j)=0.0d0
                ss(i,j)=0.0d0
            END IF
        END DO
    END DO

    DO i=1, n_atoms
        Q_tot = Q_tot+q(i)
    END DO
    if (Q_tot .ne. 0.0) then
        write(*,*) "Error: Box has total charge that is non-zero."
    end if
    close(unit=ndata)
end subroutine
