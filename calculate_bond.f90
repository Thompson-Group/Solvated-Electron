subroutine calculate_bond
    use kinds
    use common_variables
    implicit none
    ! Variable declarations
    !Loop Integers, Indices, etc.
    integer(kind=ip) :: i, b_type, b_atom_1, b_atom_2
    !Arrays for storing information, all of n_bond length.
    real(kind=dp), dimension(n_bonds):: rinst
    !Array read in from read_data.f90, stores bonding information.
    real(kind=dp) :: u_bond, r_sep
       
    !Initialize
    v_b = 0.0_dp

    !Bond Energy and Force Calculations
    IF (bond_style.eq."harmonic") THEN
        ! Do loop over number of bonds
        DO i = 1, n_bonds
            !Index calculations
            b_type  =bond_table(i,1)
            b_atom_1=bond_table(i,2)
            b_atom_2=bond_table(i,3)
            !Bond Length Calculations
            !rinst is the instantaneous bond length
            !r_sep is rinst-req
            rinst(i)=sqrt((rx(b_atom_1,b_atom_2)**2) &
                        + (ry(b_atom_1,b_atom_2)**2) &
                        + (rz(b_atom_1,b_atom_2)**2) )
            r_sep=rinst(i)-req(b_type)
            !Potential Energy Calculation
            u_bond=0.5*k_r(b_type)*r_sep**2
            v_b = v_b + u_bond
            write(*,"(3i2)") i , b_atom_1 , b_atom_2
            write(*,"(i2,3f10.5)") i , rinst(i) , u_bond , v_b

            !Force calculation 
            !Calculates the force on atom 1 as:
            !F_1 = -dU/dr * dr/dx
            !Note: du/dr = k*r_sep
            !and dr/dx = x_1/rinst(i)
            !Calculates the force on atom 2 as:
            !F_2 = -F_1

            fx_b(b_atom_1) = - k_r(b_type) * r_sep * rx(b_atom_1,b_atom_2) / rinst(i)
            fx_b(b_atom_2) = - fx_b(b_atom_1)
            fy_b(b_atom_1) = - k_r(b_type) * r_sep * ry(b_atom_1,b_atom_2) / rinst(i)
            fy_b(b_atom_2) = - fy_b(b_atom_1)
            fz_b(b_atom_1) = - k_r(b_type) * r_sep * rz(b_atom_1,b_atom_2) / rinst(i)
            fz_b(b_atom_2) = - fz_b(b_atom_1)
        END DO
    END IF


end subroutine
