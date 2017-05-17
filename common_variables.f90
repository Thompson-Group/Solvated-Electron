Module common_variables
!**************************************************************************************************************
! This module contains arrays and scalars that are common to several subroutines in this md code.
! These variables do not need to be declared in other subroutines. 
!Zeke will allocate the arrays and read in values for the variables in the read_data code
! use this module as needed in other subroutines
!Written by Mesele and Pubudu
! Hackathon Thursday January 12, 2017
!**************************************************************************************************************

  use kinds
  Implicit none

  integer(kind=ip) :: n_atoms, n_bonds, n_angles, n_dih, n_imp, coul_flag
  integer(kind=ip) :: n_a_type, n_b_type, n_angle_type, n_dih_type, n_imp_type
  real(kind=dp) :: xlo, xhi, ylo, yhi, dt, zlo, zhi, Lx, Ly, Lz, r_cut, Temp, del, alpha, v_a, v_b, v_c, v_v, v_tot
  integer(kind=ip), allocatable, dimension(:)   :: a_id, mol_id, a_type
  real(kind=dp), allocatable, dimension(:,:) :: rx, ry, rz, bond_table, angle_table, ee, ss, qq, dih_table
  real(kind=dp), allocatable, dimension(:) :: fx_v, fy_v, fz_v, fx_c, fy_c, fz_c, fx_tot, fy_tot, fz_tot
  real(kind=dp), allocatable, dimension(:) :: M, q, vx, vy, vz, ep, sig
  real(kind=dp), allocatable, dimension(:) :: x, y, z, fx_b, fy_b, fz_b, fx_a, fy_a, fz_a, k_r, req, k_ang, theta_eq
  character(len=50) :: bond_style
end module

!**************************************************************************************************************
! This module contains arrays and scalars that are common to several subroutines in this solve code.
! These variables do not need to be declared in other subroutines. 
! Written by Ward Thompson
! Hackathon Wednesday May 17, 2017
!**************************************************************************************************************

module quantum_variables

    use kinds
    implicit none
    integer(kind=ip) :: ng
    integer(kind=ip), allocatable, dimension(:) :: inxgridx, inygridx, inzgridx
    integer(kind=ip), allocatable, dimension(:) :: inxgridy, inygridy, inzgridy
    integer(kind=ip), allocatable, dimension(:) :: inxgridz, inygridz, inzgridz
    integer(kind=ip), allocatable, dimension(:) :: indx, irankyx, irankzx
    integer(kind=ip), allocatable, dimension(:) :: isx, ifx, isy, ify, isz, ifz

    real(kind=dp), allocatable, dimension(:) :: v_e
    real(kind=dp), allocatable, dimension(:,:) :: rg_e
    real(kind=dp), allocatable, dimension(:,:) :: kex, key, kez

end module quantum_variables



Module constants
!**************************************************************************************************************
! This module contains constants and the numbers used for output files.
! These variables do not need to be declared in other subroutines. 
! use this module as needed in other subroutines
!Written by Ward Thompson
! Hackathon Friday January 13, 2017
!**************************************************************************************************************

  use kinds
  Implicit none

  integer, parameter :: nxyz=20, nthermo=21, ndata=30, ninput=31, nrest=22
  real(kind=dp) :: pi = 4.0_dp*atan(1.0_dp)
  real(kind=dp), parameter :: kb=0.0019872041_dp
!  real(kind=dp), parameter :: mass_conv=2.390057361e-7_dp 
  real(kind=dp), parameter :: mass_conv=1.0_dp/4.184e-4_dp 
  real(kind=dp), parameter :: C_coul=2.40e-4_dp 

end module
