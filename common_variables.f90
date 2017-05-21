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
  real(kind=dp) :: xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz, r_cut, Temp, delta, alpha, v_a, v_b, v_c, v_v, v_q, v_tot
  integer(kind=ip), allocatable, dimension(:)   :: a_id, mol_id, a_type
  real(kind=dp), allocatable, dimension(:,:) :: rx, ry, rz, bond_table, angle_table, ee, ss, qq, dih_table
  real(kind=dp), allocatable, dimension(:) :: fx_v, fy_v, fz_v, fx_c, fy_c, fz_c, fx_tot, fy_tot, fz_tot
  real(kind=dp), allocatable, dimension(:) :: M, q, vx, vy, vz, ep, sig
  real(kind=dp), allocatable, dimension(:) :: x, y, z, fx_b, fy_b, fz_b, fx_a, fy_a, fz_a, k_r, req, k_ang, theta_eq
  real(kind=dp), allocatable, dimension(:) :: fx_q, fy_q, fz_q
  integer(kind=ip) :: n_stages
  character(len=50) :: bond_style, run_style, nvt_type
  character(len=2), allocatable, dimension(:) :: elements
  real(kind=dp) :: dt, nu
  integer(kind=ip) :: df_xyz, df_thermo, df_rest, nstep, nvt_freq
  integer(kind=ip) :: fc_cnt, qfc_cnt

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
    integer(kind=ip) :: ng, niter, actiter, Nst
    integer(kind=ip), allocatable, dimension(:) :: inxgridx, inygridx, inzgridx
    integer(kind=ip), allocatable, dimension(:) :: inxgridy, inygridy, inzgridy
    integer(kind=ip), allocatable, dimension(:) :: inxgridz, inygridz, inzgridz
    integer(kind=ip), allocatable, dimension(:) :: indxx, irankyx, irankzx
    integer(kind=ip), allocatable, dimension(:) :: isx, ifx, isy, ify, isz, ifz

    real(kind=dp), allocatable, dimension(:) :: v_e, Eigval
    real(kind=dp), allocatable, dimension(:,:) :: rg_e, Eigvec, Krylov_vectors
    real(kind=dp), allocatable, dimension(:,:) :: kex, key, kez, fg_ex, fg_ey, fg_ez
    real(kind=dp), dimension(3) :: r_e_avg
    real(kind=dp) :: xmin, xmax, del, tol, vcut, v_e_avg
    integer(kind=ip) :: nraw, eig_tol
    real(kind=dp) :: r2_e_avg

end module quantum_variables

module input_variables
    use kinds
    character(len=50), dimension(5) :: srun_style, snvt_type
    real(kind=dp), dimension(5) :: sdt, stemp, snu
    integer(kind=ip), dimension(5) :: sdf_xyz, sdf_thermo, sdf_rest, snstep, snvt_freq

end module input_variables


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
  integer, parameter :: neigs=41, nrg=42, nlanc=43
  real(kind=dp) :: pi = 4.0_dp*atan(1.0_dp)
  real(kind=dp), parameter :: kb=0.0019872041_dp
!  real(kind=dp), parameter :: mass_conv=2.390057361e-7_dp 
  real(kind=dp), parameter :: mass_conv=1.0_dp/4.184e-4_dp 
  real(kind=dp), parameter :: C_coul=2.40e-4_dp 
  real(kind=dp), parameter :: me=1.0_dp
  real(kind=dp), parameter :: angperau=0.529177249_dp
  real(kind=dp), parameter :: kcalperev=23.0605_dp
  real(kind=dp), parameter :: evperau=27.2113961_dp  
  real(kind=dp), parameter :: kcalperau=kcalperev*evperau

end module


Module pseudo_constants
!**************************************************************************************************************
! This module contains constants and the numbers used for output files.
! These variables do not need to be declared in other subroutines.
! use this module as needed in other subroutines
!Written by Ward Thompson
! Hackathon Friday May 19, 2017                                                                              
!**************************************************************************************************************  

  use kinds
  Implicit none

  ! Parameters are for the Turi-Borgis potential [JCP 117, 6186 (2002)] and currently in atomic units
!  real(kind=dp),parameter :: alpha_pol=9.7446_dp
  real(kind=dp),parameter :: alpha_pol=0.0_dp
  real(kind=dp),parameter :: qo  = -0.820_dp, qh  = 0.410_dp
  real(kind=dp),parameter :: A1o = 0.575_dp, A1h = 0.750_dp
  real(kind=dp),parameter :: B1o = 0.620_dp, B1h = 0.150_dp
  real(kind=dp),parameter :: B2o = 1.000_dp, B2h = 0.500_dp
  real(kind=dp),parameter :: B3o = 0.400_dp, B3h = 0.350_dp
  real(kind=dp),parameter :: C1o = 4.400_dp

end module


Module timings
!**************************************************************************************************************
! This module contains variables to track timings of different pieces of the code
! These variables do not need to be declared in other subroutines.
! use this module as needed in other subroutines
!Written by Ward Thompson
! May 20, 2017                                                                              
!**************************************************************************************************************  

  use kinds
  Implicit none

  real(kind=dp) :: tforces, tverlet
  real(kind=dp) :: tbonds, tang, tvdw, tcoul, tdists
  real(kind=dp) :: tgrid, tpseudo, tlanc, tqm_forces, tavgs

end module timings
