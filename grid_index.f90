!     Module to hold all the grid indices

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


!     Module to hold all the hamiltonian pieces

      module ham_data
      
         use kinds
         implicit none

      end module ham_data


!     Module to define the constants

      module constants

        use kinds
        implicit none
        real(kind=dp), parameter :: evperau=27.2113961_dp
        real(kind=dp), parameter :: kb=3.1668294e-6_dp  
        real(kind=dp), parameter :: cmperau=5.29177249e-9_dp
        real(kind=dp), parameter :: secperau=2.41889e-17_dp
        real(kind=dp), parameter :: angperau=0.529177249_dp
        real(kind=dp), parameter :: psperau=2.41889e-5_dp
        real(kind=dp), parameter :: fsperau=2.41889e-2_dp
        real(kind=dp), parameter :: cmiperau=2.195e5_dp
        real(kind=dp), parameter :: kjperkcal=4.184_dp
        real(kind=dp), parameter :: kcalperev=23.0605_dp
        real(kind=dp), parameter :: kgperau=9.109e-31_dp
        real(kind=dp), parameter :: atmperau=2.903622e8_dp
        real(kind=dp), parameter :: kcalperau=kcalperev*evperau
        real(kind=dp), parameter :: avogadro=6.02205e23_dp
        real(kind=dp), parameter :: dmperau=0.529177249e-9_dp
        real(kind=dp), parameter :: aupergmol=1822.992667_dp
        complex(kind=dp), parameter :: aye=(0.0_dp,1.0_dp)
         
      end module constants
