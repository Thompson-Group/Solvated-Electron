!***********************************************************************************
!     Subroutine to set up the quantum calculation.  
!
!  Written by A. Katiyar, P. Wimalasiri, and O. Mesele
!***********************************************************************************
   subroutine qmsetup

     use kinds
     use common_variables
     use quantum_variables
     use constants

     implicit none
     real(kind=dp) , dimension(ng) :: try

     call grid

     allocate(Eigvec(niter,niter))
     allocate(Eigval(niter))
     allocate(Krylov_vectors(0:niter,ng))
     allocate(indxx(ng))

     call kinetic
     call rankyx
     call rankzx
     call looplims
     call graphgrid
     try = 1.0_dp/(sqrt(real(ng)))

     call planczos(try)
     !call direct_diag

   end subroutine qmsetup


! 
!  Subroutine to check for quantum mechanical stages and allocate any needed arrays

   subroutine qm_allocation

     use common_variables
     use quantum_variables
     use input_variables

     implicit none
     integer(kind=ip) :: istage, nfull
     logical :: quantum

     ! First check for quantum stages

     quantum = .false.
     do istage = 1, n_stages
        if(trim(srun_style(istage)).eq.'qm_nve' .or. trim(srun_style(istage)).eq.'qm_nvt' ) then
           quantum = .true.
        endif
     enddo

     ! If there is a quantum stage (or stages) allocate arrays

     if (quantum) then
        nfull = nraw**3
        allocate(fg_ex(nfull,n_atoms)); allocate(fg_ey(nfull,n_atoms)); allocate(fg_ez(nfull,n_atoms))
        allocate(v_e(nfull)); allocate(rg_e(nfull,3))
        allocate(inxgridx(nfull)); allocate(inygridx(nfull)); allocate(inzgridx(nfull))
        allocate(inxgridy(nfull)); allocate(inygridy(nfull)); allocate(inzgridy(nfull))
        allocate(inxgridz(nfull)); allocate(inygridz(nfull)); allocate(inzgridz(nfull))
        allocate(indxx(nfull)); allocate(irankyx(nfull)); allocate(irankzx(nfull))
        allocate(isx(nfull)); allocate(ifx(nfull)); allocate(isy(nfull)); allocate(ify(nfull))
        allocate(isz(nfull)); allocate(ifz(nfull))
        allocate(kex(nraw,nraw)); allocate(key(nraw,nraw)); allocate(kez(nraw,nraw))
     endif
     
   end subroutine qm_allocation
