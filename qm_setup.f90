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
     integer(kind=ip) :: i
     real(kind=dp) , dimension(ng) :: try

     call grid

     call kinetic
     call rankyx
     call rankzx
     call looplims
!     call graphgrid
     write(6,*) ' Made it through grid indexing'
     try = 1.0_dp/(sqrt(real(ng)))

!     write(6,*) ' About to stop'
!     stop
!     write(6,*) ' Did not stop!'
     call planczos(try)
     !call direct_diag
     write(6,*) ' actiter = ',actiter
     do i = 1, eig_tol
        write(6,*) i, Eigval(i)*evperau
     enddo

     call qm_forces

!     call qm_force_check
!     stop

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
        allocate(Eigvec(niter,niter))
        allocate(Eigval(niter))
        allocate(Krylov_vectors(0:niter,nfull))
        allocate(indxx(nfull))
     endif
     
   end subroutine qm_allocation
