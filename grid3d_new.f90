!***********************************************************************************
!     Subroutine to set up the raw 3-D grid and calculate the potential 
!        The grid is truncated by a potential cuttoff.  
!        This routine has x, y, z all with the same raw dimensions.
!  Written by A. Katiyar, P. Wimalasiri, and O. Mesele
!***********************************************************************************
      subroutine grid

      use common_variables
      use quantum_variables
      use constants
      
       implicit none
!     --- Local Scalars and array ---
       integer(kind=ip) :: i, j, k, ii
       integer(kind=ip) :: ngx, ngy, ngz
       real(kind=dp), dimension(:), allocatable :: zg, xg, yg
!    ------semi - global array
       real(kind=dp), dimension(3) :: re
       real(kind=dp), dimension(n_atoms) :: fxtmp, fytmp, fztmp    
       real(kind=dp) :: vtmp
 
       allocate(zg(nraw)); allocate(yg(nraw)); allocate(xg(nraw))



      del=(xmax-xmin)/real(nraw)

!      write(nout,'(A,F12.5,A)') ' DVR dx = ',del,' angstroms'

      do i = 1, nraw
         xg(i) = xmin + dble(i-1)*del
     !    wx(i) = del
      enddo

      do i = 1, nraw
         yg(i) = xmin + dble(i-1)*del
      !   wy(i) = del
      enddo

      do i = 1, nraw
         zg(i) = xmin + dble(i-1)*del
       !  wz(i) = del
      enddo

!     Truncate the grid with a potential cuttoff Vcut with 
!       x running fastest, then y, then z.

      ngx = 0
      do k = 1, nraw
         re(3) = zg(k)
         do j = 1, nraw
            re(2) = yg(j)
            do i = 1, nraw
               re(1) = xg(i)
               call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)
               if(vtmp.le.vcut) then
!                 write(33,'(I4,4F12.5)') ngx+1, (re(ii),ii=1,3), vtmp
                  ngx = ngx + 1
                  fg_ex(ngx,:) = fxtmp(:)
                  fg_ey(ngx,:) = fytmp(:)
                  fg_ez(ngx,:) = fztmp(:)
                  v_e(ngx) = vtmp
                  rg_e(ngx,1) = xg(i)
                  rg_e(ngx,2) = yg(j)
                  rg_e(ngx,3) = zg(k)
                  inxgridx(ngx) = i
                  inygridx(ngx) = j
                  inzgridx(ngx) = k
               endif  ! cuttoff if
            enddo ! loop over x
         enddo ! loop over y
      enddo ! loop over z
      ng = ngx

!     Truncate the grid with a potential cuttoff Vcut with 
!       y running fastest, then x, then z.

      ngy = 0
      do k = 1, nraw
         re(3) = zg(k)
         do j = 1, nraw
            re(1) = xg(j)
            do i = 1, nraw
               re(2) = yg(i)
               call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)
               if(vtmp.le.vcut) then
                  ngy = ngy + 1
                  fg_ex(ngy,:) = fxtmp(:)
                  fg_ey(ngy,:) = fytmp(:)
                  fg_ez(ngy,:) = fztmp(:)
                  inygridy(ngy) = i
                  inxgridy(ngy) = j
                  inzgridy(ngy) = k
               endif  ! cuttoff if
            enddo ! loop over y
         enddo ! loop over x
      enddo ! loop over z

!     Truncate the grid with a potential cuttoff Vcut with 
!       z running fastest, then x, then y.

      ngz = 0
      do k = 1, nraw
         re(2) = yg(k)
         do j = 1, nraw
            re(1) = xg(j)
            do i = 1, nraw
               re(3) = zg(i)
               call pseudo_e_tb(re,vtmp,fxtmp,fytmp,fztmp)
               if(vtmp.le.vcut) then
                  ngz = ngz + 1
                  fg_ex(ngz,:) = fxtmp(:)
                  fg_ey(ngz,:) = fytmp(:)
                  fg_ez(ngz,:) = fztmp(:)
                  inzgridz(ngz) = i
                  inxgridz(ngz) = j
                  inygridz(ngz) = k
               endif  ! cuttoff if
            enddo ! loop over z
         enddo ! loop over x
      enddo ! loop over y
      
      !deallocate local arrays
      deallocate(zg); deallocate(yg); deallocate(xg)

      write(6,*) ' In grid: ng = ',ng,' ngy = ',ngy,' ngz = ',ngz
      v_e = v_e/kcalperau
      
      end subroutine grid
