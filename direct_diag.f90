!     Subroutine to Calculate the Eigenvalues and Eigenvectors of
!       the Hamiltonian Using Direct Diagonalization
!
      subroutine direct_diag(Eigdir,Vecdir)
      
      use kinds
      use quantum_variables
      implicit none

      integer(kind=ip) :: n, m, i, j, info
      double precision, dimension(ng,ng) :: h
      double precision, dimension(ng) :: vect, vectp1
      double precision, dimension(ng) :: Eigdir
      double precision, dimension(ng,ng) :: Vecdir
      double precision, dimension(3*ng) :: work

!     Construct the full Hamiltonian by matrix multiplication

      do i = 1, ng
         vect = 0d0
         vect(i) = 1d0

         call hvec3d(vect, vectp1)

         do j = 1, ng
            h(i,j) = vectp1(j)
         enddo
      enddo

!     Call the lapack routine to diagonalize

      call dsyev('V', 'L', ng, h, nd, Eigdir, work, 3*nd, info)

      do i = 1, 20
         write(65,*) i, Eigdir(i)
      enddo

      do i = 1, ng
         do j = 1, ng
            Vecdir(i,j) = h(i,j)
         enddo
      enddo

      end subroutine direct_diag
