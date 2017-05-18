!     Subroutine to Calculate the Eigenvalues and Eigenvectors of
!       the Hamiltonian Using a Lanczos Algorithm.
!
      subroutine planczos(try)
      
!need to add actiter, niter, tol, n_eigtol, Eigval, Eigvec, Krylov_vectors to quantum_variables      
      use quantum_variables
      implicit none

      integer(kind=ip) :: n, m, i, j, info
      real(kind=dp), dimension(ng) :: vect, vectp1, try
      real(kind=dp) ::  y, len, dummy, tol, eigval_old, error
      real(kind=dp), dimension(niter) :: Eigvals
      real(kind=dp), dimension(0:niter,0:niter) :: h
      real(kind=dp), dimension(niter,niter) :: Areal
      real(kind=dp), dimension(3*niter) :: work

      eigval_old = 0.0_dp
      
      ! Set the 0th Krylov vector to the initial guess
      do i = 1, ng 
         Krylov_vectors(0,i) = try(i)
      end do
      
      ! Normalize the 0th Krylov vector
      len = 0.0_dp
      do i = 1, ng
         len = len + Krylov_vectors(0,i)*Krylov_vectors(0,i)
      end do
      
      len = dsqrt(len)

      do i = 1, ng
         Krylov_vectors(0,i) = Krylov_vectors(0,i) / len
      end do
      
      ! Generate the Krylov space by successive applications of the Hamiltonian
      do n = 1, niter
         
         do i = 1, ng
            vect(i) = Krylov_vectors(n-1,i)
         end do

         call hvec3d(vect, vectp1)
         
         do i = 1, ng
            Krylov_vectors(n,i) = vectp1(i)
         end do
         
         ! Form the Hamiltonian matrix in the Krylov basis
         do m = 0, n-1
            
            h(m, n-1) = 0.0_dp

            do i = 1, ng
               h(m,n-1) = h(m,n-1) + Krylov_vectors(m,i)*Krylov_vectors(n,i)
            end do
            
            h(n-1, m) = h(m, n-1)
            
         end do
         
!     Orthogonalize vectors
         
         do m = 0, n-1
            
            y = 0.0_dp
            
            do i = 1, ng
               y = y + Krylov_vectors(m,i)*Krylov_vectors(n,i)
            end do
            
            do i = 1, ng
               Krylov_vectors(n,i) = Krylov_vectors(n,i) - y * Krylov_vectors(m,i)
            end do
            
         end do
         
!     Normalize vector
         
         len = 0.0_dp
         do i = 1, ng
            len = len + Krylov_vectors(n,i)*Krylov_vectors(n,i)
         end do
         
         len = dsqrt(len)
         
         do i = 1, ng
            Krylov_vectors(n,i) = Krylov_vectors(n,i) / len
         end do
         
         ! Diagonalize the Krylov space Hamiltonian 

         do i = 0, n - 1
            do j = 0, n - 1
               Areal(j+1,i+1) = h(j,i)
            end do
         end do
         
         if (n.eq.niter.or.n.gt.n_eigtol) then
            call dsyev('V', 'L', n, Areal, maxiter, Eigvals,work, 3*maxiter, info)
            do i = n, 1, -1
               Eigval(i) = Eigvals(i)
            enddo
            do i = 1, n
               do j = 1, n
                  Eigvec(i,j) = Areal(i,j)
               enddo
            enddo
            error = dabs((Eigval(n_eigtol)-eigval_old)/Eigval(n_eigtol))
            if(error.lt.tol) goto 10
            eigval_old = Eigval(n_eigtol)
         end if
      end do
            
 10   continue
      
      actiter = n

      end subroutine planczos
