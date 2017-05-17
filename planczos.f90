!     Subroutine to Calculate the Eigenvalues and Eigenvectors of
!       the Hamiltonian Using a Lanczos Algorithm.
!
      subroutine planczos(ng,niter,actiter,tol,n_eigtol,try,
     1    Eigval,Eigvec,Krylov_vectors)
      
      use sizes
      use ham_data
      implicit none
      integer :: ng, niter, n, m, i, j, info
      integer :: actiter, n_eigtol
      double precision, dimension(0:maxiter,nd) :: Krylov_vectors
      double precision, dimension(0:maxiter-1,0:maxiter) :: h
      double precision, dimension(nd) :: vect, vectp1, try
      double precision ::  y, len, dummy, tol, eigval_old, error
      double precision, dimension(maxiter) :: Eigvals, Eigval
      double precision, dimension(maxiter,maxiter) :: Areal, Eigvec
      double precision, dimension(3*maxiter) :: work
!      external          hvec3d, dsyev

      eigval_old = 0d0
      
      do i = 1, ng 
         Krylov_vectors(0,i) = try(i)
      end do
      
      len = 0.0d0
      do i = 1, ng
         len = len + Krylov_vectors(0,i)*Krylov_vectors(0,i)
      end do
      
      len = dsqrt(len)

      do i = 1, ng
         Krylov_vectors(0,i) = Krylov_vectors(0,i) / len
      end do
      
      do n = 1, niter
         
         do i = 1, ng
            vect(i) = Krylov_vectors(n-1,i)
         end do

         call hvec3d(vect, vectp1, ng)
         
         do i = 1, ng
            Krylov_vectors(n,i) = vectp1(i)
         end do
         
         do m = 0, n-1
            
!     Calculate matrix            
            
            h(m, n-1) = 0.0d0

            do i = 1, ng
               h(m,n-1) = h(m,n-1) + Krylov_vectors(m,i)*
     &            Krylov_vectors(n,i)
            end do
            
            h(n-1, m) = h(m, n-1)
            
         end do
         
!     Orthogonalize vectors
         
         do m = 0, n-1
            
            y = 0.0d0
            
            do i = 1, ng
               y = y + Krylov_vectors(m,i) * 
     &            Krylov_vectors(n,i)
            end do
            
            do i = 1, ng
               Krylov_vectors(n,i) = Krylov_vectors(n,i) - 
     &            y * Krylov_vectors(m,i)
            end do
            
         end do
         
!     Normalize vector
         
         len = 0.0d0
         do i = 1, ng
            len = len + Krylov_vectors(n,i)*Krylov_vectors(n,i)
         end do
         
         len = dsqrt(len)
         
         do i = 1, ng
            Krylov_vectors(n,i) = Krylov_vectors(n,i) / len
         end do
         
         do i = 0, n-1
            do j = 0, n-1
               Areal(j+1,i+1) = h(j,i)
            end do
         end do
         
         if (n.eq.niter.or.n.gt.n_eigtol) then
            call dsyev('V', 'L', n, Areal, maxiter, Eigvals,
     1           work, 3*maxiter, info)
!            write(45,*) n, Eigvals(1)
!            write(46,*) n, Eigvals(2)
!            write(47,*) n, Eigvals(3)
!            write(48,*) n, Eigvals(4)
            do i = n,1,-1
!               write(45,*) i,Eigvals(i)
               Eigval(i) = Eigvals(i)
            enddo
            do i = 1,n
               do j = 1,n
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
!      do i = 1, 20
!         write(45,*) i, Eigval(i) 
!      enddo
!      write(nout,*) ' Number of Lanczos Iterations = ',actiter

      return
      end subroutine planczos
