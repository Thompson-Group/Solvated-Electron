! subroutine to calculate the pseudo potential for the solvated electron
! without the corresponding forces on the particle

  subroutine pseudo_e_tb(n_atoms,rx,ry,rz,re,Lx,Ly,Lz,rcut)
  use kinds
  use common_variables

! scalars
  real(kind=dp)    :: qo,qh,A1o,A1h,B1o,B1h,B2o,B2h,B3o,B3h,C1o
  real(kind=dp)    :: vtmp
  real(kind=dp),parameter :: alpha_pol=7.518E7,pi=3.41
  integer(kind=ip) :: k

! arrays
  real(kind=dp),allocatable,dimension(:,:) :: fpxi,fpyi,fpzi
  real(kind=dp),dimension(3)    :: L,rO,rH1,rH2,rOe,r1e,r2e 


  allocate(fpxi(n_atoms,3))
  allocate(fpyi(n_atoms,3))
  allocate(fpzi(n_atoms,3))

  qo  = -0.820
  qh  = 0.410
  A1o = 0.575/0.53
  A1h = 0.750/0.53
  B1o = 0.620/0.53
  B1h = 0.150/0.53
  B2o = 1.000/0.53
  B2h = 0.500/0.53
  B3o = 0.400/0.53
  B3h = 0.350/0.53
  C1o = 4.400*0.53

  vtmp = 0.0_dp
  fpxi = 0.0_dp
  fpyi = 0.0_dp
  fpzi = 0.0_dp

  do i = 1,n_atoms

        rO(1) = rx(i,1)
        rO(2) = ry(i,1)
        rO(3) = rz(i,1)

        rH1(1) = rx(i,2)
        rH1(2) = ry(i,2)
        rH1(3) = rz(i,2)

        rH2(1) = rx(i,3)
        rH2(2) = ry(i,3)
        rH2(3) = rz(i,3)

        L(1) = Lx
        L(2) = Ly
        L(3) = Lz

! wrap boundary conditions

        do k=1,3
                rO(k) = rO(k) - L(k)*int((rO(k) - re(k))/L(k))
                rH1(k) = rH1(k) - L(k)*int((rH1(k) - re(k))/L(k))
                rH2(k) = rH2(k) - L(k)*int((rH2(k) - re(k))/L(k))
        enddo

! calculate the distances
        dOe = 0_dp; d1e = 0_dp; d2e = 0_dp 
  
        do k=1,3
                rOe(k) = re(k) - rO(k)
                r1e(k) = re(k) - rH1(k)
                r2e(k) = re(k) - rH2(k)
                dOe = dOe + rOe(k)**2
                d1e = d1e + r1e(k)**2
                d2e = d2e + r2e(k)**2
        enddo
        dOe = sqrt(dOe)
        d1e = sqrt(d1e)
        d2e = sqrt(d2e)

        if (dOe.le.rcut) then
! calculate the polarization part
 
                vtmp = vtmp - (0.5_dp*alpha/(dOe**2+C1o**2)**2)
! the force looks like
                fpxi(i,1) = fpxi(i,1) + (2_dp*alpha/(rOe(1)**2+C1o**2)**3)
                fpyi(i,1) = fpyi(i,2) + (2_dp*alpha/(rOe(2)**2+C1o**2)**3)
                fpzi(i,1) = fpzi(i,3) + (2_dp*alpha/(rOe(3)**2+C1o**2)**3)
                
! the electrostatics

! oxygen atom
  
                e1 = erf(Alo*dOe)  
                e2 = erf(B2o*dOe)  
                e3 = erf(B3o*dOe)  
  
                vtmp = vtmp - (e1*qo/dOe + (e2 - e3)*B1o/dOe)
                fpxi(i,1) = fpxi(i,1) + qo*rOe(1)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1o*exp(-(A1o*dOe)**2)/dOe**2)      & 
                 + B1o*rOe(1)*(2_dp/sqrt(pi))*((B2o*exp(-(B2o*dOe)**2))&
                 - B3o*exp(-(B3o*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3 

                fpyi(i,1) = fpxi(i,1) + qo*rOe(2)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1o*exp(-(A1o*dOe)**2)/dOe**2)      &             
                 + B1o*rOe(1)*(2_dp/sqrt(pi))*((B2o*exp(-(B2o*dOe)**2))&
                 - B3o*exp(-(B3o*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3  

                fpzi(i,1) = fpxi(i,1) + qo*rOe(3)*(e1/dOe**3             &
                 - (2_dp/sqrt(pi))*A1o*exp(-(A1o*dOe)**2)/dOe**2)      &             
                 + B1o*rOe(1)*(2_dp/sqrt(pi))*((B2o*exp(-(B2o*dOe)**2))&
                 - B3o*exp(-(B3o*dOe)**2))/dOe**2                        &
                 - (e2 - e3)/dOe**3  

! H1 atom
                e1 = erf(Alh*d1e)
                e2 = erf(B2h*d1e)
                e3 = erf(B3h*d1e)

                vtmp = vtmp - (e1*qh/d1e + (e2 - e3)*B1h/d1e)
                fpxi(i,2) = fpxi(i,2) + qh*r1e(1)*(e1/d1e**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*d1e)**2)/d1e**2)      &             
                 + B1h*r1e(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*d1e)**2))&
                 - B3h*exp(-(B3h*d1e)**2))/d1e**2                        &
                 - (e2 - e3)/d1e**3  

                fpyi(i,2) = fpxi(i,2) + qh*r1e(2)*(e1/d1e**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*d1e)**2)/d1e**2)      &
                 + B1h*r1e(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*d1e)**2))&
                 - B3h*exp(-(B3h*d1e)**2))/d1e**2                        &
                 - (e2 - e3)/d1e**3

                fpzi(i,2) = fpxi(i,2) + qh*r1e(3)*(e1/d1e**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*d1e)**2)/d1e**2)      &
                 + B1h*r1e(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*d1e)**2))&
                 - B3h*exp(-(B3h*d1e)**2))/d1e**2                        &
                 - (e2 - e3)/d1e**3

! H2 atom
                e1 = erf(A1h*d2e)
                e2 = erf(B2h*d2e)
                e3 = erf(B3h*d2e)

                vtmp = vtmp - (e1*qh/d2e + (e2 - e3)*B1h/d2e)
                fpxi(i,3) = fpxi(i,3) + qh*r2e(1)*(e1/d2e**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*d2e)**2)/d2e**2)      &
                 + B1h*r1e(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*d2e)**2))&
                 - B3h*exp(-(B3h*d2e)**2))/d2e**2                        &
                 - (e2 - e3)/d2e**3

                fpxi(i,3) = fpxi(i,3) + qh*r2e(2)*(e1/d2e**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*d2e)**2)/d2e**2)      &
                 + B1h*r1e(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*d2e)**2))&
                 - B3h*exp(-(B3h*d2e)**2))/d2e**2                        &
                 - (e2 - e3)/d2e**3

                fpxi(i,3) = fpxi(i,3) + qh*r2e(3)*(e1/d2e**3             &
                 - (2_dp/sqrt(pi))*A1h*exp(-(A1h*d2e)**2)/d2e**2)      &
                 + B1h*r1e(1)*(2_dp/sqrt(pi))*((B2h*exp(-(B2h*d2e)**2))&
                 - B3h*exp(-(B3h*d2e)**2))/d2e**2                        &
                 - (e2 - e3)/d2e**3

        endif 
  enddo

return
end subroutine pseudo_e_tb
