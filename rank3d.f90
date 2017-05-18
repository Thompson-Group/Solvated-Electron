!
!
      subroutine rankyx
!
!
      use quantum_variables
      implicit   none

      integer(kind=ip) :: j, indxj
      logical, external :: comparyx

      call indexx(comparyx)

      do j = 1, ng
         indxj = indxx(j)
         irankyx(indxj) = j
      enddo

      end subroutine rankyx
 
!
!
      subroutine rankzx
!
!
      use quantum_variables
      implicit   none

      integer(kind=ip) :: j, indxj
      logical, external :: comparzx

      call indexx(comparzx)

      do j = 1, ng
         indxj = indxx(j)
         irankzx(indxj) = j
      enddo

      end subroutine rankzx

!
!
      logical function comparyx(i,j)
!
!
      use quantum_variables
      implicit   none

      integer(kind=ip) :: i,j

      if (inzgridy(i).gt.inzgridy(j)) then
         comparyx = .true.
      elseif (inzgridy(i).lt.inzgridy(j)) then
         comparyx = .false.
      else
         
         if (inygridy(i).gt.inygridy(j)) then
            comparyx = .true.
         elseif (inygridy(i).lt.inygridy(j)) then
            comparyx = .false.
         else

            comparyx = inxgridy(i).gt.inxgridy(j)
         endif
      endif
      
      end function comparyx

!
!
      logical function comparzx(i,j)
!
!
      use quantum_variables
      implicit   none

      integer(kind=ip) :: i,j

      if (inzgridz(i).gt.inzgridz(j)) then
         comparzx = .true.
      elseif (inzgridz(i).lt.inzgridz(j)) then
         comparzx = .false.
      else

         if (inygridz(i).gt.inygridz(j)) then
            comparzx = .true.
         elseif (inygridz(i).lt.inygridz(j)) then
            comparzx = .false.
         else
            
            comparzx = inxgridz(i).gt.inxgridz(j)
         endif
      endif

      end function comparzx

!
!
      subroutine looplims
!
!
      use quantum_variables
      implicit   none

      integer(kind=ip) :: n,nrev
      integer(kind=ip) :: lastix,lastiy,lastiz,adjixn,adjiyn,adjizn
      integer(kind=ip) :: ixn,iyn,izn,rowfin,rowst
      logical diffseg

!       Run the x grid forward to determine isx.

      do n = 1, ng
         ixn = inxgridx(n)
         iyn = inygridx(n)
         izn = inzgridx(n)
         if(n.eq.1) then
            isx(1) = 0
            rowst = 1
         else
            adjixn = ixn - 1
            diffseg = izn.ne.lastiz.or.iyn.ne.lastiy.or.adjixn.ne.lastix
            if(diffseg) then
               isx(n) = 0
               rowst = n
            else
               isx(n) = rowst - n
            endif
         endif
         lastix = ixn
         lastiy = iyn
         lastiz = izn
      enddo

!     Run the x grid backward to determine ifx.

      do nrev = 1, ng
         n = ng - nrev + 1
         ixn = inxgridx(n)
         iyn = inygridx(n)
         izn = inzgridx(n)
         if(n.eq.ng) then
            ifx(ng) = 0
            rowfin = ng
         else
            adjixn = ixn + 1
            diffseg = izn.ne.lastiz.or.iyn.ne.lastiy.or.adjixn.ne.lastix
            if(diffseg) then
               ifx(n) = 0
               rowfin = n
            else
               ifx(n) = rowfin - n
            endif
         endif
         lastix = ixn
         lastiy = iyn
         lastiz = izn
      enddo

!     Run the y grid forward to determine isy
     
      do n = 1, ng
         ixn = inxgridy(n)
         iyn = inygridy(n)
         izn = inzgridy(n)
         if(n.eq.1) then
            isy(1) = 0
            rowst = 1
         else
            adjiyn = iyn - 1
            diffseg = izn.ne.lastiz.or.ixn.ne.lastix.or.adjiyn.ne.lastiy
            if(diffseg) then
               isy(n) = 0
               rowst = n
            else
               isy(n) = rowst - n
            endif
         endif
         lastix = ixn
         lastiy = iyn
         lastiz = izn
      enddo

!     Run the y grid backward to determine ify.

      do nrev = 1, ng
         n = ng - nrev + 1
         ixn = inxgridy(n)
         iyn = inygridy(n)
         izn = inzgridy(n)
         if(n.eq.ng) then
            ify(ng) = 0
            rowfin = ng
         else
            adjiyn = iyn + 1
            diffseg = izn.ne.lastiz.or.ixn.ne.lastix.or.adjiyn.ne.lastiy
            if(diffseg) then
               ify(n) = 0
               rowfin = n
            else
               ify(n) = rowfin - n
            endif
         endif
         lastix = ixn
         lastiy = iyn
         lastiz = izn
      enddo

!     Run the z grid forward to find isz.

      do n = 1, ng
         ixn = inxgridz(n)
         iyn = inygridz(n)
         izn = inzgridz(n)
         if(n.eq.1) then
            isz(1) = 0
            rowst = 1
         else
            adjizn = izn - 1
            diffseg = ixn.ne.lastix.or.iyn.ne.lastiy.or.adjizn.ne.lastiz
            if(diffseg) then
               isz(n) = 0
               rowst = n
            else
               isz(n) = rowst - n
            endif
         endif
         lastix = ixn
         lastiy = iyn
         lastiz = izn
      enddo

!     Run the z grid backward to find ifz.

      do nrev = 1, ng
         n = ng - nrev + 1
         ixn = inxgridz(n)
         iyn = inygridz(n)
         izn = inzgridz(n)
         if(n.eq.ng) then
            ifz(ng) = 0
            rowfin = ng
         else
            adjizn = izn + 1
            diffseg = ixn.ne.lastix.or.iyn.ne.lastiy.or.adjizn.ne.lastiz
            if(diffseg) then
               ifz(n) = 0
               rowfin = n
            else
               ifz(n) = rowfin - n
            endif
         endif
         lastix = ixn
         lastiy = iyn
         lastiz = izn
      enddo

      end subroutine looplims

!
!
      subroutine graphgrid
!
!
      use quantum_variables
      implicit none

      integer(kind=ip) :: n

      open(68,file='grid_properties.dat')
      write(68,*) ' Loop for applying x-kinetic energy'
      do n = 1,ng
         write(68,*) ' x grid pt# = ',n,' x-index  = ',inxgridx(n)
         write(68,*) ' x grid pt# = ',n,' y-index  = ',inygridx(n)
         write(68,*) ' x grid pt# = ',n,' z-index  = ',inzgridx(n)
         write(68,*) ' x grid pt# = ',n,' loop-st  = ',isx(n)
         write(68,*) ' x grid pt# = ',n,' loop-fn  = ',ifx(n)
      enddo
      write(68,*) ' Loop for applying y-kinetic energy'
      do n = 1,ng
         write(68,*) ' y grid pt# = ',n,' x-index  = ',inxgridy(n)
         write(68,*) ' y grid pt# = ',n,' y-index  = ',inygridy(n)
         write(68,*) ' y grid pt# = ',n,' z-index  = ',inzgridy(n)
         write(68,*) ' y grid pt# = ',n,' loop-st  = ',isy(n)
         write(68,*) ' y grid pt# = ',n,' loop-fn  = ',ify(n)
      enddo
      write(68,*) ' Loop for applying z-kinetic energy'
      do n = 1,ng
         write(68,*) ' z grid pt# = ',n,' x-index  = ',inxgridz(n)
         write(68,*) ' z grid pt# = ',n,' y-index  = ',inygridz(n)
         write(68,*) ' z grid pt# = ',n,' z-index  = ',inzgridz(n)
         write(68,*) ' z grid pt# = ',n,' loop-st  = ',isz(n)
         write(68,*) ' z grid pt# = ',n,' loop-fn  = ',ifz(n)
      enddo
      write(68,*) ' Ranking of the y-grid points in the x-grid'
      do n = 1,ng
         write(68,*) ' y-grid pt# = ',n,' is x-grid pt# = ',irankyx(n) 
      enddo
      write(68,*) ' Ranking of the z-grid points in the x-grid'
      do n = 1,ng
        write(68,*) ' z-grid pt# = ',n,' is x-grid pt# = ',irankzx(n)
      enddo
      close(68)

      end subroutine graphgrid


!
!
      subroutine hvec3d(c,d)
!
!
      use quantum_variables
      implicit none

!     --- Global Arrays ---
      real(kind=dp), dimension(ng) :: c, d
!     --- Local Scalars ---     
      integer(kind=ip) :: i,ipp,ippp, j,jp,jpp, m,mp,mpp, n,np
      real(kind=dp) :: vntz, vntx, vnty, vnvc

!     [T_x + V(x,y,z)] act on c.
      
      do n = 1, ng
         vnvc = v_e(n)*c(n)
         j = inxgridx(n)
         vntx = 0.0_dp
         do jp = isx(n), ifx(n)
            jpp = j + jp
            np  = n + jp
            vntx = vntx + kex(jpp,j)*c(np)
         enddo
         d(n) = vntx + vnvc
      enddo

!     T_y acts on c and add result from above.

      do m = 1, ng
         i = inygridy(m)
         n = irankyx(m)
         vnty = 0.0_dp
         do ipp = isy(m), ify(m)
            ippp = i + ipp
            mp  = m + ipp
            np  = irankyx(mp)
            vnty = vnty + key(ippp,i)*c(np)
         enddo
         d(n) = d(n) + vnty
      enddo

!     T_z acts on c and add result from above.

      do m = 1, ng
         i = inzgridz(m)
         n = irankzx(m)
         vntz = 0.0_dp
         do ipp = isz(m), ifz(m)
            ippp = i + ipp
            mp  = m + ipp
            np  = irankzx(mp)
            vntz = vntz + kez(ippp,i)*c(np)
         enddo
         d(n) = d(n) + vntz
      enddo

      end subroutine hvec3d





