      subroutine indexx(compar)

      use quantum_variables
      implicit none

      integer(kind=ip), parameter :: Mtmp=7, NSTACK=50
      integer(kind=ip), dimension(NSTACK) :: istack
      integer(kind=ip) :: i, indxi, indxt, ir, itemp, j, jstack, k, l
      logical, external :: compar

      do j = 1, ng
        indxx(j) = j
      enddo

      jstack = 0
      l = 1
      ir = ng
 1    if(ir-l.lt.Mtmp)then
         do 13 j=l+1,ir
            indxt=indxx(j)
            do 12 i=j-1,l,-1
               indxi = indxx(i)
               if(compar(indxt,indxi)) goto 2
               indxx(i+1)=indxi
 12         continue
            i=l-1
 2          indxx(i+1)=indxt
 13      continue
         if(jstack.eq.0) goto 40
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indxx(k)
         indxx(k)=indxx(l+1)
         indxx(l+1)=itemp
         if(compar(indxx(l),indxx(ir))) then
            itemp=indxx(l)
            indxx(l)=indxx(ir)
            indxx(ir)=itemp
         endif
         if(compar(indxx(l+1),indxx(ir))) then
            itemp=indxx(l+1)
            indxx(l+1)=indxx(ir)
            indxx(ir)=itemp
         endif
         if(compar(indxx(l),indxx(l+1))) then
            itemp=indxx(l)
            indxx(l)=indxx(l+1)
            indxx(l+1)=itemp
         endif
         i=l+1
         j=ir
         indxt=indxx(l+1)
 3       continue
         i=i+1
         if(compar(indxt,indxx(i))) goto 3
 4       continue
         j=j-1
         if(compar(indxx(j),indxt)) goto 4
         if(j.lt.i) goto 5
         itemp=indxx(i)
         indxx(i)=indxx(j)
         indxx(j)=itemp
         goto 3
 5       indxx(l+1)=indxx(j)
         indxx(j)=indxt
         jstack=jstack+2
         if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1

 40   continue

      return
      end subroutine indexx



