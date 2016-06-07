!===================================================================
      subroutine Create_NeiETable
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer :: i,j     
      real(kind(0.0d0)) :: EMax

      NeiETable=0d0
      if(NTotal .eq. 1) return
      
      do i=1,maxMol
        if(.not. isActive(i)) then
          cycle
        endif
        EMax = -1d40
        do j=1,maxMol
          if(.not. isActive(j)) then
            cycle       
          endif
          if(NeighborList(j,i)) then
            if(i .ne. j) then          
              if(ETable(j) .gt. EMax) then
                EMax = ETable(j)
              endif
            endif
           endif
        enddo
        NeiETable(i) = EMax
      enddo
       
       
      end subroutine
!=================================================================================
      subroutine Insert_NewNeiETable(nType,PairList,dE,newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer, intent(in) :: nType
      real(kind(0.0d0)), intent(in) :: PairList(:)
      real(kind(0.0d0)), intent(inout) :: dE(:), newNeiTable(:)

      integer :: i,j, nIndx
      real(kind(0.0d0)) :: EMax
       
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      
!      write(35,*) "PairList"
!      do i = 1,maxMol
!        write(35,*) i, ETable(i) + dE(i), PairList(i)
!      enddo
      
      newNeiTable = 0d0
      do i=1,maxMol
       if(isActive(i) .eqv. .false.) then
         if(i .ne. nIndx) then
            cycle
         else
           EMax = -1d40               
           do j = 1, maxMol
             if(.not. isActive(j)) cycle           
             if(PairList(j) .le. Eng_Critr(typeList(i),typeList(j))) then
               if(ETable(j) + dE(j) .gt. EMax) then
                 EMax = ETable(j) + dE(j)
               endif
             endif
           enddo            
         endif
       else
         EMax = -1d40      
         do j=1,maxMol
           if(isActive(j) .eqv. .false.) then
             if(j .ne. nIndx) then
               cycle
             else
               if(PairList(i) .le. Eng_Critr(typeList(i),typeList(j))) then
                 if(dE(j) .gt. EMax) then
                   EMax = dE(j)
                 endif
               endif
             endif
           else
             if(NeighborList(j,i)) then
               if(i .ne. j) then         
                 if(ETable(j) + dE(j) .gt. EMax) then
                   EMax = ETable(j) + dE(j)
                 endif
               endif
              endif
           endif        
         enddo
       endif
       newNeiTable(i) = EMax
      enddo
      
!      write(35,*) "NeighborTable"
!      do i = 1, maxMol
!        write(35,*) i, newNeiTable(i)
!      enddo
       
       
      end subroutine      
!=================================================================================
      subroutine Insert_NewNeiETable_Distance(nType, PairList, dE, newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer, intent(in) :: nType
      real(kind(0.0d0)), intent(in) :: PairList(:)
      real(kind(0.0d0)), intent(inout) :: dE(:), newNeiTable(:)

      integer :: i,j, nIndx
      real(kind(0.0d0)) :: EMax
       
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      
!      write(35,*) "PairList"
!      do i = 1,maxMol
!        write(35,*) i, ETable(i) + dE(i), PairList(i)
!      enddo
      
      
      do i=1,maxMol
       newNeiTable(i) = 0d0
       EMax = -1d40
       if(isActive(i) .eqv. .false.) then
         if(i .ne. nIndx) then
            cycle
         else
           do j = 1, maxMol
             if(.not. isActive(j)) cycle
             if(PairList(j) .le. Dist_Critr_sq) then
               if(ETable(j) + dE(j) .gt. EMax) then
                 EMax = ETable(j) + dE(j)
               endif
             endif
           enddo            
         endif
       else
         do j=1,maxMol
           if(isActive(j) .eqv. .false.) then
             if(j .ne. nIndx) then
               cycle
             else
               if(PairList(i) .le. Dist_Critr_sq) then
                 if(dE(j) .gt. EMax) then
                   EMax = dE(j)
                 endif
               endif
             endif
           else
             if(NeighborList(i,j)) then
               if(i .ne. j) then         
                 if(ETable(j) + dE(j) .gt. EMax) then
                   EMax = ETable(j) + dE(j)
                 endif
               endif
              endif
           endif        
         enddo
       endif
       newNeiTable(i) = EMax
      enddo
      
!      write(35,*) "NeighborTable"
!      do i = 1, maxMol
!        write(35,*) i, newNeiTable(i)
!      enddo
       
       
      end subroutine         
!=================================================================================