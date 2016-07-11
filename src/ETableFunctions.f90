      module NeighborTable
      contains
!===================================================================
      subroutine Create_NeiETable(biasArray)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      real(dp), intent(in) :: biasArray(:)
      integer :: i,j     
      integer :: iType, jType, iLowIndx, jLowIndx
      real(dp) :: EMax, ETab
      real(dp) :: biasOld, biasNew

      NeiETable=0d0
!      return
 
      if(NTotal .eq. 1) return
      iLowIndx = 0      
      do iType = 1,nMolTypes
        do i = ilowIndx+1, ilowIndx+NPART(iType)
          EMax = -1d40
          jLowIndx = 0
          do jType = 1,nMolTypes
            do j = jlowIndx+1,jlowIndx+NPART(jType)
              if(NeighborList(j,i)) then
                ETab = ETable(j) - biasArray(jType)/beta
                if(ETab .gt. EMax) then
                 EMax = ETab
               endif
              endif 
            enddo
            jLowIndx = jLowIndx + NMAX(jType)
          enddo    
          NeiETable(i) = EMax       
!          write(*,*) i, NeiETable(i)
        enddo
        iLowIndx = iLowIndx + NMAX(iType)
      enddo



!      do i=1,maxMol
!        if(.not. isActive(i)) then
!          cycle
!        endif
!        EMax = -1d40
!        do j=1,maxMol
!          if(.not. isActive(j)) then
!            cycle       
!          endif
!          if(NeighborList(j,i)) then
!             if(ETable(j) .gt. EMax) then
!               EMax = ETable(j)
!             endif
!           endif
!        enddo
!        NeiETable(i) = EMax
!      enddo
       
      end subroutine
!=================================================================================
      subroutine Insert_NewNeiETable(nType, PairList, dE, biasArray,newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer, intent(in) :: nType
      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)
      real(dp), intent(in) :: biasArray(:)

      integer :: i,j, nIndx
      integer :: iType, jType, iLowIndx, jLowIndx
      real(dp) :: EMax, ETab
       
      
      
      newNeiTable = 0d0
!      return
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx

!      write(2,*)
!      write(2,*) "NPART:", NPART
      iLowIndx = 0      
      do iType = 1,nMolTypes
        do i = ilowIndx+1,ilowIndx+NPART(iType)
          EMax = -1d40
          jLowIndx = 0
          do jType = 1,nMolTypes
            do j = jlowIndx+1,jlowIndx+NPART(jType)
!              write(2,*) i,j
              if(NeighborList(j,i)) then
                ETab = ETable(j) + dE(j) - biasArray(jType)/beta
                if(ETab .gt. EMax) then
                  EMax = ETab
                endif
              endif 
            enddo
            jLowIndx = jLowIndx + NMAX(jType)
          enddo    
!          write(2,*) i, nIndx
	  if(PairList(i) .le. Eng_Critr(iType, nType)) then
            ETab = ETable(nIndx) + dE(nIndx) - biasArray(nType)/beta
            if(ETab .gt. EMax) then
              EMax = ETab
            endif		
	  endif
          newNeiTable(i) = EMax       
        enddo
        iLowIndx = iLowIndx + NMAX(iType)
      enddo

      jLowIndx = 0
      EMax = -1d40
      do jType = 1,nMolTypes
        do j = jlowIndx+1,jlowIndx+NPART(jType)
!	  write(2,*) nindx, j
          if(PairList(j) .le. Eng_Critr(jType, nType)) then
            ETab = ETable(j) + dE(j) - biasArray(jType)/beta
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif 
        enddo
        jLowIndx = jLowIndx + NMAX(jType)
      enddo    
      newNeiTable(nIndx) = EMax
     

      end subroutine      
!=================================================================================
      subroutine Insert_NewNeiETable_Distance(nType, PairList, dE, newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer, intent(in) :: nType
      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)

      integer :: i,j, nIndx
      real(dp) :: EMax
       
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
      end module
