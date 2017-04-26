      module NeighborTable
      contains
!===================================================================
      subroutine Create_NeiETable(nType)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer :: i,j, nType
      integer :: iType, jType, iLowIndx, jLowIndx
      real(dp) :: EMax, ETab
      real(dp) :: biasOld, biasNew

      NeiETable = 0E0
      neiCount = 0
      if(NTotal .eq. 1) then
        return
      endif
      iLowIndx = 0      

      jLowIndx = 0
      do jType = 1, nType-1
        jLowIndx = jLowIndx + NMAX(jType)
      enddo

      do iType = 1,nMolTypes
        do i = ilowIndx+1, ilowIndx+NPART(iType)
          EMax = -huge(dp)
          do j = jlowIndx+1, jlowIndx+NPART(nType)
            if( NeighborList(j,i) ) then
              ETab = ETable(j)
              neiCount(i) = neiCount(i) + 1
              if(ETab .gt. EMax) then
                EMax = ETab
              endif
            endif 
          enddo
          NeiETable(i) = EMax       
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
      subroutine Insert_NewNeiETable(nType, PairList, dE, newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      implicit none
      integer, intent(in) :: nType
      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)
!      real(dp), intent(in) :: biasArray(:)

      integer :: i,j, nIndx
      integer :: iType, jType, iLowIndx, jLowIndx
      real(dp) :: EMax, ETab
       
      
      neiCount = 0
      newNeiTable = 0E0
!      return
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx

      jLowIndx = 0
      do jType = 1, nType-1
        jLowIndx = jLowIndx + NMAX(jType)
      enddo

      iLowIndx = 0
      do iType = 1,nMolTypes
        do i = ilowIndx+1,ilowIndx+NPART(iType)
          EMax = -huge(dp)
          do j = jlowIndx+1, jlowIndx+NPART(nType)
            if(NeighborList(j,i)) then
              ETab = ETable(j) + dE(j)
              neiCount(i) = neiCount(i) + 1
              if(ETab .gt. EMax) then
                EMax = ETab
              endif
            endif 
          enddo
          if(PairList(i) .le. Eng_Critr(nType, iType)) then
            ETab = ETable(nIndx) + dE(nIndx)
            neiCount(i) = neiCount(i) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif		
          endif
          newNeiTable(i) = EMax  
        enddo    
        iLowIndx = iLowIndx + NMAX(iType)
      enddo

      EMax = -huge(dp)
      do jType = 1,nMolTypes
        do j = jlowIndx+1,jlowIndx+NPART(nType)
          if(PairList(j) .le. Eng_Critr(jType, nType)) then
            ETab = ETable(j) + dE(j)
            neiCount(nIndx) = neiCount(nIndx) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif 
        enddo
      enddo    
      newNeiTable(nIndx) = EMax
     
!      if( all(neiCount .eq. 0) )then
!        write(*,*) "ERROR!"
!        do i = 1, maxMol
!          write(*,*) i, PairList(i), neiCount(i)
!        enddo
!      endif
      end subroutine      
!=================================================================================
      subroutine Insert_NewNeiETable_Distance_Old(nType, dE, newNeiTable)
      use Coords
      use EnergyTables
      use PairStorage
      use SimParameters
      implicit none
      integer, intent(in) :: nType
!      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)

      integer :: i,j, nIndx
      integer :: iType, iMol
      integer :: jType, jMol
      integer :: globIndxN, globIndx1, globIndx2
      real(dp) :: EMax
       
      neiCount = 0
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      globIndxN = molArray(nType)%mol(NPART(nType)+1)%globalIndx(1)
      do i = 1, maxMol
        newNeiTable(i) = 0E0

        EMax = -huge(dp)
        if( .not. isActive(i) ) then
          if( i .ne. nIndx ) then
             cycle
          else
            do j = 1, maxMol
              if(.not. isActive(j)) then
                cycle
              endif
              jType = typeList(j)
              jMol = subIndxList(j)
              globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
              if(rPairNew(globIndxN,globIndx2)%p%r_sq .le. Dist_Critr_sq) then
!                neiCount(i) = neiCount(i) + 1
                if(ETable(j) + dE(j) .gt. EMax) then
                  EMax = ETable(j) + dE(j)
                endif
              endif
            enddo            
          endif
        else
          iType = typeList(i)
          iMol = subIndxList(i)
          globIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
          do j = 1, maxMol
            if( .not. isActive(j) ) then
              if(j .ne. nIndx) then
                cycle
              else
                if(rPairNew(globIndxN,globIndx1)%p%r_sq .le. Dist_Critr_sq) then
                  if(dE(j) .gt. EMax) then
!                     neiCount(nIndx) = neiCount(nIndx) + 1
                     EMax = dE(j)
                   endif
                endif
              endif
            else
              if( NeighborList(i,j) ) then
                if( i .ne. j ) then         
                   if( ETable(j) + dE(j) .gt. EMax ) then
                    EMax = ETable(j) + dE(j)
                  endif
                endif
               endif
            endif        
          enddo
        endif
        newNeiTable(i) = EMax
      enddo
      
      write(2,*) "NeighborTable"
      do i = 1, maxMol
        write(2,*) i, newNeiTable(i), neiCount(i)
      enddo
      flush(2)
       
      end subroutine    

!=================================================================================
      subroutine Insert_NewNeiETable_Distance(nType, dE, newNeiTable)
      use EnergyTables
      use SimParameters
      use Coords
      use PairStorage, only: rPairNew
      implicit none
      integer, intent(in) :: nType
!      real(dp), intent(in) :: PairList(:)
      real(dp), intent(inout) :: dE(:), newNeiTable(:)
!      real(dp), intent(in) :: biasArray(:)

      integer :: i,j, nIndx
      integer :: iType, jType, iMol, jMol, iLowIndx, jLowIndx
      integer :: globIndxN, globIndx1, globIndx2
      real(dp) :: EMax, ETab
       
      
      neiCount = 0
      newNeiTable = 0E0
!      return
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      globIndxN = molArray(nType)%mol(NPART(nType)+1)%globalIndx(1)

      jLowIndx = 0
      do jType = 1, nType-1
        jLowIndx = jLowIndx + NMAX(jType)
      enddo

      iLowIndx = 0
      do iType = 1,nMolTypes
        do i = ilowIndx+1,ilowIndx+NPART(iType)
          iMol = subIndxList(i)
          globIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
          EMax = -huge(dp)
          do j = jlowIndx+1, jlowIndx+NPART(nType)
            if(NeighborList(j,i)) then
              ETab = ETable(j) + dE(j)
              neiCount(i) = neiCount(i) + 1
              if(ETab .gt. EMax) then
                EMax = ETab
              endif
            endif 
          enddo
          if(rPairNew(globIndxN,globIndx1)%p%r_sq .le. Dist_Critr_sq) then
            ETab = ETable(nIndx) + dE(nIndx)
            neiCount(i) = neiCount(i) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif		
          endif
          newNeiTable(i) = EMax  
        enddo    
        iLowIndx = iLowIndx + NMAX(iType)
      enddo


!      Calcuate the neiTable value for the newly inserted particle.

      EMax = -huge(dp)
      do jType = 1,nMolTypes
        do j = jlowIndx+1,jlowIndx+NPART(nType)
          jMol = subIndxList(j)
          globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
          if(rPairNew(globIndxN, globIndx2)%p%r_sq .le. Dist_Critr_sq) then
            ETab = ETable(j) + dE(j)
            neiCount(nIndx) = neiCount(nIndx) + 1
            if(ETab .gt. EMax) then
              EMax = ETab
            endif
          endif 
        enddo
      enddo    
      newNeiTable(nIndx) = EMax
     
!      if( all(neiCount .eq. 0) )then
!        write(*,*) "ERROR!"
!        do i = 1, maxMol
!          write(*,*) i, PairList(i), neiCount(i)
!        enddo
!      endif

!      write(2,*) "NeighborTable"
!      do i = 1, maxMol
!        write(2,*) i, newNeiTable(i), neiCount(i)
!      enddo
      end subroutine          
!=================================================================================
      end module
