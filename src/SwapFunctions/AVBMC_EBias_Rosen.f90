!****** Energy Biased Aggregation-Volume-Bias Monte Carlo (AVBMC) algorithm *******
!   This file contains the nessisary functions to impliment the energy biased swap
!   move for cluster simulations. 
!===================================================================================            
      module AVBMC_Module
      contains
!===================================================================================            
      subroutine AVBMC(E_T, acc_x, atmp_x)
      use SimParameters
      use Constants       
      implicit none
      real(dp), intent(inout) :: E_T, atmp_x, acc_x
      real(dp) :: grnd
        
      atmp_x = atmp_x + 1d0
      if(grnd() .lt. 0.5d0) then
        call AVBMC_EBias_Rosen_In(E_T, acc_x)     
      else
        call AVBMC_EBias_Rosen_Out(E_T, acc_x)    
      endif

      end subroutine
!===================================================================================
      subroutine AVBMC_EBias_Rosen_In(E_T, acc_x)
      use AcceptRates
      use AVBMC_RejectionVar
      use AVBMC_CBMC
      use CBMC_Utility
      use CBMC_Variables
      use Constants  
      use Coords
      use SimParameters
      use ForceField
!      use E_Interface
      use EnergyPointers, only: SwapIn_ECalc, Update_SubEnergies, Quick_Nei_ECalc
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions      
      use EnergyCriteria
      use DistanceCriteria
      use InterEnergy_LJ_Electro
      use EnergyTables
      use NeighborTable

      implicit none

!      interface      
!        subroutine Insert_NewNeiETable_Distance(nType,PairList,dE,newNeiTable)
!          implicit none
!          integer, intent(in) :: nType
!          real(dp), intent(in) :: PairList(:)
!          real(dp), intent(inout) :: dE(:), newNeiTable(:)
!        end subroutine
!      end interface      
      
!      interface      
!        subroutine Insert_NewNeiETable(nType,PairList,dE,newNeiTable)
!          implicit none
!          integer, intent(in) :: nType
!          real(dp), intent(in) :: PairList(:)
!          real(dp), intent(inout) :: dE(:), newNeiTable(:)
!        end subroutine
!      end interface
      
      real(dp), intent(inout) :: E_T      
      real(dp), intent(inout) :: acc_x
      logical :: rejMove     
      logical :: isIncluded(1:maxMol)
      integer :: NDiff(1:nMolTypes)
      integer :: i, iType, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nType, nIndx, bIndx
      integer :: atmType1, atmType2      
      real(dp) :: grnd
      real(dp) :: dx, dy, dz, r
      real(dp) :: genProbRatio, rosenRatio
      real(dp) :: E_Inter, E_Intra, bias_Diff
      real(dp) :: biasOld, biasNew
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      real(dp) :: newNeiETable(1:maxMol)      
      real(dp) :: x1, y1, z1
      real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out
      real(dp) :: rmin_ij
      real(dp) :: biasArray(1:nMolTypes)
  
      if(NTotal .eq. maxMol) return

!      Choose the type of molecule to be inserted      
      if(nMolTypes .eq. 1) then
        nType = 1
      else
        nType = floor(nMolTypes*grnd() + 1d0)
      endif
      atmpSwapIn(nType) = atmpSwapIn(nType) + 1d0
      atmpInSize(NTotal) = atmpInSize(NTotal) + 1d0
      if(NPART(nType) .eq. NMAX(nTYPE)) then
         return
      endif


      call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx      

!      Generate the configuration for the newly inserted molecule
      nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx      
      call Rosen_CreateSubset(nTarget, isIncluded)

      select case(regrowType(nType))
      case(0)
        call Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case(1)
        call Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)   
      case(2)
        call StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)   
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType)
!"
        stop
      end select        
      if(rejMove) then
        totalRej = totalRej + 1d0
        ovrlapRej = ovrlapRej + 1d0
        return
      endif      

!      call DEBUG_Output_NewConfig

!      Perform a check to see if the cluster criteria is statisfied or not.
      if(.not. distCriteria) then
        rejMove = .false.
!        call QuickNei_ECalc_Inter(nTargType, nTargMol, rejMove)     
        call Quick_Nei_ECalc(nTargType, nTargMol, rejMove)     
        if(rejMove) then
          totalRej = totalRej + 1d0
          critriaRej = critriaRej + 1d0
          clusterCritRej = clusterCritRej + 1d0
          return
        endif  
      endif

!      Calculate the Energy Difference Associated with the move
      E_Inter = 0d0
      E_Intra = 0d0
!      call SwapIn_EnergyCalc(E_Inter, E_Intra, PairList, dETable, rejMove) 
      call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove) 
      if(rejMove) then
        totalRej = totalRej + 1d0
        ovrlapRej = ovrlapRej + 1d0
        return
      endif
  
      biasArray = 0d0
      NDiff = 0
      NDiff(nType) = +1
      bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
      biasOld = NBias(bIndx)      
      do iType = 1, nMolTypes
        NDiff = 0
        NDiff(nType) = +1
        NDiff(iType) = NDiff(iType) - 1       
        if(NPART(iType) - 1 .lt. NMIN(iType)) cycle
        bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
        biasNew = NBias(bIndx)      
        bias_diff = biasNew - biasOld
        biasArray(iType) = bias_diff
      enddo   

!     Determine the reverse probability of this move.
      if(distCriteria) then
        call Insert_NewNeiETable_Distance(nType, PairList, dETable,  newNeiETable)  
      else
        call Insert_NewNeiETable(nType, PairList, dETable, biasArray, newNeiETable)      
      endif
      call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
      call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, biasArray, ProbSel_Out)
      
!     Calculate the umbrella sampling bias.
      NDiff = 0 
      NDiff(nType) = 1
      bIndx = getBiasIndex(NPart,NMAX)
      biasOld = NBias(bIndx)
      bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
      biasNew = NBias(bIndx)      
      bias_diff = biasNew - biasOld

!     Calculate acceptance probability and determine if the move is accepted or not          
      genProbRatio = (ProbTarg_Out * ProbSel_Out * avbmc_vol * dble(nMolTypes) * gas_dens(nType)) / (ProbTarg_In * rosenRatio)

      if( genProbRatio * exp(-beta*E_Inter + bias_diff) .gt. grnd() ) then
         acptSwapIn(nType) = acptSwapIn(nType) + 1d0        
         acptInSize(NTotal) = acptInSize(NTotal) + 1d0         
         do i=1,nAtoms(nType)      
           molArray(nType)%mol(NPART(nType)+1)%x(i) = newMol%x(i)
           molArray(nType)%mol(NPART(nType)+1)%y(i) = newMol%y(i)
           molArray(nType)%mol(NPART(nType)+1)%z(i) = newMol%z(i)
         enddo
         E_T = E_T + E_Inter + E_Intra
         acc_x = acc_x + 1d0
         isActive(nIndx) = .true.
!         nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
         if(distCriteria) then
           call NeighborUpdate_Distance(PairList,nIndx)        
         else
           call NeighborUpdate(PairList, nIndx)
         endif  
         NTotal = NTotal + 1
         ETable = ETable + dETable         
         NPART(nType) = NPART(nType) + 1 
         call Update_SubEnergies
       else
         totalRej = totalRej + 1d0
         dbalRej = dbalRej + 1d0
       endif
       end subroutine
!===================================================================================            
      subroutine AVBMC_EBias_Rosen_Out(E_T, acc_x)
      use AVBMC_CBMC
      use SimParameters
      use Constants
!      use E_Interface
      use EnergyPointers, only: SwapOut_ECalc, Update_SubEnergies
      use Coords
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use AcceptRates
      use UmbrellaFunctions
      use CBMC_Variables
      use NeighborTable
      implicit none
      
      real(dp), intent(inout) :: E_T      
      real(dp), intent(inout) :: acc_x
      
      logical :: rejMove  
      integer :: nTarget, nIndx, bIndx, iType
      integer :: nSel,nType, nMol,nTargMol,nTargType
      integer :: NDiff(1:nMolTypes)      
      real(dp) :: grnd
      real(dp) :: genProbRatio
      real(dp) :: bias_diff       
      real(dp) :: biasOld, biasNew      
      real(dp) :: E_Inter, E_Intra
      real(dp) :: dETable(1:maxMol)
      real(dp) :: ProbTargOut, ProbSel, ProbTargIn
      real(dp) :: rx, ry, rz, dist, rosenRatio
      real(dp) :: biasArray(1:nMolTypes)
      
      if(NTotal .eq. 1) return
      
!     Pick a Random Target Particle to Delete   
      biasArray = 0d0 
      bIndx = getBiasIndex(NPart,NMAX)
      biasOld = NBias(bIndx)
      do iType = 1, nMolTypes
        if(NPART(iType) .eq. 0) cycle
        NDiff = 0
        NDiff(iType) = -1
        if(NPART(iType) - 1 .lt. NMIN(iType)) cycle
        bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
        biasNew = NBias(bIndx)      
        bias_diff = biasNew - biasOld
        biasArray(iType) = bias_diff
      enddo      

      call Create_NeiETable(biasArray)
      call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut)
      call EBias_Remove_ChooseNeighbor(nTarget, biasArray, nSel, ProbSel)
      nType = typeList(nSel)
      nMol = subIndxList(nSel)
     
      atmpSwapOut(nType) = atmpSwapOut(nType) + 1d0      
      if(NPART(nType) .eq. NMIN(nType)) then
        return
      endif
      
      bias_diff = biasArray(nType)
!     Check to see that the appropriate atoms are within the insertion distance
!     in order to ensure the move is reversible. If not reject the move since
!     the reverse probility is equal to 0. 
      if(.not. distCriteria) then
        rx = molArray(nTargType)%mol(nTargMol)%x(1) - molArray(nType)%mol(nMol)%x(1)
        ry = molArray(nTargType)%mol(nTargMol)%y(1) - molArray(nType)%mol(nMol)%y(1)
        rz = molArray(nTargType)%mol(nTargMol)%z(1) - molArray(nType)%mol(nMol)%z(1)
        dist = rx*rx + ry*ry + rz*rz
        if(dist .gt. Dist_Critr_sq) then
          return
        endif
      endif
        
!     Check to see if the deletion of the particle will break the cluster
      rejMove=.false.
      call SwapOut_EnergyCriteria(nSel, rejMove)
      if(rejMove) then
        clusterCritRej = clusterCritRej + 1d0
        return
      endif


      select case(regrowType(nType))
      case(0)
        call Ridgid_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case(1)
        call Simple_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case(2)
        call StraightChain_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", nType
        stop
      end select 
!      Calculate the Energy Difference Associated with the move.
      E_Inter = 0d0
      E_Intra = 0d0
!      call SwapOut_EnergyCalc(E_Inter, E_Intra, nType, nMol, dETable)
      call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable)
      call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)
      genProbRatio = (ProbTargIn * rosenRatio) / (ProbTargOut * ProbSel * dble(nMolTypes) * avbmc_vol * gas_dens(nType))


!
!      Calculate Acceptance and determine if the move is accepted or not         
      if( genProbRatio * exp(-beta*E_Inter + bias_diff) .gt. grnd() ) then
         acptSwapOut(nType) = acptSwapOut(nType) + 1d0      
         molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
         E_T = E_T + E_Inter + E_Intra
         nIndx = molArray(nType)%mol(nMol)%indx
         call NeighborUpdate_Delete(nIndx)
         isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.         
         ETable = ETable - dETable
         ETable(nIndx) = ETable( molArray(nType)%mol(NPART(nType))%indx )
         ETable( molArray(nType)%mol(NPART(nType))%indx ) = 0d0
         NPART(nType) = NPART(nType) - 1
         NTotal = NTotal - 1         
         acc_x = acc_x + 1d0 
         call Update_SubEnergies
!         call DEBUG_Output_NeighborList
       endif
       end subroutine
!=================================================================================    
      subroutine EBias_Insert_ChooseTarget(nInsType, nTarget, nTargType, nMol, ProbSel)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer , intent(in) :: nInsType
      integer, intent(out) :: nTarget, nTargType, nMol
      real(dp), intent(out) :: ProbSel
      
      integer :: i, iType
      integer :: cnt(1:nMolTypes)
      real(dp) :: avgE(1:nMolTypes)
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
      real(dp) :: ranNum, sumInt   
        
      ProbTable = 0d0
      avgE = 0d0
      cnt = 0 
      do i = 1, maxMol
        if(isActive(i)) then
           iType = typeList(i)
           avgE(iType) = avgE(iType) + ETable(i)
           cnt(iType) = cnt(iType) + 1
        endif
      enddo

      do iType = 1, nMolTypes
        if(cnt(iType) .ne. 0) then
          avgE(iType) = avgE(iType)/dble(cnt(iType))
        endif
      enddo
      
      do i = 1, maxMol
        if(isActive(i)) then
          iType = typeList(i)        
          ProbTable(i) = exp(biasAlpha(nInsType,iType)*(ETable(i)-avgE(iType)))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm*grnd()

      sumInt = 0d0
      nTarget = 0
      do while(sumInt .lt. ranNum)
        nTarget = nTarget + 1
        sumInt = sumInt + ProbTable(nTarget)
      enddo
      
      nTargType = typeList(nTarget) 
      ProbSel = ProbTable(nTarget)/norm
      
      nMol = 0
      do i = 1, nTargType-1
        nMol = nMol + NMAX(i)
      enddo
      nMol = nTarget - nMol
      
      end subroutine
!=================================================================================    
      subroutine EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbRev)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget,nType
      real(dp), intent(in) :: newNeiETable(:)
      real(dp), intent(out) :: ProbRev
      
      integer :: i, nIndx
      integer :: cnt(1:nMolTypes)
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: norm, EMax
  
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      ProbTable = 0d0
!      EMax = newNeiETable(nTarget)
      EMax = maxval(newNeiETable)
      cnt = 0 
      do i = 1, maxMol
        if(isActive(i) .eqv. .true.) then
          ProbTable(i) = exp(beta*(newNeiETable(i)-Emax))
        elseif(i .eq. nIndx) then
          ProbTable(i) = exp(beta*(newNeiETable(i)-Emax))
        endif
      enddo
!      ProbTable(nIndx) = exp(beta*(newNeiETable(nIndx)-Emax))
      
      
      norm = sum(ProbTable)
      ProbRev = ProbTable(nTarget)/norm

      
      end subroutine
!=================================================================================    
      subroutine EBias_Insert_ReverseProbSel(nTarget, nType, dE, biasArray, ProbRev)
      use SimParameters  
      use Coords
      use EnergyTables
      use UmbrellaFunctions
      implicit none
      integer, intent(in) :: nTarget, nType
      real(dp), intent(in) :: dE(:)
      real(dp), intent(out) :: ProbRev
      real(dp) :: biasArray(:)
      
      integer :: i, iType, nIndx, bIndx
      integer :: NDiff(1:nMolTypes)
      real(dp) :: norm, bias_diff, biasOld, biasNew       



      
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      norm = 0d0
      do i = 1, maxMol
        if(NeighborList(nTarget,i)) then
          iType = typeList(i)        
          norm = norm + exp(beta*(ETable(i)+dE(i)-dE(nIndx)) + biasArray(iType))
!          norm = norm + exp(beta*(ETable(i)+dE(i)-dE(nIndx)))
        endif
      enddo
      norm = norm + exp(biasArray(nType))
      ProbRev = exp(biasArray(nType))/norm

!      norm = norm + 1d0
!      ProbRev = 1d0/norm      
      
      end subroutine
!=================================================================================    
      subroutine EBias_Remove_ChooseTarget(nTarget, nType, nMol, ProbTarget)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(out) :: nTarget, nType, nMol
      real(dp), intent(out) :: ProbTarget
      
      integer :: i
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
      real(dp) :: ranNum, sumInt   
        
      ProbTable = 0d0
      do i = 1, maxMol
        if(isActive(i)) then
          ProbTable(i) = exp(beta * NeiETable(i))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm * grnd()

      sumInt = ProbTable(1)
      nTarget = 1
      do while(sumInt .lt. ranNum)
        nTarget = nTarget + 1
        sumInt = sumInt + ProbTable(nTarget)
      enddo
      
      nType = typeList(nTarget) 
      ProbTarget = ProbTable(nTarget)/norm
      
      nMol = 0
      do i = 1, nType-1
        nMol = nMol + NMAX(i)
      enddo
      nMol = nTarget - nMol
      
      end subroutine
!=================================================================================    
      subroutine EBias_Remove_ChooseNeighbor(nTarget, biasArray, nSel, ProbTarget)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget
      real(dp), intent(in) :: biasArray(:)       
      integer, intent(out) ::  nSel
      real(dp), intent(out) :: ProbTarget
      
      integer :: i, iType
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
      real(dp) :: ranNum, sumInt   

      
      ProbTable = 0d0
      do i = 1, maxMol
        if(NeighborList(nTarget,i)) then
          iType = typeList(i)
          ProbTable(i) = exp(beta * ETable(i) + biasArray(iType))
!          ProbTable(i) = exp(beta * ETable(i))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm * grnd()
      sumInt = ProbTable(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbTable(nSel)
      enddo
      
      ProbTarget = ProbTable(nSel) / norm
      
      end subroutine 
!=================================================================================    
      subroutine EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dE, ProbSel)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nSel, nType
      integer, intent(in) :: nTarget
      real(dp), intent(in) :: dE(:)
      real(dp), intent(out) :: ProbSel

      
      integer :: i, iType
      integer :: cnt(1:nMolTypes)

      real(dp) :: avgE(1:nMolTypes)
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm       
        
      if(NTotal .eq. 2) then
         ProbSel=1d0
         return      
      endif
      
      ProbTable = 0d0
      avgE = 0d0
      cnt = 0 
      do i = 1, maxMol
        if(isActive(i)) then
          if(i .ne. nSel) then        
            iType = typeList(i)
            avgE(iType) = avgE(iType) + ETable(i) - dE(i)
            cnt(iType) = cnt(iType) + 1

          endif
        endif
      enddo

      do iType = 1, nMolTypes
        if(cnt(iType) .ne. 0) then
          avgE(iType) = avgE(iType)/dble(cnt(iType))
        endif
      enddo
      
      do i = 1, maxMol
        if(isActive(i)) then
          if(i .ne. nSel) then    
            iType = typeList(i)        
            ProbTable(i) = exp(biasAlpha(nType,iType)*((ETable(i)-dE(i))-avgE(iType)))
          endif
        endif
      enddo
    
      norm = sum(ProbTable)
      ProbSel = ProbTable(nTarget)/norm

      end subroutine

!=================================================================================    
      end module
