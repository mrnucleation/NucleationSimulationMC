!****** Energy Biased Aggregation-Volume-Bias Monte Carlo (AVBMC) algorithm *******
!   This file contains the nessisary functions to impliment the energy biased swap
!   move for cluster simulations. 
!===================================================================================            
      module Exchange_Module
      contains
!===================================================================================
      subroutine Exchange(E_T, acc_x)
      use AVBMC_RejectionVar
      use SimParameters
      use Constants  
      use ForceField
      use E_Interface
      use Coords
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions      
      use EnergyCriteria
      use DistanceCriteria
      use InterEnergy_LJ_Electro
      use EnergyTables
      use AcceptRates
      use CBMC_Variables
      use NeighborTable
      implicit none
      
      real(dp), intent(inout) :: E_T      
      real(dp), intent(inout) :: acc_x
      logical :: rejMove     
      integer :: NDiff(1:nMolTypes)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
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

      NDiff = 0 
      NDiff(nType) = 1
      call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx      

!      Generate the configuration for the newly inserted molecule
      nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx      
      select case(regrowType(nType))
      case(0)
        call Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, rosenRatio, rejMove)
      case(1)
        call Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, rosenRatio, rejMove)   
      case(2)
        call StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, rosenRatio, rejMove)   
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType)
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
        call QuickNei_ECalc_Inter(nTargType, nTargMol, rejMove)     
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
      call SwapIn_EnergyCalc(E_Inter, E_Intra, PairList, dETable, rejMove) 
      if(rejMove) then
        totalRej = totalRej + 1d0
        ovrlapRej = ovrlapRej + 1d0
        return
      endif
  

!     Determine the reverse probability of this move.
      if(distCriteria) then
        call Insert_NewNeiETable_Distance(nType, PairList, dETable, newNeiETable)  
      else
        call Insert_NewNeiETable(nType, PairList, dETable, newNeiETable)      
      endif
      call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
      call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)
      
!     Calculate the umbrella sampling bias.
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
         nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
         isActive(molArray(nType)%mol(NPART(nType)+1)%indx) = .true.

         if(distCriteria) then
           call NeighborUpdate_Distance(PairList,nIndx)        
         else
           call NeighborUpdate(PairList, nIndx)
         endif  
         NTotal = NTotal + 1
         ETable = ETable + dETable         
         NPART(nType) = NPART(nType) + 1 
!         call Create_NeiETable
         call Update_SubEnergies
       else
         totalRej = totalRej + 1d0
         dbalRej = dbalRej + 1d0
       endif
       end subroutine
!=================================================================================    
      end module
