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
      logical rejMove     
      integer :: NDiff(1:nMolTypes)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nType1, nType2, nIndx, nIndx2,  bIndx
      integer :: atmType1, atmType2      
      real(dp) :: grnd
      real(dp) :: dx, dy, dz, r
      real(dp) :: genProbRatio, rosenRatio_in, rosenRatio_out
      real(dp) :: E_Inter, E_Intra, bias_Diff
      real(dp) :: biasOld, biasNew
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      real(dp) :: newNeiETable(1:maxMol)      
      real(dp) :: x1, y1, z1
      real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out
      real(dp) :: rmin_ij

      nTarget = floor(NTotal*grnd() + 1d0)
      call Get_MolIndex(nMove, NPart, nTargType, nTargMol)
      call Uniform_ChooseNeighbor(nTarget, nMove, nNei)

      nType2 = nType1
      do while(nType2 .eq. nType1) 
        nType2 = floor(nMolTypes*grnd() + 1d0)
      enddo

         
      NDiff = 0 
      NDiff(nType1) = -1
      NDiff(nType2) = 1

      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx      

!      Generate the configuration for the newly inserted molecule
      nIndx2 = MolArray(nType)%mol(NPART(nType) + 1)%indx      
      call Rosen_CreateSubset(nTarget, isIncluded)
      isIncluded(nIndx) = .false.
      select case(regrowType(nType2))
      case(0)
        call Ridgid_RosenConfigGen(nType2, nIndx2, nTarget, nTargType, isIncluded, rosenRatio_in, rejMove)
      case(1)
        call Simple_RosenConfigGen(nType2, nIndx2, nTarget, nTargType, isIncluded, rosenRatio_in, rejMove)   
      case(2)
        call StraightChain_RosenConfigGen(nType2, nIndx2, nTarget, nTargType, isIncluded, rosenRatio_in, rejMove)   
      case default
        write(*,*) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType2)
        stop
      end select        
      if(rejMove) then
        return
      endif          

!      Perform a check to see if the cluster criteria is statisfied or not.
      if(.not. distCriteria) then
        rejMove = .false.
        call QuickNei_ECalc_Inter(nTargType, nTargMol, rejMove)     
        if(rejMove) then
          return
        endif  
      endif
      
!      call DEBUG_Output_NewConfig


      select case(regrowType(nType1))
      case(0)
         call Ridgid_RosenConfigGen_Reverse(nType1, nMol, nTarget, nTargType, rosenRatio_out)
      case(1)
         call Simple_RosenConfigGen_Reverse(nType1, nMol, nTarget, nTargType, rosenRatio_out)
      case(2)
         call StraightChain_RosenConfigGen_Reverse(nType1, nMol, nTarget, nTargType, rosenRatio_out)
      case default
         write(*,*) "Error! EBias can not regrow a molecule of regrow type:", nType1
         stop
      end select 




!      Calculate the Energy Difference Associated with the move
      E_Inter = 0d0
      E_Intra = 0d0
      call SwapIn_EnergyCalc(E_Inter, E_Intra, PairList, dETable, rejMove) 
      if(rejMove) then
        return
      endif
      
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
         isActive(molArray(nType)%mol(NPART(nType)+1)%indx) = .true.
         nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
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
       endif

       end subroutine
!=================================================================================	  
      subroutine Uniform_ChooseNeighbor(nTarget, nSel, nNei)
      use SimParameters  
      use Coords
      implicit none
      integer, intent(in) :: nTarget
      integer, intent(out) :: nSel, nNei
      integer :: i,j
      integer :: ListCur(1:60)
      real(kind(0.0d0)) :: grnd	  
	  
        
      nNei = 0
      ListCur = 0
      do j = 1, maxMol
        if( NeighborList(nTarget,j) ) then
          if( j .ne. nTarget ) then             
            nNei = nNei + 1
            ListCur(nNei) = j
          endif
        endif
      enddo

      j = floor(nNei*grnd() + 1d0)	  
      nSel = ListCur(j)
      
      if(nSel .eq. 0) then
        write(35,*) "ERROR"
        write(35,*) "NPART:", NPART
        do i = 1, maxMol
          write(35,*) (NeighborList(i,j),j=1,maxMol)      
        enddo
      endif
      
      end subroutine
!=================================================================================    
      end module
