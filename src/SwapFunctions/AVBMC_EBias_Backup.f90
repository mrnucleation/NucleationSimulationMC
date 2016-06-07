!****** Energy Biased Aggregation-Volume-Bias Monte Carlo (AVBMC) algorithm *******
!   This file contains the nessisary functions to impliment the energy biased swap
!   move for cluster simulations. 
!===================================================================================            
      subroutine AVBMC(E_T, acc_x, atmp_x)
      use SimParameters
      use Constants       
      implicit none
      real(kind(0.0d0)), intent(inout) :: E_T, atmp_x, acc_x
      real(kind(0.0d0)) :: grnd
        
      atmp_x = atmp_x + 1d0
      if(grnd() .lt. 0.5d0) then
        call AVBMC_EBias_In(E_T, acc_x)     
      else
        call AVBMC_EBias_Out(E_T, acc_x)    
      endif

      end subroutine
!===================================================================================
      subroutine AVBMC_EBias_In(E_T, acc_x)
      use SimParameters
      use Constants  
      use ForceField
      use E_Interface
      use Coords
      use UmbrellaFunctions
      use ForceField      
      use IndexingFunctions      
      use EnergyCriteria      
      use InterEnergy_LJ_Electro
      use EnergyTables
      use AcceptRates      
      implicit none

      interface
        subroutine EBias_Insert_ReverseProbTarget(nTarget,nType,newNeiETable, ProbRev)
          implicit none        
          integer, intent(in) :: nTarget,nType
          real(kind(0.0d0)), intent(in) :: newNeiETable(:)
          real(kind(0.0d0)), intent(out) :: ProbRev      
        end subroutine
      end interface

      interface
        subroutine EBias_Insert_ReverseProbSel(nTarget, nType, dE, ProbRev)
         implicit none
         integer, intent(in) :: nTarget, nType
         real(kind(0.0d0)), intent(in) :: dE(:)
         real(kind(0.0d0)), intent(out) :: ProbRev
        end subroutine
      end interface
      
      interface      
        subroutine Insert_NewNeiETable(nType,PairList,dE,newNeiTable)
          implicit none
          integer, intent(in) :: nType
          real(kind(0.0d0)), intent(in) :: PairList(:)
          real(kind(0.0d0)), intent(inout) :: dE(:), newNeiTable(:)
        end subroutine
      end interface
      
      real(kind(0.0d0)), intent(inout) :: E_T      
      real(kind(0.0d0)), intent(inout) :: acc_x
      logical rejMove     
      integer :: NDiff(1:nMolTypes)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nType, nIndx, bIndx
      real(kind(0.0d0)) :: grnd
      real(kind(0.0d0)) :: dx, dy, dz, r
      real(kind(0.0d0)) :: genProbRatio
      real(kind(0.0d0)) :: E_Diff, bias_Diff
      real(kind(0.0d0)) :: biasOld, biasNew
      real(kind(0.0d0)) :: PairList(1:maxMol)
      real(kind(0.0d0)) :: dETable(1:maxMol)
      real(kind(0.0d0)) :: newNeiETable(1:maxMol)      
      real(kind(0.0d0)) :: x1, y1, z1
      real(kind(0.0d0)) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out
      
      if(NTotal .eq. maxMol) return

!      Choose the type of molecule to be inserted      
      if(nMolTypes .eq. 1) then
        nType = 1
      else
        nType = floor(nMolTypes*grnd() + 1d0)
      endif
      atmpSwapIn(nType) = atmpSwapIn(nType) + 1d0
      
      if(NPART(nType) .eq. NMAX(nTYPE)) then
         return
      endif

      NDiff = 0 
      NDiff(nType) = 1
      call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx      

!     Uniformly generate a random distance
      r = Dist_Critr * grnd()**(1d0/3d0)
      if(r .lt. global_r_min) return

      call Generate_UnitSphere(dx, dy, dz)
      dx = r * dx
      dy = r * dy
      dz = r * dz      

      call Ridgid_ConfigGen(nType)

      x1 = molArray(nTargType)%mol(nTargMol)%x(1) + dx - newMol%x(1)
      y1 = molArray(nTargType)%mol(nTargMol)%y(1) + dy - newMol%y(1)
      z1 = molArray(nTargType)%mol(nTargMol)%z(1) + dz - newMol%z(1)

      do i=1,nAtoms(nType)
        newMol%x(i) = newMol%x(i) + x1
        newMol%y(i) = newMol%y(i) + y1
        newMol%z(i) = newMol%z(i) + z1
      enddo

      rejMove = .false.
      call QuickNei_ECalc_Inter(nTargType, nTargMol, rejMove)     
      if(rejMove) then
        clusterCritRej = clusterCritRej + 1d0
        return
      endif  

!      Calculate the Energy Difference Associated with the move
      E_Diff=0d0
      call SwapIn_EnergyCalc(E_Diff, PairList, dETable, rejMove) 
      if(rejMove) then
        return
      endif

!      call SwapIn_EnergyCriteria(nType,PairList,rejMove)
!      if(rejMove) then
!        return
!      endif  
      
!     Determine the reverse probability of this move.
      call Insert_NewNeiETable(nType, PairList, dETable, newNeiETable)
      call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
      call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)
      
!     Calculate the umbrella sampling bias.
      bIndx = getBiasIndex(NPart,NMAX)
      biasOld = NBias(bIndx)
      bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
      biasNew = NBias(bIndx)      
      bias_diff = biasNew - biasOld


!     Calculate acceptance probability and determine if the move is accepted or not          
      genProbRatio = (ProbTarg_Out * ProbSel_Out * avbmc_vol * dble(nMolTypes) * gas_dens(nType)) / ProbTarg_In
!     !write(35,*) genProbRatio
!     !write(35,*) ProbTarg_Out, ProbSel_Out, ProbTarg_In
!     !write(35,*) E_Diff, bias_diff
!     !write(35,*) "----------------------------"      
      if( genProbRatio * exp(-beta*E_Diff + bias_diff) .gt. grnd() ) then
         acptSwapIn(nType) = acptSwapIn(nType) + 1d0           
         do i=1,nAtoms(nType)      
           molArray(nType)%mol(NPART(nType)+1)%x(i) = newMol%x(i)
           molArray(nType)%mol(NPART(nType)+1)%y(i) = newMol%y(i)
           molArray(nType)%mol(NPART(nType)+1)%z(i) = newMol%z(i)
         enddo
         E_T = E_T + E_Diff
         acc_x = acc_x + 1d0
         isActive(molArray(nType)%mol(NPART(nType)+1)%indx) = .true.
         nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
         call NeighborUpdate(PairList, nIndx)
         NTotal = NTotal + 1
         ETable = ETable + dETable         
         NPART(nType) = NPART(nType) + 1 
         call Create_NeiETable         
       endif
       end subroutine
!===================================================================================            
      subroutine AVBMC_EBias_Out(E_T, acc_x)
      use SimParameters
      use Constants
      use E_Interface
      use Coords
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions
      use EnergyCriteria
      use EnergyTables
      use AcceptRates
      implicit none

      interface
        subroutine EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dE, ProbSel)
         implicit none
         integer, intent(in) :: nSel, nType
         integer, intent(in) :: nTarget
         real(kind(0.0d0)), intent(in) :: dE(:)
         real(kind(0.0d0)), intent(out) :: ProbSel
        end subroutine
      end interface      
      
      real(kind(0.0d0)), intent(inout) :: E_T      
      real(kind(0.0d0)), intent(inout) :: acc_x
      
      logical :: rejMove  
      integer :: nTarget, nIndx, bIndx
      integer :: nSel,nType, nMol,nTargMol,nTargType
      integer :: NDiff(1:nMolTypes)      
      real(kind(0.0d0)) :: grnd
      real(kind(0.0d0)) :: genProbRatio
      real(kind(0.0d0)) :: bias_diff       
      real(kind(0.0d0)) :: biasOld, biasNew      
      real(kind(0.0d0)) :: E_Diff
      real(kind(0.0d0)) :: dETable(1:maxMol)
      real(kind(0.0d0)) :: ProbTargOut, ProbSel, ProbTargIn
      real(kind(0.0d0)) :: rx, ry, rz, dist
      
      if(NTotal .eq. 1) return
      
!     Pick a Random Target Particle to Delete   
!      nTarget = floor(NTotal*grnd() + 1d0)
!      call Get_MolIndex(nTarget, NPART, nTargType, nTargMol)
!      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx 

      call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut)
      call EBias_Remove_ChooseNeighbor(nTarget, nSel, ProbSel)
      nType = typeList(nSel)
      nMol = subIndxList(nSel)
     
      atmpSwapOut(nType) = atmpSwapOut(nType) + 1d0      
      if(NPART(nType) .eq. NMIN(nType)) then
        return
      endif
      
      NDiff = 0
      NDiff(nType) = -1
      
!     Check to see that the appropriate atoms are within the insertion distance
!     in order to ensure the move is reversible. If not reject the move since
!     the reverse probility is equal to 0. 
      rx = molArray(nTargType)%mol(nTargMol)%x(1)-molArray(nType)%mol(nMol)%x(1)
      ry = molArray(nTargType)%mol(nTargMol)%y(1)-molArray(nType)%mol(nMol)%y(1)
      rz = molArray(nTargType)%mol(nTargMol)%z(1)-molArray(nType)%mol(nMol)%z(1)
      dist = rx**2 + ry**2 + rz**2
      if(dist .gt. Dist_Critr_sq) then
        return
      endif
        
!     Check to see if the deletion of the particle will break the cluster
      rejMove=.false.
      call SwapOut_EnergyCriteria(nSel, rejMove)
      if(rejMove) then
        clusterCritRej = clusterCritRej + 1d0
        return
      endif

!      Calculate the Energy Difference Associated with the move and also Calculate
!      the Q6 value used for the biased sampling.        
      E_Diff=0d0
      call SwapOut_EnergyCalc(nType, nMol, dETable, E_Diff)

      bIndx = getBiasIndex(NPart,NMAX)
      biasOld = NBias(bIndx)
      bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
      biasNew = NBias(bIndx)      
      bias_diff = biasNew - biasOld

      call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType,dETable,ProbTargIn)

  
      
      genProbRatio = ProbTargIn / (ProbTargOut * ProbSel * dble(nMolTypes) * avbmc_vol * gas_dens(nType))
      
!      Calculate Acceptance and determine if the move is accepted or not         
      if( genProbRatio * exp(-beta*E_Diff + bias_diff) .gt. grnd() ) then
         acptSwapOut(nType) = acptSwapOut(nType) + 1d0      
         molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
         E_T = E_T + E_Diff         
         nIndx = molArray(nType)%mol(nMol)%indx
         call NeighborUpdate_Delete(nIndx)
         isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.         
         acc_x = acc_x + 1d0 
         ETable = ETable - dETable
         ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
         ETable(molArray(nType)%mol(NPART(nType))%indx) = 0d0
         NPART(nType) = NPART(nType) - 1
         NTotal = NTotal - 1         
         call Create_NeiETable         
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
      real(kind(0.0d0)), intent(out) :: ProbSel
      
      integer :: i, iType
      integer :: cnt(1:nMolTypes)
      real(kind(0.0d0)) :: avgE(1:nMolTypes)
      real(kind(0.0d0)) :: ProbTable(1:maxMol)
      real(kind(0.0d0)) :: grnd, norm       
      real(kind(0.0d0)) :: ranNum, sumInt   
        
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
      subroutine EBias_Insert_ReverseProbTarget(nTarget,nType,newNeiETable, ProbRev)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget,nType
      real(kind(0.0d0)), intent(in) :: newNeiETable(:)
      real(kind(0.0d0)), intent(out) :: ProbRev
      
      integer :: i, nIndx
      integer :: cnt(1:nMolTypes)
      real(kind(0.0d0)) :: ProbTable(1:maxMol)
      real(kind(0.0d0)) :: norm, EMax
  
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      ProbTable = 0d0
      EMax = newNeiETable(nTarget)
      cnt = 0 
      do i = 1, maxMol
        if(isActive(i) .or. (i .eq. nIndx)) then
          ProbTable(i) = exp(beta*(newNeiETable(i)-Emax))
        endif
      enddo
!      ProbTable(nIndx) = exp(beta*(newNeiETable(nIndx)-Emax))
      
      
      norm = sum(ProbTable)
      ProbRev = ProbTable(nTarget)/norm

      
      end subroutine
!=================================================================================    
      subroutine EBias_Insert_ReverseProbSel(nTarget, nType, dE, ProbRev)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget, nType
      real(kind(0.0d0)), intent(in) :: dE(:)
      real(kind(0.0d0)), intent(out) :: ProbRev
      
      integer :: i, nIndx
      integer :: cnt(1:nMolTypes)
      real(kind(0.0d0)) :: norm       
  
      cnt = 0 
      nIndx = molArray(nType)%mol(NPART(nType)+1)%indx
      norm = 0d0
      do i = 1, maxMol
        if(NeighborList(nTarget,i)) then
          norm = norm + exp(beta*(ETable(i)+dE(i)-dE(nIndx)))
        endif
      enddo
      norm = norm + 1d0
      ProbRev = 1d0/norm
      
      end subroutine
!=================================================================================    
      subroutine EBias_Remove_ChooseTarget(nTarget, nType, nMol, ProbTarget)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(out) :: nTarget, nType, nMol
      real(kind(0.0d0)), intent(out) :: ProbTarget
      
      integer :: i
      real(kind(0.0d0)) :: ProbTable(1:maxMol)
      real(kind(0.0d0)) :: grnd, norm       
      real(kind(0.0d0)) :: ranNum, sumInt   
        
      ProbTable = 0d0
      do i = 1, maxMol
        if(isActive(i)) then
          ProbTable(i) = exp(beta * NeiETable(i))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm * grnd()

      sumInt = 0d0
      nTarget = 0
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
      subroutine EBias_Remove_ChooseNeighbor(nTarget, nSel, ProbTarget)
      use SimParameters  
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget      
      integer, intent(out) ::  nSel
      real(kind(0.0d0)), intent(out) :: ProbTarget
      
      integer :: i
      real(kind(0.0d0)) :: ProbTable(1:maxMol)
      real(kind(0.0d0)) :: grnd, norm       
      real(kind(0.0d0)) :: ranNum, sumInt   
        
      ProbTable = 0d0
      do i = 1, maxMol
        if(NeighborList(nTarget,i)) then
          ProbTable(i) = exp(beta * ETable(i))
        endif
      enddo

      norm = sum(ProbTable)
      ranNum = norm * grnd()

      sumInt = 0d0
      nSel = 0
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
      real(kind(0.0d0)), intent(in) :: dE(:)
      real(kind(0.0d0)), intent(out) :: ProbSel

      
      integer :: i, iType
      integer :: cnt(1:nMolTypes)

      real(kind(0.0d0)) :: avgE(1:nMolTypes)
      real(kind(0.0d0)) :: ProbTable(1:maxMol)
      real(kind(0.0d0)) :: grnd, norm       
        
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
