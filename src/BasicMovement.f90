!===========================================================================================
      module SimpleMCMoves_Module
      contains
!===========================================================================================
      subroutine SingleAtom_Translation(E_T,acc_x,atmp_x)
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies 
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables
      use CBMC_Variables      
      use PairStorage, only: UpdateDistArray
      implicit none
      
      real(dp),intent(inout) :: E_T,acc_x,atmp_x      
      real(dp) max_distx

      logical, parameter :: useIntra(1:4) = [.true., .true., .true., .true.]
      
      logical rejMove      
      integer nType,nMol,nIndx,nMove, nAtom
      real(dp) :: grnd 
      real(dp) :: dx,dy,dz      
      real(dp) :: E_Diff,E_Inter, E_Intra
      type (displacement) :: disp(1:1)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      

      max_distx = 0.1E0
      rejMove = .false.
      atmp_x = atmp_x + 1E0
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)

      call Get_MolIndex(nMove, NPart, nType, nMol)
!      if(nType .ne. 1) then
!        return
!      endif
      nIndx = MolArray(nType)%mol(nMol)%indx
      if(regrowType(nType) .eq. 0) return
      
      nAtom = floor(nAtoms(nType)*grnd() + 1E0)
      

!     Generate a random translational displacement             
      dx = max_distx * (2E0*grnd() - 1E0)
      dy = max_distx * (2E0*grnd() - 1E0)
      dz = max_distx * (2E0*grnd() - 1E0)

!     Construct the Displacement Vectors for each atom in the molecule that was chosen.'
      disp(1)%molType = int(nType,atomIntType)
      disp(1)%molIndx = int(nMol,atomIntType)
      disp(1)%atmIndx = int(nAtom,atomIntType)
        
      disp(1)%x_old => MolArray(nType)%mol(nMol)%x(nAtom)
      disp(1)%y_old => MolArray(nType)%mol(nMol)%y(nAtom)
      disp(1)%z_old => MolArray(nType)%mol(nMol)%z(nAtom)

      disp(1)%x_new = disp(1)%x_old + dx
      disp(1)%y_new = disp(1)%y_old + dy
      disp(1)%z_new = disp(1)%z_old + dz

      
!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
      E_Inter = 0E0
      E_Intra = 0E0
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:1), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:1), PairList, dETable, useIntra, rejMove)
      if(rejMove) return
      E_Diff = E_Inter + E_Intra

!     Calculate Acceptance and determine if the move is accepted or not     
      if(E_Diff .le. 0E0) then
        disp(1)%x_old = disp(1)%x_new
        disp(1)%y_old = disp(1)%y_new
        disp(1)%z_old = disp(1)%z_new
        E_T = E_T + E_Diff
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        if(distCriteria) then
          if(disp(1)%atmIndx .eq. 1) then
            call NeighborUpdate_Distance(PairList,nIndx)        
          endif
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies
      elseif(exp(-beta*E_Diff) .gt. grnd()) then
        disp(1)%x_old = disp(1)%x_new
        disp(1)%y_old = disp(1)%y_new
        disp(1)%z_old = disp(1)%z_new
        E_T = E_T + E_Diff
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        call UpdateDistArray
        if(distCriteria) then
          if(disp(1)%atmIndx .eq. 1) then
            call NeighborUpdate_Distance(PairList,nIndx)        
          endif
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      endif

      	  
      end subroutine
!===========================================================================================
      subroutine Translation(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables
      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none
      
      real(dp),intent(inout) :: E_T, acc_x, atmp_x      
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: iAtom, nType, nMol, nIndx, nMove
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: grnd
      real(dp) :: dx, dy, dz      
      real(dp) :: E_Inter, E_Intra
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
      if(NTotal .eq. 1) return
      
      rejMove = .false.
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)

      atmp_x = atmp_x + 1E0
      atmpTrans(nType) = atmpTrans(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx
!     Generate a random translational displacement             
      dx = max_dist(nType) * (2E0*grnd() - 1E0)
      dy = max_dist(nType) * (2E0*grnd() - 1E0)
      dz = max_dist(nType) * (2E0*grnd() - 1E0)

!     Construct the Displacement Vectors for each atom in the molecule that was chosen.
      do iAtom=1,nAtoms(nType)
        disp(iAtom)%molType = int(nType, atomIntType)
        disp(iAtom)%molIndx = int(nMol, atomIntType)
        disp(iAtom)%atmIndx = int(iAtom, atomIntType)
        
        disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
        disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
        disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
        
        disp(iAtom)%x_new = disp(iAtom)%x_old + dx
        disp(iAtom)%y_new = disp(iAtom)%y_old + dy
        disp(iAtom)%z_new = disp(iAtom)%z_old + dz        
      enddo
      
!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
      E_Inter = 0E0
      E_Intra = 0E0
      !call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) then
        return
      endif
      
      biasDiff = 0E0
!      write(*,*) useUmbrella
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
        

!     Calculate Acceptance and determine if the move is accepted or not     
      if(biasEnergy .le. 0E0) then
        do iAtom=1,nAtoms(nType)      
          disp(iAtom)%x_old = disp(iAtom)%x_new
          disp(iAtom)%y_old = disp(iAtom)%y_new
          disp(iAtom)%z_old = disp(iAtom)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x + 1E0 
        acptTrans(nType) = acptTrans(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList,nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif    
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      elseif(exp(-biasEnergy) .gt. grnd()) then
        do iAtom=1,nAtoms(nType)      
          disp(iAtom)%x_old = disp(iAtom)%x_new
          disp(iAtom)%y_old = disp(iAtom)%y_new
          disp(iAtom)%z_old = disp(iAtom)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        acptTrans(nType) = acptTrans(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList, nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      endif

      	  
      end subroutine
!===========================================================================================
      subroutine Rotation(E_T,acc_x,atmp_x)
      use SimParameters    
      implicit none
      real(dp), intent(inout) :: atmp_x,acc_x,E_T
      real(dp) :: ran_num, grnd


      if(NTotal .eq. 1) then
!       acc_x=acc_x+1E0
       return            
      endif      
      

      ran_num = grnd()
      if(ran_num .lt. 1E0/3E0) then
         call Rot_xy(E_T, acc_x, atmp_x)      
      elseif(ran_num .lt. 2E0/3E0) then
         call Rot_xz(E_T, acc_x, atmp_x)  
      else
         call Rot_yz(E_T, acc_x, atmp_x)  
      endif      
      
      end subroutine   

!=======================================================      
      subroutine Rot_xy(E_T,acc_x,atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies 
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables
      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T,acc_x,atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: i,nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd   
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term,s_term
      real(dp) :: x_scale, y_scale
      real(dp) :: xcm,ycm,angle
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if(nAtoms(nType) .eq. 1) then
        return
      endif
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array 
      do i=1,nAtoms(nType)
        disp(i)%molType = int(nType,2)
        disp(i)%molIndx = int(nMol,2)
        disp(i)%atmIndx = int(i,2)
        
        disp(i)%x_old => MolArray(nType)%mol(nMol)%x(i)
        disp(i)%y_old => MolArray(nType)%mol(nMol)%y(i)
        disp(i)%z_old => MolArray(nType)%mol(nMol)%z(i)
      enddo
      
!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd()-1E0)
      c_term = cos(angle)
      s_term = sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion. 
      xcm=0E0
      ycm=0E0
      do i=1,nAtoms(nType)
        atmType = atomArray(nType,i)
        xcm = xcm + atomData(atmType)%mass*disp(i)%x_old
        ycm = ycm + atomData(atmType)%mass*disp(i)%y_old
      enddo

      xcm = xcm/totalMass(nType)   
      ycm = ycm/totalMass(nType)

!     Generate a random translational displacement      
      do i=1,nAtoms(nType)
        disp(i)%z_new = disp(i)%z_old
        x_scale = disp(i)%x_old - xcm
        y_scale = disp(i)%y_old - ycm
        disp(i)%x_new = c_term*x_scale - s_term*y_scale + xcm
        disp(i)%y_new = s_term*x_scale + c_term*y_scale + ycm
      enddo

!     Calculate the Energy Difference Associated with the move   
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) return

      biasDiff = 0E0
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
      if(biasEnergy .le. 0E0) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x+1E0
        acptRot(nType) = acptRot(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList, nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      elseif(exp(-biasEnergy) .gt. grnd()) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x+1E0
        acptRot(nType) = acptRot(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList, nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      endif
      end subroutine
!=======================================================      
      subroutine Rot_xz(E_T,acc_x,atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables      
      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T,acc_x,atmp_x

      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: i,nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: angle
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd   
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term,s_term
      real(dp) :: x_scale, z_scale
      real(dp) :: xcm,zcm
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if(nAtoms(nType) .eq. 1) then
        return
      endif
      atmp_x = atmp_x+1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array 
      do i=1,nAtoms(nType)
        disp(i)%molType = int(nType,atomIntType)
        disp(i)%molIndx = int(nMol,atomIntType)
        disp(i)%atmIndx = int(i,atomIntType)
        
        disp(i)%x_old => MolArray(nType)%mol(nMol)%x(i)
        disp(i)%y_old => MolArray(nType)%mol(nMol)%y(i)
        disp(i)%z_old => MolArray(nType)%mol(nMol)%z(i)
      enddo
      
!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType) * (2E0 * grnd() - 1E0)
      c_term=cos(angle)
      s_term=sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion. 
      xcm=0E0
      zcm=0E0
      do i=1,nAtoms(nType)
        atmType = atomArray(nType,i)
        xcm = xcm + atomData(atmType)%mass*disp(i)%x_old
        zcm = zcm + atomData(atmType)%mass*disp(i)%z_old
      enddo

      xcm = xcm/totalMass(nType)    
      zcm = zcm/totalMass(nType)

!     Generate a random translational displacement      
      do i=1,nAtoms(nType)
        disp(i)%y_new = disp(i)%y_old
        x_scale = disp(i)%x_old - xcm
        z_scale = disp(i)%z_old - zcm
        disp(i)%x_new = c_term*x_scale - s_term*z_scale + xcm
        disp(i)%z_new = s_term*x_scale + c_term*z_scale + zcm
      enddo

!     Calculate the Energy Difference Associated with the move   
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) return

      biasDiff = 0E0
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
      if(biasEnergy .le. 0E0) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        acptRot(nType) = acptRot(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList,nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call Update_SubEnergies        
        call UpdateDistArray
      elseif(exp(-biasEnergy) .gt. grnd()) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        acptRot(nType) = acptRot(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList, nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      endif

      end subroutine
!=======================================================      
      subroutine Rot_yz(E_T,acc_x,atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none
      
      real(dp), intent(inout) :: E_T,acc_x,atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: i,nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: angle
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd   
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term,s_term
      real(dp) :: y_scale, z_scale
      real(dp) :: ycm,zcm
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if(nAtoms(nType) .eq. 1) then
        return
      endif
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array 
      do i=1,nAtoms(nType)
        disp(i)%molType = int(nType,atomIntType)
        disp(i)%molIndx = int(nMol,atomIntType)
        disp(i)%atmIndx = int(i,atomIntType)
        
        disp(i)%x_old => MolArray(nType)%mol(nMol)%x(i)
        disp(i)%y_old => MolArray(nType)%mol(nMol)%y(i)
        disp(i)%z_old => MolArray(nType)%mol(nMol)%z(i)
      enddo
      
!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd()-1E0)
      c_term=cos(angle)
      s_term=sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion. 
      ycm=0E0
      zcm=0E0
      do i=1,nAtoms(nType)
        atmType = atomArray(nType,i)
        ycm = ycm + atomData(atmType)%mass*disp(i)%y_old
        zcm = zcm + atomData(atmType)%mass*disp(i)%z_old
      enddo

      ycm = ycm/totalMass(nType)
      zcm = zcm/totalMass(nType)

!     Generate a random translational displacement      
      do i=1,nAtoms(nType)
        disp(i)%x_new = disp(i)%x_old
        y_scale = disp(i)%y_old - ycm
        z_scale = disp(i)%z_old - zcm
        disp(i)%y_new = c_term*y_scale - s_term*z_scale + ycm
        disp(i)%z_new = s_term*y_scale + c_term*z_scale + zcm
      enddo

!     Calculate the Energy Difference Associated with the move   
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) return

      biasDiff = 0E0
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
      if(biasEnergy .le. 0E0) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x+1E0
        acptRot(nType) = acptRot(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList,nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      elseif(exp(-biasEnergy) .gt. grnd()) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter
        ETable = ETable + dETable
        acc_x = acc_x+1E0
        acptRot(nType) = acptRot(nType) + 1E0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList,nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
!        call Create_NeiETable
        call UpdateDistArray
        call Update_SubEnergies        
      endif
      end subroutine

!===========================================================================================
      end module
