!=============================================================
      subroutine CBMC(E_T, acc_x, atmp_x)  
      use SimParameters
      use CBMC_Variables
      use Coords
      use E_Interface
      use EnergyCriteria
      use EnergyTables
      use Forcefield
      use DistanceCriteria      
      use IndexingFunctions
      implicit none
      
      real(kind(0.0d0)),intent(inout) :: E_T, acc_x, atmp_x     
    
      logical, parameter :: useIntra(1:4) = [.true., .true., .true., .true.]

      logical :: rejMove      
      logical :: regrown(1:maxAtoms)
      integer :: i, j, nType, nMol, nIndx, nMove, nPath, nAtom
      integer :: iAtom, nAtomLoc
      logical :: regrowDirection
      real(kind(0.0d0)) :: grnd, ranNum, sumInt
      real(kind(0.0d0)) :: dx,dy,dz      
      real(kind(0.0d0)) :: E_Diff
      type (displacement) :: disp(1:maxAtoms)
      real(kind(0.0d0)) :: PairList(1:maxMol)
      real(kind(0.0d0)) :: dETable(1:maxMol)
      real(kind(0.0d0)) :: rosenProb_New, rosenProb_Old, rosenRatio      
      real(kind(0.0d0)) :: E_Inter, E_Intra

      rejMove = .false.
      atmp_x = atmp_x+1d0
!     Randomly Select a Molecule from the cluster and obtain its molecule type and index
      ranNum = grnd()
      sumInt = probTypeCBMC(1)
      nType = 1
      do while(sumInt .lt. ranNum)
        nType = nType + 1
        sumInt = probTypeCBMC(nType)
      enddo
      if(NPART(nType) .eq. 0d0) then
        return
      endif
      nMol = floor(NPART(nType)*grnd() + 1d0)
      nIndx = molArray(nType)%mol(nMol)%indx
      

      select case(regrowType(nType))
      case(1, 2)
        regrown = .true.
!        Begin by selecting an atom in the chain to serve as the regrowth point        
        nAtom = floor(grnd()*pathArray(nType)%pathMax(1) + 1d0)
!        Choose the direction to regrow the chain        
        if(grnd() .gt. 0.5d0) then
          regrowDirection = .true.
        else
          regrowDirection = .false.
        endif
!        If regrowDirection = .true. then regrowth is performed in the positive chain direction, otherwise
!        it is done in the negative direction.  In this step the
        do i =1, pathArray(nType)%pathMax(1)
          if(pathArray(nType)%path(1, i) .eq. nAtom) then
            nAtomLoc = i
          endif
        enddo
        if(regrowDirection) then
         do i = nAtomLoc, pathArray(nType)%pathMax(1)
           iAtom = pathArray(nType)%path(1, i)
           regrown(iAtom) = .false. 
         enddo
        else
          do i = 1, nAtomLoc
            iAtom = pathArray(nType)%path(1, i)
            regrown(iAtom) = .false. 
          enddo
        endif
!        In the event that all atoms have been selected for regrowth, randomly pick one chain's terminal
!        atoms to remain grown.
        if( all(regrown .eqv. .false. ) ) then
          regrown(nAtom) = .true.
        endif
        call StraightChain_Partial_ConfigGen(nType, nMol, regrown(1:maxatoms), regrowDirection, nAtom, rosenProb_New, rejMove)
        if(rejMove) then
          return
        endif
        call StraightChain_Partial_ConfigGen_Reverse(nType, nMol, regrown(1:maxatoms), regrowDirection, nAtom, rosenProb_Old)
      case default                   
        write(*,*) "ERROR! CBMC is not implimented for this regrowth type! :", regrowType(nType)
        stop
      end select
      

      do i=1,nAtoms(nType)
        if(.not. regrown(i)) then
          disp(i)%molType = int(nType,2)
          disp(i)%molIndx = int(nMol,2)
          disp(i)%atmIndx = int(i,2)
        
          disp(i)%x_old => MolArray(nType)%mol(nMol)%x(i)
          disp(i)%y_old => MolArray(nType)%mol(nMol)%y(i)
          disp(i)%z_old => MolArray(nType)%mol(nMol)%z(i)
        
          disp(i)%x_new = newMol%x(i)
          disp(i)%y_new = newMol%y(i)
          disp(i)%z_new = newMol%z(i)
        else
          disp(i)%molType = int(nType,2)
          disp(i)%molIndx = int(nMol,2)
          disp(i)%atmIndx = int(i,2)
        
          disp(i)%x_old => MolArray(nType)%mol(nMol)%x(i)
          disp(i)%y_old => MolArray(nType)%mol(nMol)%y(i)
          disp(i)%z_old => MolArray(nType)%mol(nMol)%z(i)
        
          disp(i)%x_new = disp(i)%x_old
          disp(i)%y_new = disp(i)%y_old
          disp(i)%z_new = disp(i)%z_old
        endif
      enddo

      E_Inter = 0d0
      E_Intra = 0d0
      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, .true., useIntra, rejMove)
      if(rejMove) then
        return      
      endif

      rosenRatio = rosenProb_Old/rosenProb_New
!      write(2,*) E_Inter, E_Intra
      if(rosenRatio*exp(-beta*E_Inter) .gt. grnd()) then
        do i=1,nAtoms(nType)      
          disp(i)%x_old = disp(i)%x_new
          disp(i)%y_old = disp(i)%y_new
          disp(i)%z_old = disp(i)%z_new
        enddo
        E_T = E_T + E_Inter + E_Intra
        ETable = ETable + dETable
        acc_x = acc_x + 1d0
        if(distCriteria) then
          call NeighborUpdate_Distance(PairList, nIndx)        
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
        call Update_SubEnergies        
      endif
      
      end subroutine

