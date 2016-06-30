      module CBMC_Utility
      contains
!=======================================================================
!     The purpose of this subroutine is to create a list of molecules
!     that will be included in the intermolecular component of the Rosenbluth
!     weight.  
      subroutine Rosen_CreateSubset(nTarget, included)
      use Coords
      use ForceField
      use Constants
      use SimParameters
      implicit none
      integer, intent(in) :: nTarget
      logical, intent(out) :: included(1:maxMol)
      
      logical:: inc_firstPass(1:maxMol)
      integer :: jIndx, kIndx
      
      included = .false.
      inc_firstPass = .false.
      
!      Loop 1: Adds the neighbors of the target molecule to the list      
      do jIndx = 1, maxMol
         if(.not. isActive(jIndx)) then
           cycle
         endif
         if(NeighborList(nTarget, jIndx)) then
            inc_firstPass(jIndx) = .true.
            included(jIndx) = .true.         
         endif
      enddo


!      Loop 2: Adds the neighbors of the molecules added in the previous step to the list      
!      do jIndx = 1, maxMol
!       if(.not. isActive(jIndx)) then
!         cycle
!       endif
!       if(inc_firstPass(jIndx) .eqv. .true.) then
!         do kIndx = 1, maxMol
!           if(.not. isActive(kIndx)) then
!             cycle
!           endif           
!           if(NeighborList(jIndx, kIndx)) then
!             included(kIndx) = .true.         
!           endif
!         enddo
!       endif
!      enddo

      included(nTarget) = .true.
      
      
      end subroutine 
!=======================================================================
!     This function is similar to the Rosen_CreateSubset function.  The
!     only difference being this is intended for calculating the reverse probability. Thus it
!     will excluse nIndx from the subset.
      subroutine Rosen_CreateSubset_Reverse(nTarget, nIndx ,included)
      use Coords
      use ForceField
      use Constants
      use SimParameters
      implicit none
      integer, intent(in) :: nTarget, nIndx
      logical, intent(out) :: included(1:maxMol)

      logical:: inc_firstPass(1:maxMol)
      integer :: jIndx, kIndx
      
      included = .false.
      inc_firstPass = .false.
      
!      Loop 1: Adds the neighbors of the target molecule to the list      
      do jIndx = 1, maxMol
         if(NeighborList(nTarget, jIndx)) then
            if(.not. isActive(jIndx)) then
              cycle
            endif   
            if(jIndx .eq. nIndx) then
              cycle
            endif                   
            inc_firstPass(jIndx) = .true.
            included(jIndx) = .true.         
         endif
      enddo


!      Loop 2: Adds the neighbors of the molecules added in the previous step to the list      
!      do jIndx = 1, maxMol
!       if(inc_firstPass(jIndx) .eqv. .true.) then
!         if(.not. isActive(jIndx)) then
!           cycle
!         endif       
!         do kIndx = 1, maxMol
!           if(NeighborList(jIndx, kIndx)) then
!             if(.not. isActive(kIndx)) then
!               cycle
!             endif           
!             if(kIndx .eq. nIndx) then
!               cycle
!             endif              
!             included(kIndx) = .true.         
!           endif
!         enddo
!       endif
!      enddo
      
      included(nTarget) = .true.
      
      end subroutine 
!=======================================================================
      subroutine FindAtomsFromPath(nType, regrown, pathNum, Atm4, Atm1, Atm2, Atm3)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer, intent(in) :: Atm4, pathNum, nType
      logical, intent(in) :: regrown(1:maxAtoms)
      integer, intent(out) ::  Atm1, Atm2, Atm3

      integer :: i, atm4Pos
      integer :: atm4plus1, atm4minus1

      do i = 1, pathArray(nType)%pathMax(pathNum)
        if(pathArray(nType)%path(pathNum,i) .eq. Atm4 ) then
          atm4Pos = i
          exit
        endif
      enddo
  
      atm4plus1 = atm4Pos + 1
      atm4minus1 = atm4Pos - 1

      if(atm4minus1 .lt. 1) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos+1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos+2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos+3)
        return
      endif

      if(atm4plus1 .gt. pathArray(nType)%pathMax(pathNum)) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos-1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos-2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos-3)
        return
      endif

      if( regrown(pathArray(nType)%path(pathNum,atm4plus1)) ) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos+1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos+2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos+3)
        return
      endif


      if(.not. regrown(pathArray(nType)%path(pathNum,atm4minus1))) then
        Atm3 = pathArray(nType)%path(pathNum,atm4Pos-1)
        Atm2 = pathArray(nType)%path(pathNum,atm4Pos-2)
        Atm1 = pathArray(nType)%path(pathNum,atm4Pos-3)
        return
      endif

      write(*,*) "ERROR in FindAtomsFromPath subroutine! Could not find the torsional members"
      write(*,*) "from the stored pathway for atom:", Atm4
      stop
      
      end subroutine 
!=======================================================================
      subroutine FindNextBondFromPath(nType, regrown, pathNum, Atm2, Atm1)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer, intent(in) :: Atm2, pathNum, nType
      logical, intent(in) :: regrown(1:maxAtoms)
      integer, intent(out) ::  Atm1

      integer :: i, atm2Pos
      integer :: atm2plus1, atm2minus1

      do i = 1, pathArray(nType)%pathMax(pathNum)
        if(pathArray(nType)%path(pathNum,i) .eq. Atm2 ) then
          atm2Pos = i
          exit
        endif
      enddo
  
      atm2plus1 = atm2Pos + 1
      atm2minus1 = atm2Pos - 1

      if(atm2minus1 .lt. 1) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos+1)
        return
      endif

      if(atm2plus1 .gt. pathArray(nType)%pathMax(pathNum)) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos-1)
        return
      endif

      if( regrown(pathArray(nType)%path(pathNum,atm2plus1)) ) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos+1)
        return
      endif


      if(.not. regrown(pathArray(nType)%path(pathNum,atm2minus1))) then
        Atm1 = pathArray(nType)%path(pathNum,atm2Pos-1)
        return
      endif

      write(*,*) "ERROR in FindNextBondFromPath subroutine! Could not find the bond members"
      write(*,*) "from the stored pathway for atom:", Atm2
      stop
      
      end subroutine 
!=======================================================================
      subroutine FindNextAngleFromPath(nType, regrown, pathNum, Atm3, Atm1, Atm2)
      use Forcefield
      use SimParameters
      use Coords
      use CBMC_Variables
      implicit none
      integer, intent(in) :: Atm3, pathNum, nType
      logical, intent(in) :: regrown(1:maxAtoms)
      integer, intent(out) ::  Atm1, Atm2

      integer :: i, atm3Pos
      integer :: atm3plus1, atm3minus1

      do i = 1, pathArray(nType)%pathMax(pathNum)
        if(pathArray(nType)%path(pathNum,i) .eq. Atm3 ) then
          atm3Pos = i
          exit
        endif
      enddo
  
      atm3plus1 = atm3Pos + 1
      atm3minus1 = atm3Pos - 1

      if(atm3minus1 .lt. 1) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos+1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos+2)
        return
      endif

      if(atm3plus1 .gt. pathArray(nType)%pathMax(pathNum)) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos-1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos-2)
        return
      endif

      if( regrown(pathArray(nType)%path(pathNum,atm3plus1)) ) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos+1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos+2)
        return
      endif


      if(.not. regrown(pathArray(nType)%path(pathNum,atm3minus1))) then
        Atm2 = pathArray(nType)%path(pathNum,atm3Pos-1)
        Atm1 = pathArray(nType)%path(pathNum,atm3Pos-2)
        return
      endif

      write(*,*) "ERROR in FindNextAngleFromPath subroutine! Could not find the torsional members"
      write(*,*) "from the stored pathway for atom:", Atm3
      stop
      
      end subroutine 
!=======================================================================
      subroutine ChooseRosenTrial(atomNumber, nType, trialPos, regrown, E_Trial, overlap, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use CBMC_Variables
      implicit none
      logical, intent(in) :: overlap(:)
      integer, intent(in) :: atomNumber, nType
      type(SimpleAtomCoords), intent(in) :: trialPos(:)
      logical, intent(inout) :: regrown(:)
      real(dp), intent(inout) :: E_Trial(:)
      real(dp), intent(inout) :: rosenRatio

 
      logical,intent(out) :: rejMove
      integer :: nSel, iRosen
      real(dp) :: E_Min, ProbRosen(1:maxRosenTrial), rosenNorm, grnd
      real(dp) :: ranNum, sumInt

      rejMove = .false.
      E_Min = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
         ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Min))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)/rosenNorm
      regrown(atomNumber) = .true.
      newMol%x(atomNumber) = trialPos(nSel)%x 
      newMol%y(atomNumber) = trialPos(nSel)%y 
      newMol%z(atomNumber) = trialPos(nSel)%z

      end subroutine
!=======================================================================
      end module
