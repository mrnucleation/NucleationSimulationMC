 !========================================================  
      module ForceFieldInput

     
      interface 
        subroutine MCMoveSub
          use VarPrecision
          implicit none
        end subroutine
      end interface 
 
      type MoveArray
        procedure(MCMoveSub), pointer, nopass :: commonFunction => NULL()
      end type

!========================================================  
      contains
!========================================================  
      subroutine SetForcefieldType(potenType)
      use AnalysisMain, only: ReadAnalysisInput
      use ForceField, only: ForceFieldName
      use EnergyPointers
      use SwapBoundary
      use ParallelVar, only: nout
      use UmbrellaSamplingNew, only: ReadInput_Umbrella
      use WHAM_Functions
      implicit none
      character(len=10) :: potenType

      select case(trim(adjustl(potenType)))
      case("LJ_Q")
        write(nout,*) "Forcefield Type: Standard Lennard-Jones w/ Eletrostatic"
        ForceFieldName = "LJ_Q"
        call ReadForcefield_LJ_Q
        Detailed_ECalc => Detailed_EnergyCalc_LJ_Q
        Shift_ECalc => Shift_EnergyCalc_LJ_Q
        SwapIn_ECalc => SwapIn_EnergyCalc_LJ_Q
        SwapOut_ECalc => SwapOut_EnergyCalc_LJ_Q
        Rosen_Mol_New => Rosen_BoltzWeight_Molecule_New
        Rosen_Mol_Old => Rosen_BoltzWeight_Molecule_Old
        Quick_Nei_ECalc => QuickNei_ECalc_Inter_LJ_Q
        boundaryFunction => Bound_MaxMin
        commonFunction => Allocate_Common_Variables_LJ_Q
      case("Pedone")
        write(nout,*) "Forcefield Type: Pedone"
        ForceFieldName = "Pedone"
        call ReadForcefield_Pedone
        Detailed_ECalc => Detailed_EnergyCalc_Pedone
        Shift_ECalc => Shift_EnergyCalc_Pedone
        SwapIn_ECalc => SwapIn_EnergyCalc_Pedone
        SwapOut_ECalc => SwapOut_EnergyCalc_Pedone
        Rosen_Mol_New => Rosen_BoltzWeight_Pedone_New
        Rosen_Mol_Old => Rosen_BoltzWeight_Pedone_Old
        Quick_Nei_ECalc => QuickNei_ECalc_Inter_Pedone
        boundaryFunction => Bound_PedoneChargeBalance
        commonFunction => Allocate_Common_Variables_Pedone
      case("custompairwise")
      case default
        stop "Unknown potential type given in forcefield input"
      end select
 

      end subroutine


!========================================================            
      subroutine ScriptForcefield(lineStore)
      use VarPrecision
      use SimParameters
      use CBMC_Variables
      use Coords
      use EnergyTables
      use Units
      implicit none
      character(len=100), intent(in) :: lineStore(:)

      character(len=30) :: dummy, command, command2, stringValue
      logical :: logicValue
      integer :: iLine, nLines
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      if(nMolTypes .eq. 0) then
        write(*,*) "ERROR! Forcefield input has been called before the number of molecule types have been defined"
        stop
      endif

      nAtomTypes = 0
      nBondTypes = 0
      nAngleTypes = 0
      nTorsionalTypes = 0

      do iLine = 1, nLines
        if(lineBuffer .gt. 0) then
          lineBuffer = lineBuffer - 1
          cycle
        endif
        call GetCommand(lineStore(iLine), command, lineStat)    
        if(lineStat .gt. 0) then
          cycle
        endif 
        select case( trim(adjustl(command)) )
          case("define")
            call FindCommandBlock(iLine, lineStore, "end_define", lineBuffer)
            call DefineForcefield(lineStore(iLine:iLine+lineBuffer)
          case("create")
            call FindCommandBlock(iLine, lineStore, "end_create", lineBuffer)

          case default
            write(*,*) "ERROR! Invalid command in Forcefield file on line:", iLine
            write(*,*) lineStore(iLine)
            stop "INPUT ERROR!"
        end select

      enddo



     
      end subroutine   
!========================================================            
      subroutine DefineForcefield(lineStore)
      use VarPrecision
      use SimParameters
      use CBMC_Variables
      use Coords
      use EnergyTables
      use Units
      implicit none
      character(len=100), intent(in) :: lineStore(:)
      character(len=30) :: dummy, defType, stringValue
      logical :: logicValue
      integer :: iLine, nLines
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      read(lineStore(1),*) dummy, defType, intValue
      select case( trim(adjustl(defType)) )
        case("atomtypes")
          if( allocated(atomData) ) then
            write(*,*) "ERROR! AtomTypes already defined in the forcefield file"
            stop
          endif
          nAtomTypes = intValue
          ALLOCATE (atomData(1:nAtomTypes), STAT = AllocateStatus)
          ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
          ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
          ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus) 
        case("bondtypes")
          nBondTypes = intValue
          ALLOCATE (bondData(1:nBondTypes), STAT = AllocateStatus)
        case("angletypes")
          nAngleTypes = intValue
          ALLOCATE (bendData(1:nAngleTypes), STAT = AllocateStatus)
        case("torsiontypes")
          nTorsionalTypes = intValue
          ALLOCATE (torsData(1:nTorsionalTypes), STAT = AllocateStatus)
        case("rmin")
          if(nAtomTypes .eq. 0) then
            write(*,*) "ERROR! RMin is called before the number of atom types"
            write(*,*) "have been defined"
            stop
          endif
          
        case default 
          write(*,*) "ERROR! Invalid command in Forcefield file on line:", iLine
          write(*,*) lineStore(iLine)
          stop "INPUT ERROR!"
      end select
     
      end subroutine 
!========================================================
      subroutine FindCommandBlock(iLine, lineStore, endCommand ,lineBuffer)
      implicit none
      integer, intent(in) :: iLine
      character(len=100), intent(in) :: lineStore(:)      
      character(len=*), intent(in) :: endCommand
      integer, intent(out) :: lineBuffer
      logical :: found
      integer :: i, lineStat, nLines
      character(len=35) :: dummy 


      dummy = " "
      nLines = size(lineStore)
      found = .false.
      do i = iLine + 1, nLines
        call GetCommand(lineStore(i), dummy, lineStat)
        if( trim(adjustl(dummy)) .eq. trim(adjustl(endCommand)) ) then
          lineBuffer = i - iLine
          found = .true.
          exit
        endif
      enddo

      if(.not. found) then
        write(*,*) "ERROR! A command block was opened in the input script, but no closing END statement found!"
        stop
      endif

      end subroutine
!========================================================            
!     This subrotuine searches a given input line for the first command. 
      subroutine GetCommand(line, command, lineStat)
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: line
      character(len=25), intent(out) :: command
      integer, intent(out) :: lineStat
      integer :: i, sizeLine, lowerLim, upperLim

      sizeLine = len(line)
      lineStat = 0
      i = 1
!      Find the first non-blank character in the string
      do while(i .le. sizeLine)
        if(ichar(line(i:i)) .ne. ichar(' ')) then
           !If a non-blank character is found, check first to see if it is the comment character.
          if(ichar(line(i:i)) .eq. ichar('#')) then
            lineStat = 1
            return
          else
            exit
          endif
        endif
        i = i + 1
      enddo
!      If no characters are found the line is empty, 
      if(i .ge. sizeLine) then
        lineStat = 1
        return
      endif
      lowerLim = i

!      
      do while(i .le. sizeLine)
        if(line(i:i) .eq. " ") then
          exit
        endif
        i = i + 1
      enddo
      upperLim = i

      command = line(lowerLim:upperLim)
     
      end subroutine
!================================================================ 
      subroutine Allocate_Common_Variables_LJ_Q
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      
      ALLOCATE (atomData(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (bondData(1:nBondTypes), STAT = AllocateStatus)
      ALLOCATE (bendData(1:nAngleTypes), STAT = AllocateStatus)
      ALLOCATE (torsData(1:nTorsionalTypes), STAT = AllocateStatus)
      ALLOCATE (impropData(1:nImproperTypes), STAT = AllocateStatus)

      ALLOCATE (r_min(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT = AllocateStatus) 

      ALLOCATE (ep_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes), STAT = AllocateStatus)
     
      ep_tab = 0d0
      sig_tab = 0d0       
      q_tab = 0d0
      r_min_tab = 0d0

      ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)
      
      ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (nIntraNonBond(1:nMolTypes), STAT = AllocateStatus)      
      ALLOCATE (nBonds(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (nAngles(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (nTorsional(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (nImproper(1:nMolTypes), STAT = AllocateStatus)

      ALLOCATE (acptTrans(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (acptRot(1:nMolTypes),   STAT = AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes),  STAT = AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (acptSwapOut(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (acptInSize(1:maxMol), STAT = AllocateStatus)
      ALLOCATE (atmpInSize(1:maxMol), STAT = AllocateStatus)


      ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus) 
      
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
      end subroutine
!================================================================================
      end module
