      module ForceFieldInput
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

      character(len=30) :: dummy, command, stringValue
      logical :: logicValue
      integer :: iLine, nLines
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      do iLine = 1, nLines

        call GetCommand(lineStore(iLine), command, lineStat)    
        if(lineStat .gt. 0) then
          cycle
        endif 
        select case( trim(adjustl(command)) )
          case("define")
            read(line,*) dummy, command, realValue
          case("create")
            read(line,*) dummy, command, realValue
          case default
            write(*,*) "ERROR! Invalid command in Forcefield file on line:", iLine
            write(*,*) lineStore(iLine)
            stop "INPUT ERROR!"
        end select

      enddo



     
      end subroutine   
!========================================================            
      subroutine DefineForcefield(line)
      use VarPrecision
      use SimParameters
      use CBMC_Variables
      use Coords
      use EnergyTables
      use Units
      implicit none
      character(len=100), intent(in) :: line
      character(len=30) :: dummy, command, stringValue
      logical :: logicValue
      integer :: iLine, nLines
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      do iLine = 1, nLines

        read(line,*) dummy, command, intValue
        select case( trim(adjustl(command)) )
          case("atomtypes")

          case default
            write(*,*) "ERROR! Invalid command in Forcefield file on line:", iLine
            write(*,*) lineStore(iLine)
            stop "INPUT ERROR!"
        end select

      enddo



     
      end subroutine 
!========================================================
      subroutine FindForceBlock(iLine, lineStore, endCommand ,lineBuffer)
      implicit none
      integer, intent(in) :: iLine
      character(len=100), intent(in) :: lineStore(:)      
      character(len=20), intent(in) :: endCommand
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
!================================================================================
      end module
