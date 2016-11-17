      module ForceFieldInput
!========================================================  
      contains
!========================================================  
      subroutine ScriptForcefield
      use AnalysisMain, only: ReadAnalysisInput
      use ForceField, only: ForceFieldName
      use EnergyPointers
      use SwapBoundary
      use ParallelVar, only: nout
      use UmbrellaSamplingNew, only: ReadInput_Umbrella
      use WHAM_Functions
      implicit none
      character(len=15) :: labelField 
      character(len=10) :: potenType

      open(unit=55, file="input_forcefield.dat", status='OLD')
      read(55,*) 
      read(55,*) labelField, potenType
      write(35,*) "-------------------------"
      write(35,*) labelField, potenType
      write(35,*) "-------------------------"

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
 
      read(54,*)
      read(54,*)
      read(54,*)
      call ReadAnalysisInput(54)
      read(54,*)
      read(54,*)
      read(54,*)
      call ReadInput_Umbrella(54)
      close(54)                
      if(useWHAM) then
        call WHAM_Initialize
      endif

      end subroutine
!================================================================================
      end module
