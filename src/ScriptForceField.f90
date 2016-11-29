 !========================================================  
      module ForceFieldInput
      use VarPrecision
      use Units
     
      abstract interface 
        subroutine CommonSub
          use VarPrecision
          implicit none
        end subroutine
      end interface 
 
      abstract interface 
        subroutine interSub(lineStore)
          implicit none
          character(len=100), intent(in) :: lineStore(:)
        end subroutine
      end interface 




      private
      logical :: fieldTypeSet = .false.
      procedure (CommonSub), pointer :: commonFunction => NULL()
      procedure (interSub), pointer :: interFunction => NULL()
      real(dp) :: convEng
      real(dp) :: convDist
      real(dp) :: convAng


      public :: SetForcefieldType, ScriptForcefield, fieldTypeSet

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
      character(len=10) :: lenDefault, distDefault, angDefaults


      lenDefault = "kb"
      distDefault = "ang"
      angDefaults = "deg"

      convEng = FindEngUnit(lenDefault)
      convDist = FindLengthUnit(distDefault)
      convAng = FindAngularUnit(angDefaults)

      select case(trim(adjustl(potenType)))
      case("lj_q")
        write(nout,*) "Forcefield Type: Standard Lennard-Jones w/ Eletrostatic"
        ForceFieldName = "LJ_Q"
!        call ReadForcefield_LJ_Q
        Detailed_ECalc => Detailed_EnergyCalc_LJ_Q
        Shift_ECalc => Shift_EnergyCalc_LJ_Q
        SwapIn_ECalc => SwapIn_EnergyCalc_LJ_Q
        SwapOut_ECalc => SwapOut_EnergyCalc_LJ_Q
        Rosen_Mol_New => Rosen_BoltzWeight_Molecule_New
        Rosen_Mol_Old => Rosen_BoltzWeight_Molecule_Old
        Quick_Nei_ECalc => QuickNei_ECalc_Inter_LJ_Q
        boundaryFunction => Bound_MaxMin
        commonFunction => Allocate_LJ_Q
        interFunction => Read_LJ_Q
      case("pedone")
        write(nout,*) "Forcefield Type: Pedone"
        ForceFieldName = "Pedone"
!        call ReadForcefield_Pedone
        Detailed_ECalc => Detailed_EnergyCalc_Pedone
        Shift_ECalc => Shift_EnergyCalc_Pedone
        SwapIn_ECalc => SwapIn_EnergyCalc_Pedone
        SwapOut_ECalc => SwapOut_EnergyCalc_Pedone
        Rosen_Mol_New => Rosen_BoltzWeight_Pedone_New
        Rosen_Mol_Old => Rosen_BoltzWeight_Pedone_Old
        Quick_Nei_ECalc => QuickNei_ECalc_Inter_Pedone
        boundaryFunction => Bound_PedoneChargeBalance
!        commonFunction => Allocate_Pedone
!        interFunction => Read_Pedone
      case("custompairwise")
      case default
        stop "Unknown potential type given in forcefield input"
      end select
 

      end subroutine


!========================================================            
      subroutine ScriptForcefield(lineStore)
      use CBMC_Variables
      use Coords
      use EnergyTables
      use ForceField
      use SimParameters
      use Units
      use VarPrecision
      implicit none
      character(len=100), intent(in) :: lineStore(:)

      character(len=25) :: dummy, command, command2, stringValue
      logical :: logicValue
      integer :: iLine, nLines, lineBuffer, lineStat
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

      lineBuffer = 0
      nLines = size(lineStore)
      do iLine = 1, nLines
        if(lineBuffer .gt. 0) then
          lineBuffer = lineBuffer - 1
          cycle
        endif
        call GetCommand(lineStore(iLine), command, lineStat)    
        if(lineStat .gt. 0) then
          cycle
        endif 
        call LowerCaseLine(command)
        select case( adjustl(trim(command)) )
          case("define")
            call FindCommandBlock(iLine, lineStore, "end_define", lineBuffer)
            call DefineForcefield(lineStore(iLine:iLine+lineBuffer) )
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
      use Coords
      use ForceField
      use SimParameters
      use Units
      use VarPrecision
      implicit none
      character(len=100), intent(in) :: lineStore(:)
      character(len=30) :: dummy, defType, stringValue
      logical :: logicValue
      integer :: i, j, iLine, nLines, nParam
      integer :: intValue, AllocateStat
      real(dp) :: realValue

      nLines = size(lineStore)
      read(lineStore(1),*) dummy, defType
      call LowerCaseLine(defType)
      select case( trim(adjustl(defType)) )
        case("atomtypes")
          if( allocated(atomData) ) then
            write(*,*) "ERROR! AtomTypes already defined in the forcefield file"
            stop
          endif
          read(lineStore(1),*) dummy, defType, intValue
          nAtomTypes = intValue
          ALLOCATE (atomData(1:nAtomTypes), STAT = AllocateStat)
          write(*,*) "Here!"
          call commonFunction
          write(*,*) "Here!"
          call interFunction(lineStore)
          write(*,*) "Here!"
        case("bondtypes")
          read(lineStore(1),*) dummy, defType, intValue
          nBondTypes = intValue
          ALLOCATE (bondData(1:nBondTypes), STAT = AllocateStat)
          i = 0
          do iLine = 2, nLines-1
            i = i + 1
            read(lineStore(iLine),*) bondData(i)%bondName, bondData(i)%k_eq, bondData(i)%r_eq
            bondData(i)%k_eq = bondData(i)%k_eq * convEng
            bondData(i)%r_eq = bondData(i)%r_eq * convDist
            if(echoInput) then
              write(35,*) bondData(i)%bondName, bondData(i)%k_eq, bondData(i)%r_eq    
            endif        
          enddo          
        case("angletypes")
          read(lineStore(1),*) dummy, defType, intValue
          nAngleTypes = intValue
          ALLOCATE (bendData(1:nAngleTypes), STAT = AllocateStat)
          i = 0
          do iLine = 2, nLines-1
            i = i + 1
            read(lineStore(iLine),*) bendData(i)%angleName, bendData(i)%k_eq, bendData(i)%ang_eq
            if(echoInput) then
              write(35,*) bendData(i)%angleName, bendData(i)%k_eq, bendData(i)%ang_eq 
            endif          
            bendData(i)%k_eq = bendData(i)%k_eq * convEng
            bendData(i)%ang_eq = bendData(i)%ang_eq * convAng
          enddo   
        case("torsiontypes")
          read(lineStore(1),*) dummy, defType, intValue
          nTorsionalTypes = intValue
          ALLOCATE (torsData(1:nTorsionalTypes), STAT = AllocateStat)
          i = 0
          do iLine = 2, nLines-1
            i = i + 1
            read(lineStore(iLine),*) torsData(i)%torsName, nParam
            allocate(torsData(i)%a(1:nParam),STAT = AllocateStat)
            read(lineStore(iLine),*) torsData(i)%torsName, torsData(i)%nPara, (torsData(i)%a(j),j=1,nParam)
            if(echoInput) then
              write(35,*) torsData(i)%torsName, torsData(i)%nPara, (torsData(i)%a(j),j=1,nParam)
            endif         
            do j=1,nParam
              torsData(i)%a(j) = torsData(i)%a(j) * convEng
            enddo   
          enddo        
        case("rmin")
          if(nAtomTypes .eq. 0) then
            write(*,*) "ERROR! RMin is called before the number of atom types"
            write(*,*) "have been defined"
            stop
          endif
          call SetRMin(lineStore)
        case default 
          write(*,*) "ERROR! Invalid command in Forcefield file on line:"
          write(*,*) lineStore(1)
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
      subroutine SetRMin(lineStore)
      use VarPrecision
      use ForceField
      use ForceFieldFunctions
      use SimParameters
      implicit none
      character(len=100), intent(in) :: lineStore(:)  

      logical :: custom
      integer :: i, j, iLine, nLines
      integer :: indx1, indx2
      character(len=30) :: dummy, dummy2, mixingRule
      procedure (MixRule), pointer :: rmin_func => null()   
      real(dp) :: curVal


      nLines = size(lineStore)
      read(lineStore(1),*) dummy, dummy2, mixingRule

      custom = .false.
      select case(trim(adjustl(mixingRule)))
        case("geometric")
          rmin_func => GeoMean_MixingFunc
        case("average")
           rmin_func => Mean_MixingFunc
        case("min")
           rmin_func => Min_MixingFunc
        case("max")
           rmin_func => Max_MixingFunc
        case("custom")
           custom = .true. 
        case default
          stop "Error! Invalid R_Min Mixing Rule Type"
      end select


      r_min = 0E0_dp
      r_min_tab = 0E0_dp
      if(custom) then           
        do iLine = 2, nLines-1
          read(lineStore(iLine), *) indx1, indx2, curVal 
          if(echoInput) then
            write(35,*) indx1, indx2, curVal
            flush(35)
           endif
           r_min_tab(indx1, indx2) = curVal
           r_min_tab(indx2, indx1) = curVal
        enddo
        do i = 1,nAtomTypes
          do j = 1,nAtomTypes
            r_min_tab(i,j) = r_min_tab(i,j)**2
          enddo
        enddo

      else

        do iLine = 2, nLines-1
          read(lineStore(iLine), *) indx1, curVal 
          r_min(indx1) = curVal
        enddo
        do i = 1,nAtomTypes
          do j = i,nAtomTypes
            r_min_tab(i,j) = rmin_func(r_min(i), r_min(j))**2
            r_min_tab(j,i) = r_min_tab(i,j)
          enddo
        enddo

      endif

      write(35,*) "Rmin Table:"
      do i = 1, nAtomTypes
        write(35,*) (sqrt(r_min_tab(i,j)), j= 1, nAtomTypes)
      enddo
     
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
      command = " "
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
      subroutine Allocate_LJ_Q
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      

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

      ALLOCATE (max_dist(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT = AllocateStatus) 
      
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
      end subroutine
!===================================================================================
      subroutine Read_LJ_Q(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      implicit none
      character(len=100), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines-1
        i = i + 1
        read(lineStore(iLine),*) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
                                 atomData(i)%sig, atomData(i)%q, atomData(i)%mass
      enddo

!      Generate the look up tables for the inter molecular interactions
      if(echoInput) then
        write(35,*) "---------------------------------------------"
        write(35,*) "Interaction Table"
        write(35,*) " i "," j ", " eps ", " sig ", " q "
      endif     
      do i = 1,nAtomTypes
        do j = i,nAtomTypes
          ep_tab(i,j)  = 4d0*ep_func(atomData(i)%ep, atomData(j)%ep)          
          sig_tab(i,j) = sig_func(atomData(i)%sig, atomData(j)%sig)**2
          q_tab(i,j) = atomData(i)%q * atomData(j)%q * 1.671d5

          ep_tab(j,i) = ep_tab(i,j)
          sig_tab(j,i) = sig_tab(i,j)
          q_tab(j,i) = q_tab(i,j)
          if(echoInput) then
            write(35,*) i,j, ep_tab(i,j)/4d0, sqrt(sig_tab(i,j)), q_tab(i,j) 
          endif          
        enddo
      enddo
      if(echoInput) then
        write(35,*) "---------------------------------------------"
        flush(35)
      endif


      end subroutine

!================================================================================
      end module
