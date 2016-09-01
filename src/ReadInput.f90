!========================================================            
      subroutine ReadParameters(seed,outFreq_Traj,outFreq_Screen,outFreq_GCD,screenEcho)
      use CBMC_Variables
      use Constants
      use Coords
      use EnergyTables
      use ForceField
      use SimParameters
      use Units
      use MoveTypeModule
      use ParallelVar
      use WHAM_Module
      use UmbrellaFunctions
      use VarPrecision
      implicit none
      
      logical, intent(OUT)  :: screenEcho
      integer, intent(OUT) :: seed
!      integer(kind=8), intent(OUT) :: ncycle,nmoves
      integer, intent(OUT)  :: outFreq_Traj,outFreq_Screen,outFreq_GCD
      
      integer :: i,j
      logical :: useBias   
      integer :: AllocateStatus            
      character(len=30) :: labelField 
      character(len=30) :: moveName_temp
      character(len=30) :: fileName
      real(dp) :: dummyCycle
      real(dp) :: norm

!      integer, allocatable :: NArray(:)

      open(unit=54,file="input_Parameters.dat",status='OLD')       

!        Comment Space     
      read(54,*)
      read(54,*)
      read(54,*)

!      Number of Cycles    
      read(54,*) labelField, seed
      read(54,*) labelField, dummyCycle
       NCYCLE=nint(dummyCycle)
      read(54,*) labelField, dummyCycle
       ncycle2=nint(dummyCycle)
      read(54,*) labelField, nMolTypes
   


      allocate( NPART(1:nMolTypes),STAT = AllocateStatus )    
      allocate( NPART_new(1:nMolTypes),STAT = AllocateStatus )   
      allocate( NMIN(1:nMolTypes),STAT = AllocateStatus )     
      allocate( NMAX(1:nMolTypes),STAT = AllocateStatus )     
      allocate( gas_dens(1:nMolTypes),STAT = AllocateStatus )      
      allocate( nRosenTrials(1:nMolTypes),STAT = AllocateStatus )       
      
      read(54,*) labelField, (NMIN(j),j=1,nMolTypes) 
      read(54,*) labelField, (NMAX(j),j=1,nMolTypes)
      maxMol = sum(NMAX)      
      do j=1,nMolTypes
        if(NMIN(j) .gt. NMAX(j) ) then
          write(nout,*) "ERROR! NMin is greater than NMax in input file"
          write(nout,*) "Molecule Type:", j
          write(nout,*) "NMin:", NMIN(j)
          write(nout,*) "NMax:", NMAX(j)
          stop "Input File Error! See Report for details"          
        endif        
      enddo


      
      
      read(54,*) labelField, temperature 
       beta = 1d0/temperature
      read(54,*) labelField,screenEcho
      read(54,*) labelField, useBias
      if(useBias) then
        backspace(54)
        read(54,*) labelField, useBias, fileName
!        call AllocateUmbrellaBias(fileName)
      else 
!        call BlankUmbrellaBias
      endif
      read(54,*) labelField, multipleInput
      read(54,*) labelField, (nRosenTrials(j),j=1,nMolTypes) 
      maxRosenTrial = maxval(nRosenTrials)
      allocate(rosenTrial(1:maxRosenTrial))

      read(54,*) labelField, Dist_Critr
      read(54,*) labelField, softCutoff

!       Comment Space      
      read(54,*)
      read(54,*)
      read(54,*)           
!      Monte Carlo Move Parameters   
      read(54,*) labelField, (gas_dens(j),j=1,nMolTypes)    

!       Comment Space             
      read(54,*)
      read(54,*)
      read(54,*)  
      read(54,*) labelField, outFreq_Screen  
      read(54,*) labelField, outFreq_Traj
      call GCD_Calc(outFreq_Screen,outFreq_Traj,outFreq_GCD)

      read(54,*) labelField, outputEngUnits
      outputEConv = FindEngUnit(outputEngUnits)
      read(54,*) labelField, outputLenUnits   
      outputLenConv = FindLengthUnit(outputLenUnits)     
!       Comment Space      
      read(54,*)
      read(54,*)
      read(54,*) 

       ! This block allows the user to define which types of Monte Carlo moves will be used
       ! during the simulation.
      norm = 0d0
      read(54,*) labelField, nMoveTypes
      if(nMoveTypes .le. 0) then
        write(*,*) "ERROR! The user has specified an invalid number of Monte Carlo moves"
        write(*,*) "Please specify at least one valid Monte Carlo move to continue"
        write(*,*) labelField, nMoveTypes
        stop
      endif

      allocate(mcMoveArray(1:nMoveTypes), STAT = AllocateStatus)
      allocate(moveProbability(1:nMoveTypes), STAT = AllocateStatus)
      allocate(movesAccepted(1:nMoveTypes), STAT = AllocateStatus)
      allocate(movesAttempt(1:nMoveTypes), STAT = AllocateStatus)
      allocate(accptRate(1:nMoveTypes), STAT = AllocateStatus)
      allocate(moveName(1:nMoveTypes), STAT = AllocateStatus)
      norm = 0d0
      avbmcUsed = .false.
      cbmcUsed = .false.
      do i = 1, nMoveTypes
        read(54,*) moveName_temp, moveProbability(i)
        norm = norm + moveProbability(i)
        select case( trim(adjustl(moveName_temp)) )
        case("translation")
          mcMoveArray(i) % moveFunction => Translation
          moveName(i) = "Translation"
        case("rotation")
          mcMoveArray(i) % moveFunction => Rotation
          moveName(i) = "Rotation"
        case("avbmc")
          mcMoveArray(i) % moveFunction => AVBMC
          moveName(i) = "AVBMC"
          avbmcUsed = .true.
        case("cbmc")
          mcMoveArray(i) % moveFunction => CBMC
          moveName(i) = "CBMC"
          cbmcUsed = .true.
        case("exchange")
          mcMoveArray(i) % moveFunction => Exchange
          moveName(i) = "Exchange"
        case("singleatom_translation")
          mcMoveArray(i) % moveFunction => SingleAtom_Translation
          moveName(i) = "Single Atom Translation"
        case default
          write(*,*) "ERROR! Invalid move type specified in input file"
          write(*,*) moveName, moveProbability(i)
          stop
        end select
!        moveName(i) = moveName_temp
      enddo

      do i =1, nMoveTypes
        moveProbability(i) = moveProbability(i)/norm
      enddo
      if(nMoveTypes .gt. 1) then
        do i = 2, nMoveTypes
          moveProbability(i) = moveProbability(i) + moveProbability(i-1)
        enddo
      endif

!       Comment Space      
      read(54,*)
      read(54,*)
      read(54,*) 

      allocate(refSizeNumbers(1:nMolTypes),STAT = AllocateStatus)

      read(54,*) labelField, useWHAM
      read(54,*) labelField, (refSizeNumbers(j), j=1,nMolTypes)
!      refBin = getBiasIndex(refSizeNumbers, NMAX)
      refBin = 1
      read(54,*) labelField, intervalWHAM
      read(54,*) labelField, tolLimit
      read(54,*) labelField, maxSelfConsist
      read(54,*) labelField, whamEstInterval
      read(54,*) labelField, equilInterval
      if(useWham) then
        nWhamItter = ceiling(dble(ncycle)/dble(intervalWHAM))
!        call WHAM_Initialize
      endif

      
      read(54,*)
      read(54,*)
      read(54,*)
!     Cluster Critiera Parameters
      allocate( Eng_Critr(1:nMolTypes,1:nMolTypes), stat = AllocateStatus) 
      read(54,*)
      do i=1,nMolTypes      
        read(54,*) (Eng_Critr(i,j),j=1,nMolTypes)
      enddo
      read(54,*)      
      
      allocate(biasAlpha(1:nMolTypes,1:nMolTypes), stat = AllocateStatus)      
      do i = 1,nMolTypes
         read(54,*) (biasAlpha(i,j),j=1,nMolTypes)     
         do j = 1, nMolTypes
            biasAlpha(i,j) = biasAlpha(i,j)/temperature
         enddo         
      enddo      


      end subroutine
!========================================================  
      subroutine ReadForcefield
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
!========================================================  
!     This subroutine reads the user specified forcefield.  The first half of this
!     code deals with defining the types of atoms,bonds, etc. to be used in the simulation.
!     This is done to such that degenerate angles only need to be defined once. 
!     The second half of this routine deals with defining the molecule configuration using
!     the atoms, angles, etc. defined in the first section of the routine.
      subroutine ReadForcefield_LJ_Q
      use CBMC_Variables
      use Constants
      use Coords
      use ForceField
      use ForceFieldPara_LJ_Q
      use ForceFieldFunctions
      use PairStorage, only:createDistArrays, SetStorageFlags
      use SimParameters
      use Units
      implicit none
      logical :: custom
      integer i,j,h,AllocateStatus,nParam
      integer :: nAtomsMax, nBondsMax,nAnglesMax    
      integer :: nTorsMax, nImpropMax, nNonBondMax
      character(len=15) :: labelField 
      character(len=10) :: mixingRule
      character(len=10) :: unitsEnergy,unitsDistance,unitsAngular
      real(dp) :: convEng, convDist, convAng
      procedure (MixRule), pointer :: ep_func => null()
      procedure (MixRule), pointer :: sig_func => null()      
      procedure (MixRule), pointer :: rmin_func => null()   

!     Open forcefield file      
!      open(unit=55,file="input_forcefield.dat",status='OLD')
!     Blank Space for Comments
      read(55,*) 
      read(55,*)
      read(55,*)   
       
!      Specifies the number of atoms, bonds, etc. types to be used in the simulation.
      read(55,*) labelField, nAtomTypes
      read(55,*) labelField, nBondTypes
      read(55,*) labelField, nAngleTypes
      read(55,*) labelField, nTorsionalTypes
      read(55,*) labelField, nImproperTypes
       
      call Allocate_Common_Variables_LJ_Q

!     Blank Space for Comments
      read(55,*) 
      read(55,*)

!     Epsilon Mixing Rule
      read(55,*) labelField, mixingRule 
      select case(mixingRule)
        case("geometric")
           ep_func => GeoMean_MixingFunc
!           write(6,*) labelField
        case("average")
           ep_func => Mean_MixingFunc
        case("min")
           ep_func => Min_MixingFunc
        case("max")
           ep_func => Max_MixingFunc
!        case("fromfile")
          !Not Yet Implimented
        case default
          stop "Error! Invalid Epsilon Mixing Rule Type"
      end select

!     Sigma Mixing Rule    
      read(55,*) labelField, mixingRule 
!       write(6,*) labelField, mixingRule      
      select case(mixingRule)
        case("geometric")
           sig_func => GeoMean_MixingFunc
        case("average")
           sig_func => Mean_MixingFunc
        case("min")
           sig_func => Min_MixingFunc
        case("max")
           sig_func => Max_MixingFunc
!        case("fromfile")
          !Not Yet Implimented
        case default
          stop "Error! Invalid Sigma Mixing Rule Type"
      end select
      
      read(55,*) 
      read(55,*)      
      read(55,*) labelField, unitsEnergy
      read(55,*) labelField, unitsDistance
      read(55,*) labelField, unitsAngular
      convEng = FindEngUnit(unitsEnergy)
      convDist = FindLengthUnit(unitsDistance)
      convAng = FindAngularUnit(unitsAngular)
       
!      Begin reading the intermolecular forcefield from the file
      read(55,*)
      read(55,*)      
      do i = 1,nAtomTypes
        read(55,*) labelField, atomData(i)%Symb, atomData(i)%ep, atomData(i)%sig, atomData(i)%q, atomData(i)%mass
        if(echoInput) then
          write(35,*) labelField, atomData(i)%Symb, atomData(i)%ep, atomData(i)%sig, atomData(i)%q, atomData(i)%mass
        endif
        atomData(i)%ep = atomData(i)%ep * convEng
        atomData(i)%sig = atomData(i)%sig * convDist   
!        r_min_sq(i) = r_min(i)*r_min(i)

      enddo
      
!      Generate the look up tables for the inter molecular interactions
      if(echoInput) then
        write(35,*) "---------------------------------------------"
        write(35,*) "Interaction Table"
        write(35,*) " i "," j ", " eps ", " sig ", " q "
      endif     
      do i = 1,nAtomTypes
        do j = i,nAtomTypes
!          ep_tab(i,j)  = ep_func(atomData(i)%ep, atomData(j)%ep)
          ep_tab(i,j)  = 4d0*ep_func(atomData(i)%ep, atomData(j)%ep)          
          sig_tab(i,j) = sig_func(atomData(i)%sig, atomData(j)%sig)**2
          q_tab(i,j) = atomData(i)%q * atomData(j)%q * 1.671d5
!          r_min_tab(i,j) = ( (r_min(i)+r_min(j))/2d0 )**2
!          r_min_tab(i,j) = min(r_min(i),r_min(j))**2
!           r_min_tab(i,j) = rmin_func(r_min(i), r_min(j))**2

          ep_tab(j,i) = ep_tab(i,j)
          sig_tab(j,i) = sig_tab(i,j)
          q_tab(j,i) = q_tab(i,j)
          if(echoInput) then
            write(35,*) i,j, ep_tab(i,j)/4d0, sqrt(sig_tab(i,j)), q_tab(i,j) 
          endif          
!          r_min_tab(j,i) = r_min_tab(i,j)
        enddo
      enddo
      if(echoInput) then
        write(35,*) "---------------------------------------------"
        flush(35)
      endif
      
!     R_Min Mixing Rule    
      read(55,*)
      read(55,*) labelField, mixingRule 
!       write(6,*) labelField, mixingRule      
      custom = .false.
      select case(mixingRule)
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

      if(custom) then           
        do i = 1,nAtomTypes
          read(55,*) (r_min_tab(i,j), j=1,nAtomTypes)
          if(echoInput) then
            write(35,*) (r_min_tab(i,j), j=1,nAtomTypes)
            flush(35)
           endif
        enddo
        do i = 1,nAtomTypes
          do j = 1,nAtomTypes
            r_min_tab(i,j) = r_min_tab(i,j)**2
          enddo
        enddo
      else
        do i = 1,nAtomTypes
          read(55,*) r_min(i)
        enddo
        do i = 1,nAtomTypes
          do j = i,nAtomTypes
            r_min_tab(i,j) = rmin_func(r_min(i), r_min(j))**2
            r_min_tab(j,i) = r_min_tab(i,j)
          enddo
        enddo
      endif

      global_r_min = minval(r_min)
 
      write(35,*) "Rmin Table:"
      do i = 1, nAtomTypes
        write(35,*) (sqrt(r_min_tab(i,j)), j= 1, nAtomTypes)
      enddo

 
!     Begin reading the Bond type definitions from file
      read(55,*)
      read(55,*)      
      do i=1,nBondTypes
        read(55,*) labelField, bondData(i)%k_eq, bondData(i)%r_eq
        bondData(i)%k_eq = bondData(i)%k_eq * convEng
        bondData(i)%r_eq = bondData(i)%r_eq * convDist
        if(echoInput) then
          write(35,*) labelField, bondData(i)%k_eq, bondData(i)%r_eq    
        endif        
      enddo
!     Begin reading the Bending Angle type definitions from file
      read(55,*)
      read(55,*)      
      do i=1,nAngleTypes
        read(55,*) labelField, bendData(i)%k_eq, bendData(i)%ang_eq
        if(echoInput) then
          write(35,*) labelField, bendData(i)%k_eq, bendData(i)%ang_eq 
        endif          
        bendData(i)%k_eq = bendData(i)%k_eq * convEng
        bendData(i)%ang_eq = bendData(i)%ang_eq * convAng
      enddo      

!     Begin reading the Torsional Angle type definitions from file
      read(55,*)
      read(55,*)      
      do i=1,nTorsionalTypes
        read(55,*) labelField, nParam
        backspace(55)
        allocate(torsData(i)%a(1:nParam),STAT = AllocateStatus)
        read(55,*) labelField, torsData(i)%nPara, (torsData(i)%a(j),j=1,nParam)
        if(echoInput) then
          write(35,*) labelField, torsData(i)%nPara, (torsData(i)%a(j),j=1,nParam)
        endif         
        do j=1,nParam
          torsData(i)%a(j) = torsData(i)%a(j) * convEng
        enddo   
      enddo          

!     Begin reading the Improper Angle type definitions from file
      read(55,*)
      read(55,*)      
      do i=1,nImproperTypes
        read(55,*) labelField, nParam
        backspace(55)
        allocate(impropData(i)%a(1:nParam), STAT = AllocateStatus)
        read(55,*) labelField, impropData(i)%nPara, (impropData(i)%a(j),j=1,nParam)
        do j=1,nParam
          impropData(i)%a(j) = impropData(i)%a(j) * convEng
        enddo
      enddo                

!     Blank space for Comments
      read(55,*) 
       
!     This section constructs the internal molecular configuration 
!     using the previously defined atom, bond, bend, etc. types.
      do i=1,nMolTypes
        read(55,*)
        read(55,*)       
        read(55,*) labelField, nAtoms(i)
        read(55,*) labelField, nIntraNonBond(i)        
        read(55,*) labelField, nBonds(i)
        read(55,*) labelField, nAngles(i)
        read(55,*) labelField, nTorsional(i)
        read(55,*) labelField, nImproper(i)       
      enddo

      nAtomsMax = maxval(nAtoms)
      nNonBondMax = maxval(nIntraNonBond)      
      nBondsMax = maxval(nBonds)      
      nAnglesMax = maxval(nAngles)
      nTorsMax = maxval(nTorsional)
      nImpropMax= maxval(nImproper)
  
      ALLOCATE (atomArray(1:nMolTypes,1:nAtomsMax), STAT = AllocateStatus)
      ALLOCATE (nonBondArray(1:nMolTypes,1:nNonBondMax), STAT = AllocateStatus)     
      ALLOCATE (bondArray(1:nMolTypes,1:nBondsMax), STAT = AllocateStatus)
      ALLOCATE (bendArray(1:nMolTypes,1:nAnglesMax), STAT = AllocateStatus)
      ALLOCATE (torsArray(1:nMolTypes,1:nTorsMax), STAT = AllocateStatus)
      ALLOCATE (impropArray(1:nMolTypes,1:nImpropMax), STAT = AllocateStatus)
  
!      ALLOCATE (atomUseByBond(1:nMolTypes,1:nAtomsMax), STAT = AllocateStatus)
!      ALLOCATE (atomUseByBend(1:nMolTypes,1:nAtomsMax), STAT = AllocateStatus)
!      ALLOCATE (atomUseByTors(1:nMolTypes,1:nAtomsMax), STAT = AllocateStatus)
!      ALLOCATE (atomUseByImprop(1:nMolTypes,1:nAtomsMax), STAT = AllocateStatus)

      do i = 1, maxRosenTrial
        allocate(rosenTrial(i)%x(1:nAtomsMax))      
        allocate(rosenTrial(i)%y(1:nAtomsMax)) 
        allocate(rosenTrial(i)%z(1:nAtomsMax))         
      enddo
     
      do i=1,nMolTypes
        !Atom Definition for each molecular species
        read(55,*)       
        read(55,*)
        read(55,*)
        do j=1,nAtoms(i)  
          read(55,*) labelField, atomArray(i,j)
          if(echoInput) then
             write(35,*) labelField, atomArray(i,j)
           endif           
        enddo
!        Intra Nonbonded (1-5) interaction definition for each molecular species
        read(55,*)
        read(55,*)
        do j=1,nIntraNonBond(i) 
          read(55,*) labelField, (nonBondArray(i,j)%nonMembr(h),h=1,2)
          if(echoInput) then
            write(35,*) labelField, (nonBondArray(i,j)%nonMembr(h),h=1,2)
          endif                
        enddo        
        !Bond Definition for each molecular species
        read(55,*)
        read(55,*)
        do j=1,nBonds(i)  
          read(55,*) labelField, bondArray(i,j)%bondType, (bondArray(i,j)%bondMembr(h),h=1,2)
          if(echoInput) then
            write(35,*) labelField, bondArray(i,j)%bondType, (bondArray(i,j)%bondMembr(h),h=1,2)
          endif                 
        enddo
        !Bending Angle Definition for each molecular species          
        read(55,*)
        read(55,*)
        do j=1,nAngles(i)  
          read(55,*) labelField, bendArray(i,j)%bendType, (bendArray(i,j)%bendMembr(h),h=1,3)
          if(echoInput) then
            write(35,*) labelField, bendArray(i,j)%bendType, (bendArray(i,j)%bendMembr(h),h=1,3)
          endif          
        enddo
        !Torsion Angle Definition for each molecular species          
        read(55,*)
        read(55,*)
        do j=1,nTorsional(i)  
          read(55,*) labelField, torsArray(i,j)%TorsType, (torsArray(i,j)%torsMembr(h),h=1,4)
          if(echoInput) then
            write(35,*) labelField, torsArray(i,j)%TorsType, (torsArray(i,j)%torsMembr(h),h=1,4)
          endif       
        enddo
        !Improper Angle Definition for each molecular species         
        read(55,*)
        read(55,*)
        do j=1,nImproper(i)  
          read(55,*) labelField, impropArray(i,j)%ImpropType, (impropArray(i,j)%impropMembr(h),h=1,4)
          if(echoInput) then
            write(35,*) labelField, impropArray(i,j)%ImpropType, (impropArray(i,j)%impropMembr(h),h=1,4)
          endif        
        enddo
      enddo

     
      vmdAtoms = 0
      do i = 1,nMolTypes
        vmdAtoms = vmdAtoms + (NMAX(i)*nAtoms(i))
      enddo
      
      totalMass = 0d0
      do i = 1, nMolTypes
        do j = 1,nAtoms(i)
          totalMass(i) = totalMass(i) + atomData(atomArray(i,j))%mass
        enddo
      enddo
      
      call AllocateCoordinateArrays
      call CreateJointArray  
      call SetStorageFlags(q_tab)
      call IntegrateBendAngleProb

      close(55)

 
      
      end subroutine
!========================================================  
!     This subroutine reads the user specified forcefield.  The first half of this
!     code deals with defining the types of atoms,bonds, etc. to be used in the simulation.
!     This is done to such that degenerate angles only need to be defined once. 
!     The second half of this routine deals with defining the molecule configuration using
!     the atoms, angles, etc. defined in the first section of the routine.
      subroutine ReadForcefield_Pedone
      use SimParameters
      use Constants
      use Coords
      use ForceField
      use ForceFieldPara_Pedone
      use ForceFieldFunctions
      use Units
      use CBMC_Variables
      use PairStorage, only: CreateDistArrays, SetStorageFlags
      implicit none
      logical :: custom
      integer :: i,j,h,AllocateStatus,nParam
      integer :: nAtomsMax, nBondsMax,nAnglesMax    
      integer :: nTorsMax, nImpropMax, nNonBondMax
      character(len=15) :: labelField 
      character(len=10) :: mixingRule
      character(len=10) :: unitsEnergy,unitsDistance,unitsAngular
      real(dp) :: convEng, convDist, convAng
      procedure (MixRule), pointer :: rmin_func => null()   

!     Open forcefield file      
!      open(unit=55,file="input_forcefield.dat",status='OLD')
!     Blank Space for Comments

      call Allocate_Common_Variables_Pedone

      nAtoms = 1

      nAtomTypes = nMolTypes

      read(55,*) 
      read(55,*)
      read(55,*)   
       
!      Specifies the number of atoms, bonds, etc. types to be used in the simulation.



      read(55,*) labelField, unitsEnergy
      read(55,*) labelField, unitsDistance
      read(55,*) labelField, unitsAngular
      convEng = FindEngUnit(unitsEnergy)
      convDist = FindLengthUnit(unitsDistance)
      convAng = FindAngularUnit(unitsAngular)

      read(55,*)
      read(55,*)      
      do i = 1, nMolTypes
        read(55,*) labelField, pedoneData(i)%Symb, pedoneData(i)%repul, pedoneData(i)%rEq, &
                    pedoneData(i)%alpha, pedoneData(i)%delta, pedoneData(i)%q, pedoneData(i)%mass
        if(echoInput) then
          write(35,*) labelField, pedoneData(i)%Symb, pedoneData(i)%repul, pedoneData(i)%rEq, &
                    pedoneData(i)%alpha, pedoneData(i)%delta, pedoneData(i)%q, pedoneData(i)%mass
        endif
        pedoneData(i)%repul = pedoneData(i)%repul * convEng
        pedoneData(i)%delta = pedoneData(i)%delta * convEng
        pedoneData(i)%rEq = pedoneData(i)%rEq * convDist   
!        r_min_sq(i) = r_min(i)*r_min(i)

      enddo

      q_tab = 0d0
      do i = 1,nAtomTypes
        do j = i,nAtomTypes
          q_tab(i,j) = pedoneData(i)%q * pedoneData(j)%q * 1.671d5
          q_tab(j,i) = q_tab(i,j)
          if(echoInput) then
             write(35,*) i, j, q_tab(i,j) 
           endif          
        enddo
      enddo
 
      repul_tab = 0d0
      alpha_Tab = 0d0
      rEq_tab = 0d0
      D_Tab = 0d0
      do i = 1,nAtomTypes
        repul_tab(i,1) = pedoneData(i)%repul
        alpha_Tab(i,1) = pedoneData(i)%alpha
        rEq_tab(i,1) = pedoneData(i)%rEq
        D_Tab(i,1) = pedoneData(i)%delta

        repul_tab(1,i) = repul_tab(i,1)
        alpha_Tab(1,i) = alpha_Tab(i,1)
        rEq_tab(1,i) = rEq_tab(i,1)
        D_Tab(1,i) = D_Tab(i,1)
         
      enddo

   
        if(echoInput) then
          do i = 1, nAtomTypes 
            do j = 1, nAtomTypes 
              write(35,*) i, j, repul_tab(i,j), alpha_Tab(i,j), rEq_tab(i,j), D_Tab(i,j), q_tab(i,j)
            enddo
          enddo
        endif 

!     R_Min Mixing Rule    
      read(55,*)
      read(55,*) labelField, mixingRule 
!       write(6,*) labelField, mixingRule      
      custom = .false.
      select case(mixingRule)
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

      if(custom) then           
        do i = 1,nAtomTypes
          read(55,*) (r_min_tab(i,j), j=1,nAtomTypes)
        enddo
        do i = 1,nAtomTypes
          do j = 1,nAtomTypes
            r_min_tab(i,j) = r_min_tab(i,j)**2
          enddo
        enddo
      else
        do i = 1,nAtomTypes
          read(55,*) r_min(i)
        enddo
        do i = 1,nAtomTypes
          do j = i,nAtomTypes
            r_min_tab(i,j) = rmin_func(r_min(i), r_min(j))**2
            r_min_tab(j,i) = r_min_tab(i,j)
          enddo
        enddo
      endif

      do i = 1, maxRosenTrial
        allocate(rosenTrial(i)%x(1:1))      
        allocate(rosenTrial(i)%y(1:1)) 
        allocate(rosenTrial(i)%z(1:1))         
      enddo

      ALLOCATE (atomArray(1:nMolTypes,1:1), STAT = AllocateStatus)
      do i = 1, nMolTypes
        atomArray(i,1) = i
        atomData(i)%Symb = pedoneData(i)%Symb
      enddo
 
      vmdAtoms = 0
      do i = 1,nMolTypes
        vmdAtoms = vmdAtoms + (NMAX(i)*nAtoms(i))
      enddo
      
      totalMass = 0d0
      do i = 1, nMolTypes
        do j = 1,nAtoms(i)
          totalMass(i) = totalMass(i) + atomData(atomArray(i,j))%mass
        enddo
      enddo

      ALLOCATE (bornRad(1:nMolTypes), STAT = AllocateStatus)
      read(55,*)
      read(55,*) labelField, implcSolvent 
      do i = 1, nAtomTypes
        read(55,*) labelField, bornRad(i)
      enddo

      call AllocateCoordinateArrays
      call CreateJointArray  
      call SetStorageFlags(q_tab)
!      call createDistArrays
      flush(35)

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
!================================================================ 
      subroutine Allocate_Common_Variables_Pedone
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      
      ALLOCATE (pedoneData(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (atomData(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (nAtoms(1:nMolTypes), STAT = AllocateStatus)

      nAtoms = 1

      ALLOCATE (r_min(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_sq(1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (r_min_tab(1:nMolTypes, 1:nMolTypes), STAT = AllocateStatus) 

      ALLOCATE (alpha_Tab(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (D_Tab(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (repul_tab(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (rEq_tab(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus)
      ALLOCATE (q_tab(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus)

      ALLOCATE (totalMass(1:nMolTypes), STAT = AllocateStatus)
     
      repul_tab = 0d0
      D_Tab = 0d0       
      alpha_Tab = 0d0   
      rEq_tab = 0d0
      q_tab = 0d0
      r_min_tab = 0d0



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

