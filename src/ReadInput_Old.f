!========================================================            
      subroutine ReadParameters(seed,ncycle,nmoves,outFreq_Traj,
     &                        outFreq_Screen,outFreq_GCD,screenEcho)
      use CBMC_Variables
      use Constants
      use Coords
      use EnergyTables
      use ForceField
      use SimParameters
      use Units
      use ParallelVar
      use WHAM_Module
      use UmbrellaFunctions
      implicit none
      
      logical, intent(OUT)  :: screenEcho
      integer, intent(OUT) :: seed
      integer(kind=8), intent(OUT) :: ncycle,nmoves
      integer, intent(OUT)  :: outFreq_Traj,outFreq_Screen,outFreq_GCD
      
      integer :: i,j
      logical :: useBias   
      integer :: AllocateStatus            
      character(len=30) :: labelField 
      character(len=30) :: fileName
      real(kind(0.0d0)) :: dummyCycle

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
       nMoves=nint(dummyCycle)
      read(54,*) labelField, nMolTypes
   


      allocate( NPART(1:nMolTypes),STAT = AllocateStatus )    
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
!        write(6,*) labelField, useBias, fileName
        call AllocateUmbrellaBias(fileName)
      else 
        call BlankUmbrellaBias
      endif

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
!      read(54,*) labelField, r_min
!       r_min_sq=r_min**2
      read(54,*) labelField, (gas_dens(j),j=1,nMolTypes)    
       
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

      allocate(refSizeNumbers(1:nMolTypes),STAT = AllocateStatus)

      read(54,*) labelField, useWHAM
      read(54,*) labelField, (refSizeNumbers(j), j=1,nMolTypes)
      refBin = getBiasIndex(refSizeNumbers, NMAX)
      read(54,*) labelField, intervalWHAM
      read(54,*) labelField, tolLimit
      read(54,*) labelField, maxSelfConsist
      read(54,*) labelField, whamEstInterval
      read(54,*) labelField, equilInterval
      if(useWham) then
        nWhamItter = ceiling(dble(ncycle)/dble(intervalWHAM))
        call WHAM_Initialize
      endif

      
      read(54,*)
      read(54,*)
      read(54,*)
!     Cluster Critiera Parameters
      allocate( Eng_Critr(1:nMolTypes,1:nMolTypes),
     &         stat = AllocateStatus) 
      read(54,*)
      do i=1,nMolTypes      
        read(54,*) (Eng_Critr(i,j),j=1,nMolTypes)
      enddo
      read(54,*)      
      
      allocate(biasAlpha(1:nMolTypes,1:nMolTypes),
     &         stat = AllocateStatus)      
      do i = 1,nMolTypes
         read(54,*) (biasAlpha(i,j),j=1,nMolTypes)     
         do j = 1, nMolTypes
            biasAlpha(i,j) = biasAlpha(i,j)/temperature
         enddo         
      enddo      


      
      close(54)       
      end            
!========================================================  
!     This subroutine reads the user specified forcefield.  The first half of this
!     code deals with defining the types of atoms,bonds, etc. to be used in the simulation.
!     This is done to such that degenerate angles only need to be defined once. 
!     The second half of this routine deals with defining the molecule configuration using
!     the atoms, angles, etc. defined in the first section of the routine.
      subroutine ReadForcefield
      use SimParameters
      use Constants
      use Coords
      use ForceField
      use ForceFieldFunctions
      use Units
      use CBMC_Variables
      implicit none
      logical :: custom
      integer i,j,h,AllocateStatus,nParam
      integer :: nAtomsMax, nBondsMax,nAnglesMax    
      integer :: nTorsMax, nImpropMax, nNonBondMax
      character(len=15) :: labelField 
      character(len=10) :: mixingRule
      character(len=10) :: unitsEnergy,unitsDistance,unitsAngular
      real(kind(0.0d0)) :: convEng, convDist, convAng
      procedure (MixRule), pointer :: ep_func => null()
      procedure (MixRule), pointer :: sig_func => null()      
      procedure (MixRule), pointer :: rmin_func => null()   

!     Open forcefield file      
      open(unit=54,file="input_forcefield.dat",status='OLD')
!     Blank Space for Comments
      read(54,*) 
      read(54,*)
      read(54,*)   
       
!      Specifies the number of atoms, bonds, etc. types to be used in the simulation.
      read(54,*) labelField, nAtomTypes
      read(54,*) labelField, nBondTypes
      read(54,*) labelField, nAngleTypes
      read(54,*) labelField, nTorsionalTypes
      read(54,*) labelField, nImproperTypes


       
      call Allocate_Common_Variables

!     Blank Space for Comments
      read(54,*) 
      read(54,*)

!     Epsilon Mixing Rule
      read(54,*) labelField, mixingRule 
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
      read(54,*) labelField, mixingRule 
C       write(6,*) labelField, mixingRule      
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



       
      read(54,*) 
      read(54,*)      
      read(54,*) labelField, unitsEnergy
      read(54,*) labelField, unitsDistance
      read(54,*) labelField, unitsAngular
      convEng = FindEngUnit(unitsEnergy)
      convDist = FindLengthUnit(unitsDistance)
      convAng = FindAngularUnit(unitsAngular)
       
!      Begin reading the intermolecular forcefield from the file
      read(54,*)
      read(54,*)      
      do i = 1,nAtomTypes
        read(54,*) labelField,atomData(i)%Symb, atomData(i)%ep, 
     &        atomData(i)%sig, atomData(i)%q, atomData(i)%mass
        if(echoInput) then
          write(35,*) labelField,atomData(i)%Symb, atomData(i)%ep, 
     &        atomData(i)%sig, atomData(i)%q, atomData(i)%mass
        endif
        atomData(i)%ep = atomData(i)%ep * convEng
        atomData(i)%sig = atomData(i)%sig * convDist   
        r_min_sq(i) = r_min(i)*r_min(i)

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
             write(35,*) i,j, ep_tab(i,j)/4d0, sqrt(sig_tab(i,j)),
     &                   q_tab(i,j) 
           endif          
!          r_min_tab(j,i) = r_min_tab(i,j)
        enddo
      enddo
      if(echoInput) then
        write(35,*) "---------------------------------------------"
        flush(35)
      endif
      
!     R_Min Mixing Rule    
      read(54,*)
      read(54,*) labelField, mixingRule 
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
          read(54,*) (r_min_tab(i,j), j=1,nAtomTypes)
        enddo
        do i = 1,nAtomTypes
          do j = 1,nAtomTypes
            r_min_tab(i,j) = r_min_tab(i,j)**2
          enddo
        enddo
      else
        do i = 1,nAtomTypes
          read(54,*) r_min(i)
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
      read(54,*)
      read(54,*)      
      do i=1,nBondTypes
        read(54,*) labelField, bondData(i)%k_eq, bondData(i)%r_eq
        bondData(i)%k_eq = bondData(i)%k_eq * convEng
        bondData(i)%r_eq = bondData(i)%r_eq * convDist
        if(echoInput) then
          write(35,*) labelField, bondData(i)%k_eq, bondData(i)%r_eq    
        endif        
      enddo
!     Begin reading the Bending Angle type definitions from file
      read(54,*)
      read(54,*)      
      do i=1,nAngleTypes
        read(54,*) labelField, bendData(i)%k_eq, bendData(i)%ang_eq
        if(echoInput) then
          write(35,*) labelField, bendData(i)%k_eq, bendData(i)%ang_eq 
        endif          
        bendData(i)%k_eq = bendData(i)%k_eq * convEng
        bendData(i)%ang_eq = bendData(i)%ang_eq * convAng
      enddo      

!     Begin reading the Torsional Angle type definitions from file
      read(54,*)
      read(54,*)      
      do i=1,nTorsionalTypes
        read(54,*) labelField, nParam
          backspace(54)
        allocate(torsData(i)%a(1:nParam),
     &            STAT = AllocateStatus)
        read(54,*) labelField, torsData(i)%nPara,
     &               (torsData(i)%a(j),j=1,nParam)
        if(echoInput) then
          write(35,*) labelField, torsData(i)%nPara,
     &               (torsData(i)%a(j),j=1,nParam)
        endif         
        do j=1,nParam
          torsData(i)%a(j) = torsData(i)%a(j) * convEng
        enddo   
      enddo          

!     Begin reading the Improper Angle type definitions from file
      read(54,*)
      read(54,*)      
      do i=1,nImproperTypes
        read(54,*) labelField, nParam
          backspace(54)
        allocate(impropData(i)%a(1:nParam),
     &            STAT = AllocateStatus)
        read(54,*) labelField, impropData(i)%nPara,
     &               (impropData(i)%a(j),j=1,nParam)
        do j=1,nParam
          impropData(i)%a(j) = impropData(i)%a(j) * convEng
        enddo
      enddo                

!     Blank space for Comments
      read(54,*) 
       
!     This section constructs the internal molecular configuration 
!     using the previously defined atom, bond, bend, etc. types.
      do i=1,nMolTypes
        read(54,*)
        read(54,*)       
        read(54,*) labelField, nAtoms(i)
        read(54,*) labelField, nIntraNonBond(i)        
        read(54,*) labelField, nBonds(i)
        read(54,*) labelField, nAngles(i)
        read(54,*) labelField, nTorsional(i)
        read(54,*) labelField, nImproper(i)       
      enddo

      nAtomsMax = maxval(nAtoms)
      nNonBondMax = maxval(nIntraNonBond)      
      nBondsMax = maxval(nBonds)      
      nAnglesMax = maxval(nAngles)
      nTorsMax = maxval(nTorsional)
      nImpropMax= maxval(nImproper)
  
      ALLOCATE (atomArray(1:nMolTypes,1:nAtomsMax),
     &        STAT = AllocateStatus)
      ALLOCATE (nonBondArray(1:nMolTypes,1:nNonBondMax),
     &        STAT = AllocateStatus)     
      ALLOCATE (bondArray(1:nMolTypes,1:nBondsMax),
     &        STAT = AllocateStatus)
      ALLOCATE (bendArray(1:nMolTypes,1:nAnglesMax),
     &        STAT = AllocateStatus)
      ALLOCATE (torsArray(1:nMolTypes,1:nTorsMax),
     &        STAT = AllocateStatus)
      ALLOCATE (impropArray(1:nMolTypes,1:nImpropMax),
     &        STAT = AllocateStatus)
  
      ALLOCATE (atomUseByBond(1:nMolTypes,1:nAtomsMax),
     &        STAT = AllocateStatus)
      ALLOCATE (atomUseByBend(1:nMolTypes,1:nAtomsMax),
     &        STAT = AllocateStatus)
      ALLOCATE (atomUseByTors(1:nMolTypes,1:nAtomsMax),
     &        STAT = AllocateStatus)
      ALLOCATE (atomUseByImprop(1:nMolTypes,1:nAtomsMax),
     &        STAT = AllocateStatus)

      do i = 1, maxRosenTrial
        allocate(rosenTrial(i)%x(1:nAtomsMax))      
        allocate(rosenTrial(i)%y(1:nAtomsMax)) 
        allocate(rosenTrial(i)%z(1:nAtomsMax))         
      enddo
     
      do i=1,nMolTypes
        !Atom Definition for each molecular species
        read(54,*)       
        read(54,*)
        read(54,*)
        do j=1,nAtoms(i)  
          read(54,*) labelField, atomArray(i,j)
          if(echoInput) then
             write(35,*) labelField, atomArray(i,j)
           endif           
!          write(6,*) labelField, atomArray(i,j)
        enddo
!       !Intra Nonbonded (1-5) interaction definition for each molecular species
        read(54,*)
        read(54,*)
        do j=1,nIntraNonBond(i) 
          read(54,*) labelField, (nonBondArray(i,j)%nonMembr(h),h=1,2)
          if(echoInput) then
           write(35,*) labelField, (nonBondArray(i,j)%nonMembr(h),h=1,2)
          endif                
!          write(6,*) labelField, (nonBondArray(i,j)%nonMembr(h),h=1,2)
        enddo        
        !Bond Definition for each molecular species
        read(54,*)
        read(54,*)
        do j=1,nBonds(i)  
          read(54,*) labelField, bondArray(i,j)%bondType,
     &                  (bondArray(i,j)%bondMembr(h),h=1,2)
          if(echoInput) then
            write(35,*) labelField, bondArray(i,j)%bondType,
     &                  (bondArray(i,j)%bondMembr(h),h=1,2)
          endif                 
     
        enddo
        !Bending Angle Definition for each molecular species          
        read(54,*)
        read(54,*)
        do j=1,nAngles(i)  
          read(54,*) labelField, bendArray(i,j)%bendType,
     &                  (bendArray(i,j)%bendMembr(h),h=1,3)
          if(echoInput) then
            write(35,*) labelField, bendArray(i,j)%bendType,
     &                  (bendArray(i,j)%bendMembr(h),h=1,3)
          endif          
        enddo
        !Torsion Angle Definition for each molecular species          
        read(54,*)
        read(54,*)
        do j=1,nTorsional(i)  
          read(54,*) labelField, torsArray(i,j)%TorsType,
     &                  (torsArray(i,j)%torsMembr(h),h=1,4)
          if(echoInput) then
            write(35,*) labelField, torsArray(i,j)%TorsType,
     &                  (torsArray(i,j)%torsMembr(h),h=1,4)
          endif       
        enddo
        !Improper Angle Definition for each molecular species         
        read(54,*)
        read(54,*)
        do j=1,nImproper(i)  

          read(54,*) labelField, impropArray(i,j)%ImpropType,
     &                  (impropArray(i,j)%impropMembr(h),h=1,4)
          if(echoInput) then
            write(35,*) labelField, impropArray(i,j)%ImpropType,
     &                  (impropArray(i,j)%impropMembr(h),h=1,4)
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
      
      close(54)

 
      
      end subroutine

!================================================================ 
      subroutine Allocate_Common_Variables
      use SimParameters
      use ForceField
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      
      ALLOCATE (atomData(1:nAtomTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (bondData(1:nBondTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (bendData(1:nAngleTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (torsData(1:nTorsionalTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (impropData(1:nImproperTypes),
     &        STAT = AllocateStatus)

      ALLOCATE (r_min(1:nAtomTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes),
     &        STAT = AllocateStatus) 

      ALLOCATE (ep_tab(1:nAtomTypes,1:nAtomTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes),
     &        STAT = AllocateStatus)
     
      ep_tab = 0d0
      sig_tab = 0d0       
      q_tab = 0d0
      r_min_tab = 0d0

      ALLOCATE (totalMass(1:nMolTypes),
     &        STAT = AllocateStatus)
      
      ALLOCATE (nAtoms(1:nMolTypes),
     &          STAT = AllocateStatus)
      ALLOCATE (nIntraNonBond(1:nMolTypes),
     &        STAT = AllocateStatus)      
      ALLOCATE (nBonds(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (nAngles(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (nTorsional(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (nImproper(1:nMolTypes),
     &        STAT = AllocateStatus)


      ALLOCATE (acptTrans(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (acptRot(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (acptSwapOut(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (acptInSize(1:maxMol),
     &        STAT = AllocateStatus)
      ALLOCATE (atmpInSize(1:maxMol),
     &        STAT = AllocateStatus)


      ALLOCATE (max_dist(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes),
     &        STAT = AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes),
     &        STAT = AllocateStatus) 
      
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  
      end    
