!========================================================            
      subroutine ReadParameters(seed,ncycle,nmoves,max_dist,max_rot,
     &	                  outFreq_Traj,outFreq_Screen,outFreq_GCD,
     &	                  screenEcho)
      use SimParameters
      use Constants
      use ForceField
      use Units
      implicit none
      
      logical, intent(INOUT)  :: screenEcho
      integer, intent(INOUT) :: seed
      integer(kind=8), intent(INOUT) :: ncycle,nmoves
      integer, intent(INOUT)  :: outFreq_Traj,outFreq_Screen,outFreq_GCD
      real(kind(0.0d0)), intent(INOUT) :: max_dist,max_rot
      
      integer :: i,j
      logical :: useBias	  
      integer :: AllocateStatus	        
      character(len=30) :: labelField 
      character(len=30) :: fileName
      real(kind(0.0d0)) :: dummyCycle

      open(unit=54,file="input_Parameters.dat",status='OLD')	  

!	    Comment Space	  
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
	  
      allocate( NPART(1:nMolTypes) )	  
      allocate( NMIN(1:nMolTypes) )	  
      allocate( NMAX(1:nMolTypes) )
      allocate( gas_dens(1:nMolTypes) )
	  
      read(54,*) labelField, (NMIN(j),j=1,nMolTypes)	  
      read(54,*) labelField, (NMAX(j),j=1,nMolTypes)
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
!       Comment Space	  
      read(54,*)
      read(54,*)
      read(54,*)       	  
!      Monte Carlo Move Parameters	  
      read(54,*) labelField, r_min
       r_min_sq=r_min**2
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
      use ForceField
      use ForceFieldFunctions
      use Units
      implicit none
      integer i,j,h,AllocateStatus,nParam
      integer :: nAtomsMax, nBondsMax,nAnglesMax	  
      integer :: nTorsMax, nImpropMax
      character(len=15) :: labelField 
      character(len=15) :: mixingRule, units1,units2
      character(len=15) :: unitsEnergy,unitsDistance,unitsAngular
      real(kind(0.0d0)) :: convEng, convDist, convAng
      procedure (MixRule), pointer :: ep_func => null()
      procedure (MixRule), pointer :: sig_func => null()      
	  
      open(unit=54,file="input_forcefield.dat",status='OLD')! Open forcefield file
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

C       write(6,*) labelField, nAtomTypes
C       write(6,*) labelField, nBondTypes
C       write(6,*) labelField, nAngleTypes
C       write(6,*) labelField, nTorsionalTypes
C       write(6,*) labelField, nImproperTypes	  
	  
      call Allocate_Common_Variables

!     Blank Space for Comments
      read(54,*) 
      read(54,*)

!     Epsilon Mixing Rule
      read(54,*) labelField, mixingRule 
C       write(6,*) labelField, mixingRule 	  
      select case(mixingRule)
        case("geometric")
           ep_func => GeoMean_MixingFunc
        case("average")
           ep_func => Mean_MixingFunc
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

C       write(6,*) labelField, unitsEnergy, convEng
C       write(6,*) labelField, unitsDistance, convDist
C       write(6,*) labelField, unitsAngular, convAng
	  
!      Begin reading the intermolecular forcefield from the file
      read(54,*)
      read(54,*)      
      do i = 1,nAtomTypes
        read(54,*) atomData(i)%Symb, atomData(i)%ep, atomData(i)%sig, 
     &	           atomData(i)%q, atomData(i)%mass
C         write(6,*) atomData(i)%Symb, atomData(i)%ep, atomData(i)%sig, 
C      &	           atomData(i)%q, atomData(i)%mass	 
        atomData(i)%ep = atomData(i)%ep * convEng
        atomData(i)%sig = atomData(i)%sig * convDist		
      enddo

!      Generate the look up tables for the inter molecular interactions
      do i = 1,nAtomTypes
        do j = i,nAtomTypes
          ep_tab(i,j)  = ep_func(atomData(i)%ep, atomData(j)%ep)
          sig_tab(i,j) = sig_func(atomData(i)%sig, atomData(j)%sig)**2
          q_tab(i,j) = atomData(i)%q * atomData(j)%q * 1.671d5
		  
          ep_tab(j,i) = ep_tab(i,j)
          sig_tab(j,i) = sig_tab(i,j)
          q_tab(j,i) = q_tab(i,j)
        enddo
      enddo
 
	  
!     Begin reading the Bond type definitions from file
      read(54,*)
      read(54,*)      
      do i=1,nBondTypes
        read(54,*) labelField, bondData(i)%k_eq, bondData(i)%r_eq
!        write(6,*) labelField, bondData(i)%k_eq, bondData(i)%r_eq		
        bondData(i)%k_eq = bondData(i)%k_eq * convEng / convDist**2
        bondData(i)%r_eq = bondData(i)%r_eq * convDist
      enddo
!     Begin reading the Bending Angle type definitions from file
      read(54,*)
      read(54,*)      
      do i=1,nAngleTypes
        read(54,*) labelField, bendData(i)%k_eq, bendData(i)%ang_eq
!        write(6,*) labelField, bendData(i)%k_eq, bendData(i)%ang_eq		
        bendData(i)%k_eq = bendData(i)%k_eq * convEng / convAng**2
        bendData(i)%ang_eq = bendData(i)%ang_eq * convAng
      enddo      

!     Begin reading the Torsional Angle type definitions from file
      read(54,*)
      read(54,*)      
      do i=1,nTorsionalTypes
        read(54,*) labelField, nParam
     	backspace(54)
        allocate(torsData(i)%a(1:nParam),
     &	        STAT = AllocateStatus)
        read(54,*) labelField, torsData(i)%nPara,
     &	           (torsData(i)%a(j),j=1,nParam)
!        write(6,*) labelField, torsData(i)%nPara,
!     &	           (torsData(i)%a(j),j=1,nParam)	 
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
     &	        STAT = AllocateStatus)
        read(54,*) labelField, impropData(i)%nPara,
     &	           (impropData(i)%a(j),j=1,nParam)
!        write(6,*) labelField, impropData(i)%nPara,
!     &	           (impropData(i)%a(j),j=1,nParam)	 
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
        read(54,*) labelField, nBonds(i)
        read(54,*) labelField, nAngles(i)
        read(54,*) labelField, nTorsional(i)
        read(54,*) labelField, nImproper(i)		
      enddo

      nAtomsMax = maxval(nAtoms)
      nBondsMax = maxval(nBonds)
      nAnglesMax = maxval(nAngles)
      nTorsMax = maxval(nTorsional)
      nImpropMax= maxval(nImproper)
	  
      ALLOCATE (atomArray(1:nMolTypes,1:nAtomsMax),
     &	        STAT = AllocateStatus)
      ALLOCATE (bondArray(1:nMolTypes,1:nBondsMax),
     &	        STAT = AllocateStatus)
      ALLOCATE (bendArray(1:nMolTypes,1:nAnglesMax),
     &	        STAT = AllocateStatus)
      ALLOCATE (torsArray(1:nMolTypes,1:nTorsMax),
     &	        STAT = AllocateStatus)
      ALLOCATE (impropArray(1:nMolTypes,1:nImpropMax),
     &	        STAT = AllocateStatus)
	  

C       write(6,*) (nAtoms(i),i=1,nMolTypes)
C       write(6,*) (nBonds(i),i=1,nMolTypes)
C       write(6,*) (nAngles(i),i=1,nMolTypes)
C       write(6,*) (nTorsional(i),i=1,nMolTypes)
C       write(6,*) (nImproper(i),i=1,nMolTypes)
	  
      do i=1,nMolTypes
        !Atom Definition for each molecular species
        read(54,*)		
        read(54,*)
        read(54,*)
        do j=1,nAtoms(i)  
          read(54,*) labelField, atomArray(i,j)
        enddo
        !Bond Definition for each molecular species
        read(54,*)
        read(54,*)
        do j=1,nBonds(i)  
          read(54,*) labelField, bondArray(i,j)%bondType,
     &                  (bondArray(i,j)%bondMembr(h),h=1,2)
        enddo
        !Bending Angle Definition for each molecular species		
        read(54,*)
        read(54,*)
        do j=1,nAngles(i)  
          read(54,*) labelField, bendArray(i,j)%bendType,
     &                  (bendArray(i,j)%bendMembr(h),h=1,3)
        enddo
        !Torsion Angle Definition for each molecular species		
        read(54,*)
        read(54,*)		
        do j=1,nTorsional(i)  
          read(54,*) labelField, torsArray(i,j)%TorsType,
     &                  (torsArray(i,j)%torsMembr(h),h=1,4)
        enddo
        !Improper Angle Definition for each molecular species		
        read(54,*)
        read(54,*)		
        do j=1,nImproper(i)  
          read(54,*) labelField, impropArray(i,j)%TorsType,
     &                  (impropArray(i,j)%torsMembr(h),h=1,4)
        enddo		
      enddo

      
      close(54)
      end subroutine

!================================================================ 
      subroutine Allocate_Common_Variables
      use SimParameters
      use ForceField	  
      implicit none
      integer :: AllocateStatus
      !
      ALLOCATE (atomData(1:nAtomTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (bondData(1:nBondTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (bendData(1:nAngleTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (torsData(1:nTorsionalTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (impropData(1:nImproperTypes),
     &	        STAT = AllocateStatus)

	 
      ALLOCATE (ep_tab(1:nAtomTypes,1:nAtomTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (sig_tab(1:nAtomTypes,1:nAtomTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes,1:nAtomTypes),
     &	        STAT = AllocateStatus)
	 
      ALLOCATE (nAtoms(1:nMolTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (nBonds(1:nMolTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (nAngles(1:nMolTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (nTorsional(1:nMolTypes),
     &	        STAT = AllocateStatus)
      ALLOCATE (nImproper(1:nMolTypes),
     &	        STAT = AllocateStatus)


      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	  
      end    
