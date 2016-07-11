
!======================================================
      module CoordinateTypes
      use VarPrecision

      type FloatPointer
        real(dp),pointer :: pnt
      end type

      type SimpleMolCoords
        real(dp),allocatable :: x(:), y(:), z(:)
      end type

      type SimpleAtomCoords
        real(dp) :: x, y, z
      end type

      type Molecule
        integer(kind=2) :: indx
!        integer :: indx        
        real(dp),allocatable :: x(:), y(:), z(:)
      end type

      type MolArrayType
        type(Molecule),allocatable :: mol(:)
      end type

      type MolPointers
        integer(kind=2) :: molType,molIndx
!        integer :: molType,molIndx        
        type(FloatPointer),allocatable :: x(:), y(:), z(:)
      end type

      type Displacement
        integer(kind=2) :: molType, atmIndx, molIndx
!        integer :: molType,molIndx,atmIndx        
        real(dp) :: x_new, y_new, z_new
        real(dp),pointer :: x_old, y_old, z_old
      end type

      type TrialCoordinates
        integer(kind=2) :: molType, molIndx
        real(dp), allocatable :: x(:), y(:), z(:)
      end type

      end module

!==============================================================
      module AcceptRates
      use VarPrecision
      real(dp), allocatable :: acptTrans(:), acptRot(:)
      real(dp), allocatable :: atmpTrans(:), atmpRot(:)      
      
      real(dp), allocatable :: acptSwapIn(:), atmpSwapIn(:)
      real(dp), allocatable :: acptSwapOut(:), atmpSwapOut(:)

      real(dp), allocatable :: acptInSize(:), atmpInSize(:)

      real(dp) :: angGen_accpt, angGen_atmp
      real(dp) :: dihedGen_accpt, dihedGen_atmp
      real(dp) :: distGen_accpt, distGen_atmp
      real(dp) :: clusterCritRej
      
      end module  
!==============================================================
      module AVBMC_RejectionVar
      use CoordinateTypes
        
      real(dp) :: totalRej
      real(dp) :: ovrlapRej, dbalRej, critriaRej, boundaryRej

      real(dp) :: totalRej_out
      real(dp) :: boundaryRej_out
      
      end module
!======================================================
      module ParallelVar
        integer myid, p_size, ierror, tag, nout
      end module         
!======================================================
      module IonBias
      use VarPrecision
        logical, parameter :: biasIon = .false.
        integer, parameter :: rBins = 300
        real(dp) :: ionDistance, dR
        real(dp),allocatable :: NR_Hist(:,:), NR_Bias(:,:)  
        
      end module         

!==============================================================
      module Coords
      use CoordinateTypes
        
      type(MolArrayType), allocatable, target :: MolArray(:)
      type(MolPointers), allocatable :: JointArray(:) 
      
      type(TrialCoordinates) :: newMol, newMol2
      type(SimpleMolCoords),allocatable :: rosenTrial(:)
      
      logical, allocatable :: NeighborList(:,:)
      integer, allocatable :: NumNei(:)   
!      real(dp), allocatable :: PairTable(:,:)
      integer, allocatable :: typeList(:)
      integer, allocatable :: subIndxList(:)
      
      type(SimpleMolCoords), allocatable :: gasConfig(:)
      
      end module
    
!==============================================================
      module CBMC_Variables
      use VarPrecision

      type BondNumber
        integer(kind=2),allocatable :: atom(:)
      end type

      type Pathing
        integer :: nTerminal, nLinker,nHub      
        integer(kind=2),allocatable :: termAtoms(:)
        integer(kind=2),allocatable :: hubAtoms(:) 
        integer :: nPaths
        integer(kind=2),allocatable :: path(:,:)        
        integer(kind=2),allocatable :: pathMax(:)
      end type      
 
 
      integer :: maxRosenTrial
      integer, allocatable :: nRosenTrials(:)      
      integer, allocatable :: regrowType(:)
      integer, allocatable :: regrowOrder(:,:)
      real(dp), allocatable :: probTypeCBMC(:)
      type(BondNumber),allocatable :: topolArray(:)
      type(Pathing),allocatable :: pathArray(:)
      
      integer, allocatable :: usedByPath(:,:)
      integer, allocatable :: atomPathIndex(:,:)
      
      end module


!==============================================================
!      module CBMC_Pointers
!      
!      interface 
!        subroutine growFunc(nType)
!          integer, intent(in) :: nType
!        end subroutine
!      end interface
!      
!      type GrowPointer
!        procedure(growFunc), pointer :: pnt
!      end type
!      
!      type(GrowPointer), allocatable :: growthFunc(:)
!      
!      end module
!==============================================================
      module EnergyTables
      use VarPrecision

      real(dp), allocatable :: ETable(:)
      real(dp), allocatable :: NeiETable(:)
      real(dp), allocatable :: biasAlpha(:,:)
      
      real(dp) :: E_Inter_T, E_Inter_Diff
      real(dp) :: E_NBond_T, E_NBond_Diff      
      real(dp) :: E_Stretch_T, E_Strch_Diff
      real(dp) :: E_Bend_T, E_Bend_Diff
      real(dp) :: E_Torsion_T, E_Tors_Diff
      
      end module      

 !================================================================ 
      module IndexingFunctions
      contains
 
!     ----------------------------------------------------
      subroutine Get_MolIndex(nMove,NPart,molType,molIndx)
       implicit none
       integer,intent(in) :: nMove, NPart(:)
       integer,intent(inout) :: molIndx, molType       
       integer :: i, sizeN,curLimit
       
       sizeN = size(NPART)
!       curLimit = NPART(1)
       curLimit = 0
       molIndx = 1

       molType = 1
       do i=1,sizeN
         curLimit = curLimit + NPART(i)  
         if(nMove .le. curLimit) then       
           molIndx = nMove - curLimit + NPART(i)  
           molType = i
           return
         endif
       enddo

      write(35,*) "NotFound", nMove, molType,molIndx
      end subroutine
!     ----------------------------------------------------
      pure subroutine Get_SubIndex(nIndx,nType,NMAX)
       implicit none
       integer,intent(in) :: nIndx, NMAX(:)
       integer,intent(out) :: nType   
       integer :: i, sizeN,curLimit
       
       sizeN = size(NMAX)
       curLimit = 0
       do i = 1, nType - 1
         curLimit = curLimit + NMAX(i)
       enddo
       nType = nIndx - curLimit
       
      end subroutine
      
!     ----------------------------------------------------      
      pure integer function Get_MolType(Indx,NMAX)
       implicit none
       integer,intent(in) :: Indx
       integer,intent(in) :: NMAX(:)
       integer :: i, sizeN,curLimit
       
       sizeN = size(NMAX)
       curLimit = 0
       Get_MolType = 0
       do i=1,sizeN
         curLimit = curLimit + NMAX(i)                
         if(Indx .le. curLimit) then       
           Get_MolType = i
           return
         endif
       enddo
       
      end function
!     ----------------------------------------------------
      
      end module

!==============================================================
      module Histogram
      use VarPrecision
      integer, parameter :: nBin_Hist = 1000
      real(dp) :: d_ang      
      real(dp) :: HistAngle(0:nBin_Hist)

      real(dp) :: d_r 
      real(dp) :: HistDist(0:nBin_Hist)
      
      end module      
!==============================================================
      module SimParameters
      use VarPrecision
      logical, parameter :: echoInput = .false.      
      logical, parameter :: distCriteria = .false.


      integer(kind=8) :: ncycle,ncycle2       
      logical :: ridgidMolecules, multipleInput
      integer(kind=2) :: nMolTypes
      integer,allocatable :: NPart(:),NMin(:),NMax(:)
      logical,allocatable :: isActive(:)
      integer :: NTotal
      integer :: maxMol, maxAtoms, vmdAtoms
      real(dp) :: avbmc_vol
      
      integer :: umbrellaLimit
      real(dp),allocatable :: NHist(:), NBias(:), E_Avg(:)
      real(dp),allocatable :: gas_dens(:)
      real(dp),allocatable :: max_dist(:), max_dist_single(:)
      real(dp),allocatable :: max_rot(:)
        
      real(dp) :: temperature, beta
!      real(dp) :: r_min, r_min_sq
      real(dp) :: global_r_min

      real(dp) :: softCutoff

!     Cluster Criteria block      
      real(dp) :: Dist_Critr
      real(dp) :: Dist_Critr_sq
!      real(dp), allocatable :: Dist_Critr
!      real(dp), allocatable :: Dist_Critr_sq      
      real(dp) :: NeighRej
      real(dp),allocatable :: Eng_Critr(:,:)
      real(dp) :: ECritMax
      
      character(len=10) :: outputEngUnits,outputLenUnits
      real(dp) :: outputEConv,outputLenConv
        
      end module
!================================================================ 
      module UmbrellaFunctions

      contains

      pure integer function getBiasIndex(NPart,NMAX)
       implicit none
       integer,intent(in) :: NPart(:), NMAX(:)
       integer :: i, sizeN
       integer :: curIndx,maxIndx
        
       sizeN = size(NPART)
!       curIndx = 0
       curIndx = NPart(sizeN)
       maxIndx = NMAX(sizeN) + 1 
       do i=1,sizeN-1
           curIndx = curIndx + maxIndx*NPart(sizeN-i)
           maxIndx = maxIndx * (NMAX(sizeN-i) + 1)
       enddo
       getBiasIndex = curIndx + 1
       
      end function
!     ---------------------------------------------------------------
      pure integer function getNewBiasIndex(NPart,NMAX, increment)
       implicit none
       integer,intent(in) :: NPart(:), NMAX(:), increment(:)
       integer :: i, sizeN
       integer :: curIndx,maxIndx
        
       sizeN = size(NPART)
!       curIndx = 0
       curIndx = NPart(sizeN) + increment(sizeN)
       maxIndx = NMAX(sizeN) + 1 
       do i=1,sizeN-1
           curIndx = curIndx + maxIndx*(NPart(sizeN-i) + increment(sizeN-i))
           maxIndx = maxIndx * (NMAX(sizeN-i) + 1)
       enddo
       getNewBiasIndex = curIndx + 1
       
      end function      
    
!     ---------------------------------------------------------------      
      end module   
!======================================================
      module WHAM_Module
      use VarPrecision

      logical, parameter :: WHAM_ExtensiveOutput = .false.

      integer, allocatable :: refSizeNumbers(:)
      integer :: refBin
      logical :: useWHAM
      integer :: intervalWHAM
      integer :: maxSelfConsist
      real(dp) :: tolLimit
      integer :: equilInterval
      integer :: whamEstInterval
!      real(dp) :: unsampledIncrmnt = 1.5d0

      integer :: nWhamItter, nCurWhamItter

      real(dp),allocatable :: WHAM_Numerator(:)
      real(dp),allocatable :: WHAM_Denominator(:,:)
      real(dp),allocatable :: HistStorage(:)
      real(dp),allocatable :: FreeEnergyEst(:)
      real(dp),allocatable :: BiasStorage(:,:)
      real(dp),allocatable :: NewBias(:)
      real(dp),allocatable :: ProbArray(:)      
      real(dp),allocatable :: TempHist(:)
        
      end module  
!================================================================   
