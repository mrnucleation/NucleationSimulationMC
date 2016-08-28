!==========================================================================
    module UmbrellaSamplingNew
    use VarPrecision
    implicit none
    private

    type BiasVariablePointer
      integer :: varType
      integer, pointer :: intVar
      real(dp), pointer :: realVar
    end type

    interface
      subroutine UDispFunc(disp)
        use CoordinateTypes
        type(Displacement), intent(in) :: disp(:)
      end subroutine
    end interface

    type DispUmbrellaArray
      procedure(UDispFunc), pointer, nopass :: func
    end type

    type SwapUmbrellaArray
      procedure(), pointer, nopass :: func
    end type

    logical :: useUmbrella
    integer :: nBiasVariables, umbrellaLimit
    integer :: curUIndx
    integer, allocatable :: curVarIndx
    integer, allocatable :: binIndx(:)
    integer, allocatable :: varMax(:), varMin(:)
    integer, allocatable :: indexCoeff(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UBinSize(:)
    type(BiasVariablePointer), allocatable :: biasvar(:)
    type(BiasVariablePointer), allocatable :: biasvarnew(:)

!    integer :: nDispFunc, nSwapInFunc, nSwapOutFunc
!    type(DispUmbrellaArray), allocatable :: DispUmbrella(:)
!    type(SwapUmbrellaArray), allocatable :: SwapInUmbrella(:)
!    type(SwapUmbrellaArray), allocatable :: SwapOutUmbrella(:)

    public :: AllocateUmbrellaVariables, ReadInput_Umbrella, AllocateUmbrellaArray

!==========================================================================================
    contains
!==========================================================================================
    subroutine AllocateUmbrellaVariables
    implicit none
    integer :: AllocateStatus
        
    allocate( VarMin(1:nBiasVariables), STAT = AllocateStatus )
    allocate( VarMax(1:nBiasVariables), STAT = AllocateStatus )
    allocate( UBinSize(1:nBiasVariables), STAT = AllocateStatus )
    allocate( biasvar(1:nBiasVariables), STAT = AllocateStatus )
    allocate( biasvarnew(1:nBiasVariables), STAT = AllocateStatus )
    allocate( binIndx(1:nBiasVariables), STAT = AllocateStatus )
    allocate( indexCoeff(1:nBiasVariables), STAT = AllocateStatus )

!    allocate( DispUmbrella(1:nBiasVariables), STAT = AllocateStatus )
!    allocate( SwapInUmbrella(1:nBiasVariables), STAT = AllocateStatus )
!    allocate( SwapOutUmbrella(1:nBiasVariables), STAT = AllocateStatus )

    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    end subroutine
!==========================================================================================
    subroutine ReadInput_Umbrella(fileUnit)
      use MiscelaniousVars
      use SimpleDistPair, only: nDistPair, pairArrayIndx
      use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes
      implicit none
      integer, intent(in) :: fileUnit
      integer :: iUmbrella, AllocateStatus
      integer :: indxVar
      real(dp) :: binSize
      character(len=30) :: labelField 
      character(len=30) :: umbrellaName

      read(fileUnit, *) labelField, nBiasVariables
      if(nBiasVariables .lt. 0) then
        write(*,*) "ERROR! The user has specified an invalid number of Umbrella Sampling Variables"
        write(*,*) labelField, nBiasVariables
        stop
      endif
      if(nBiasVariables .eq. 0) then
        useUmbrella = .false.
      else
        useUmbrella = .true.
      endif

!      nDispFunc = 0
!      nSwapInFunc = 0
!      nSwapOutFunc = 0
      do iUmbrella = 1, nBiasVariables
        read(fileUnit, *) umbrellaName
        select case( trim(adjustl(umbrellaName)) )
        case("clustersize")
          backspace(fileUnit)
          read(fileUnit, *) labelField, indxVar
          if(indxVar .le. 0) then
            write(*,*) "Error! An invalid molecule type has been chosen!"
            write(*,*) "Defined Mol Types:", nMolTypes
            write(*,*) "Chosen Mol Type:", indxVar
            stop
          endif
          if(indxVar .gt. nMolTypes) then
            write(*,*) "Error! An invalid molecule type has been chosen!"
            write(*,*) "Defined Mol Types:", nMolTypes
            write(*,*) "Chosen Mol Type:", indxVar
            stop
          endif
          biasvar(iUmbrella) % varType = 1
          biasvar(iUmbrella) % intVar => NPart(indxVar)
          biasvarnew(iUmbrella) % varType = 1
          biasvarnew(iUmbrella) % intVar => NPart_New(indxVar)
          varMax(iUmbrella) = NMAX(indxVar)
          varMin(iUmbrella) = NMin(indxVar)
          UBinSize(iUmbrella) = 1E0
        case("pairdist")
          indxVar = 0
          backspace(fileUnit)
          read(fileUnit, *) labelField, indxVar, varMin(iUmbrella), varMax(iUmbrella), binSize
          if(nDistPair .le. 0) then
            write(*,*) "Error! An invalid distance variable has been chosen!"
            write(*,*) "Defined Distance Pairs:", nDistPair
            write(*,*) "Chosen Distance Pair:", indxVar
            stop
          endif
          if(indxVar .gt. nDistPair) then
            write(*,*) "Error! An invalid distance variable has been chosen!"
            write(*,*) "Defined Distance Pairs:", nDistPair
            write(*,*) "Chosen Distance Pair:", indxVar
            stop
          endif
          biasvar(iUmbrella) % varType = 2
          biasvar(iUmbrella) % realVar => miscCoord(pairArrayIndx(indxVar))
          biasvarnew(iUmbrella) % varType = 2
          biasvarnew(iUmbrella) % realVar => miscCoord_New(pairArrayIndx(indxVar))
          UBinSize(iUmbrella) = binSize
        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) umbrellaName
          stop
        end select
      enddo


    end subroutine
!==========================================================================================
    subroutine AllocateUmbrellaArray
    use UmbrellaFunctions
    use SimParameters
    use WHAM_Module
    implicit none
    integer :: i, j
    integer :: AllocateStatus
        
     umbrellaLimit = 1
     do i = 1, nBiasVariables 
       umbrellaLimit = umbrellaLimit * (VarMax(i) - VarMin(i) + 1)
     enddo
        
     allocate(UBias(1:umbrellaLimit), STAT = AllocateStatus)
     allocate(UHist(1:umbrellaLimit), STAT = AllocateStatus)

     indexCoeff(1) = 1
     do i = 2, nBiasVariables 
       indexCoeff(i) = 1
       do j = 1, i-1
         indexCoeff(i) = indexCoeff(i) + indexCoeff(j) * (VarMax(j) - VarMin(j))
       enddo
     enddo      

     if(useWHAM) then
       call WHAM_Initialize
     endif
      
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
     end subroutine
!==========================================================================
     function getBiasIndex() result(biasIndx)
     integer :: biasIndx
     integer :: iBias
      

     do iBias = 1, nBiasVariables
       if(biasvar(iBias) % varType .eq. 1) then
         binIndx(iBias) = nint( biasvar(iBias) % intVar * UBinSize(iBias) )
       elseif(biasvar(iBias) % varType .eq. 2) then
         binIndx(iBias) = nint( biasvar(iBias) % realVar * UBinSize(iBias) )
       endif
     enddo


     biasIndx = 1
     do iBias = 1, nBiasVariables
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - VarMin(iBias) )
     enddo
     

     end function
!==========================================================================
     function getNewBiasIndex() result(biasIndx)
!     real(dp), intent(in) :: newArray(:)
     integer :: biasIndx
     integer :: iBias
      

     do iBias = 1, nBiasVariables
       if(biasvarnew(iBias) % varType .eq. 1) then
         binIndx(iBias) = nint( biasvarnew(iBias)%intVar * UBinSize(iBias) )
       elseif(biasvar(iBias) % varType .eq. 2) then
         binIndx(iBias) = nint( biasvarnew(iBias)%realVar * UBinSize(iBias) )
       endif
     enddo


     biasIndx = 1
     do iBias = 1, nBiasVariables
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - VarMin(iBias) )
     enddo
     

     end function
!==========================================================================
    end module
!==========================================================================


