!==========================================================================
    module UmbrellaSamplingNew
    use VarPrecision
    implicit none
    private

    type BiasVariablePointer
      real(dp), pointer :: var
    end type

    type DispUmbrellaArray
      procedure(), pointer, nopass :: func
    end type

    logical :: useUmbrella
    integer :: nBiasVariables, umbrellaLimit
    integer, allocatable :: binIndx(:)
    integer, allocatable :: varMax(:), varMin(:)
    integer, allocatable :: indexCoeff(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: binSize(:)
    type(BiasVariablePointer), allocatable :: biasvar(:)
    type(BiasVariablePointer), allocatable :: biasvarnew(:)

    type(DispUmbrellaArray), allocatable :: DispUmbrella(:)
    type(DispUmbrellaArray), allocatable :: SwapInUmbrella(:)
    type(DispUmbrellaArray), allocatable :: SwapOutUmbrella(:)

    public :: AllocateUmbrellaBias

!==========================================================================
    contains
!==========================================================================================
    subroutine AllocateUmbrellaVariables
    implicit none
    integer :: AllocateStatus
        
    allocate( VarMin(1:nBiasVariables), STAT = AllocateStatus )
    allocate( VarMax(1:nBiasVariables), STAT = AllocateStatus )
    allocate( binSize(1:nBiasVariables), STAT = AllocateStatus )
    allocate( biasvar(1:nBiasVariables), STAT = AllocateStatus )
    allocate( biasvarnew(1:nBiasVariables), STAT = AllocateStatus )
    allocate( binIndx(1:nBiasVariables), STAT = AllocateStatus )
    allocate( indexCoeff(1:nBiasVariables), STAT = AllocateStatus )

    allocate( DispUmbrella(1:nBiasVariables), STAT = AllocateStatus )
    allocate( SwapInUmbrella(1:nBiasVariables), STAT = AllocateStatus )
    allocate( SwapOutUmbrella(1:nBiasVariables), STAT = AllocateStatus )

    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
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
!==========================================================================================
     subroutine ReadUmbrellaInput(fileUnit)
     use UmbrellaFunctions
     use SimParameters
     use WHAM_Module
     implicit none
     integer, intent(in) :: fileUnit
     integer :: i,j
     integer :: AllocateStatus
     integer,allocatable :: arrayIndx(:)   
     integer :: curIndx
     real(dp) :: curValue, defaultVal
        
     umbrellaLimit = 1
     do i = 1, nMolTypes
       umbrellaLimit = umbrellaLimit*(VarMax(i) - VarMin(i) + 1)
     enddo
        
     allocate(UBias(1:umbrellaLimit), STAT = AllocateStatus)
     allocate(UHist(1:umbrellaLimit), STAT = AllocateStatus)
       
     NHist = 0E0
     NBias = 0E0
     open( unit = 25, file = trim( adjustl(fileName) ), status='OLD')
!      write(35,*) ""
!      do i = 1,umbrellaLimit
     do i = 1, nint(1d7)
        read(25,*,end=30) (arrayIndx(j), j=1,nMolTypes), curValue
!          write(35,*) (arrayIndx(j), j=1,nMolTypes), curValue
        curIndx = getBiasIndex(arrayIndx,NMAX)
        if(curIndx .le. umbrellaLimit) then
          if(curIndx .gt. 0) then
            NBias(curIndx) = curValue
          endif
        endif
     enddo
30   close(25)

     deallocate(arrayIndx)

     if(useWHAM) then
       call WHAM_Initialize
     endif
       
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
     end subroutine
!==========================================================================
     function getBiasIndex result(biasIndx)
     integer :: biasIndx
     integer :: iBias
      

     do iBias = 1, nBiasVariables
       binIndx(iBias) = nint( biasvar(iBias) * binSize(iBias) )
     enddo


     biasIndx = 1
     do iBias = 1, nBiasVariables
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - VarMin(iBias) )
     enddo
     

     end function
!==========================================================================
     function getNewBiasIndex(newArray) result(biasIndx)
     real(dp), intent(in) :: newArray(:)
     integer :: biasIndx
     integer :: iBias
      

     do iBias = 1, nBiasVariables
       binIndx(iBias) = nint( newArray(iBias) * binSize(iBias) )
     enddo


     biasIndx = 1
     do iBias = 1, nBiasVariables
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - VarMin(iBias) )
     enddo
     

     end function
!==========================================================================
    end module
!==========================================================================


