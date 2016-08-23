!==========================================================================
    module UmbrellaSamplingNew
    use VarPrecision
    implicit none
    private

    type BiasVariablePointer
      real(dp), pointer :: var
    end type

    type BiasIndex
      integer, allocatable :: binNumber(:)
      real(dp), allocatable :: binValue(:)
    end type
   

    logical :: useUmbrella
    integer :: nBiasVariables, umbrellaLimit
    integer, allocatable :: VarMax(:), VarMin(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: binSize(:)
    type(BiasVariablePointer), allocatable :: biasvar(:)
    type(BiasVariablePointer), allocatable :: biasvarnew(:)

    public :: AllocateUmbrellaBias

!==========================================================================
    contains
!==========================================================================================
    subroutine AllocateUmbrellaVariables
    use UmbrellaFunctions
    use SimParameters
    use WHAM_Module
    implicit none
    integer :: i,j
    integer :: AllocateStatus
    integer :: curIndx
    real(dp) :: curValue, defaultVal
        
    allocate(VarMin(1:nBiasVariables), STAT = AllocateStatus)
    allocate(VarMax(1:nBiasVariables), STAT = AllocateStatus)
    allocate(binSize(1:nBiasVariables), STAT = AllocateStatus)
    allocate(biasvar(1:nBiasVariables), STAT = AllocateStatus)
    allocate(biasvarnew(1:nBiasVariables), STAT = AllocateStatus)
     
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    end subroutine
!==========================================================================================
    subroutine AllocateUmbrellaArray(fileName)
    use UmbrellaFunctions
    use SimParameters
    use WHAM_Module
    implicit none
    integer :: i,j
    integer :: AllocateStatus
    integer,allocatable :: arrayIndx(:)   
    integer :: curIndx
    character(len=30) :: fileName   
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
      do i = 1,umbrellaLimit
!      do i = 1,nint(1d7)
         read(25,*,end=30) (arrayIndx(j), j=1,nMolTypes), curValue
!         write(35,*) (arrayIndx(j), j=1,nMolTypes), curValue
         curIndx = getBiasIndex(arrayIndx,NMAX)
         if(curIndx .le. umbrellaLimit) then
           if(curIndx .gt. 0) then
             NBias(curIndx) = curValue
           endif
         endif
      enddo
30    close(25)

      deallocate(arrayIndx)

      if(useWHAM) then
        call WHAM_Initialize
      endif
       
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    end subroutine

!==========================================================================
    end module
!==========================================================================

