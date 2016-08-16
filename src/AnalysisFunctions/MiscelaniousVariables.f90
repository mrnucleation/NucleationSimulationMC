!===========================================================
!   The purpose of this module is to provide a set of general purpose
!   memory variables that can be used in the creation of analysis funcitons
!   such as radial distribution functions.  
    module MiscelaniousVars
    use VarPrecision

    type Histograms
      character(len=10) :: histName
      integer :: nBins
      real(dp) :: binSize, varMax
      real(dp), allocatable :: binCount(:)
      real(dp), allocatable :: binValue(:)      
    end type
    
    integer :: nGeometry = 0
    integer :: nHistArrays = 0
    real(dp), allocatable  :: miscCoord(:)
    real(dp), allocatable  :: miscCoord_New(:)
    type(Histograms), allocatable :: miscHist(:)

    contains
!===========================================================
    subroutine ReserveSpace_Coord(nReserve, startIndx, endIndx)
    implicit none
    integer, intent(in) :: nReserve
    integer, intent(out) :: startIndx, endIndx

    startIndx = nReserve + 1
    nHistArrays = nHistArrays + nReserve
    endIndx = nHistArrays

    end subroutine

!===========================================================
    subroutine ReserveSpace_Histograms(nReserve, startIndx, endIndx)
    implicit none
    integer, intent(in) :: nReserve
    integer, intent(out) :: startIndx, endIndx

    startIndx = nReserve + 1
    nGeometry = nGeometry + nReserve
    endIndx = nGeometry

    end subroutine
!===========================================================
    subroutine AllocateMiscArrays
    implicit none
    integer :: AllocationStat

    if(nGeometry .ne. 0) then
      allocate(miscCoord(1:nGeometry), Stat = AllocationStat)
      allocate(miscCoord_New(1:nGeometry), Stat = AllocationStat)
    endif

    if(nHistArrays .ne. 0) then
      allocate(miscHist(1:nHistArrays), Stat = AllocationStat)
    endif

    end subroutine
!===========================================================
    end module
!===========================================================
