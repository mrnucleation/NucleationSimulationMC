!===========================================================
!   The purpose of this module is to provide a set of general purpose
!   memory variables that can be used in the creation of analysis funcitons
!   such as radial distribution functions.  
    module MiscelaniousVars
    use VarPrecision

    type Geometry
      character(len=10) :: varName
      real(dp), target :: varValue     
    end type

    type Histograms
      character(len=10) :: histName
      integer :: nBins
      integer :: binSize
      real(dp), allocatable :: binIndex(:)
      real(dp), allocatable :: binValue(:)      
    end type
    
    integer :: nGeometry
    integer :: nHistArrays
    type(Geometry), allocatable  :: miscCoord(:)
    type(Geometry), allocatable  :: miscCoord_New(:)
    type(Histograms), allocatable :: miscHist(:)

    contains

    end module

!===========================================================
