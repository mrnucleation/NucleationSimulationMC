!======================================================
     module AnalysisMain
     use RadialDistribution
     use SimpleDistPair
     use MiscelaniousVars

     interface 
       subroutine TrialFunction(disp)
        use CoordinateTypes
        type(Displacement), intent(in) :: disp(:)
       end subroutine
     end interface

     type TrialFunctionArray
       procedure(TrialFunction), pointer, nopass :: func
     end type

     type AnalysisFunctionArray
       procedure(), pointer, nopass :: func
     end type

     private
     logical :: useAnalysis
     integer :: nPostMove, nMidMove, nOutput, nTrialVar
     type(TrialFunctionArray), allocatable :: TrialArray(:)
     type(AnalysisFunctionArray), allocatable :: postMoveArray(:)
     type(AnalysisFunctionArray), allocatable :: outputArray(:)

     public :: useAnalysis, DummyAnalysisTest2, PostMoveAnalysis, OutputAnalysis

!======================================================
     contains
!======================================================
     subroutine DummyAnalysisTest
     implicit none

     useAnalysis = .true.
     nPostMove = 1
     nOutput = 1
     allocate(postMoveArray(1:nPostMove))
     allocate(outputArray(1:nOutput))
     postMoveArray(1)%func => Calc_RadialDist
     outputArray(1)%func => Output_RadialDist
     nRadialDist = 1
     call Initialize_RadialDist
     call AllocateArrays
     miscHist(1)%binSize = 0.01E0
     miscHist(1)%sizeInv = 1E0/miscHist(1)%binSize
     miscHist(1)%nBins = 1000
     miscHist(1)%fileName = "RadialDistribution.txt"

!     miscHist(2)%binSize = 0.01E0
!     miscHist(2)%sizeInv = 1E0/miscHist(1)%binSize
!     miscHist(2)%nBins = 1000
!     miscHist(2)%fileName = "RadialDistribution2.txt"
     call AllocateHistBins
     call SetRadialParameters(1, 1, 1, 1, 1)
!     call SetRadialParameters(2, 1, 1, 1, 2)


     end subroutine
!======================================================
     subroutine DummyAnalysisTest2
     implicit none

     useAnalysis = .true.
     nPostMove = 1
     nOutput = 0
     allocate(postMoveArray(1:nPostMove))
     postMoveArray(1)%func => CalcDistPairs
     nDistPair = 1
     call Initialize_DistPair
     call AllocateArrays
     call SetPairVariables(1, 1, 1, 1, 1, 2, 1)

     end subroutine
!======================================================
!     subroutine ReadAnalysisInput(fileUnit)
!     implicit none
!     integer, intent(in) :: fileUnit
!     integer :: 
!
!
!     end subroutine
!======================================================
     subroutine TrialPositionAnalysis(disp)
     use CoordinateTypes
     implicit none
     type(Displacement), intent(in) :: disp(:)
     integer :: iTrialVar

     do iTrialVar = 1, nTrialVar
       call TrialArray(iTrialVar)%func(disp)
     enddo

     end subroutine
!======================================================
     subroutine PostMoveAnalysis
     implicit none
     integer :: iPostMove

     do iPostMove = 1, nPostMove
       call postMoveArray(iPostMove)%func
     enddo

     end subroutine
!======================================================
     subroutine OutputAnalysis
     use ParallelVar
     implicit none
     integer :: iOutput

     do iOutput = 1, nOutput
       call outputArray(iOutput)%func
     enddo

     end subroutine
!======================================================
     end module
!======================================================
