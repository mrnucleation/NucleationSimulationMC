!======================================================
     module AnalysisMain
     use RadialDistribution
     use SimpleDistPair
     use MiscelaniousVars

     type AnalysisFunctionArray
       procedure(), pointer, nopass :: func
     end type

     private
     logical :: useAnalysis
     integer :: nPostMove, nMidMove, nOutput
     type(AnalysisFunctionArray), allocatable :: postMoveArray(:)
     type(AnalysisFunctionArray), allocatable :: outputArray(:)

     public :: useAnalysis, DummyAnalysisTest, PostMoveAnalysis, OutputAnalysis

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
     call AllocateHistBins
     radType1(1) = 1
     radType2(1) = 1
     radAtom1(1) = 1
     radAtom2(1) = 1


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
