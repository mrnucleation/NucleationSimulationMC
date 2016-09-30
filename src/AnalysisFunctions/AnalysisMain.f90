!======================================================
     module AnalysisMain
     use RadialDensity
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
     integer :: nAnalysisVar
     integer :: nPostMove, nOutput, nTrialVar
     type(TrialFunctionArray), allocatable :: TrialArray(:)
     type(AnalysisFunctionArray), allocatable :: postMoveArray(:)
     type(AnalysisFunctionArray), allocatable :: outputArray(:)
     type(TrialFunctionArray), allocatable :: UpdateArray(:)

     public :: useAnalysis, ReadAnalysisInput, DummyAnalysisTest2, PostMoveAnalysis, OutputAnalysis

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
     call AllocateHistBins
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
     call AllocateMiscArrays
     call SetPairVariables(1, 1, 1, 1, 1, 2, 1)

     end subroutine
!======================================================
     subroutine ReadAnalysisInput(fileUnit)
      use MiscelaniousVars
      use SimpleDistPair, only: nDistPair, pairArrayIndx
      use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes
      use Q6Functions, only: CalcQ6, Initialize_q6, useQ6, q6Dist, q6DistSq
      implicit none
      integer, intent(in) :: fileUnit
      integer :: iAnalysis, AllocateStatus
      integer :: indxVar, nBins
      integer :: iRadial, iDistPair, iRadDens, iPostMove, iOutput
      integer :: type1, type2, mol1, mol2, atom1, atom2
      real(dp) :: binSize, realVar
      character(len=200), allocatable :: inputLines(:)
      character(len=30) :: labelField 
      character(len=30) :: analysisName, fileName

      read(fileUnit, *) labelField, nAnalysisVar
      if(nAnalysisVar .lt. 0) then
        write(*,*) "ERROR! The user has specified an invalid number of Analysis Variables"
        write(*,*) labelField, nAnalysisVar
        stop
      endif
      if(nAnalysisVar .eq. 0) then
        useAnalysis = .false.
        return
      else
        useAnalysis = .true.
      endif

      allocate( inputLines(1:nAnalysisVar), stat = AllocateStatus ) 
!      Due to the nature how the analysis arrays are structured, it is nessiary
!      to first gather the lines
      do iAnalysis = 1, nAnalysisVar
        read(fileUnit, "(A)") inputLines(iAnalysis)
!        write(*,*) inputLines(iAnalysis)
      enddo


!      Begin by counting how many entries are required for each function array.
!      These will be used in the next step to allocate the arrays. 
      nPostMove = 0
      nOutput = 0
      nDistPair = 0 
      nRadialDist = 0
      nRadialDens = 0
      do iAnalysis = 1, nAnalysisVar
        read(inputLines(iAnalysis), *) analysisName

        select case( trim(adjustl(analysisName)) )
        case("radialdistribution")
          nRadialDist = nRadialDist + 1
          if(nRadialDist .eq. 1) then
            nPostMove = nPostMove + 1
            nOutput = nOutput + 1
          endif
        case("pairdist")
          nDistPair = nDistPair + 1
          if(nDistPair .eq. 1) then
            nPostMove = nPostMove + 1
          endif
        case("q6")
          if(useQ6 .eqv. .false.) then
            useQ6 = .true.
            nPostMove = nPostMove + 1
          else
            stop "ERROR! The Q6 analysis function has been defined more than once in the input script."
          endif
        case("radialdensity")
          nRadialDens = nRadialDens + 1
          if(nRadialDens .eq. 1) then
            nPostMove = nPostMove + 1
            nOutput = nOutput + 1
          endif
        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) analysisName
          stop
        end select
      enddo

!      Now that we know how much memory is required, allocate all relevent arrays
      if(nPostMove .ne. 0) then
        allocate( postMoveArray(1:nPostMove) )
      endif
      if(nOutput .ne. 0) then
        allocate( outputArray(1:nOutput) )
      endif 
      call Initialize_RadialDist
      call Initialize_DistPair
      call Initialize_RadialDens
      call Initialize_Q6
      call AllocateMiscArrays

      iOutPut = 0
      iPostMove = 0
      iRadial = 0
      iDistPair = 0
      iRadDens = 0
      do iAnalysis = 1, nAnalysisVar
        read(inputLines(iAnalysis), *) analysisName

        select case( trim(adjustl(analysisName)) )
        case("radialdistribution")
          iRadial = iRadial + 1
          read(inputLines(iAnalysis), *)  analysisName, type1, atom1, type2, atom2, binSize, nBins, fileName

          call SetRadialParameters(iRadial, type1, type2, atom1, atom2)
          call SetRadialHist(iRadial, binSize, nBins, fileName)

          if(iRadial .eq. 1) then
            iPostMove = iPostMove + 1
            postMoveArray(iPostMove)%func => Calc_RadialDist
            iOutPut = iOutPut + 1
            outputArray(iOutPut)%func => Output_RadialDist
          endif 
        case("pairdist")
          iDistPair = iDistPair + 1
          read(inputLines(iAnalysis), *)  analysisName, type1, mol1, atom1, type2, mol2, atom2
          call SetPairVariables(iDistPair, Type1, Mol1, Atom1, Type2, Mol2, Atom2)
          if(iDistPair .eq. 1) then
            iPostMove = iPostMove + 1
            postMoveArray(iPostMove)%func => CalcDistPairs
          endif 
        case("q6")
          read(inputLines(iAnalysis), *)  analysisName, realVar
          iPostMove = iPostMove + 1
          q6Dist = realVar
          q6DistSq = q6Dist*q6Dist
          postMoveArray(iPostMove)%func => CalcQ6
        case("radialdensity")
          read(inputLines(iAnalysis), *)  analysisName, type1, binSize, nBins, fileName
          iPostMove = iPostMove + 1
          iRadDens = iRadDens + 1
          iOutPut = iOutPut + 1
          call SetDensityParameters(iRadDens, type1)
          call SetDensityHist(iRadDens, binSize, nBins, fileName)

          postMoveArray(iPostMove)%func => Calc_RadialDensity
          outputArray(iOutPut)%func => Output_RadialDensity
        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) analysisName
          stop
        end select
      enddo

      call AllocateHistBins

      deallocate( inputLines ) 

     end subroutine
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
