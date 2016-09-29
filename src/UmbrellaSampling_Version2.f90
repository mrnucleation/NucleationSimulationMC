!==========================================================================
    module UmbrellaSamplingNew
    use VarPrecision
    use SimpleDistPair
    implicit none
    private

    type BiasVariablePointer
      integer :: varType
      integer, pointer :: intVar
      real(dp), pointer :: realVar
    end type

    abstract interface
      subroutine UDispFunc(disp)
        use CoordinateTypes
        implicit none
        type(Displacement), intent(in) :: disp(:)
      end subroutine
    end interface

    abstract interface
      subroutine USwapOutFunc(nType, nMol)
        use CoordinateTypes
        implicit none
        integer, intent(in) :: nType, nMol
      end subroutine
    end interface

    type DispUmbrellaArray
      procedure(UDispFunc), pointer, nopass :: func
    end type

    type SwapInUmbrellaArray
      procedure(), pointer, nopass :: func
    end type

    type SwapOutUmbrellaArray
      procedure(USwapOutFunc), pointer, nopass :: func
    end type

    logical :: useUmbrella, UScreenOut
    integer :: nBiasVariables, umbrellaLimit
    integer :: curUIndx
    integer, allocatable :: curVarIndx
    integer, allocatable :: binIndx(:)
    integer, allocatable :: varMax(:), varMin(:)
    integer, allocatable :: binMax(:), binMin(:)
    integer, allocatable :: indexCoeff(:)
    integer, allocatable :: UArray(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UBinSize(:)
    real(dp), allocatable :: varValues(:)
    character(len=20), allocatable :: inputFile
    character(len=10), allocatable :: outputFormat(:)
    character(len=100) :: screenFormat
    type(BiasVariablePointer), allocatable :: biasvar(:)
    type(BiasVariablePointer), allocatable :: biasvarnew(:)

    integer :: nDispFunc, nSwapInFunc, nSwapOutFunc
    type(DispUmbrellaArray), allocatable :: DispUmbrella(:)
    type(SwapInUmbrellaArray), allocatable :: SwapInUmbrella(:)
    type(SwapOutUmbrellaArray), allocatable :: SwapOutUmbrella(:)

    public :: AllocateUmbrellaVariables, ReadInput_Umbrella, AllocateUmbrellaArray, UmbrellaHistAdd
    public :: useUmbrella, OutputUmbrellaHist, GetUmbrellaBias_Disp, findVarValues, getBiasIndex
    public :: nBiasVariables, umbrellaLimit, UBias, UHist, UBinSize, outputFormat, curUIndx
    public :: GetUmbrellaBias_SwapIn, GetUmbrellaBias_SwapOut, ScreenOutputUmbrella, screenFormat
!==========================================================================================
    contains
!==========================================================================================
    subroutine AllocateUmbrellaVariables
    implicit none
    integer :: AllocateStatus
        
    allocate( binMin(1:nBiasVariables), STAT = AllocateStatus )
    allocate( binMax(1:nBiasVariables), STAT = AllocateStatus )
    allocate( UBinSize(1:nBiasVariables), STAT = AllocateStatus )
    allocate( biasvar(1:nBiasVariables), STAT = AllocateStatus )
    allocate( biasvarnew(1:nBiasVariables), STAT = AllocateStatus )
    allocate( binIndx(1:nBiasVariables), STAT = AllocateStatus )
    allocate( indexCoeff(1:nBiasVariables), STAT = AllocateStatus )
    allocate( outputFormat(1:nBiasVariables), STAT = AllocateStatus )

    allocate( DispUmbrella(1:nDispFunc), STAT = AllocateStatus )
    allocate( SwapInUmbrella(1:nSwapInFunc), STAT = AllocateStatus )
    allocate( SwapOutUmbrella(1:nSwapOutFunc), STAT = AllocateStatus )

    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    end subroutine
!==========================================================================================
    subroutine ReadInput_Umbrella(fileUnit)
      use MiscelaniousVars
      use SimpleDistPair, only: nDistPair, pairArrayIndx, CalcDistPairs_New
      use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes, echoInput
      use Q6Functions, only: q6ArrayIndx, CalcQ6_Disp, CalcQ6_SwapIn, CalcQ6_SwapOut
      use ParallelVar, only: nout
      use WHAM_Module, only: refBin, refSizeNumbers
      implicit none
      integer, intent(in) :: fileUnit
      integer :: iUmbrella, AllocateStatus
      integer :: indxVar, stat
      integer :: iDisp, iSwapIn, iSwapOut
      real(dp) :: binSize, valMax, valMin
      real(dp), allocatable :: refVals(:)
      character(len=100), allocatable :: inputLines(:)
      character(len=100) :: refLine
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
        write(nout,*) "Umbrella Sampling Used? = ", useUmbrella
        return
      else
        useUmbrella = .true.
        write(nout,*) "Umbrella Sampling Used? = ", useUmbrella
        write(nout,*) "Number of Umbrella Variables:", nBiasVariables
      endif

      read(fileUnit, "(A)") refLine
      allocate( inputLines(1:nBiasVariables) )
      do iUmbrella = 1, nBiasVariables
        read(fileUnit, "(A)") inputLines(iUmbrella)
        if(echoinput) then
          write(35, "(A)") inputLines(iUmbrella)
        endif 
      enddo

     allocate( refVals(1:nBiasVariables) )

     nDispFunc = 0
     nSwapInFunc = 0
     nSwapOutFunc = 0
     do iUmbrella = 1, nBiasVariables
        read(inputLines(iUmbrella), *) umbrellaName
        select case( trim(adjustl(umbrellaName)) )
        case("clustersize")
!          nSwapInFunc = nSwapInFunc + 1
!          nSwapOutFunc = nSwapOutFunc + 1
        case("totalclustersize")
!          nSwapInFunc = nSwapInFunc + 1
!          nSwapOutFunc = nSwapOutFunc + 1
        case("pairdist")
          nDispFunc = nDispFunc + 1
          nSwapInFunc = nSwapInFunc + 1
          nSwapOutFunc = nSwapOutFunc + 1
        case("q6")
          nDispFunc = nDispFunc + 1
          nSwapInFunc = nSwapInFunc + 1
          nSwapOutFunc = nSwapOutFunc + 1
        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) umbrellaName
          stop
        end select
      enddo



      call AllocateUmbrellaVariables

      iDisp = 0
      iSwapIn = 0
      iSwapOut = 0
      do iUmbrella = 1, nBiasVariables
        read(inputLines(iUmbrella), *) umbrellaName
        select case( trim(adjustl(umbrellaName)) )
        case("clustersize")
          read(inputLines(iUmbrella), *) labelField, indxVar
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
          binMax(iUmbrella) = NMAX(indxVar)
          binMin(iUmbrella) = NMin(indxVar)
          UBinSize(iUmbrella) = 1E0
          outputFormat(iUmbrella) = "2x,F5.1,"
          
!          iSwapIn = iSwapIn + 1
!          iSwapOut = iSwapOut + 1
        case("clustersize")
          read(inputLines(iUmbrella), *) labelField
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
          biasvar(iUmbrella) % intVar => NTotal
          biasvarnew(iUmbrella) % varType = 1
          biasvarnew(iUmbrella) % intVar => NTotal_New
          binMax(iUmbrella) = maxMol
          binMin(iUmbrella) = 1
          UBinSize(iUmbrella) = 1E0
          outputFormat(iUmbrella) = "2x,F5.1,"
          
!          iSwapIn = iSwapIn + 1
!          iSwapOut = iSwapOut + 1
        case("pairdist")
          indxVar = 0
          read(inputLines(iUmbrella), *) labelField, indxVar, valMin, valMax , binSize
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
          binMin(iUmbrella) = nint(valMin / binSize)
          binMax(iUmbrella) = nint(valMax / binSize)
          outputFormat(iUmbrella) = "2x,F12.6,"

          iDisp = iDisp + 1
          DispUmbrella(iDisp) % func => CalcDistPairs_New
          iSwapIn = iSwapIn + 1
          SwapInUmbrella(iSwapIn) % func => CalcDistPairs_SwapIn
          iSwapOut = iSwapOut + 1
          SwapOutUmbrella(iSwapOut) % func => CalcDistPairs_SwapOut

        case("q6")
          read(inputLines(iUmbrella), *) labelField, valMin, valMax, binSize
          biasvar(iUmbrella) % varType = 2
          biasvar(iUmbrella) % realVar => miscCoord(q6ArrayIndx)
          biasvarnew(iUmbrella) % varType = 2
          biasvarnew(iUmbrella) % realVar => miscCoord_New(q6ArrayIndx)
          UBinSize(iUmbrella) = binSize
          binMin(iUmbrella) = nint(valMin / binSize)
          binMax(iUmbrella) = nint(valMax / binSize)
          outputFormat(iUmbrella) = "2x,F12.6,"

          iDisp = iDisp + 1
          DispUmbrella(iDisp) % func => CalcQ6_Disp
          iSwapIn = iSwapIn + 1
          SwapInUmbrella(iSwapIn) % func => CalcQ6_SwapIn
          iSwapOut = iSwapOut + 1
          SwapOutUmbrella(iSwapOut) % func => CalcQ6_SwapOut
!          write(*,*) labelField, valMin, valMax, binSize
        case default
          write(*,*) "ERROR! Invalid variable type specified in input file"
          write(*,*) umbrellaName
          stop
        end select
      enddo
      write(screenFormat,*) "(", (trim(adjustl(outputFormat(iUmbrella))), iUmbrella=1,nBiasVariables), "1x)"

      allocate(varValues(1:nBiasVariables))
      allocate(UArray(1:nBiasVariables))
      deallocate(inputLines)


      call AllocateUmbrellaArray
      call ReadUmbrellaInput

!      Use the input to specify the reference bin
      allocate(refSizeNumbers(1:nBiasVariables),STAT = AllocateStatus)
      read(refLine, *) labelField, (refVals(iUmbrella), iUmbrella=1,nBiasVariables)
      call getUIndexArray(refVals, refBin, stat)
      do iUmbrella = 1, nBiasVariables
        refSizeNumbers(iUmbrella) = refVals(iUmbrella)
      enddo
      if(stat .eq. 1) then
        write(*,*) "ERROR! An invalid ref bin was specified!"
        write(*,*) refLine
        stop
      endif

      deallocate(refVals)


    end subroutine
!==========================================================================================
    subroutine AllocateUmbrellaArray
    use ParallelVar
    implicit none
    integer :: i, j
    integer :: AllocateStatus
        


     indexCoeff(1) = 1
     do i = 2, nBiasVariables 
       indexCoeff(i) = 1
       do j = 1, i-1
         indexCoeff(i) = indexCoeff(i) + indexCoeff(j) * (binMax(j) - binMin(j))
       enddo
     enddo      


     umbrellaLimit = 1
     do i = 1, nBiasVariables 
       umbrellaLimit = umbrellaLimit + indexCoeff(i) * (binMax(i) - binMin(i))
     enddo

     write(nout,*) "Number of Umbrella Bins:", umbrellaLimit
       
     allocate(UBias(1:umbrellaLimit), STAT = AllocateStatus)
     allocate(UHist(1:umbrellaLimit), STAT = AllocateStatus)

     UBias = 0E0_dp
     UHist = 0E0_dp


      
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
     end subroutine
!==========================================================================================
    subroutine ReadUmbrellaInput
    implicit none
    integer :: AllocateStatus
    integer :: j, iInput, iBias, inStat, biasIndx
    real(dp), allocatable :: varValue(:)
    real(dp) :: curBias

    open(unit=80, file="TempIn.txt")
    allocate(varValue(1:nBiasVariables), STAT = AllocateStatus )
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    UBias = 0E0_dp
    do iInput = 1, nint(1d7)
      read(80, *, IOSTAT=inStat) (varValue(j), j=1,nBiasVariables), curBias
      if(inStat .lt. 0) then
        exit
      endif

      call getUIndexArray(varValue, biasIndx, inStat) 
!      write(*,*) varValue, biasIndx
      if(inStat .eq. 1) then
        cycle
      endif

      UBias(biasIndx) = curBias
    enddo

    deallocate(varValue)

    close(80)

    end subroutine
!==========================================================================================
     subroutine UmbrellaHistAdd
     implicit none

     curUIndx = getBiasIndex()
     UHist(curUIndx) = UHist(curUIndx) + 1E0_dp
 
     end subroutine
!==========================================================================================
     subroutine GetUmbrellaBias_Disp(disp, biasDiff, rejMove)
     use CoordinateTypes
     use SimParameters, only: NPART, NPART_new
     implicit none
     logical, intent(out) :: rejMove
     type(Displacement), intent(in) :: disp(:)
     real(dp), intent(out) :: biasDiff
     integer :: iDispFunc, newUIndx, sizeDisp
     real(dp) :: biasOld, biasNew

     if(.not. useUmbrella) then
       biasDiff = 0E0_dp
       return
     endif

     NPART_New = NPART
     NTotal_New = NTotal
     rejMove = .false.
     curUIndx = getBiasIndex()
     biasOld = UBias(curUIndx)
     if(nDispFunc .ne. 0) then
       sizeDisp = size(disp)
       do iDispFunc = 1, nDispFunc
         call DispUmbrella(iDispFunc) % func(disp(1:sizeDisp))
       enddo
     endif

     call getNewBiasIndex(newUIndx, rejMove)
     if(rejMove) then
!       write(*,*) "BLAH!"
       return
     endif
     biasNew = UBias(newUIndx)
     biasDiff = biasNew - biasOld

 
     end subroutine
!==========================================================================================
     subroutine GetUmbrellaBias_SwapIn(biasDiff, rejMove)
     use SimParameters, only: NPART, NPART_new
     implicit none
     logical, intent(out) :: rejMove
     real(dp), intent(out) :: biasDiff
     integer :: iSwapFunc, newUIndx, sizeDisp
     real(dp) :: biasOld, biasNew

     if(.not. useUmbrella) then
       biasDiff = 0E0_dp
       return
     endif
 
     rejMove = .false.
     curUIndx = getBiasIndex()
     biasOld = UBias(curUIndx)
     if(nSwapInFunc .ne. 0) then
       do iSwapFunc = 1, nSwapInFunc
         call SwapInUmbrella(iSwapFunc) % func
       enddo
     endif

     call getNewBiasIndex(newUIndx, rejMove)
     if(rejMove) then
       return
     endif
     biasNew = UBias(newUIndx)
     biasDiff = biasNew - biasOld
 
!     write(2,*) biasNew, biasOld, biasDiff

     end subroutine
!==========================================================================================
     subroutine GetUmbrellaBias_SwapOut(nType, nMol, biasDiff, rejMove)
     use SimParameters, only: NPART, NPART_new
     implicit none
     logical, intent(out) :: rejMove
     integer, intent(in) :: nType, nMol
     real(dp), intent(out) :: biasDiff
     integer :: iSwapFunc, newUIndx, sizeDisp
     real(dp) :: biasOld, biasNew

     if(.not. useUmbrella) then
       biasDiff = 0E0_dp
       return
     endif

     rejMove = .false.
     curUIndx = getBiasIndex()
     biasOld = UBias(curUIndx)
     if(nSwapOutFunc .ne. 0) then
       do iSwapFunc = 1, nSwapInFunc
         call SwapOutUmbrella(iSwapFunc) % func(nType, nMol)
       enddo
     endif

     call getNewBiasIndex(newUIndx, rejMove)
     if(rejMove) then
       return
     endif
     biasNew = UBias(newUIndx)
     biasDiff = biasNew - biasOld

 
     end subroutine
!==========================================================================================
    subroutine ScreenOutputUmbrella
    use ParallelVar
    implicit none
    integer :: iBias, j

    do iBias = 1, nBiasVariables
      if(biasvar(iBias)%varType .eq. 1) then
        varValues(iBias) = real(biasvar(iBias)%intVar,dp)
      else
        varValues(iBias) = biasvar(iBias)%realVar
      endif
    enddo
    write(nout,screenFormat) (varValues(j), j=1,nBiasVariables)

    end subroutine
!==========================================================================================
    subroutine OutputUmbrellaHist
    implicit none
    integer :: iUmbrella, iBias
!    integer, allocatable :: UArray(:)
    character(len = 100) :: outputString


!    allocate(UArray(1:nBiasVariables))
!    allocate(varValues(1:nBiasVariables))

    write(outputString, *) "(", (trim(outputFormat(iBias)), iBias =1,nBiasVariables), "2x, F18.1)"
    open(unit=60, file="TemporaryHist.txt")
    do iUmbrella = 1, umbrellaLimit
      call findVarValues(iUmbrella, UArray)
      do iBias = 1, nBiasVariables
        varValues(iBias) = real( UArray(iBias), dp) * UBinSize(iBias)
      enddo
      if(UHist(iUmbrella) .ne. 0E0_dp) then
        write(60,outputString) (varValues(iBias), iBias =1,nBiasVariables), UHist(iUmbrella)
      endif
    enddo 
    flush(60)
    close(60)

!    deallocate(UArray)
!    deallocate(varValues)

    end subroutine

!==========================================================================
     function getBiasIndex() result(biasIndx)
     integer :: biasIndx
     integer :: iBias
      

     do iBias = 1, nBiasVariables
       if(biasvar(iBias) % varType .eq. 1) then
         binIndx(iBias) = floor( biasvar(iBias) % intVar / UBinSize(iBias) )
!         write(*,*) iBias, biasvar(iBias)%intVar, binIndx(iBias)
       elseif(biasvar(iBias) % varType .eq. 2) then
         binIndx(iBias) = floor( biasvar(iBias) % realVar / UBinSize(iBias) )
!         write(*,*) iBias, biasvar(iBias)%realvar, binIndx(iBias)
       endif
     enddo


     biasIndx = 1
     do iBias = 1, nBiasVariables
!       write(*,*) biasIndx,indexCoeff(iBias), binIndx(iBias), binMin(iBias) 
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - binMin(iBias) )
     enddo

!     write(*,*) biasIndx
 
     end function
!==========================================================================
     subroutine getNewBiasIndex(biasIndx, rejMove)
     logical, intent(out) :: rejMove
     integer, intent(out) :: biasIndx
     integer :: iBias
      
     rejMove = .false.
     do iBias = 1, nBiasVariables
       if(biasvarnew(iBias) % varType .eq. 1) then
         binIndx(iBias) = biasvarnew(iBias)%intVar
       elseif(biasvar(iBias) % varType .eq. 2) then
         binIndx(iBias) = floor( biasvarnew(iBias)%realVar / UBinSize(iBias) )
       endif
       if(binIndx(iBias) .lt. binMin(iBias) ) then
         rejMove = .true.
         return
       endif
       if(binIndx(iBias) .gt. binMax(iBias) ) then
         rejMove = .true.
         return
       endif

     enddo

     biasIndx = 1
     do iBias = 1, nBiasVariables
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - binMin(iBias) )
     enddo
     if(biasIndx .gt. umbrellaLimit ) then
       rejMove = .true.
       return
     endif

     end subroutine
!==========================================================================
     subroutine getUIndexArray(varArray, biasIndx, stat) 
     real(dp), intent(in) :: varArray(:)
     integer, intent(out) :: biasIndx, stat
     integer :: iBias
      
     stat = 0
     do iBias = 1, nBiasVariables
!       binIndx(iBias) = floor( varArray(iBias) / UBinSize(iBias) + 1E-8 )
       binIndx(iBias) = nint( varArray(iBias) / UBinSize(iBias) )
!       write(*,*) nint( varArray(iBias) / UBinSize(iBias) ), floor( varArray(iBias) / UBinSize(iBias) + 1E-7 )
       if(binIndx(iBias) .gt. binMax(iBias)) then
         stat = 1
         return
       endif
       if(binIndx(iBias) .lt. binMin(iBias)) then
         stat = 1
         return
       endif
     enddo


     biasIndx = 1
     do iBias = 1, nBiasVariables
       biasIndx = biasIndx + indexCoeff(iBias) * ( binIndx(iBias) - binMin(iBias) )
     enddo
     

     end subroutine
!===========================================================================
     subroutine findVarValues(UIndx, UArray)
     implicit none
     integer, intent(in) :: UIndx
     integer, intent(inout) :: UArray(:)
     integer :: i, iBias
     integer :: remainder, curVal

     remainder = UIndx - 1
     do i = 1, nBiasVariables
       iBias = nBiasVariables - i + 1
       curVal = int( real(remainder, dp)/real(indexCoeff(iBias),dp) )
       UArray(iBias) = curVal + binMin(iBias)
       remainder = remainder - curVal * indexCoeff(iBias)
     enddo
    
     end subroutine
!==========================================================================
    end module
!==========================================================================


