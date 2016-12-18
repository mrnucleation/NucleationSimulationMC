!=========================================================================
      module WHAM_Functions
      use WHAM_Module
      implicit none

!=========================================================================
      contains
!=========================================================================
      subroutine WHAM_Initialize
      use ParallelVar
!      use SimParameters
      use WHAM_Module
      use UmbrellaSamplingNew
      implicit none
      integer :: AllocateStatus

      if(myid .eq. 0) then
        allocate(WHAM_Numerator(1:umbrellaLimit), STAT = AllocateStatus)
        allocate(WHAM_Denominator(1:umbrellaLimit,1:nWhamItter+1), STAT = AllocateStatus)
        allocate(HistStorage(1:umbrellaLimit), STAT = AllocateStatus)
        allocate(BiasStorage(1:umbrellaLimit,1:nWhamItter+1), STAT = AllocateStatus)
        allocate(FreeEnergyEst(1:umbrellaLimit), STAT = AllocateStatus)
        allocate(ProbArray(1:umbrellaLimit), STAT = AllocateStatus)

!        allocate(NewBias(1:umbrellaLimit), STAT = AllocateStatus)
      
        write(nout,*) "Allocated WHAM Variables"
      
  
        WHAM_Numerator = 0E0
        WHAM_Denominator = 0E0
        HistStorage = 0E0
        BiasStorage = 0E0
        nCurWhamItter = 1
!        tolLimit = 1d-2
        open(unit = 96, file="WHAM_TempHist.incomp")
        open(unit = 97, file="WHAM_Potential.incomp")
        open(unit = 98, file="WHAM_Mid_DG.incomp")
      endif

      allocate(TempHist(1:umbrellaLimit), STAT = AllocateStatus)      
      allocate(NewBias(1:umbrellaLimit), STAT = AllocateStatus)
      NewBias = 0E0
      TempHist = 0E0
      end subroutine

!=========================================================================
!     This subroutine periodically adjusts the Umbrella Sampling Bias
!     by collecting histogram data from across 
      subroutine WHAM_AdjustHist
      use ParallelVar     
      use UmbrellaSamplingNew
      use WHAM_Module
      implicit none
      include 'mpif.h' 

      integer :: arraySize, i, j, cnt, maxbin, maxbin2
      real(dp) :: norm, ratio, maxBias, denomSum
      real(dp) :: F_Estimate(1:nWhamItter),F_Old(1:nWhamItter), fSum     
      real(dp) :: tol, refBias, perError, probNorm


      write(nout,*) "Halting for WHAM"
      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
      
!      This block condences the histogram data from all the different processors
!      into one collective array on the root (myid = 0) processor.        
      arraySize = size(UHist)     
      if(myid .eq. 0) then
        TempHist = 0E0
      endif
      call MPI_REDUCE(UHist, TempHist, arraySize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)       

!      write(*,*) "UHIST!"
      if(myid .eq. 0) then
!        This block calculates the terms needed to 
        norm = sum(TempHist)
        do i = 1, umbrellaLimit
          BiasStorage(i,nCurWhamItter) = UBias(i) 
          if(TempHist(i) .ne. 0E0) then
            WHAM_Numerator(i) = WHAM_Numerator(i) + TempHist(i)*TempHist(i)/norm
            WHAM_Denominator(i, nCurWhamItter) = TempHist(i)*exp(UBias(i))
!            WHAM_Numerator(i) = WHAM_Numerator(i) + TempHist(i)
!            WHAM_Denominator(i, nCurWhamItter) = norm*exp(UBias(i))
            HistStorage(i) = HistStorage(i) + TempHist(i)
          endif
        enddo 

!        This block solves for the free energy terms required by WHAM.  This is done
!        self-consistently.
        ProbArray = 0E0
        F_Estimate = 0E0
        tol = tolLimit + 1E0
        cnt = 0
        do while(tol .gt. tolLimit)
          cnt = cnt + 1
!          Infinite Loop protection
          if(cnt .gt. maxSelfConsist) then 
            write(35,*) "Self Consistent Limit Hit"
            exit
          endif

          ProbArray = 0E0
          do j = 1, nCurWhamItter
            F_Old(j) = F_Estimate(j)
          enddo
!            If bin #i has been sampled at any point in the simulation, estimate the unbiased probability
!            based on the current guess value for F
          do i = 1, umbrellaLimit
            if(WHAM_Numerator(i) .ne. 0E0) then
              denomSum = 0E0
!              maxBias = minval(F_Estimate)
              do j = 1, nCurWhamItter
                if(WHAM_Denominator(i,j) .gt. 0E0) then
                  denomSum = denomSum + WHAM_Denominator(i,j)*exp(-F_Estimate(j))
                endif
              enddo
              if(denomSum .ne. 0E0) then
                ProbArray(i) = WHAM_Numerator(i)/denomSum
              endif
            else
              ProbArray(i) = 0E0
            endif
          enddo 

          norm = sum(ProbArray)
          do i = 1, umbrellaLimit
            ProbArray(i) = ProbArray(i)/norm
          enddo 
!          Once all the unbiased probabilities have been estimated, use these unbiased estimates
!          to calculate a new estimate for F
          do j = 1, nCurWhamItter
            fSum = 0E0
!            maxBias = maxval(BiasStorage(:,j))
            do i = 1, umbrellaLimit
              if(ProbArray(i) .ne. 0E0) then
                fSum = fSum + ProbArray(i)*exp(BiasStorage(i,j))
!                fSum = fSum + ProbArray(i)*exp(BiasStorage(i,j) - maxBias)
              endif
            enddo
!            F_Estimate(j) = log(fSum)
            F_Estimate(j) = log(fSum)
            F_Estimate(j) = (F_Estimate(j) + F_Old(j))*0.5E0
         enddo 
!         Calculate the average change in F from the previous estimate and determine 
!         if there has been a significant change to the F values.
         tol = 0E0
         do j = 1, nCurWhamItter
           tol = tol + abs(F_Estimate(j) - F_Old(j))
         enddo
       enddo

!        Using the new estimates for the unbiased probability, calculate the free energy of nucleation
!        and modify the umbrella sampling bias to
        NewBias = 0E0
        maxbin = maxloc(HistStorage,1)
        if(mod(nCurWhamItter,whamEstInterval) .eq. 0) then
!          maxbin = maxloc(HistStorage,1)
          do i = 1, umbrellaLimit
            if(ProbArray(i) .gt. 0E0) then
              FreeEnergyEst(i) = -log(ProbArray(i)/ProbArray(maxbin))
              NewBias(i) = FreeEnergyEst(i)
            endif
          enddo
          maxBias = -huge(dp)
          do i = 1, umbrellaLimit
            if(ProbArray(i) .gt. 0E0) then
              if(maxBias .lt. FreeEnergyEst(i)) then
                maxBias = FreeEnergyEst(i)
              endif
            endif
          enddo
          do i = 1, umbrellaLimit
            if(ProbArray(i) .le. 0E0) then
              NewBias(i) = maxBias + 1E0
              FreeEnergyEst(i) = maxBias + 1E0
            endif
          enddo
!          call WHAM_CurveSmoothing(NewBias, HistStorage)
        else
          maxbin2 = maxloc(TempHist,1)
          do i = 1, umbrellaLimit
            if(ProbArray(i) .gt. 0E0) then
              FreeEnergyEst(i) = -log(ProbArray(i)/ProbArray(maxbin))
            endif
            if(TempHist(i) .ge. 1E0) then
              NewBias(i) = UBias(i) - UBias(maxbin2) - log(TempHist(i)/TempHist(maxbin2))
            endif
          enddo
!          maxBias = UBias(maxbin2)
          do i = 1, umbrellaLimit
            if(TempHist(i) .lt. 1E0) then
              NewBias(i) = UBias(i) - UBias(maxbin2) + log(TempHist(maxbin2))
            endif
          enddo
        endif
!        Rescale the pontential such that the reference free energy is set to 0
        refBias = NewBias(refBin)
        do i = 1, umbrellaLimit
          NewBias(i) = NewBias(i) - refBias
        enddo
        refBias = FreeEnergyEst(refBin)
        do i = 1, umbrellaLimit
          FreeEnergyEst(i) = FreeEnergyEst(i) - refBias
        enddo
        

      endif
!     End of processor 0 only block

      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
!      Distribute the new free energy estimate to all threads so that they can continue the simulation
!      with the new free energy. 
      arraySize = size(NewBias)      
      call MPI_BCast(NewBias, arraySize, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 

      do i = 1, umbrellaLimit
!        write(*,*) i, UBias(i), Newbias(i)
        UBias(i) = NewBias(i)
        UHist(i) = 0E0
      enddo 
      if(myid .eq. 0) then
        call WHAM_MidSimOutput
      endif

      nCurWhamItter = nCurWhamItter + 1
       
      end subroutine
!==================================================================================
      subroutine WHAM_MidSimOutput
      use UmbrellaSamplingNew
      use ParallelVar
      use SwapBoundary
      use WHAM_Module
      implicit none
      logical :: goToNext
      integer :: arraySize
      integer :: i,j
      integer :: UArray(1:nBiasVariables)
      real(dp) :: refBias
      character(len = 100) :: outputString

      write(outputString, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.10)"

!        This block exports the calculated free energy to a file
      rewind(96)
      do i = 1, umbrellaLimit
        if(HistStorage(i) .ne. 0E0 ) then
          call findVarValues(i, UArray)
          write(96, outputString) (UArray(j)*UBinSize(j),j=1,nBiasVariables), HistStorage(i)
        endif
      enddo
      flush(96)

!        This block exports the current umbrella bias
      rewind(97)
      do i = 1, umbrellaLimit
        call findVarValues(i, UArray)
        write(97, outputString) (UArray(j)*UBinSize(j),j=1,nBiasVariables), UBias(i)
      enddo
      flush(97)
 

!        This block exports the calculated free energy to a file
      rewind(98)
      do i = 1, umbrellaLimit
        if(ProbArray(i) .ne. 0E0 ) then
          call findVarValues(i, UArray)
          write(98, outputString) (UArray(j)*UBinSize(j),j=1,nBiasVariables), FreeEnergyEst(i)
        endif
      enddo
      flush(98)



           
      end subroutine
!==================================================================================
      subroutine WHAM_Finalize
      use UmbrellaSamplingNew
      use ParallelVar
      use SwapBoundary
      use WHAM_Module
      implicit none
      include 'mpif.h' 
      logical :: goToNext
      integer :: arraySize
      integer :: i,j
      integer :: UArray(1:nBiasVariables)
      real(dp) :: norm, ratio, maxBias, refBias, probNorm
      character(len = 100) :: outputString, outputString2

      write(outputString, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.10)"
      write(outputString2, *) "(", (trim(outputFormat(j)), j =1,nBiasVariables), "2x, F18.1)"

!      call WHAM_AdjustHist

      if(myid .eq. 0) then
!        This block exports the calculated free energy to a file
        open(unit = 92, file="WHAM_DG_Output.txt")
        refBias = FreeEnergyEst(refBin)
        do i = 1, umbrellaLimit
          if(ProbArray(i) .gt. 0E0) then
            call findVarValues(i, UArray)
            write(92, outputString) (UArray(j)*UBinSize(j),j=1,nBiasVariables), FreeEnergyEst(i)
          endif
        enddo
        close(92)

!        This block exports the final bias potential to a text file. 
        open(unit = 92, file="WHAM_FinalBias.txt")
!        refBias = UBias(refBin)
        do i = 1, umbrellaLimit
          call findVarValues(i, UArray)
          write(92, outputString) (UArray(j)*UBinSize(j),j=1,nBiasVariables), UBias(i)
        enddo
        close(92)


!        This block exports the calculated probabilities to a file
        open(unit = 92, file="WHAM_Probabilities.txt")
        probNorm = sum(ProbArray)
        do i = 1, umbrellaLimit
          if(ProbArray(i) .gt. 0E0) then
            call findVarValues(i, UArray)
            write(92, outputString) (UArray(j)*UBinSize(j),j=1,nBiasVariables), ProbArray(i)/probNorm
          endif
        enddo
        close(92)



        open(unit=36, file="WHAM_OverallHist.txt")
        do i = 1, umbrellaLimit
          if(HistStorage(i) .gt. 0E0) then
            call findVarValues(i, UArray)
            write(36, outputString2) (UArray(j)*UBinSize(j),j=1,nBiasVariables), HistStorage(i)
          endif
        enddo
        close(36)

        if(WHAM_ExtensiveOutput) then
          open(unit=36, file="WHAM_Hist_Numerator.txt")
          do i = 1, umbrellaLimit
            if(HistStorage(i) .gt. 0E0) then
              call findVarValues(i, UArray)
              write(36, *) (UArray(j)*UBinSize(j),j=1,nBiasVariables), WHAM_Numerator(i)
            endif
          enddo
          close(36)

          open(unit=36, file="WHAM_Hist_Denominator.txt")
          do i = 1, umbrellaLimit
            if(HistStorage(i) .gt. 0E0) then
              call findVarValues(i, UArray)
              write(36, *) (UArray(j)*UBinSize(j),j=1,nBiasVariables), (WHAM_Denominator(i,j), j=1,nCurWhamItter)
            endif
          enddo
          close(36)
        endif
      endif

             
      end subroutine
!=========================================================================     
      end module
