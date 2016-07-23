!==========================================================================================
      subroutine IonHistAdd(E_T)
      use IonBias
      use SimParameters
      implicit none
      integer :: i, bin, sizeN
      integer :: curIndx,maxIndx
      real(dp),intent(inout) :: E_T
      
      sizeN = size(NPART)
!       curIndx = 0
      curIndx = NPart(sizeN)
      maxIndx = NMAX(sizeN) + 1 
      do i=1,sizeN-1
        curIndx = curIndx + maxIndx*NPart(sizeN-i)
        maxIndx = maxIndx * (NMAX(sizeN-i) + 1)
      enddo
      bin = floor(ionDistance*dR)
      NHist(curIndx+1, bin) = NHist(curIndx+1, bin) + 1d0
      E_Avg(curIndx+1) = E_Avg(curIndx+1) + E_T
      
      end subroutine
  !========================================================    
!     This subrotuine collects the cluster size Histogram from across each processor
!     and collapses them into a single output file for data analysis.  
      subroutine CollapseIon
      use IonBias
      use SimParameters
      use ParallelVar  
      implicit none
      include 'mpif.h'  
      integer :: i, j, arraySize, lowerBound(1:2), upperBound(1:2)
      real(dp), allocatable :: Hist_temp(:)


      lowerBound = lbound(NR_Hist)
      upperBound = ubound(NR_Hist)
      arraySize = size(NR_Hist)
!      write(35,*) "Array Bounds:", lowerBound, upperBound
!      flush(35)
      allocate( Hist_Temp(lowerBound(1):upperBound(1), lowerBound(2):upperBound(2)))

      Hist_temp = 0d0
      
      call Single_Ion_HistOutput
      write(nout,*) "Compressing Histogram Data..."           
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      call MPI_REDUCE(NR_Hist, Hist_temp, arraySize, 
     &           MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,
     &           ierror)  

      write(nout,*) "Finished Compressing Histogram Data..."      
      if(myid .eq. 0) then
        do i = lowerBound(1), upperBound(1)
          do j = lowerBound(2), upperBound(2)
            NR_Hist(i,j) = Hist_temp(i,j)
          enddo
        enddo
      endif

 
      end subroutine
!==================================================================================
      subroutine Ion_HistOutput
      use IonBias
      use SimParameters
      use ParallelVar
      implicit none
      integer :: i,j
      integer :: NArray(1:nMolTypes)
            
      open(unit = 92, file="IonPair_Histogram.txt")

      NArray = 0
!      NArray(nMolTypes) = NArray(nMolTypes) + 1      
      do i = 1, umbrellaLimit
        if(nMolTypes .gt. 1) then
          do j = 1,nMolTypes-1
            if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1))then
              NArray(nMolTypes - j + 1) = 0
              NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
            endif
          enddo
         endif
!       write(nout,*) NArray
        do j = 0, rBin
          if(NR_Hist(i,j) .ne. 0d0) then
            write(92, *) (NArray(j),j=1,nMolTypes), j*dR, NR_Hist(i)    
          endif       
        enddo

        NArray(nMolTypes) = NArray(nMolTypes) + 1
      enddo
      close(92)
       
      end subroutine
!==================================================================================
      subroutine Single_Ion_HistOutput
      use SimParameters
      use ParallelVar
      implicit none
      integer :: i,j
      integer :: NArray(1:nMolTypes)
     
      NArray = 0
      write(35,*)"---------------------------------"
      write(35,*)"Cluster Size Histogram:"
!      NArray(nMolTypes) = NArray(nMolTypes) + 1      
      do i = 1, umbrellaLimit
       do j = 1,nMolTypes-1
         if(NArray(nMolTypes - j + 1) .gt. NMAX(nMolTypes - j + 1)) then
          NArray(nMolTypes - j + 1) = 0
          NArray(nMolTypes - j) = NArray(nMolTypes - j) + 1          
         endif
       enddo
!       write(nout,*) NArray
       if(NHist(i) .ne. 0d0) then
         write(35,*) (NArray(j),j=1,nMolTypes), NHist(i)    
       endif       
       NArray(nMolTypes) = NArray(nMolTypes) + 1       
      enddo   
      write(35,*)"---------------------------------"       
!      flush(35)
      end subroutine      
