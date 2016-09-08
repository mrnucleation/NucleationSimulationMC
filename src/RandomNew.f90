!=======================================================
      real(dp) function grnd()
      use VarPrecision
      implicit none
      real(dp) :: r

      call RANDOM_NUMBER(r)
      grnd = r

     
      end function
!=======================================================
      subroutine sgrnd(seed)
      implicit none
      integer,intent(inout) :: seed
      integer :: i,n
      integer, allocatable :: tempSeed(:)
      
      call RANDOM_SEED(size=n)      
      allocate(tempSeed(1:n))
      tempSeed = seed + 37 * (/ (i - 1, i = 1, n) /)

      call RANDOM_SEED(put=tempSeed)
     
      deallocate(tempSeed)
     
      end subroutine      
!=======================================================
      subroutine SelectBin_IntTable(integralTable, nBin)
      use VarPrecision
      implicit none
      real(dp), intent(in) :: integralTable(:)
      integer, intent(out) :: nBin
      integer :: iBin, binMin, binMax, binMid
      integer :: lowIndx, highIndx
      real(dp) ::  grnd, ranNum

      binMin = LBOUND(integralTable,1)
      binMax = size(integralTable) + binMin - 1
      ranNum = grnd()

!      write(*,*) binMin, binMax
 
      do while(abs(binMax - binMin) > 50) 
        binMid = nint(0.2d0*binMax + 0.8d0*binMin)
        if(integralTable(binMid) .gt. ranNum) then
          binMax = binMid
        endif
        if(integralTable(binMid) .lt. ranNum) then
          binMin = binMid
        endif
      enddo


      do iBin = binMin, binMax
        if(integralTable(iBin) .ge. ranNum) then
          nBin = iBin
          return
        endif
      enddo



      end subroutine  
!=======================================================

