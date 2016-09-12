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

