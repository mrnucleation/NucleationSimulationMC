!================================================================================
      module SwapBoundary

      interface 
        function boundFunction(NPART, NDiff) result(rejMove)
          logical :: rejMove
          integer, intent(in) :: NPART(:), NDiff(:)
        end function
      end interface

      procedure(boundFunction), pointer :: BoundaryFunction => null()



!================================================================================
      contains
!================================================================================
      function Bound_MaxMin(NPART, NDiff) result(rejMove)
      use SimParameters, only: NMAX, NMIN
      logical :: rejMove
      integer, intent(in) :: NPART(:), NDiff(:)
      integer :: i, nTypes

      nTypes = size(NPART)
      rejMove = .false.
      do i = 1, nTypes
        if( (NPART(i)+NDiff(i)) .gt. NMAX(i)  ) then
          rejMove = .true.
          return
        endif
        if( (NPART(i)+NDiff(i)) .lt. NMIN(i)  ) then
          rejMove = .true.
          return
        endif
      enddo

      end function bound_MaxMin
!================================================================================
      function Bound_PedoneChargeBalance(NPART, NDiff) result(rejMove)
      use SimParameters, only: NMAX, NMIN
      use ForceFieldPara_Pedone
      logical :: rejMove
      integer, intent(in) :: NPART(:), NDiff(:)
      integer :: i, nTypes
      real(dp) :: charge
      real(dp) :: q

      nTypes = size(NPART)
      rejMove = .false.
      do i = 1, nTypes
        if( (NPART(i)+NDiff(i)) .gt. NMAX(i)  ) then
          rejMove = .true.
          return
        endif
        if( (NPART(i)+NDiff(i)) .lt. NMIN(i)  ) then
          rejMove = .true.
          return
        endif
      enddo

      charge = 0d0
      do i = 1, nTypes
        q = pedoneData(i)%q
        charge = charge + q*(NPART(i)+NDiff(i))
      enddo


      if(abs(charge) .gt. 4d0) then
        rejMove = .true.
        return
      endif

      end function bound_PedoneChargeBalance

!================================================================================
      end module
!================================================================================
