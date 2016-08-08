!================================================================================
      module SwapBoundary

      interface 
        function boundFunction(NDiff) result(rejMove)
        logical :: rejMove
        integer, intent(in) :: NDiff(:)
        end function
      end interface

      procedure(boundFunction), pointer :: boundaryFunction => null()



!================================================================================
      contains
!================================================================================
      function Bound_MaxMin(NDiff) result(rejMove)
      use SimParameters
      logical :: rejMove
      integer, intent(in) :: NDiff(:)
      integer :: i

      rejMove = .false.
      do i = 1, nMolTypes
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
      function Bound_PedoneChargeBalance(NDiff) result(rejMove)
      use SimParameters
      use ForceFieldPara_Pedone
      logical :: rejMove
      integer, intent(in) :: NDiff(:)
      integer :: i
      real(dp) :: charge
      real(dp) :: q

      rejMove = .false.
      do i = 1, nMolTypes
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
      do i = 1, nMolTypes
        q = pedoneData(i)%q
        charge = charge + q*(NPART(i)+NDiff(i))
      enddo

      if(abs(charge) .gt. 4d0) then
        rejMove = .true.
      endif

      end function bound_PedoneChargeBalance

!================================================================================
      end module
!================================================================================
