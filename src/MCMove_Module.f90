!======================================================
      module MoveTypeModule
      use AVBMC_Module, only: AVBMC
      use CBMC_Module, only: CBMC
      use Exchange_Module, only: Exchange
      use SimpleMCMoves_Module, only: Translation, Rotation, SingleAtom_Translation
      use VarPrecision

      abstract interface 
        subroutine MCMoveSub(E_T, acc_x, atmp_x)
          use VarPrecision
          implicit none
          real(dp), intent(inout) :: E_T, acc_x, atmp_x
        end subroutine
      end interface 
 
      type MoveArray
        procedure(MCMoveSub), pointer, nopass :: moveFunction => NULL()
      end type

      integer :: nMoveTypes
      type(MoveArray), allocatable :: mcMoveArray(:)
      real(dp), allocatable :: moveProbability(:)
      real(dp), allocatable :: movesAccepted(:), movesAttempt(:)
      real(dp), allocatable, target :: accptRate(:)
      character(len=35), allocatable :: moveName(:)

      logical :: avbmcUsed, cbmcUsed
      contains
!======================================================
      subroutine CalcAcceptanceRates
      implicit none
      integer :: iMoves

      do iMoves = 1, nMoveTypes
        accptRate(iMoves) = movesAccepted(iMoves)/movesAttempt(iMoves)
      enddo

      end subroutine
!======================================================
      subroutine ReadInput_MCMove(fileUnit)
      implicit none
      integer, intent(in) :: fileUnit
      integer :: iMoves, AllocateStatus
      real(dp) :: norm
      character(len=30) :: labelField 
      character(len=30) :: moveName_temp

      norm = 0d0
      read(fileUnit, *) labelField, nMoveTypes
      if(nMoveTypes .le. 0) then
        write(*,*) "ERROR! The user has specified an invalid number of Monte Carlo moves"
        write(*,*) "Please specify at least one valid Monte Carlo move to continue"
        write(*,*) labelField, nMoveTypes
        stop
      endif

      allocate(mcMoveArray(1:nMoveTypes), STAT = AllocateStatus)
      allocate(moveProbability(1:nMoveTypes), STAT = AllocateStatus)
      allocate(movesAccepted(1:nMoveTypes), STAT = AllocateStatus)
      allocate(movesAttempt(1:nMoveTypes), STAT = AllocateStatus)
      allocate(accptRate(1:nMoveTypes), STAT = AllocateStatus)
      allocate(moveName(1:nMoveTypes), STAT = AllocateStatus)
      norm = 0d0
      avbmcUsed = .false.
      cbmcUsed = .false.
      do iMoves = 1, nMoveTypes
        read(fileUnit, *) moveName_temp, moveProbability(iMoves)
        norm = norm + moveProbability(iMoves)
        select case( trim(adjustl(moveName_temp)) )
        case("translation")
          mcMoveArray(iMoves) % moveFunction => Translation
          moveName(iMoves) = "Translation"
        case("rotation")
          mcMoveArray(iMoves) % moveFunction => Rotation
          moveName(iMoves) = "Rotation"
        case("avbmc")
          mcMoveArray(iMoves) % moveFunction => AVBMC
          moveName(iMoves) = "AVBMC"
          avbmcUsed = .true.
        case("cbmc")
          mcMoveArray(iMoves) % moveFunction => CBMC
          moveName(iMoves) = "CBMC"
          cbmcUsed = .true.
        case("exchange")
          mcMoveArray(iMoves) % moveFunction => Exchange
          moveName(iMoves) = "Exchange"
        case("singleatom_translation")
          mcMoveArray(iMoves) % moveFunction => SingleAtom_Translation
          moveName(iMoves) = "Single Atom Translation"
        case default
          write(*,*) "ERROR! Invalid move type specified in input file"
          write(*,*) moveName, moveProbability(iMoves)
          stop
        end select
!        moveName(i) = moveName_temp
      enddo

      do iMoves =1, nMoveTypes
        moveProbability(iMoves) = moveProbability(iMoves)/norm
      enddo
      if(nMoveTypes .gt. 1) then
        do iMoves = 2, nMoveTypes
          moveProbability(iMoves) = moveProbability(iMoves) + moveProbability(iMoves-1)
        enddo
      endif

      end subroutine
!======================================================
      end module
!======================================================
