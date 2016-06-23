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
      character(len=35), allocatable :: moveName(:)

      logical :: avbmcUsed, cbmcUsed
      end module
!======================================================
