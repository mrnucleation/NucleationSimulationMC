!================================================================================================================
     module EnergyPointers
     use E_Interface
     interface 
       subroutine DetailedInterface(E_T,rejMove)
         use VarPrecision
         implicit none
         logical , intent(inout) :: rejMove
         real(dp), intent(inout) :: E_T
       end subroutine
     end interface 

     interface 
       subroutine ShiftInterface(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, useInter)
         use CoordinateTypes
         use VarPrecision
         implicit none
         logical, intent(in), optional :: useInter
         logical, intent(in) :: useIntra(1:4)
         type(Displacement), intent(in) :: disp(:)
         logical, intent(inout) :: rejMove
         real(dp), intent(inout) :: dETable(:)
         real(dp), intent(out) :: E_Intra, E_Inter
         real(dp), intent(InOut) :: PairList(:)
       end subroutine
     end interface 

     interface 
       subroutine SwapInInterface(E_Inter, E_Intra, PairList, dETable, rejMove, useInter)
         use VarPrecision
         implicit none
         logical, intent(out) :: rejMove
         logical, intent(in), optional :: useInter
         real(dp), intent(out) :: E_Inter, E_Intra
         real(dp), intent(inout) :: PairList(:), dETable(:)  
       end subroutine
     end interface 

     interface 
       subroutine SwapOutInterface(E_Inter, E_Intra, nType, nMol, dETable, useInter)
         use VarPrecision
         implicit none
         logical, intent(in), optional :: useInter
         real(dp), intent(out) :: E_Inter, E_Intra      
         integer, intent(in) :: nType, nMol
         real(dp), intent(inout) :: dETable(:)   
       end subroutine
     end interface 
     
     procedure(DetailedInterface), pointer :: Detailed_ECalc => NULL()
     procedure(ShiftInterface), pointer  :: Shift_ECalc => NULL()
     procedure(SwapInInterface), pointer :: SwapIn_ECalc => NULL()
     procedure(SwapOutInterface), pointer :: SwapOut_ECalc => NULL()

     end module
!================================================================================================================
