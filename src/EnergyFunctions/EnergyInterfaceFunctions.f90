!****************************************************************************************

!****************************************************************************************
      module E_Interface
      use CoordinateTypes
!=============================================================================      
      ! interface
         ! subroutine TotalE(E_T)
           ! real(dp), intent(inout) :: E_T
         ! end subroutine

         ! subroutine ShiftE(E_Trial,disp)
           ! real(dp), intent(inout) :: E_Trial
           ! type(Displacement), intent(in) :: disp(:)
         ! end subroutine          
      ! end interface  
        
       
      ! procedure(TotalE), pointer :: TotalECalc => NULL()
      ! procedure(TotalE), pointer :: TotalECalc => NULL()      
!=============================================================================
      contains
!=============================================================================      
      subroutine Detailed_EnergyCalc(E_T,rejMove)
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use BondStretchFunctions
      use BendingFunctions
      use TorsionalFunctions
      use ImproperAngleFunctions
      use EnergyCriteria
      use DistanceCriteria      
      use SimParameters
      implicit none
      
      logical , intent(inout) :: rejMove
      integer :: i,j
      real(dp), intent(inout) :: E_T
      real(dp) :: PairList(1:maxMol,1:maxMol)
      
      E_T = 0d0
      call Detailed_ECalc_Inter(E_T,PairList)

      if(distCriteria) then
        call Detailed_DistanceCriteria(PairList,rejMove)
      else
        call Detailed_EnergyCriteria(PairList,rejMove)      
      endif
      
      write(35,*) "Pairlist:"
      do i = 1, maxMol
        if(isActive(i)) then
          do j = 1, maxMol
            if(isActive(j)) then
              write(35,*) i, j, PairList(i,j)            
            endif
          enddo
        endif
      enddo
      write(35,*)
      
      call Detailed_ECalc_IntraNonBonded(E_T)
      call Detailed_ECalc_BondStretch(E_T)
      call Detailed_ECalc_Bending(E_T)
      call Detailed_ECalc_Torsional(E_T)      
!      call Detailed_ECalc_Improper(E_T)
      
      end subroutine
!=============================================================================      
!     This function contains the energy and cluster criteria functions for any move
!     where a molecule is moved within a cluster.  This function takes a set of displacement
!     vectors (the "disp" variable) and returns the change in energy for both the Intra- and
!     Inter-molecular components.  This function can be used for for moves that any number of
!     atoms in a given molecule, but it can not be used if more than one molecule changes
!     in a given move.  
      subroutine Shift_EnergyCalc(E_Inter, E_Intra, disp, PairList, dETable, useInter, useIntra, rejMove)
      use SimParameters, only: distCriteria, beta, softcutoff, NTotal
      use BendingFunctions      
      use BondStretchFunctions
      use CBMC_Variables   
      use Coords      
      use DistanceCriteria            
      use EnergyCriteria
      use EnergyTables      
      use ImproperAngleFunctions      
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use TorsionalFunctions
      implicit none
      
      logical, intent(in) :: useInter
      logical, intent(in) :: useIntra(1:4)
      type(Displacement), intent(in) :: disp(:)
      logical, intent(inout) :: rejMove
      real(dp), intent(inout) :: dETable(:)
      real(dp), intent(out) :: E_Intra, E_Inter
      real(dp), intent(InOut) :: PairList(:)
      
      integer :: nIndx, nDisp
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      nDisp = size(disp)
      rejMove = .false.
      E_Inter = 0d0
      E_Intra = 0d0
      E_NonBond = 0d0
      E_Stretch = 0d0
      E_Bend = 0d0
      E_Torsion = 0d0
      E_Improper = 0d0
      dETable = 0d0

!     Begin by calculating the intermolecular potential. If any atoms overlap the move will be rejected
!     immediately.      
      if(useInter) then
        if(NTotal .gt. 1) then
          dETable = 0d0
          PairList = 0d0
          call Shift_ECalc_Inter(E_Inter,disp, PairList, dETable, rejMove)
          if(rejMove) then
            return
          endif

          if(E_Inter*beta .gt. softCutOff) then
            rejMove = .true.
            return
          endif
        endif

!        Using the data collected from the intermolecular function, check to see that the new position
!        satisfies the distance criteria.      
        nIndx = MolArray( disp(1)%molType )%mol( disp(1)%molIndx )%indx
        if(distCriteria) then
          if(any(disp(1:nDisp)%atmIndx .eq. 1)) then
            call Shift_DistanceCriteria(PairList, nIndx, rejMove)      
          endif
        else      
          call Shift_EnergyCriteria(PairList, nIndx, rejMove)
        endif        
        if(rejMove) then
          return
        endif
      endif
      
!      This block contains the calculations for all Intramolecular interactions.  For moves
!      that do not change the internal configuration of a molecule (Translation for example)
!      the variable "useIntra" is set to false which skips these calculations. 
      if(regrowType(disp(1)%molType) .ne. 0) then
         if(useIntra(1))  then
           call Shift_ECalc_IntraNonBonded(E_NonBond, disp)    
         endif
         if(useIntra(2))  then
           call Shift_ECalc_BondStretch(E_Stretch, disp)
         endif
         if(useIntra(3)) then
           call Shift_ECalc_Bending(E_Bend, disp)
         endif
         if(useIntra(4))  then
           call Shift_ECalc_Torsional(E_Torsion, disp)
         endif
!           call Shift_ECalc_Improper(E_Improper, disp)      
!         write(2,*) useIntra
!         write(2,*) E_Stretch, E_Bend
         E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
                  
      endif
!      Store the change to
      if(useInter) then
        E_Inter_Diff = E_Inter
      endif
      if(any(useIntra .eqv. .true.)) then
        E_NBond_Diff = E_NonBond
        E_Strch_Diff = E_Stretch
        E_Bend_Diff = E_Bend
        E_Tors_Diff = E_Torsion
      endif
      
      end subroutine
!=============================================================================      
      subroutine SwapIn_EnergyCalc(E_Inter, E_Intra, PairList, dETable, rejMove)
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use BondStretchFunctions
      use BendingFunctions
      use TorsionalFunctions
      use ImproperAngleFunctions
      use EnergyCriteria
      use EnergyTables         
      use Coords
      use CBMC_Variables      
      implicit none
      
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Inter, E_Intra
      real(dp), intent(inout) :: PairList(:), dETable(:)      
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      
      rejMove = .false.      
      
      E_Inter = 0d0
      E_Intra = 0d0
      E_NonBond = 0d0
      E_Stretch = 0d0
      E_Bend = 0d0      
      E_Torsion = 0d0      
      E_Improper = 0d0            
      PairList = 0d0
      dETable = 0d0
!           
      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove)
      if(rejMove) then
        return
      endif
      
!     This block contains the calculations for all Intramolecular interactions.  
      if(regrowType(newMol%molType) .ne. 0) then      
        call NewMol_ECalc_IntraNonBonded(E_NonBond)    
        call NewMol_ECalc_BondStretch(E_Stretch)
        call NewMol_ECalc_Bending(E_Bend)
        call NewMol_ECalc_Torsional(E_Torsion)
!        call NewMol_ECalc_Improper(E_Improper)      
        E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
      endif          

      E_Inter_Diff = 0d0
      E_NBond_Diff = 0d0
      E_Strch_Diff = 0d0
      E_Bend_Diff = 0d0
      E_Tors_Diff = 0d0
      
      E_Inter_Diff = E_Inter
      E_NBond_Diff = E_NonBond
      E_Strch_Diff = E_Stretch
      E_Bend_Diff = E_Bend
      E_Tors_Diff = E_Torsion

      
      end subroutine
!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.  
      subroutine SwapOut_EnergyCalc(E_Inter, E_Intra, nType, nMol, dETable)
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use BondStretchFunctions
      use BendingFunctions
      use TorsionalFunctions
      use ImproperAngleFunctions
      use EnergyCriteria
      use EnergyTables         
      use Coords
      use CBMC_Variables
      implicit none
      
      real(dp), intent(out) :: E_Inter, E_Intra      
      integer, intent(in) :: nType, nMol
      real(dp), intent(inout) :: dETable(:)      
      
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper
      
      E_Inter = 0d0
      E_Intra = 0d0
      E_NonBond = 0d0
      E_Stretch = 0d0
      E_Bend = 0d0      
      E_Torsion = 0d0      
      E_Improper = 0d0            
      dETable = 0d0

      call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter)
      E_Inter = -E_Inter
      
!     This block contains the calculations for all Intramolecular interactions.  
      if(regrowType(nType) .ne. 0) then        
        call Mol_ECalc_IntraNonBonded(nType, nMol, E_NonBond)  
        call Mol_ECalc_BondStretch(nType, nMol, E_Stretch)
        call Mol_ECalc_Bending(nType, nMol, E_Bend)
        call Mol_ECalc_Torsional(nType, nMol, E_Torsion)
!        call Mol_ECalc_Improper(nType, nMol, E_Improper)      
         E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
         E_Intra = -E_Intra
      endif
  
      E_Inter_Diff = 0d0
      E_NBond_Diff = 0d0
      E_Strch_Diff = 0d0
      E_Bend_Diff = 0d0
      E_Tors_Diff = 0d0

      E_Inter_Diff = E_Inter
      E_NBond_Diff = -E_NonBond
      E_Strch_Diff = -E_Stretch
      E_Bend_Diff = -E_Bend
      E_Tors_Diff = -E_Torsion

  
      end subroutine      

!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.  
      subroutine Update_SubEnergies
      use EnergyTables
      implicit none
      
      E_Inter_T = E_Inter_T + E_Inter_Diff
      E_NBond_T = E_NBond_T + E_NBond_Diff
      E_Stretch_T = E_Stretch_T + E_Strch_Diff
      E_Bend_T = E_Bend_T + E_Bend_Diff
      E_Torsion_T = E_Torsion_T + E_Tors_Diff      
      
      end subroutine
!=============================================================================      
      
      end module E_Interface
      
      
