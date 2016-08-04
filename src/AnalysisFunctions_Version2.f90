!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies
      module SimpleDistPair
        type PairData
          integer :: type1, mol1, atom1
          integer :: type2, mol2, atom2
        end type

        integer :: nDistPair
        integer, allocatable :: distPairIndx(:)
        type(PairData), allocatable :: distPairData(:)
  
        contains
     !--------------------------------------------------------------------------------
        subroutine Initialize_DistPair
           use MiscelaniousVars
           implicit none 
           integer :: AllocationStatus

           allocate(distPairIndx(1:nDistPair), STAT = AllocationStatus)
           allocate(distPairData(1:nDistPair), STAT = AllocationStatus)
        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs
          use MiscelaniousVars
          use Coords
          implicit none 
          integer :: iDistPair
          integer :: type1, mol1, atom1
          integer :: type2, mol2, atom2
          real(dp) :: rx, ry, rz, r

          do iDistPair =1, nDistPair
            type1 = distPairData(iDistPair)%type1
            type2 = distPairData(iDistPair)%type2
            mol1 = distPairData(iDistPair)%mol1
            mol2 = distPairData(iDistPair)%mol2
            atom1 = distPairData(iDistPair)%atom1
            atom2 = distPairData(iDistPair)%atom2
            rx = MolArray(type1)%mol(mol1)%x(atom1) - MolArray(type2)%mol(mol2)%x(atom2)
            ry = MolArray(type1)%mol(mol1)%y(atom1) - MolArray(type2)%mol(mol2)%y(atom2)
            rz = MolArray(type1)%mol(mol1)%z(atom1) - MolArray(type2)%mol(mol2)%z(atom2)
            r = rx*rx + ry*ry + rz*rz
            r = dsqrt(r)
            miscCoord(distPairIndx(iDistPair))%varValue = r
          enddo


        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs_New(disp)
          use MiscelaniousVars
          use Coords
          implicit none 
          type(Displacement), intent(in) :: disp(:)
          integer :: nDisp, iDistPair
          integer :: type1, mol1, atom1
          integer :: type2, mol2, atom2
          integer :: dispIndx, indx1, indx2
          real(dp) :: rx, ry, rz, r

          dispIndx = disp(1)%molIndx

          nDisp = size(disp)

          do iDistPair = 1, nDistPair
            if(dispIndx .ne. iIndx) then
              if(dispIndx .ne. jIndx) then
                miscCoord_New(distPairIndx(iDistPair))%varValue = r  
              endif
            endif
            type1 = distPairData(iDistPair)%type1
            type2 = distPairData(iDistPair)%type2
            mol1 = distPairData(iDistPair)%mol1
            mol2 = distPairData(iDistPair)%mol2
            atom1 = distPairData(iDistPair)%atom1
            atom2 = distPairData(iDistPair)%atom2
            rx = MolArray(type1)%mol(mol1)%x(atom1) - MolArray(type2)%mol(mol2)%x(atom2)
            ry = MolArray(type1)%mol(mol1)%y(atom1) - MolArray(type2)%mol(mol2)%y(atom2)
            rz = MolArray(type1)%mol(mol1)%z(atom1) - MolArray(type2)%mol(mol2)%z(atom2)
            r = rx*rx + ry*ry + rz*rz
            r = dsqrt(r)
            miscCoord_New(distPairIndx(iDistPair))%varValue = r            
          enddo


        end subroutine

     !--------------------------------------------------------------------------------
      end module
!====================================================================================== 



!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies

      module RadialDist
        integer :: nRadialDist
        integer, allocatable :: radHistIndx(:)
        integer, allocatable :: atomPairIndex(:)
        contains
   !--------------------------------------------------------------------------------
        subroutine Initialize_RadialDist
        use MiscelaniousVars
        implicit none 
        integer :: iRadial
        integer :: nBins, AllocationStatus

        do iRadial = 1, nRadialDist
          nBins = miscHist( radhistIndx(i) )%nBins
          allocate( miscHist(radHistIndx(i))%binIndex(1:nBins), stat = AllocationStatus )
          allocate( miscHist(radHistIndx(i))%binValue(1:nBins), stat = AllocationStatus )
        enddo

        end subroutine
        
   !--------------------------------------------------------------------------------
!        subroutine Calc_RadialDist
!        use Coords
!        implicit none
!        integer :: iRadial
!        do iRadial = 1, nRadialDist
!         
!        enddo
!        
!      
!        end subroutine  

      end module
!======================================================================================    
