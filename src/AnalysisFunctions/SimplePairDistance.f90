!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies
      module SimpleDistPair
        use PairStorage

        integer :: nDistPair
        integer, allocatable :: pairArrayIndx(:)
        integer, allocatable :: pairGloIndx1(:), pairGloIndx2(:)
        integer, allocatable :: molIndx1(:), molIndx2(:)  

        contains
     !--------------------------------------------------------------------------------
        subroutine Initialize_DistPair
           use MiscelaniousVars
           implicit none 
           integer :: AllocationStatus
           integer :: startIndx, endIndx, iPair

           allocate(pairArrayIndx(1:nDistPair), STAT = AllocationStatus)
           allocate(pairGloIndx1(1:nDistPair), STAT = AllocationStatus)
           allocate(pairGloIndx2(1:nDistPair), STAT = AllocationStatus)
           allocate(molIndx1(1:nDistPair), STAT = AllocationStatus)
           allocate(molIndx2(1:nDistPair), STAT = AllocationStatus)
           
           call ReserveSpace_Coord(nDistPair, startIndx, endIndx)

           do iPair = 1, nDistPair
             pairArrayIndx(iPair) = startIndx + iPair - 1
           enddo
        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs
          use MiscelaniousVars
          use Coords
          implicit none 
          integer :: iDistPair
          integer :: gloIndx1, gloIndx2
          real(dp) :: r, r_sq

          do iDistPair =1, nDistPair
            gloIndx1 = pairGloIndx1(iDistPair)
            gloIndx2 = pairGloIndx2(iDistPair)
            r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq           
            r = dsqrt(r_sq)
            miscCoord(pairArrayIndx(iDistPair)) = r
          enddo


        end subroutine
     !--------------------------------------------------------------------------------
        subroutine CalcDistPairs_New(disp)
          use MiscelaniousVars
          use Coords
          implicit none 
          type(Displacement), intent(in) :: disp(:)
!          logical :: 
          integer :: sizeDisp, iDistPair, iDisp
          integer :: type1, mol1, atom1
          integer :: dispIndx, indx1, indx2
          real(dp) :: rx, ry, rz, r

          sizeDisp = size(disp)
          dispIndx = disp(1)%molIndx

          do iDistPair = 1, nDistPair
!            if(molIndx1(iDistPair) .eq. 

!            rx = MolArray(type1)%mol(mol1)%x(atom1) - MolArray(type2)%mol(mol2)%x(atom2)
!            ry = MolArray(type1)%mol(mol1)%y(atom1) - MolArray(type2)%mol(mol2)%y(atom2)
!            rz = MolArray(type1)%mol(mol1)%z(atom1) - MolArray(type2)%mol(mol2)%z(atom2)
            r = rx*rx + ry*ry + rz*rz
            r = dsqrt(r)
            miscCoord_New(pairArrayIndx(iDistPair)) = r            
          enddo


        end subroutine

     !--------------------------------------------------------------------------------
      end module
!====================================================================================== 


