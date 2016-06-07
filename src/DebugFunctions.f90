!=======================================================
      subroutine DEBUG_Output_NeighborList
      use SimParameters
      use Coords
      implicit none
      integer i,j  
      write(35,*) "-------------------------------------------"
      write(35,*) "NeighborList:"
      write(35,*) "NPART=", NPART
      do i=1,maxMol
        write(35,*) (NeighborList(i,j), j=1,maxMol)
      enddo

        
      end subroutine
!=======================================================
      subroutine DEBUG_Output_NewConfig
      use VarPrecision
      use SimParameters
      use Coords
      use Forcefield
      implicit none
      integer :: iType, iMol, iAtom, atmType
      write(35,*) "-------------------------------------------"
      write(35,*) "Trial Configuration:"
      write(35,*) "NPART=", NPART
      write(35,*)
      do iType = 1, nMolTypes
        do iMol = 1, NPart(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(35,*) atomData(atmType)%Symb, &
                        MolArray(iType)%mol(iMol)%x(iAtom), &
                        MolArray(iType)%mol(iMol)%y(iAtom), &
                        MolArray(iType)%mol(iMol)%z(iAtom)
          enddo
        enddo
      enddo

      iType = newMol%molType
      do iAtom = 1, nAtoms(iType)
        atmType = atomArray(iType,iAtom)
        write(35,*) atomData(atmType)%Symb, newMol%x(iAtom), newMol%y(iAtom), newMol%z(iAtom)
      enddo

        
      end subroutine
!==================================================================
