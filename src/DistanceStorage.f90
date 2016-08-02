!=====================================================================
    module PairStorage
    use VarPrecision

    type DistArray
      real(dp) :: r_sq
      real(dp) :: E_Pair
    end type

    type DistArrayNew
      integer :: indx1, indx2
      real(dp) :: r_sq
      real(dp) :: E_Pair
    end type

    integer, allocatable :: atomSumList(:)
    type(DistArray), allocatable :: distStorage(:,:)
    type(DistArrayNew), allocatable :: distStorage_new(:)

    contains
    !----------------------------------------------------------
     integer function getDistIndex(iType, iMol, iAtom)
     use ForceField, only: nAtoms
     use SimParameters, only: NMAX
     implicit none
     integer, intent(in) :: iType, iMol, iAtom
     integer :: indx

     indx = atomSumList(iType)
     indx = indx + nAtoms(iType)*iMol + iAtom

     end function
    !----------------------------------------------------------

    end module
!=====================================================================
