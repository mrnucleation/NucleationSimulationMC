!=====================================================================
!   This module contains data arrays whose purpose is to save the distance
!   pairs during the course of the simulation. The primary storage variable
!   is allocated as a 1D array to minimize the amount of memory that needs
!   to be used to describe each distance pair.  To facilitate easy access
!   a 2D array of pointers (rPair) is used so that the programer can quickly
!   access and modify the distance storage using matrix notation.

    module PairStorage
    use VarPrecision

     !Primary distance/energy pair storage variable defintinio
    type DistArray
      integer :: indx1, indx2
      real(dp) :: r_sq
      real(dp) :: E_Pair
    end type

     !Variable defintion for storing trial distances
    type DistArrayNew
      integer :: indx1, indx2
      real(dp) :: r_sq
      real(dp) :: E_Pair
    end type

     !Variable used to create a 2D array of pointers to make
     !accessing the distance storage easier to work with
    type DistPointer
      type(DistArray), pointer :: p      
    end type

    integer :: nPairs, nTotalAtoms, nNewDist
    type(DistArray), allocatable, target :: distStorage(:)
    type(DistPointer), allocatable :: rPair(:,:)
    type(DistArrayNew), allocatable :: newDist(:)

    contains
    !----------------------------------------------------------
     subroutine createDistArrays
      use ForceField, only: nAtoms
      use SimParameters, only: NMAX, nMolTypes, maxAtoms
      implicit none
      integer :: iType, AllocationStat
      integer :: i, j, cnt

      nTotalAtoms = 0
      do iType = 1, nMolTypes
        nTotalAtoms = nTotalAtoms + NMAX(iType)*nAtoms(iType)
      enddo

      nPairs = nint( dble(nTotalAtoms * (nTotalAtoms - 1)) / 2d0 )
      allocate(distStorage(0:nPairs), stat = AllocationStat)
      allocate(rPair(1:nTotalAtoms, 1:nTotalAtoms), stat = AllocationStat)
      allocate(newDist(1:nTotalAtoms), stat = AllocationStat) 

      cnt = 0
      do i = 1, nTotalAtoms-1
        do j = i+1, nTotalAtoms
          cnt = cnt + 1
          distStorage(cnt)%indx1 = i
          distStorage(cnt)%indx2 = j
          distStorage(cnt)%r_sq = 0d0
          distStorage(cnt)%E_Pair = 0d0
          rPair(i,j)%p => distStorage(cnt)
          rPair(j,i)%p => distStorage(cnt)
        enddo
      enddo
     
      distStorage(0)%indx1 = 0
      distStorage(0)%indx2 = 0
      distStorage(0)%r_sq = 0d0
      distStorage(0)%E_Pair = 0d0
      do i = 1, nTotalAtoms
        rPair(i,i)%p => distStorage(0)
      enddo 

     end subroutine
    !----------------------------------------------------------
     subroutine createDistArrays
      use ForceField, only: nAtoms
      use SimParameters, only: NMAX, nMolTypes, maxAtoms
      implicit none
      integer :: iType, AllocationStat
      integer :: i, j, cnt

      nTotalAtoms = 0
      do iType = 1, nMolTypes
        nTotalAtoms = nTotalAtoms + NMAX(iType)*nAtoms(iType)
      enddo

      nPairs = nint( dble(nTotalAtoms * (nTotalAtoms - 1)) / 2d0 )
      allocate(distStorage(0:nPairs), stat = AllocationStat)
      allocate(rPair(1:nTotalAtoms, 1:nTotalAtoms), stat = AllocationStat)
      allocate(newDist(1:nTotalAtoms), stat = AllocationStat) 

      cnt = 0
      do i = 1, nTotalAtoms-1
        do j = i+1, nTotalAtoms
          cnt = cnt + 1
          distStorage(cnt)%indx1 = i
          distStorage(cnt)%indx2 = j
          distStorage(cnt)%r_sq = 0d0
          distStorage(cnt)%E_Pair = 0d0
          rPair(i,j)%p => distStorage(cnt)
          rPair(j,i)%p => distStorage(cnt)
        enddo
      enddo
     
      distStorage(0)%indx1 = 0
      distStorage(0)%indx2 = 0
      distStorage(0)%r_sq = 0d0
      distStorage(0)%E_Pair = 0d0
      do i = 1, nTotalAtoms
        rPair(i,i)%p => distStorage(0)
      enddo 

     end subroutine
    !----------------------------------------------------------
     subroutine calculateAllDistPairs
      use Coords
      use ForceField, only: nAtoms
      use SimParameters, only: NMAX, nMolTypes, maxAtoms
      implicit none
      integer :: iType, AllocationStat
      integer :: i, j, cnt


     end subroutine
    !----------------------------------------------------------

    end module
!=====================================================================


