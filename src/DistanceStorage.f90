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
      integer :: arrayIndx
      integer :: indx1, indx2
      real(dp) :: r_sq
      real(dp) :: E_Pair
    end type

     !Variable defintion for storing trial distances
    type DistArrayNew
      integer :: oldIndx
      integer :: nType2, nMol2, nAtom2
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
!=====================================================================
      subroutine CreateDistArrays
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
      allocate(newDist(1:nPairs), stat = AllocationStat) 

      cnt = 0
      do i = 1, nTotalAtoms-1
        do j = i+1, nTotalAtoms
          cnt = cnt + 1
          distStorage(cnt)%arrayIndx = cnt
          distStorage(cnt)%indx1 = i
          distStorage(cnt)%indx2 = j
          distStorage(cnt)%r_sq = 0d0
          distStorage(cnt)%E_Pair = 0d0
          rPair(i,j)%p => distStorage(cnt)
          rPair(j,i)%p => distStorage(cnt)
        enddo
      enddo

      distStorage(cnt)%arrayIndx = 0
      distStorage(0)%indx1 = 0
      distStorage(0)%indx2 = 0
      distStorage(0)%r_sq = 0d0
      distStorage(0)%E_Pair = 0d0
      do i = 1, nTotalAtoms
        rPair(i,i)%p => distStorage(0)
      enddo 

     end subroutine
!=====================================================================
     subroutine CalcAllDistPairs
      use Coords
      use ForceField
      use SimParameters, only: NMAX, NPART, nMolTypes, maxAtoms
      implicit none
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: globIndx1, globIndx2, jMolMin
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij   
   
      do iType = 1,nMolTypes
        do jType = iType, nMolTypes
          do iMol=1,NPART(iType)
            do jMol = jMolMin,NPART(jType)
              do iAtom = 1,nAtoms(iType)
                atmType1 = atomArray(iType,iAtom)
                globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
                do jAtom = 1,nAtoms(jType)       
                  globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                  atmType2 = atomArray(jType, jAtom) 
                  rmin_ij = r_min_tab(atmType1,atmType2)          
                  rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom) 
                  r = rx*rx + ry*ry + rz*rz
                  if(r .lt. rmin_ij) then
                    if(iMol .ne. jMol) then
                      stop "ERROR! Overlaping atoms found in the current configuration!"
                    endif
                  endif 
                  rPair(globIndx1, globIndx2)%p%r_sq = r
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo


     end subroutine
!=====================================================================
     subroutine CalcNewDistPairs(disp, rejMove)
      use Coords
      use ForceField
      use SimParameters, only: NMAX, NPART, nMolTypes, maxAtoms
      implicit none
      type(Displacement), intent(in) :: disp(:)
      logical, intent(out) :: rejMove
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iDisp, sizeDisp, gloIndx1, gloIndx2, oldIndx
      integer :: jMolMin
      real(dp) :: rx,ry,rz,r
      real(dp) :: rmin_ij   

   
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      sizeDisp = size(disp)
      

      nNewDist = 0
      do iDisp = 1, sizeDisp
        iAtom = disp(iDisp)%atmIndx
        atmType1 = atomArray(iType, iAtom)
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            rmin_ij = r_min_tab(atmType1,atmType2)     
            do jMol = 1, NPART(jType)
              gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
              if(gloIndx1 .eq. gloIndx2 ) then
                cycle
              endif
!               Distance for the New position
              rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
!             If r_new is less than r_min reject the move.              
              if(r .lt. rmin_ij) then
                if(jMol .ne. iMol) then
                  rejMove = .true.
                  return
                endif
              endif    
              nNewDist = nNewDist + 1
              oldIndx = rPair(gloIndx1, gloIndx2)%p%arrayIndx
              newDist(nNewDist)%oldIndx = oldIndx
              newDist(nNewDist)%nType2 = jType
              newDist(nNewDist)%nMol2 = jMol
              newDist(nNewDist)%nAtom2 = jAtom
              newDist(nNewDist)%r_sq = r
            enddo
          enddo
        enddo
      enddo


     end subroutine
!=====================================================================
    end module
!=====================================================================


