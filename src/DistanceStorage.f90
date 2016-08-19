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
      logical :: storeRValue
      integer :: arrayIndx
      integer :: indx1, indx2
      real(dp) :: r_sq, r
      real(dp) :: E_Pair
    end type

     !Variable defintion for storing trial distances
    type DistArrayNew
      integer :: indx1, indx2
      real(dp) :: r_sq, r
      real(dp) :: E_Pair
    end type

     !Variable used to create a 2D array of pointers to make
     !accessing the distance storage easier to work with
    type DistPointer
      type(DistArray), pointer :: p      
    end type

    integer :: nMaxPairs, nTotalAtoms, nNewDist
    integer, allocatable, target :: oldIndxArray(:)
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

      nMaxPairs = nint( dble(nTotalAtoms * (nTotalAtoms - 1)) / 2d0 )
      allocate(distStorage(0:nMaxPairs), stat = AllocationStat)
      allocate(rPair(1:nTotalAtoms, 1:nTotalAtoms), stat = AllocationStat)
      allocate(newDist(1:nMaxPairs), stat = AllocationStat) 
      allocate(oldIndxArray(1:nMaxPairs), stat = AllocationStat) 

      cnt = 0
      do i = 1, nTotalAtoms-1
        do j = i+1, nTotalAtoms
          cnt = cnt + 1
          distStorage(cnt)%arrayIndx = cnt
          distStorage(cnt)%indx1 = i
          distStorage(cnt)%indx2 = j
          distStorage(cnt)%r_sq = 0d0
          distStorage(cnt)%r = 0d0
          distStorage(cnt)%E_Pair = 0d0
          distStorage(cnt)%storeRValue = .false.
          rPair(i,j)%p => distStorage(cnt)
          rPair(j,i)%p => distStorage(cnt)
        enddo
      enddo

!     To avoid problems with unassocaited pointers, the diagonal of the pair array pointer
!     is set to reference the 0-th value of the array.  This will hopefully make errors
!     with 
      distStorage(0)%arrayIndx = 0
      distStorage(0)%indx1 = 0
      distStorage(0)%indx2 = 0
      distStorage(0)%r_sq = 0d0
      distStorage(0)%E_Pair = 0d0
      do i = 1, nTotalAtoms
        rPair(i,i)%p => distStorage(0)
      enddo 

     end subroutine
!=====================================================================
     subroutine SetStorageFlags(q_tab)
     use Coords
     use ForceField
     use SimParameters, only: NMAX, NPART, nMolTypes, maxAtoms
     implicit none
     real(dp), intent(in) :: q_tab(:,:)
     integer :: iType,jType,iMol,jMol,iAtom,jAtom
     integer(kind=atomIntType) :: atmType1, atmType2      
     integer :: globIndx1, globIndx2 
     real(dp) :: q_ij

     do iType = 1,nMolTypes
       do jType = iType, nMolTypes
         do iMol=1,NMAX(iType)
           do jMol = 1,NMAX(jType)
             do iAtom = 1,nAtoms(iType)
               atmType1 = atomArray(iType,iAtom)
               globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
               do jAtom = 1,nAtoms(jType)       
                 globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                 atmType2 = atomArray(jType, jAtom)
                 q_ij = q_tab(atmType1, atmType2)
                 if(q_ij .ne. 0E0) then
                   rPair(globIndx1, globIndx2) % p % storeRValue = .true.
                 else
                   rPair(globIndx1, globIndx2) % p % storeRValue = .false.
                 endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

     end subroutine
!=====================================================================
     subroutine CalcAllDistPairs
      use Coords
      use ForceField
      use SimParameters, only: NMAX, NPART, nMolTypes, maxAtoms
      implicit none
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer :: iIndx, jIndx
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: globIndx1, globIndx2 
      real(dp) :: rx, ry, rz, r_sq, q_ij
      real(dp) :: rmin_ij   

      do iType = 1,nMolTypes
        do jType = iType, nMolTypes
          do iMol=1,NPART(iType)
            iIndx = molArray(iType)%mol(iMol)%indx
            do jMol = 1,NPART(jType)
              jIndx = molArray(iType)%mol(iMol)%indx
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
                  r_sq = rx*rx + ry*ry + rz*rz
                  if(r_sq .lt. rmin_ij) then
                    if(iIndx .ne. jIndx) then
                      stop "ERROR! Overlaping atoms found in the current configuration!"
                    endif
                  endif 
                  rPair(globIndx1, globIndx2) % p % r_sq = r_sq
                  rPair(globIndx1, globIndx2) % p % r = sqrt(r_sq)
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
      integer :: iIndx, jIndx
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iDisp, sizeDisp, gloIndx1, gloIndx2, oldIndx
      integer :: jMolMin
      real(dp) :: rx,ry,rz,r_sq
      real(dp) :: rmin_ij   

   
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = molArray(iType)%mol(iMol)%indx
      sizeDisp = size(disp)
      rejMove = .false.
      nNewDist = 0

      do iDisp = 1, sizeDisp
        iAtom = disp(iDisp)%atmIndx
        atmType1 = atomArray(iType, iAtom)
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1, nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            rmin_ij = r_min_tab(atmType1,atmType2)     
            do jMol = 1, NPART(jType)
              gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
              if(gloIndx1 .eq. gloIndx2 ) then
                cycle
              endif
              jIndx = molArray(jType)%mol(jMol)%indx
!               Distance for the New position
              rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz
!             If r_new is less than r_min reject the move.              
              if(r_sq .lt. rmin_ij) then
                if(iIndx .ne. jIndx) then
                  rejMove = .true.
                  return
                endif
              endif    
              nNewDist = nNewDist + 1
              oldIndxArray(nNewDist) = rPair(gloIndx1, gloIndx2)%p%arrayIndx
              newDist(nNewDist)%indx1 = gloIndx1
              newDist(nNewDist)%indx2 = gloIndx2
              newDist(nNewDist)%r_sq = r_sq
!              newDist(nNewDist)%E_Pair = 0d0
              if( rPair(gloIndx1, gloIndx2)%p%storeRValue ) then
                newDist(nNewDist)%r = sqrt(r_sq)
              endif
            enddo
          enddo
        enddo
      enddo


     end subroutine
!=====================================================================
     subroutine CalcSwapInDistPairs(rejMove)
      use Coords
      use ForceField
      use SimParameters, only: NMAX, NPART, nMolTypes, maxAtoms
      implicit none
      logical, intent(out) :: rejMove
      integer :: iType,jType,iMol,jMol,iAtom,jAtom
      integer(kind=atomIntType) :: atmType1,atmType2      
      integer :: iDisp, gloIndx1, gloIndx2, oldIndx
      integer :: jMolMin
      real(dp) :: rx,ry,rz,r_sq
      real(dp) :: rmin_ij   

   
      iType = newMol%molType
      iMol = NPART(iType)+1
      rejMove = .false.
      nNewDist = 0
      do iAtom = 1, nAtoms(iType)
        atmType1 = atomArray(iType, iAtom)
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1, nAtoms(jType)        
            atmType2 = atomArray(jType,jAtom)
            rmin_ij = r_min_tab(atmType1, atmType2)
            do jMol = 1, NPART(jType)
              gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
              rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r_sq = rx*rx + ry*ry + rz*rz
              if(r_sq .lt. rmin_ij) then
                rejMove = .true.
!                write(*,*) iType, iMol, iAtom
!                write(*,*) jType, jMol, jAtom
!                write(*,*) sqrt(r_sq), sqrt(rmin_ij)
                return
              endif
              nNewDist = nNewDist + 1
              oldIndxArray(nNewDist) = rPair(gloIndx1, gloIndx2) % p % arrayIndx
              newDist(nNewDist)%indx1 = gloIndx1
              newDist(nNewDist)%indx2 = gloIndx2
              newDist(nNewDist)%r_sq = r_sq
!              newDist(nNewDist)%E_Pair = 0d0
              if( rPair(gloIndx1, gloIndx2)%p%storeRValue ) then
                newDist(nNewDist)%r = sqrt(r_sq)
              endif
            enddo
          enddo
        enddo
      enddo


     end subroutine
!=====================================================================
     subroutine UpdateDistArray
      implicit none
      integer :: iPair
!      integer :: oldIndx
!      integer :: oldIndx


      do iPair = 1, nNewDist
        distStorage(oldIndxArray(iPair))%r_sq = newDist(iPair)%r_sq
        distStorage(oldIndxArray(iPair))%E_Pair = newDist(iPair)%E_Pair
        if( distStorage(oldIndxArray(iPair))%storeRValue ) then
          distStorage(oldIndxArray(iPair))%r = newDist(iPair)%r
        endif
      enddo


     end subroutine
!=====================================================================
     subroutine UpdateDistArray_SwapOut(nType, nMol)
      use Coords
      use Forcefield, only: nAtoms
      use SimParameters, only: NMAX, NPART, nMolTypes, maxAtoms
      implicit none
      integer, intent(in) :: nType, nMol
      integer :: iAtom, nMol2
      integer :: gloIndx1, gloIndx2, gloIndx3


      nMol2 = NPART(nType)
      if(nMol .eq. nMol2) then
        return
      endif


      do iAtom = 1, nAtoms(nType)
        gloIndx1 = molArray(nType)%mol(nMol)%globalIndx(iAtom)
        gloIndx2 = molArray(nType)%mol(nMol2)%globalIndx(iAtom)
        do gloIndx3 = 1, nTotalAtoms
          rPair(gloIndx1, gloIndx3)%p%r_sq = rPair(gloIndx2, gloIndx3)%p%r_sq
          rPair(gloIndx1, gloIndx3)%p%E_Pair = rPair(gloIndx2, gloIndx3)%p%E_Pair
          rPair(gloIndx1, gloIndx3)%p%storeRValue = rPair(gloIndx2, gloIndx3)%p%storeRValue
          if( rPair(gloIndx2, gloIndx3)%p%storeRValue ) then
            rPair(gloIndx1, gloIndx3)%p%r = rPair(gloIndx2, gloIndx3)%p%r
          endif
        enddo
      enddo


     end subroutine
!=====================================================================
    end module
!=====================================================================


