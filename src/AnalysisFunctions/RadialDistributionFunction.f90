

!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies
      module RadialDistribution
      use PairStorage
      use MiscelaniousVars
      
 
      integer :: nRadialDist
      integer, allocatable :: radHistIndx(:)
      integer, allocatable :: type1(:), atom1(:)
      integer, allocatable :: type2(:), atom2(:)
      contains

!======================================================================================    
      subroutine Initialize_RadialDist
      implicit none 
      integer :: iRadial
      integer :: AllocationStatus
      integer :: startIndx, endIndx

      allocate( type1(1:nRadialDist), stat = AllocationStatus )
      allocate( type2(1:nRadialDist), stat = AllocationStatus )
      allocate( atom1(1:nRadialDist), stat = AllocationStatus )
      allocate( atom2(1:nRadialDist), stat = AllocationStatus )
      allocate( radHistIndx(1:nRadialDist), stat = AllocationStatus )

      call ReserveSpace_Histograms(nRadialDist, startIndx, endIndx)

      do iRadial = 1, nRadialDist
        radHistIndx(iRadial) = startIndx + iRadial - 1
      enddo

      end subroutine
!======================================================================================    
      subroutine Calc_RadialDist
        use SimParameters, only: NPART
        use Coords
        implicit none
        integer :: iRadial
        integer :: bin
        integer :: nType1, nType2, nAtom1, nAtom2
        integer :: jMol, iMol
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        real(dp) :: r_sq, r

        do iRadial = 1, nRadialDist
          nType1 = type1(iRadial)
          nType2 = type2(iRadial)
          nAtom1 = atom1(iRadial) 
          nAtom2 = atom2(iRadial) 
          radialIndx = radHistIndx(iRadial)
          if(nType1 .eq. nType2) then
            do iMol = 1, NPART(nType1) - 1
              do jMol = iMol+1, NPART(nType1)
                gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
                gloIndx2 = molArray(nType2)%mol(iMol)%globalIndx(nAtom2)
                r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
                r = sqrt(r_sq)
                bin = floor(miscHist(radialIndx)%binSize * r)
                if(bin .le. miscHist(radialIndx)%nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                endif
              enddo
            enddo
          else
            do iMol = 1, NPART(nType1)
              do jMol = 1, NPART(nType2)
                gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
                gloIndx2 = molArray(nType2)%mol(iMol)%globalIndx(nAtom2)
                r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
                r = sqrt(r_sq)
                bin = floor(miscHist(radialIndx)%binSize * r)
                if(bin .le. miscHist(radialIndx)%nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                endif
              enddo
            enddo
          endif
        enddo
        
      
      end subroutine  
!======================================================================================    
      subroutine Calc_RadialDist
        use SimParameters, only: NPART
        use Coords
        implicit none
        integer :: iRadial
        integer :: bin
        integer :: nType1, nType2, nAtom1, nAtom2
        integer :: jMol, iMol
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        real(dp) :: r_sq, r

        do iRadial = 1, nRadialDist
          nType1 = type1(iRadial)
          nType2 = type2(iRadial)
          nAtom1 = atom1(iRadial) 
          nAtom2 = atom2(iRadial) 
          radialIndx = radHistIndx(iRadial)
          if(nType1 .eq. nType2) then
            do iMol = 1, NPART(nType1) - 1
              do jMol = iMol+1, NPART(nType1)
                gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
                gloIndx2 = molArray(nType2)%mol(iMol)%globalIndx(nAtom2)
                r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
                r = sqrt(r_sq)
                bin = floor(miscHist(radialIndx)%binSize * r)
                if(bin .le. miscHist(radialIndx)%nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                endif
              enddo
            enddo
          else
            do iMol = 1, NPART(nType1)
              do jMol = 1, NPART(nType2)
                gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
                gloIndx2 = molArray(nType2)%mol(iMol)%globalIndx(nAtom2)
                r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
                r = sqrt(r_sq)
                bin = floor(miscHist(radialIndx)%binSize * r)
                if(bin .le. miscHist(radialIndx)%nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                endif
              enddo
            enddo
          endif
        enddo
        
      
      end subroutine
!======================================================================================    
      subroutine Output_RadialDist
        implicit none
        integer :: iRadial, iBin
        integer :: bin
        integer :: nType1, nType2, nAtom1, nAtom2
        integer :: jMol, iMol
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        real(dp) :: r_sq, r

        do iRadial = 1, nRadialDist

        enddo
        
      
      end subroutine 
!======================================================================================    
      end module
!======================================================================================    
