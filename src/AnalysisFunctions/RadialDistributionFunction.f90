

!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies
      module RadialDistribution
      use PairStorage
      use MiscelaniousVars
      
!      private
      integer :: nRadialDist
      integer, allocatable :: radHistIndx(:)
      integer, allocatable :: radType1(:), radAtom1(:)
      integer, allocatable :: radType2(:), radAtom2(:)

      public :: Initialize_RadialDist
      public :: Calc_RadialDist
      public :: Output_RadialDist

      contains

!======================================================================================    
      subroutine Initialize_RadialDist
      implicit none 
      integer :: iRadial
      integer :: AllocationStatus
      integer :: startIndx, endIndx

      allocate( radType1(1:nRadialDist), stat = AllocationStatus )
      allocate( radType2(1:nRadialDist), stat = AllocationStatus )
      allocate( radAtom1(1:nRadialDist), stat = AllocationStatus )
      allocate( radAtom2(1:nRadialDist), stat = AllocationStatus )
      allocate( radHistIndx(1:nRadialDist), stat = AllocationStatus )

      call ReserveSpace_Histograms(nRadialDist, startIndx, endIndx)

      do iRadial = 1, nRadialDist
        radHistIndx(iRadial) = startIndx + iRadial - 1
!        write(*,*) iRadial, radHistIndx(iRadial)
      enddo

      end subroutine
!======================================================================================    
      subroutine Calc_RadialDist
        use SimParameters, only: NPART
        use Coords
        implicit none
        integer :: iRadial
        integer :: bin, nBins
        integer :: nType1, nType2, nAtom1, nAtom2
        integer :: jMol, iMol
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        real(dp) :: r_sq, r

        do iRadial = 1, nRadialDist
          nType1 = radType1(iRadial)
          nType2 = radType2(iRadial)
          nAtom1 = radAtom1(iRadial) 
          nAtom2 = radAtom2(iRadial) 
          radialIndx = radHistIndx(iRadial)
!          write(*,*) iRadial, radialIndx
          if(nType1 .eq. nType2) then
            do iMol = 1, NPART(nType1) - 1
              do jMol = iMol+1, NPART(nType1)
                gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
                gloIndx2 = molArray(nType2)%mol(jMol)%globalIndx(nAtom2)
                r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
                r = sqrt(r_sq)
                bin = floor(r * miscHist(radialIndx)%sizeInv)
                nBins = miscHist(radialIndx)%nBins
!                write(*,*) radialIndx, r, r_sq, bin
                if(bin .le. nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  miscHist(radialIndx)%binCount(nBins+1) = miscHist(radialIndx)%binCount(nBins+1) + 2d0
                endif
              enddo
            enddo
          else
            do iMol = 1, NPART(nType1)
              do jMol = 1, NPART(nType2)
                gloIndx1 = molArray(nType1)%mol(iMol)%globalIndx(nAtom1)
                gloIndx2 = molArray(nType2)%mol(jMol)%globalIndx(nAtom2)
                r_sq = rPair(gloIndx1, gloIndx2)%p%r_sq
                r = sqrt(r_sq)
                bin = floor(r * miscHist(radialIndx)%sizeInv)
                nBins = miscHist(radialIndx)%nBins
                if(bin .le. nBins) then
                  miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 2d0
                else 
                  miscHist(radialIndx)%binCount(nBins+1) = miscHist(radialIndx)%binCount(nBins+1) + 2d0
                endif
              enddo
            enddo
          endif
        enddo
        
      
      end subroutine  
!======================================================================================    
      subroutine Output_RadialDist
        use Constants
        implicit none
        integer :: iRadial, iBin
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        real(dp) :: d_bin
        real(dp) :: r, norm, rFactor

        do iRadial = 1, nRadialDist
          open(unit = 80, file = miscHist(iRadial)%fileName)
          d_bin = miscHist(iRadial)%binSize
          norm = sum(miscHist(iRadial)%binCount)
!          write(*,*) norm
          write(80,*) "Number of Bins:", miscHist(iRadial)%nBins
          write(80,*) "Bin Size:", d_bin
          write(80,*)
          do iBin = 0, miscHist(iRadial)%nBins
            r = iBin * d_bin
            rFactor = norm * 4E0/3E0 * pi * ( (r+d_bin)**3 - r**3 )
            write(80,*) r, miscHist(iRadial)%binCount(iBin)/rFactor
          enddo
          close(80)
        enddo
        
      
      end subroutine 
!======================================================================================    
      end module
!======================================================================================    
