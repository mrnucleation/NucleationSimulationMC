

!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial density distribution studies
      module RadialDensity
      use PairStorage
      use MiscelaniousVars
      
      private
      integer :: nRadialDens
      integer, allocatable :: radHistIndx(:)
      integer, allocatable :: radType(:)

      public :: nRadialDens
      public :: Initialize_RadialDens
      public :: SetDensityParameters
      public :: Calc_RadialDensity
      public :: Output_RadialDensity

      contains

!======================================================================================    
      subroutine Initialize_RadialDens
      implicit none 
      integer :: iRadial
      integer :: AllocationStatus
      integer :: startIndx, endIndx

      if(nRadialDens .eq. 0) then
        return
      endif

      allocate( radType(1:nRadialDens), stat = AllocationStatus )

      call ReserveSpace_Histograms(nRadialDens, startIndx, endIndx)

      do iRadial = 1, nRadialDens
        radHistIndx(iRadial) = startIndx + iRadial - 1
      enddo

      end subroutine
!======================================================================================    
      subroutine SetDensityParameters(iRadial, type1)
      implicit none 
      integer, intent(in) :: iRadial, type1

      radType(iRadial) = type1

      end subroutine
!======================================================================================    
      subroutine Calc_RadialDensity
        use ParallelVar, only: nout
        use SimParameters, only: NPART, nMolTypes
        use Coords
        implicit none
        integer :: iRadial
        integer :: bin, nBins
        integer :: iType, nType
        integer :: iMol
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        integer :: x_cnt, y_cnt, z_cnt
        real(dp) :: r_sq, r, rx, ry, rz
        real(dp) :: xcm, ycm, zcm

        xcm = 0d0 
        x_cnt = 0
        ycm = 0d0 
        y_cnt = 0
        zcm = 0d0 
        z_cnt = 0
        do iType = 1, nMolTypes
          do iMol = 1, NPART(iType)
            xcm = xcm + molArray(iType)%mol(iMol)%x(1)
            ycm = ycm + molArray(iType)%mol(iMol)%y(1)
            zcm = zcm + molArray(iType)%mol(iMol)%z(1) 
            x_cnt = x_cnt + 1
            y_cnt = y_cnt + 1
            z_cnt = z_cnt + 1
          enddo
        enddo

        xcm = xcm/real(x_cnt, dp)
        ycm = ycm/real(y_cnt, dp)
        zcm = zcm/real(z_cnt, dp)

        do iRadial = 1, nRadialDens
          nType = radType(iRadial)
          radialIndx = radHistIndx(iRadial)
          do iMol = 1, NPART(nType)
            rx = molArray(nType)%mol(iMol)%x(1) - xcm
            ry = molArray(nType)%mol(iMol)%y(1) - ycm
            rz = molArray(nType)%mol(iMol)%z(1) - zcm
            r_sq = rx*rx + ry*ry + rz*rz
            r = sqrt(r_sq)
            bin = floor(r * miscHist(radialIndx)%sizeInv)
            if(bin .le. miscHist(radialIndx)%nBins) then
              miscHist(radialIndx)%binCount(bin) = miscHist(radialIndx)%binCount(bin) + 1d0
            else 
              nBins = miscHist(radialIndx)%nBins
              miscHist(radialIndx)%binCount(nBins+1) = miscHist(radialIndx)%binCount(nBins+1) + 1d0
            endif            
          enddo
        enddo
      
      end subroutine  
!======================================================================================    
      subroutine Output_RadialDensity
        use Constants
        implicit none
        integer :: iRadial, iBin
        integer :: radialIndx
        integer :: gloIndx1, gloIndx2
        real(dp) :: d_bin
        real(dp) :: r, norm, rFactor

        do iRadial = 1, nRadialDens
          open(unit = 80, file = miscHist(iRadial)%fileName)
          d_bin = miscHist(iRadial)%binSize
          norm = 0E0
          do iBin = 0, miscHist(iRadial)%nBins
            norm = norm + miscHist(iRadial)%binCount(iBin)
          enddo
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
