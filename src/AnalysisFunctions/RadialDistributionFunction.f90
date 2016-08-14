

!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform radial distribution studies
      module RadialDist
      use PairStorage

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
