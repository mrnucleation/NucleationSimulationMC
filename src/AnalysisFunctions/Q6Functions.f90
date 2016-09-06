!======================================================================================    
!      This module contains the functions nessisary to initalize and 
!      perform single pair calculations.  This can be used
!      in the Umbrella Sampling algorithms or be used to study
!      di
     module Q6Functions
       use PairStorage
       use Constants
       private


       logical :: useQ6 = .false.
       integer :: q6ArrayIndx
       real(dp), parameter :: q6Constant = 4d0*pi/13d0
       real(dp) :: Q6real(0:12) , Q6Img(0:12)
       real(dp) :: dQ6real(0:12) , dQ6Img(0:12)
       real(dp) :: q6Dist


       public :: Initialize_q6, CalcQ6, useQ6, q6Dist

       contains
     !--------------------------------------------------------------------------------
       subroutine Initialize_q6
         use MiscelaniousVars
         use Coords
         use SimParameters
         implicit none 
         integer :: AllocationStatus
         integer :: startIndx, endIndx
         integer :: gloIndx1, gloIndx2
         integer :: iType, iMol, jType, jMol

         if(.not. useQ6) then
           return
         endif

         call ReserveSpace_Coord(1, startIndx, endIndx)
         q6ArrayIndx = startIndx

         do iType = 1, nMolTypes
           do iMol = 1, NMAX(iType)
             gloIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
             do jType = 1, nMolTypes
               do jMol = 1, NMAX(jType)
                 gloIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
                 rPair(gloIndx1,gloIndx2) % p % storeRParts = .true.
                 rPair(gloIndx1,gloIndx2) % p % storeRValue = .true.
               enddo
             enddo
           enddo
         enddo 


       end subroutine
     !--------------------------------------------------------------------------------
       subroutine CalcQ6
         use MiscelaniousVars
         use Coords
         use SimParameters
         implicit none 
         integer :: nNeigh, m
         integer :: iType, iMol, jType, jMol
         integer :: gloIndx1, gloIndx2
         integer :: jIndx, iIndx
         real(dp) :: rx,ry,rz,r
         real(dp) :: theta,phi,junk,sTerm,cTerm
         real(dp) :: q6Par
         real(dp) :: rPart,iPart

         nNeigh = 0
         q6par = 0E0
         Q6real = 0E0
         Q6img = 0E0	  
         if(NTotal .eq. 1) then
           q6Par = 1E0
           return	   
         endif
      
         do iType = 1, nMolTypes
           do iMol = 1, NPART(iType)
             do jType = 1, nMolTypes
               iIndx = molArray(iType)%mol(iMol)%indx
               do jMol = 1, NPART(jType)
                 jIndx = molArray(jType)%mol(jMol)%indx
                 if(iIndx .eq. jIndx) then
                   cycle
                 endif
                 gloIndx1 = molArray(iType)%mol(iMol)%globalIndx(1)
                 gloIndx2	= molArray(jType)%mol(jMol)%globalIndx(1)
                 r = rPair(gloIndx1, gloIndx2)%p%r

                 if(r .le. q6Dist) then	
                   nNeigh = nNeigh + 1
                   rx = rPair(gloIndx1, gloIndx2)%p%rx
                   ry = rPair(gloIndx1, gloIndx2)%p%ry
                   rz = rPair(gloIndx1, gloIndx2)%p%rz
                   phi = atan2(ry,rx)
                   theta = acos(rz/r)   
!                   write(2,*) rx,ry,rz
!                   write(2,*) phi, theta
                   do m = 0, 12
                     call Harmonics(theta, phi, m-6, rPart, iPart)	
                     Q6real(m) = Q6real(m) + rPart
                     Q6img(m) = Q6img(m) + iPart	
                   enddo
                 endif
               enddo
             enddo
           enddo
         enddo

         do m = 0, 12
           q6Par = q6par + Q6real(m)*Q6real(m) + Q6img(m)*Q6img(m)	  	   
         enddo
         q6par = q6par * q6Constant
         q6par = sqrt(q6par)/nNeigh
!         write(2,*) "Q6:", q6par
!         write(2,*)
 	     miscCoord(q6ArrayIndx) = q6par
      end subroutine

     !--------------------------------------------------------------------------------
      subroutine Harmonics(Theta, Phi, m, RealVal, ImgVal)
      use Constants	  
      implicit none	
      integer, intent(in) :: m
      real(dp), intent(in) :: Theta,Phi
      real(dp) :: RealVal, ImgVal, cTerm, sTerm, cTermsq, sTermsq
      real(dp), parameter :: const0 = sqrt(13d0/pi)/32d0
      real(dp), parameter :: const1 = sqrt(273d0/(2d0*pi))/16d0
      real(dp), parameter :: const2 = sqrt(1365d0/pi)/64d0
      real(dp), parameter :: const3 = sqrt(1365d0/pi)/32d0
      real(dp), parameter :: const4 = sqrt(91d0/(2d0*pi))*3d0/32d0
      real(dp), parameter :: const5 = sqrt(1001d0/pi)*3d0/32d0
      real(dp), parameter :: const6 = sqrt(3003d0/pi)/64d0
      integer :: i,n

      RealVal = 0E0
      ImgVal = 0E0
      cterm = cos(Theta)
      sterm = sin(Theta)
	  
      select case (abs(m))
      case(0)
        ctermsq = cterm * cterm
        RealVal = 231d0*ctermsq*ctermsq*ctermsq - 315d0*ctermsq*ctermsq
        RealVal = RealVal + 105d0*ctermsq - 5d0
        RealVal = RealVal * const0
      case(1) 
        ctermsq = cterm * cterm
        RealVal = 33d0*ctermsq*ctermsq*cterm
        RealVal = RealVal - 30d0*cterm*ctermsq + 5d0*cterm 
        RealVal = sterm * RealVal 
        RealVal = RealVal * const1
        if(m .gt. 0) then
	       RealVal = -RealVal
        endif
      case(2)
        ctermsq = cterm * cterm
        RealVal = sterm * sterm
        RealVal = RealVal * (33d0*ctermsq*ctermsq - 18d0*ctermsq + 1d0)
        RealVal = RealVal * const2
      case(3)	   
        RealVal = sTerm * sTerm * sTerm
        RealVal = RealVal * cterm * (11d0*cterm*cterm-3d0)
        RealVal = RealVal * const3
        if(m .gt. 0) then
	      RealVal = -RealVal
        endif
      case(4)
        sTermsq = sterm * sterm
        RealVal = sTermsq * sTermsq * (11d0*cterm*cterm - 1d0)
        RealVal = RealVal * const4
      case(5)
        sTermsq = sTerm * sTerm
        RealVal = sTerm * sTermsq * sTermsq * cterm
        RealVal = RealVal * const5
        if(m .gt. 0) then
          RealVal = -RealVal
        endif
      case(6)
        sTermsq = sTerm * sTerm
        RealVal = sTermsq * sTermsq *sTermsq
        RealVal = RealVal * const6
      end select

      ImgVal = RealVal * sin(m * Phi)
      RealVal = RealVal * cos(m * Phi)
	
	
      If(Abs(ImgVal) .lt. 1d-13) then
        ImgVal = 0d0
      endif

      If(Abs(RealVal) .lt. 1d-13) then
        RealVal = 0d0
      endif

      end subroutine
     !--------------------------------------------------------------------------------
      end module
!====================================================================================== 



