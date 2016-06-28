!==============================================================
!  This file contains a set of small functions that are either used in creating forcefield
!  tables such as mixing rules or contain functional forms that are used in the energy calculations
!  themselves. 
!==============================================================
      module ForceFieldVariableType
      use VarPrecision
!     This block defines custom variable types for forcefield
      type AtomDef
        character(len=4) :: Symb
        real(dp) :: sig, ep, q, mass
      end type
        
      type BondDef
        real(dp) :: k_eq,r_eq
      end type

      integer, parameter :: nBendHistBins = 200    
      real(dp), parameter :: bendBinWidth = 4d0*datan(1d0)/dble(nBendHistBins)    
      real(dp), parameter :: startProb = 0.05d0
      type BendAngleDef
        real(dp) :: k_eq, ang_eq

        integer :: startBin
        real(dp) :: accptConstant
        real(dp) :: Prob(1:nBendHistBins)
!        real(dp) :: Hist(0:nBendHistBins)
      end type
        
      type TorsionAngleDef
        integer(kind=2) :: nPara
        real(dp),allocatable :: a(:)
      end type

      type NonBondedIndex
!        integer :: nonMembr(1:2)
        integer(kind=2) :: nonMembr(1:2)        
      end type
      
      type BondIndex
        integer(kind=2) :: bondType
        integer(kind=2) :: bondMembr(1:2)
!        integer :: bondType
!        integer :: bondMembr(1:2)        
      end type



      type BendAngleIndex
        integer(kind=2) :: bendType   
        integer(kind=2) :: bendMembr(1:3)
!        integer :: bendType   
!        integer :: bendMembr(1:3)

      end type
        
      type TorsionAngleIndex
        integer(kind=2) :: TorsType         
        integer(kind=2) :: torsMembr(1:4)
!        integer :: TorsType         
!        integer :: torsMembr(1:4)        
      end type
      
      type ImproperAngleIndex
        integer(kind=2) :: ImpropType         
        integer(kind=2) :: impropMembr(1:4)
      end type
      end module    

!==============================================================
      module ForceField
      use ForceFieldVariableType
      use VarPrecision
      integer :: nAtomTypes
      integer :: nBondTypes     
      integer :: nAngleTypes
      integer :: nTorsionalTypes
      integer :: nImproperTypes       
        
      type(AtomDef), allocatable :: atomData(:)   
      type(BondDef), allocatable :: bondData(:)
      type(BendAngleDef), allocatable :: bendData(:)
      type(TorsionAngleDef), allocatable :: torsData(:)
      type(TorsionAngleDef), allocatable :: impropData(:)

      real(dp), allocatable :: r_min(:), r_min_sq(:)
      real(dp), allocatable :: r_min_tab(:,:)

      integer, allocatable :: nAtoms(:)
      integer, allocatable :: nIntraNonBond(:)      
      integer, allocatable :: nBonds(:)  
      integer, allocatable :: nAngles(:)
      integer, allocatable :: nTorsional(:)
      integer, allocatable :: nImproper(:)

      real(dp),allocatable :: totalMass(:)
      
      integer(kind=2), allocatable :: atomArray(:,:)
!      integer, allocatable :: atomArray(:,:)
      type(NonBondedIndex), allocatable :: nonBondArray(:,:)      
      type(BondIndex), allocatable :: bondArray(:,:)
      type(BendAngleIndex), allocatable :: bendArray(:,:)
      type(TorsionAngleIndex), allocatable :: torsArray(:,:)
      type(ImproperAngleIndex), allocatable :: impropArray(:,:)

!      character(len=16), allocatable :: atomUseByBond(:,:)
!      character(len=66), allocatable :: atomUseByBend(:,:)      
!      character(len=66), allocatable :: atomUseByTors(:,:)
!      character(len=66), allocatable :: atomUseByImprop(:,:)
      
      end module
!==============================================================
      module ForceFieldPara_LJ_Q
      use ForceFieldVariableType
      use VarPrecision

      real(dp), allocatable :: ep_tab(:,:),sig_tab(:,:)
      real(dp), allocatable :: q_tab(:,:)       
      
      end module
!==============================================================
      module ForceFieldPara_Pedone
      use ForceFieldVariableType
      use VarPrecision

      type AtomDefPedone
        character(len=4) :: Symb
        real(dp) :: repul, rEq, q, alpha, delta, mass
      end type

      type(AtomDef), allocatable :: atomData(:)
      real(dp), allocatable :: repul_tab(:,:), rEq_tab(:,:)
      real(dp), allocatable :: q_tab(:,:), alpha_Tab(:,:), D_Tab(:,:)
      
      end module
!==============================================================

      module ForceFieldFunctions
        interface
          double precision function MixRule(par1,par2)
            real(kind(0.0d0)), Intent(IN) :: par1,par2
          end function
        end interface
        interface
          double precision function TorsEng(angle,coeffArray)
           real(kind(0.0d0)), Intent(IN) :: angle
           real(kind(0.0d0)), Intent(IN) :: coeffArray(:)  
          end function
        end interface        
      contains
!        !----------------------------------------------------------------
!       !Geometric Mean Mixing rule
       double precision function GeoMean_MixingFunc(par1,par2)
         real(kind(0.0d0)), Intent(IN) :: par1,par2
         
         if((par1 .eq. 0) .or.(par2 .eq. 0)) then         
           GeoMean_MixingFunc = 0d0
         else
           GeoMean_MixingFunc = dsqrt(par1*par2)
         endif
       end function
!       !----------------------------------------------------------------
!      !Standard Mean(average) Mixing rule        
       double precision function Mean_MixingFunc(par1,par2)
         real(kind(0.0d0)), Intent(IN) :: par1,par2
         if((par1 .eq. 0d0) .or.(par2 .eq. 0d0)) then
           Mean_MixingFunc = 0d0
         else
           Mean_MixingFunc = (par1+par2)*0.5d0         
         endif

       end function
!       !----------------------------------------------------------------
!      !(Fill in Later) Mixing rule        
       double precision function MeanMax_MixingFunc(par1,par2)
         real(kind(0.0d0)), Intent(IN) :: par1,par2
         if((par1 .gt. 0d0) .and. (par2 .gt. 0d0)) then
           MeanMax_MixingFunc = abs((par1+par2)*0.5d0)
         else
           MeanMax_MixingFunc = max(abs(par1), abs(par2))
         endif

       end function
!       !----------------------------------------------------------------
!      !Minimum Mixing rule        
       double precision function Min_MixingFunc(par1,par2)
         real(kind(0.0d0)), Intent(IN) :: par1,par2
         Min_MixingFunc = min(par1,par2)
       end function
!       !----------------------------------------------------------------
!      !Maximum Mixing rule        
       double precision function Max_MixingFunc(par1,par2)
         real(kind(0.0d0)), Intent(IN) :: par1,par2
         Max_MixingFunc = max(par1,par2)
       end function
!       !----------------------------------------------------------------
!       ! Standard Harmonic Potential Energy Function
       pure double precision function Harmonic(angle,k_eq,ang_eq)
         real(kind(0.0d0)), Intent(IN) :: angle,k_eq,ang_eq
         Harmonic = 0.5d0*k_eq*(angle-ang_eq)**2
       end function
!       !----------------------------------------------------------------
!       ! Standard Harmonic Potential Energy Function that accepts a Torsional Parameter Array
       pure double precision function TorsHarmonic(angle,coeffArray)
         real(kind(0.0d0)), Intent(IN) :: angle
         real(kind(0.0d0)), Intent(IN) :: coeffArray(:)  
         
         TorsHarmonic = 0.5d0*coeffArray(1)*(angle-coeffArray(2))**2
       end function       
!       !----------------------------------------------------------------
!       ! Standard Torsional Potential using a series of Cos^n functionals
       pure double precision function Torsion_Cos_N(angle,coeffArray)
         integer :: nCoeff,i
         real(kind(0.0d0)), Intent(IN) :: angle
         real(kind(0.0d0)), Intent(IN) :: coeffArray(:)            
         real(kind(0.0d0)) :: eng
             
         nCoeff = size(coeffArray)
         eng = coeffArray(1)
         do i=2,nCoeff
           if(coeffArray(i) .ne. 0d0) then
             eng = eng + coeffArray(i)*cos(angle)**(i-1)
           endif
         enddo
         Torsion_Cos_N = eng
       end function
!       !----------------------------------------------------------------       
!       ! Torsional Potential using a series of Cos(nx) functions      
       pure double precision function Torsion_CosNx(angle,coeffArray)
         integer :: nCoeff,i
         real(kind(0.0d0)), Intent(IN) :: angle
         real(kind(0.0d0)), Intent(IN) :: coeffArray(:)            
         real(kind(0.0d0)) :: eng
             
         nCoeff = size(coeffArray)
         eng = coeffArray(1)
         do i=2,nCoeff
           if(coeffArray(i) .ne. 0d0) then
             eng = eng + coeffArray(i)*cos(dble(i-1)*angle)
           endif
         enddo
         Torsion_CosNx = eng
       end function        
!       !----------------------------------------------------------------       
!       ! Trappe Torsional Potential using a series of Cos(nx) functions      
       pure double precision function Trappe_CosNx(angle,coeffArray)
         integer :: nCoeff
         real(kind(0.0d0)), Intent(IN) :: angle
         real(kind(0.0d0)), Intent(IN) :: coeffArray(:)            
         real(kind(0.0d0)) :: eng
             
         nCoeff = size(coeffArray)
         eng = coeffArray(1)
         eng = eng + coeffArray(2)*(1d0+cos(angle))
         eng = eng + coeffArray(3)*(1d0-cos(2d0*angle))
         eng = eng + coeffArray(4)*(1d0+cos(3d0*angle))
         Trappe_CosNx = eng
       end function
      
!       !----------------------------------------------------------------       
      end module

!==============================================================
