      module AVBMC_CBMC
      contains
!=======================================================================
!     The purpose of this subroutine is generate a set of trial configurations that will
!     used for the Aggregation-Volume-Bias insertion move. This subroutine is intended to generate
!     configurations for perfectly ridgid molecules. Once a set of configurations is generated
!     one of these configurations will be randomly chosen based on the Rosenbluth weight for each
!     configuraiton. 
      subroutine Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPointers, only: Rosen_Mol_New
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType, nIndx
      logical, intent(out) :: rejMove      
      real(dp), intent(out) :: rosenRatio
      logical, intent(in) :: isIncluded(:)

      logical :: overlap(1:maxRosenTrial)
      integer :: atmType1, atmType2
      integer :: i, iRosen, nSel, nTargetMol
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: dx, dy, dz
      real(dp) :: r, x1, y1, z1      
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt, rmin_ij

      
      newMol%molType = nType      
!      call Rosen_CreateSubset(nTarget, isIncluded)
      E_Trial = 0d0
      rejMove = .false.   
      nTargetMol = subIndxList(nTarget)

      atmType1 = atomArray(nType, 1)
      atmType2 = atomArray(nTargType, 1)
      rmin_ij = r_min_tab(atmType2, atmType1)
      rmin_ij = sqrt(rmin_ij)

      do iRosen = 1, nRosenTrials(nType)
        rosenTrial(iRosen)%x(1:nAtoms(nType)) = gasConfig(nType)%x(1:nAtoms(nType))
        rosenTrial(iRosen)%y(1:nAtoms(nType)) = gasConfig(nType)%y(1:nAtoms(nType))
        rosenTrial(iRosen)%z(1:nAtoms(nType)) = gasConfig(nType)%z(1:nAtoms(nType))
      
        if(nAtoms(nType) .ne. 1) then
          x_rid_cm = gasConfig(nType)%x(1)
          y_rid_cm = gasConfig(nType)%y(1)
          z_rid_cm = gasConfig(nType)%z(1)    
        
!          Rotate xz axis
          rotang = two_pi*grnd()
          c_term = dcos(rotang)
          s_term = dsin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = rosenTrial(iRosen)%x(i) - x_rid_cm
            z_shift = rosenTrial(iRosen)%z(i) - z_rid_cm        
            rosenTrial(iRosen)%x(i) = c_term*x_shift - s_term*z_shift + x_rid_cm
            rosenTrial(iRosen)%z(i) = s_term*x_shift + c_term*z_shift + z_rid_cm
          enddo   
        
!          Rotate yz axis
          rotang = two_pi*grnd()
          c_term = dcos(rotang)
          s_term = dsin(rotang)
          do i = 1,nAtoms(nType)
            y_shift = rosenTrial(iRosen)%y(i) - y_rid_cm
            z_shift = rosenTrial(iRosen)%z(i) - z_rid_cm        
            rosenTrial(iRosen)%y(i) = c_term*y_shift - s_term*z_shift + y_rid_cm
            rosenTrial(iRosen)%z(i) = s_term*y_shift + c_term*z_shift + z_rid_cm
          enddo   
  
!          Rotate xy axis
          rotang = two_pi*grnd()
          c_term = dcos(rotang)
          s_term = dsin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = rosenTrial(iRosen)%x(i) - x_rid_cm
            y_shift = rosenTrial(iRosen)%y(i) - y_rid_cm        
            rosenTrial(iRosen)%x(i) = c_term*x_shift - s_term*y_shift + x_rid_cm
            rosenTrial(iRosen)%y(i) = s_term*x_shift + c_term*y_shift + y_rid_cm
          enddo    
        endif

        r = Dist_Critr * grnd()**(1d0/3d0)
!        if(r .lt. rmin_ij) then
!          overlap(iRosen) = .true.
!        endif

        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz      
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - rosenTrial(iRosen)%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - rosenTrial(iRosen)%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - rosenTrial(iRosen)%z(1)

        do i=1,nAtoms(nType)
          rosenTrial(iRosen)%x(i) = rosenTrial(iRosen)%x(i) + x1
          rosenTrial(iRosen)%y(i) = rosenTrial(iRosen)%y(i) + y1
          rosenTrial(iRosen)%z(i) = rosenTrial(iRosen)%z(i) + z1
        enddo 

!        call Rosen_BoltzWeight_Molecule_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
        call Rosen_Mol_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
      enddo

      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen) - E_Max))
      enddo
!      if(all(ProbRosen .eq. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      
      ranNum = grnd()*rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
         nSel = nSel+1
         sumInt = sumInt + ProbRosen(nSel)
      enddo
      
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif

!      Update the coordinates      
      newMol%x(1:nAtoms(nType)) = rosenTrial(nSel)%x(1:nAtoms(nType))
      newMol%y(1:nAtoms(nType)) = rosenTrial(nSel)%y(1:nAtoms(nType))
      newMol%z(1:nAtoms(nType)) = rosenTrial(nSel)%z(1:nAtoms(nType))
      
      rosenRatio = (ProbRosen(nSel)*dble(nRosenTrials(nType)))/rosenNorm
      
      end subroutine
!=======================================================================
!     This subroutine is used for
      subroutine Ridgid_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPointers, only: Rosen_Mol_Old
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nMol, nTarget, nTargType
      real(dp), intent(out) :: rosenRatio

      logical :: isIncluded(1:maxMol)
      logical :: overlap(1:maxRosenTrial)
      integer :: i, iRosen, nSel, nIndx, nTargetMol
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: dx, dy, dz
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: r, x1, y1, z1       

      ProbRosen = 0d0      
      E_Trial = 0d0      
      nIndx = molArray(nType)%mol(nMol)%indx
!      write(2,*) NPART
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
!      call Rosen_BoltzWeight_Molecule_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))
      call Rosen_Mol_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))
      nTargetMol = subIndxList(nTarget)
      do iRosen = 2, nRosenTrials(nType)
        newMol%x(1:nAtoms(nType)) = gasConfig(nType)%x(1:nAtoms(nType))
        newMol%y(1:nAtoms(nType)) = gasConfig(nType)%y(1:nAtoms(nType))
        newMol%z(1:nAtoms(nType)) = gasConfig(nType)%z(1:nAtoms(nType))
      
        if(nAtoms(nType) .ne. 1) then
          x_rid_cm = gasConfig(nType)%x(1)
          y_rid_cm = gasConfig(nType)%y(1)
          z_rid_cm = gasConfig(nType)%z(1)    
        
!          Rotate xz axis
          rotang = two_pi * grnd()
          c_term = dcos(rotang)
          s_term = dsin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = newMol%x(i) - x_rid_cm
            z_shift = newMol%z(i) - z_rid_cm        
            newMol%x(i) = c_term*x_shift - s_term*z_shift + x_rid_cm
            newMol%z(i) = s_term*x_shift + c_term*z_shift + z_rid_cm
          enddo   
        
!          Rotate yz axis
          rotang = two_pi*grnd()
          c_term = dcos(rotang)
          s_term = dsin(rotang)
          do i = 1,nAtoms(nType)
            y_shift = newMol%y(i) - y_rid_cm
            z_shift = newMol%z(i) - z_rid_cm        
            newMol%y(i) = c_term*y_shift - s_term*z_shift + y_rid_cm
            newMol%z(i) = s_term*y_shift + c_term*z_shift + z_rid_cm
          enddo   
  
!          Rotate xy axis
          rotang = two_pi*grnd()
          c_term = dcos(rotang)
          s_term = dsin(rotang)
          do i = 1,nAtoms(nType)
            x_shift = newMol%x(i) - x_rid_cm
            y_shift = newMol%y(i) - y_rid_cm        
            newMol%x(i) = c_term*x_shift - s_term*y_shift + x_rid_cm
            newMol%y(i) = s_term*x_shift + c_term*y_shift + y_rid_cm
          enddo    
        endif

        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz      
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - newMol%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - newMol%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - newMol%z(1)

        do i=1,nAtoms(nType)
          newMol%x(i) = newMol%x(i) + x1
          newMol%y(i) = newMol%y(i) + y1
          newMol%z(i) = newMol%z(i) + z1
        enddo 

!        call Rosen_BoltzWeight_Molecule_Old( newMol%x(:), newMol%y(:), newMol%z(:), &
!                                   nType, isIncluded, E_Trial(iRosen) )
        call Rosen_Mol_Old( newMol%x(:), newMol%y(:), newMol%z(:), &
                                   nType, isIncluded, E_Trial(iRosen) )
      enddo

      E_Max = minval(E_Trial)
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen) - E_Max))
      enddo

      rosenNorm = sum(ProbRosen)
      rosenRatio = (ProbRosen(1)*dble(nRosenTrials(nType)))/rosenNorm
      
      end subroutine
!=======================================================================
      subroutine Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPointers, only: Rosen_Mol_New
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType, nIndx
      logical, intent(out) :: rejMove      
      real(dp), intent(out):: rosenRatio
      logical, intent(in) :: isIncluded(:)
      
!      logical :: isIncluded(1:maxMol)
      logical :: overlap(1:maxRosenTrial)
      integer :: atmType1, atmType2
      integer :: i, iRosen, nSel, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: bendType, bondType  
      integer :: bendType1, bendType2, bendType3
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob
      real(dp) :: r1, r2, r3
      real(dp) :: k_bend, ang_eq, ang
      real(dp) :: ang1, ang2, dihed
      real(dp) :: x1, y1, z1, dx, dy, dz 
      real(dp) :: rmin_ij
      type(SimpleAtomCoords) :: v1, v2, v3
      
      newMol%molType = nType      
!      call Rosen_CreateSubset(nTarget, isIncluded)
      E_Trial = 0d0
      rejMove = .false.      
      nTargetMol = subIndxList(nTarget)

      atmType1 = atomArray(nType, 1)
      atmType2 = atomArray(nTargType, 1)
      rmin_ij = r_min_tab(atmType2, atmType1)
      rmin_ij = sqrt(rmin_ij)

      do iRosen = 1, nRosenTrials(nType)
      
!        Initialize the first atom coordinates to 0      
        rosenTrial(iRosen)%x = 0d0
        rosenTrial(iRosen)%y = 0d0
        rosenTrial(iRosen)%z = 0d0

!        Depending the molecule type, choose an appropriate algorithm
!        and regrow the atoms starting from atom 1. 
        select case(nAtoms(nType))
        case(1)
          continue
        case(2)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)        
          bondType = bondArray(nType,1)%bondType
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r, k_bond, r_eq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          rosenTrial(iRosen)%x(Atm2) = r * dx
          rosenTrial(iRosen)%y(Atm2) = r * dy
          rosenTrial(iRosen)%z(Atm2) = r * dz
        case(3)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(1, 3)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r, k_bond, r_eq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r*dx
          v1%y = -r*dy
          v1%z = -r*dz
          rosenTrial(iRosen)%x(atm2) = r * dx
          rosenTrial(iRosen)%y(atm2) = r * dy
          rosenTrial(iRosen)%z(atm2) = r * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r, k_bond, r_eq, Prob)
          bendType = bendArray(nType,1)%bendType
!          k_bend = bendData(bendType)%k_eq
!          ang_eq = bendData(bendType)%ang_eq
!          call GenerateBendAngle(ang, k_bend, ang_eq, Prob)
          call GenerateBendAngle(ang, bendType, Prob)
          call Generate_UnitCone(v1, r, ang, v2)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
        case(4)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r1, k_bond, r_eq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          rosenTrial(iRosen)%x(atm2) = r1 * dx
          rosenTrial(iRosen)%y(atm2) = r1 * dy
          rosenTrial(iRosen)%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r2, k_bond, r_eq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r3, k_bond, r_eq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType1)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          call GenerateTwoBranches(ang1, ang2, dihed, bendType1, bendType2, bendType3, Prob)
          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
        case default
         stop "Error! Molecule has too many atoms for a simple regrowth"
        end select

        r = Dist_Critr * grnd()**(1d0/3d0)
        if(r .lt. rmin_ij) then
          overlap(iRosen) = .true.
        endif        
        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz   
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - rosenTrial(iRosen)%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - rosenTrial(iRosen)%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - rosenTrial(iRosen)%z(1)
        
        do i=1,nAtoms(nType)
          rosenTrial(iRosen)%x(i) = rosenTrial(iRosen)%x(i) + x1 
          rosenTrial(iRosen)%y(i) = rosenTrial(iRosen)%y(i) + y1 
          rosenTrial(iRosen)%z(i) = rosenTrial(iRosen)%z(i) + z1 
        enddo 

!        call Rosen_BoltzWeight_Molecule_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
        call Rosen_Mol_New(iRosen, nType, isIncluded, E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo

!      if(all(ProbRosen .le. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
         nSel = nSel+1
         sumInt = sumInt + ProbRosen(nSel)
      enddo

      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif

!      Update the coordinates      
      newMol%x(1:nAtoms(nType)) = rosenTrial(nSel)%x(1:nAtoms(nType))
      newMol%y(1:nAtoms(nType)) = rosenTrial(nSel)%y(1:nAtoms(nType))
      newMol%z(1:nAtoms(nType)) = rosenTrial(nSel)%z(1:nAtoms(nType))
      
      rosenRatio = (ProbRosen(nSel)*dble(nRosenTrials(nType)))/rosenNorm
      
      end subroutine
!=======================================================================
      subroutine Simple_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
!      use Rosenbluth_Functions
      use EnergyPointers, only: Rosen_Mol_Old
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType, nMol
      real(dp), intent(out):: rosenRatio
      
      logical :: isIncluded(1:maxMol)

      integer :: i, iRosen, nSel, nIndx, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: bondType, bendType
      integer :: bendType1, bendType2, bendType3
      real(dp) :: E_Trial(1:maxRosenTrial)      
      real(dp) :: grnd,rotang
      real(dp) :: c_term,s_term
      real(dp) :: x_shift,y_shift,z_shift
      real(dp) :: x_rid_cm, y_rid_cm,z_rid_cm
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob
      real(dp) :: r1, r2, r3
      real(dp) :: k_bend, ang_eq, ang
      real(dp) :: ang1, ang2, dihed
      real(dp) :: x1, y1, z1, dx, dy, dz 
      type(SimpleAtomCoords) :: v1, v2, v3
      
      
      ProbRosen = 0d0      
      E_Trial = 0d0      
      nTargetMol = subIndxList(nTarget)
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
!      call Rosen_BoltzWeight_Molecule_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
!                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))      

      call Rosen_Mol_Old(molArray(nType)%mol(nMol)%x(:), molArray(nType)%mol(nMol)%y(:), &
                                 molArray(nType)%mol(nMol)%z(:), nType, isIncluded, E_Trial(1))     

      do iRosen = 2, nRosenTrials(nType)
      
        newMol%x = 0d0
        newMol%y = 0d0
        newMol%z = 0d0

!        Depending the molecule type, choose an appropriate algorithm
!        and regrow the atoms starting from atom 1. 
        select case(nAtoms(nType))
        case(1)
          continue
        case(2)
          bondType = bondArray(nType,1)%bondType
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r, k_bond, r_eq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          newMol%x(2) = r*dx
          newMol%y(2) = r*dy
          newMol%z(2) = r*dz
        case(3)
          Atm1 = pathArray(nType)%path(1,1)
          Atm2 = pathArray(nType)%path(1,2)
          Atm3 = pathArray(nType)%path(1,3)
          call FindBond(nType,Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r, k_bond, r_eq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r*dx
          v1%y = -r*dy
          v1%z = -r*dz
          newMol%x(atm2) = r*dx
          newMol%y(atm2) = r*dy
          newMol%z(atm2) = r*dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r, k_bond, r_eq, Prob)
          bendType = bendArray(nType,1)%bendType
!          k_bend = bendData(bendType)%k_eq
!          ang_eq = bendData(bendType)%ang_eq
!          call GenerateBendAngle(ang, k_bend, ang_eq, Prob)
          call GenerateBendAngle(ang, bendType, Prob)
          call Generate_UnitCone(v1, r, ang, v2)
          newMol%x(atm3) = newMol%x(atm2) + v2%x 
          newMol%y(atm3) = newMol%y(atm2) + v2%y
          newMol%z(atm3) = newMol%z(atm2) + v2%z
        case(4)
          Atm1 = pathArray(nType)%path(1, 1)
          Atm2 = pathArray(nType)%path(1, 2)
          Atm3 = pathArray(nType)%path(2, 1)
          Atm4 = pathArray(nType)%path(3, 1)
          call FindBond(nType, Atm1, Atm2, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r1, k_bond, r_eq, Prob)
          call Generate_UnitSphere(dx, dy, dz)
          v1%x = -r1*dx
          v1%y = -r1*dy
          v1%z = -r1*dz
          newMol%x(atm2) = r1 * dx
          newMol%y(atm2) = r1 * dy
          newMol%z(atm2) = r1 * dz
          call FindBond(nType, Atm2, Atm3, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r2, k_bond, r_eq, Prob)
          call FindBond(nType, Atm2, Atm4, bondType)
          k_bond = bondData(bondType)%k_eq
          r_eq = bondData(bondType)%r_eq
          call GenerateBondLength(r3, k_bond, r_eq, Prob)
          call FindAngle(nType, Atm1, Atm2, Atm3, bendType1)
          call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
          call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
          call GenerateTwoBranches(ang1, ang2, dihed, bendType1, bendType2, bendType3, Prob)
          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
          newMol%x(atm3) = newMol%x(atm2) + v2%x 
          newMol%y(atm3) = newMol%y(atm2) + v2%y
          newMol%z(atm3) = newMol%z(atm2) + v2%z
          newMol%x(atm4) = newMol%x(atm2) + v3%x 
          newMol%y(atm4) = newMol%y(atm2) + v3%y
          newMol%z(atm4) = newMol%z(atm2) + v3%z
        case default
         stop "Error! Molecule has too many atoms for a simple regrowth"
        end select
        
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        dx = r * dx
        dy = r * dy
        dz = r * dz   
        
        x1 = molArray(nTargType)%mol(nTargetMol)%x(1) + dx - newMol%x(1)
        y1 = molArray(nTargType)%mol(nTargetMol)%y(1) + dy - newMol%y(1)
        z1 = molArray(nTargType)%mol(nTargetMol)%z(1) + dz - newMol%z(1)
        
        do i=1,nAtoms(nType)
          newMol%x(i) = newMol%x(i) + x1
          newMol%y(i) = newMol%y(i) + y1
          newMol%z(i) = newMol%z(i) + z1
        enddo 

!        call Rosen_BoltzWeight_Molecule_Old( newMol%x(:), newMol%y(:), newMol%z(:), nType, isIncluded, E_Trial(iRosen) )
        call Rosen_Mol_Old( newMol%x(:), newMol%y(:), newMol%z(:), nType, isIncluded, E_Trial(iRosen) )

      enddo

      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      
      rosenRatio = (ProbRosen(1)*dble(nRosenTrials(nType)))/rosenNorm
      
      end subroutine
!=======================================================================
      subroutine StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      use SimParameters
      use Coords
      use ForceField
      use Constants
      use Rosenbluth_Functions_LJ_Q
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nTarget, nTargType
      real(dp), intent(out):: rosenRatio
      logical, intent(out) :: rejMove
      logical, intent(in) :: isIncluded(:)

!      logical :: isIncluded(1:maxMol)
      integer :: i, iRosen, iAtom, nSel, nIndx, nTargetMol
      integer :: nRegrown
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: atom1_Pos, cnt
      integer :: bondType, bendType, torsType
      logical :: overlap(1:maxRosenTrial)
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial), E_Complete
      real(dp) :: grnd,rotang
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob
      real(dp) :: k_bend, ang_eq, bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      real(dp) :: q1, q2
      type(SimpleAtomCoords) :: trialPos(1:maxRosenTrial)
      type(SimpleAtomCoords) :: v1, v2, v3
      
   
  
!      call Rosen_CreateSubset(nTarget, isIncluded)
      nTargetMol = subIndxList(nTarget)
      newMol%molType = nType 
      newMol%x = 0d0
      newMol%y = 0d0
      newMol%z = 0d0
      regrown = .false.
      rejMove = .false.
      nRegrown = 0


!      Atm1 = 1
!      Atm2 = regrowOrder(nType, 2) 
!
!      q1 = 
!      q2 = 

!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      
      E_Trial = 0d0
      E_Complete = 0d0
      overlap = .false.
      do iRosen = 1, nRosenTrials(nType)
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(1)
        trialPos(iRosen)%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(1)
        trialPos(iRosen)%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(1)
        call Rosen_BoltzWeight_Atom_New(nType, 1, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo

      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = ProbRosen(nSel)*dble(nRosenTrials(nType))/rosenNorm
      regrown(1) = .true.
      nRegrown = nRegrown + 1
      newMol%x(1) = trialPos(nSel)%x 
      newMol%y(1) = trialPos(nSel)%y 
      newMol%z(1) = trialPos(nSel)%z

!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      Atm1 = 1
      Atm2 = regrowOrder(nType, 2) 
      call FindBond(nType,Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      overlap = .false.
      E_Trial = 0d0
      do iRosen = 1, nRosenTrials(nType)
        call GenerateBondLength(r, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos(iRosen)%x = r*dx + newMol%x(Atm1) 
        trialPos(iRosen)%y = r*dy + newMol%y(Atm1)
        trialPos(iRosen)%z = r*dz + newMol%z(Atm1)
        call Rosen_BoltzWeight_Atom_New(nType, Atm2, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      if(all(ProbRosen .le. 0d0)) then
        rejMove = .true.
        return
      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)*dble(nRosenTrials(nType))/rosenNorm
      regrown(Atm2) = .true.
      nRegrown = nRegrown + 1
      newMol%x(Atm2) = trialPos(nSel)%x 
      newMol%y(Atm2) = trialPos(nSel)%y 
      newMol%z(Atm2) = trialPos(nSel)%z

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif
      v1%x = newMol%x(Atm1) - newMol%x(Atm2)
      v1%y = newMol%y(Atm1) - newMol%y(Atm2)
      v1%z = newMol%z(Atm1) - newMol%z(Atm2)
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!      k_bend = bendData(bendType)%k_eq
!      ang_eq = bendData(bendType)%ang_eq
      overlap = .false.
      do iRosen = 1, nRosenTrials(nType)
        call GenerateBondLength(r, k_bond, r_eq, Prob)
!        call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
        call GenerateBendAngle(bend_angle, bendType, Prob)
        call Generate_UnitCone(v1, r, bend_angle, v2)
        trialPos(iRosen)%x = v2%x + newMol%x(Atm2) 
        trialPos(iRosen)%y = v2%y + newMol%y(Atm2)
        trialPos(iRosen)%z = v2%z + newMol%z(Atm2)
        call Rosen_BoltzWeight_Atom_New(nType, Atm3, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
!      if(all(ProbRosen .le. 0d0)) then
!        rejMove = .true.
!        return
!      endif
      rosenNorm = sum(ProbRosen)
      ranNum = grnd() * rosenNorm
      sumInt = ProbRosen(1)
      nSel = 1
      do while(sumInt .lt. ranNum)
        nSel = nSel + 1
        sumInt = sumInt + ProbRosen(nSel)
      enddo
      if(overlap(nSel) .eqv. .true.) then
        rejMove = .true.
        return
      endif
      rosenRatio = rosenRatio*ProbRosen(nSel)*dble(nRosenTrials(nType))/rosenNorm
      regrown(Atm3) = .true.
      nRegrown = nRegrown + 1
      newMol%x(Atm3) = trialPos(nSel)%x 
      newMol%y(Atm3) = trialPos(nSel)%y 
      newMol%z(Atm3) = trialPos(nSel)%z


!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = newMol%x(Atm1) - newMol%x(Atm3)
        v1%y = newMol%y(Atm1) - newMol%y(Atm3)
        v1%z = newMol%z(Atm1) - newMol%z(Atm3)

        v2%x = newMol%x(Atm2) - newMol%x(Atm3)
        v2%y = newMol%y(Atm2) - newMol%y(Atm3)
        v2%z = newMol%z(Atm2) - newMol%z(Atm3)
        overlap = .false.
        E_Trial = 0d0
        do iRosen = 1, nRosenTrials(nType)
          call GenerateBondLength(r, k_bond, r_eq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos(iRosen)%x = v3%x + newMol%x(Atm3) 
          trialPos(iRosen)%y = v3%y + newMol%y(Atm3)
          trialPos(iRosen)%z = v3%z + newMol%z(Atm3)
          call Rosen_BoltzWeight_Atom_New(nType, Atm4, trialPos(iRosen), isIncluded,  E_Trial(iRosen), overlap(iRosen))
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
!        if(all(ProbRosen .le. 0d0)) then
!          rejMove = .true.
!          return
!        endif
        rosenNorm = sum(ProbRosen)
        ranNum = grnd() * rosenNorm
        sumInt = ProbRosen(1)
        nSel = 1
        do while(sumInt .lt. ranNum)
          nSel = nSel + 1
          sumInt = sumInt + ProbRosen(nSel)
        enddo
        if(overlap(nSel) .eqv. .true.) then
          rejMove = .true.
          return
        endif
        rosenRatio = rosenRatio*ProbRosen(nSel)*dble(nRosenTrials(nType))/rosenNorm
        regrown(Atm4) = .true.
        nRegrown = nRegrown + 1
        newMol%x(Atm4) = trialPos(nSel)%x 
        newMol%y(Atm4) = trialPos(nSel)%y 
        newMol%z(Atm4) = trialPos(nSel)%z
      enddo


      end subroutine
!=======================================================================
      subroutine StraightChain_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      use SimParameters
      use Coords
      use ForceField
      use Constants
      use Rosenbluth_Functions_LJ_Q
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nType, nMol, nTarget, nTargType
      real(dp), intent(out):: rosenRatio

      integer :: i, iRosen, iAtom, nSel, nIndx, nTargetMol
      integer :: Atm1, Atm2, Atm3, Atm4
      integer :: atom1_Pos
      integer :: bondType, bendType, torsType, cnt
      logical :: overlap
      logical :: isIncluded(1:maxMol)
      logical :: regrown(1:maxAtoms)
      real(dp) :: E_Trial(1:maxRosenTrial), E_Complete
      real(dp) :: grnd,rotang
      real(dp) :: E_Max, ProbRosen(1:maxRosenTrial), rosenNorm
      real(dp) :: ranNum, sumInt
      real(dp) :: k_bond, r_eq, r, Prob
      real(dp) :: k_bend, ang_eq, bend_angle, tors_angle
      real(dp) :: dx, dy, dz 
      type(SimpleAtomCoords) :: trialPos
      type(SimpleAtomCoords) :: v1, v2, v3
      
      ProbRosen = 0d0      
      E_Trial = 0d0      
      nTargetMol = subIndxList(nTarget)
      nIndx = molArray(nType)%mol(nMol)%indx
      call Rosen_CreateSubset_Reverse(nTarget, nIndx, isIncluded)
      regrown = .false.

!      Begin the regrowth process by choosing an insertion site for the first atom in the chain
      E_Trial = 0d0
      E_Complete = 0d0
      rosenRatio = 1d0
      call Rosen_BoltzWeight_Atom_Old(nType, nMol, 1, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        r = Dist_Critr * grnd()**(1d0/3d0)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r * dx + molArray(nTargType)%mol(nTargetMol)%x(1)
        trialPos%y = r * dy + molArray(nTargType)%mol(nTargetMol)%y(1)
        trialPos%z = r * dz + molArray(nTargType)%mol(nTargetMol)%z(1)
        call Rosen_BoltzWeight_Atom_New(nType, 1, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
      regrown(1) = .true.


!      Having inserted the first atom, when we go to insert the second atom we must account for the fact that
!      the second atom must be inserted at the correct bond distance from the first atom.
      Atm1 = 1
      Atm2 = regrowOrder(nType, 2) 
      call FindBond(nType,Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      call Rosen_BoltzWeight_Atom_Old(nType, nMol, Atm2, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        call GenerateBondLength(r, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        trialPos%x = r*dx + molArray(nType)%mol(nMol)%x(Atm1)
        trialPos%y = r*dy + molArray(nType)%mol(nMol)%y(Atm1)
        trialPos%z = r*dz + molArray(nType)%mol(nMol)%z(Atm1)
        call Rosen_BoltzWeight_Atom_New(nType, Atm2, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo
      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
      regrown(Atm2) = .true.

!       For the third atom we must begin to include the bending angle in choosing its trial positions. The first thing
!       we must consider is if the second regrown atom was a terminal atom or a linker atom. If it was a terminal atom 
!       then the bending angle must use the 1st atom as the central atom for the angle generation.  However if it was a link in the chain then the 2nd 
!       atom regrown should be used as the central atom.
      if(topolArray(nType)%atom(Atm2) .eq. 2) then
        Atm1 = regrowOrder(nType, 1) 
        Atm2 = regrowOrder(nType, 2) 
        Atm3 = regrowOrder(nType, 3) 
      else
        Atm1 = regrowOrder(nType, 2) 
        Atm2 = regrowOrder(nType, 1) 
        Atm3 = regrowOrder(nType, 3) 
      endif

      v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm2)
      v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm2)
      v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm2)
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!      k_bend = bendData(bendType)%k_eq
!      ang_eq = bendData(bendType)%ang_eq

      call Rosen_BoltzWeight_Atom_Old(nType, nMol, Atm3, isIncluded,  E_Trial(1))
      do iRosen = 2, nRosenTrials(nType)
        call GenerateBondLength(r, k_bond, r_eq, Prob)
!        call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
        call GenerateBendAngle(bend_angle, bendType, Prob)
        call Generate_UnitCone(v1, r, bend_angle, v2)
        trialPos%x = v2%x + molArray(nType)%mol(nMol)%x(Atm2)
        trialPos%y = v2%y + molArray(nType)%mol(nMol)%y(Atm2)
        trialPos%z = v2%z + molArray(nType)%mol(nMol)%z(Atm2)
        call Rosen_BoltzWeight_Atom_New(nType, Atm3, trialPos, isIncluded,  E_Trial(iRosen), overlap)
      enddo

      E_Max = minval(E_Trial)
      ProbRosen = 0d0
      do iRosen = 1, nRosenTrials(nType)
        ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
      enddo
      rosenNorm = sum(ProbRosen)
      rosenRatio = rosenRatio*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
      regrown(Atm3) = .true.

!      Now that three atoms have been regrown, all remaining atoms must take the torsional angles into account.
      cnt = 3
      do while(any(regrown .eqv. .false.))
        cnt = cnt + 1
        atm4 = regrowOrder(nType, cnt) 
        call FindAtomsFromPath(nType, regrown, 1, Atm4, Atm1, Atm2, Atm3)
        call FindBond(nType, Atm3, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call FindAngle(nType, Atm2, Atm3, Atm4, bendType)
!        k_bend = bendData(bendType)%k_eq
!        ang_eq = bendData(bendType)%ang_eq
        call FindTorsion(nType, Atm1, Atm2, Atm3, Atm4, torsType)
        v1%x = molArray(nType)%mol(nMol)%x(Atm1) - molArray(nType)%mol(nMol)%x(Atm3)
        v1%y = molArray(nType)%mol(nMol)%y(Atm1) - molArray(nType)%mol(nMol)%y(Atm3)
        v1%z = molArray(nType)%mol(nMol)%z(Atm1) - molArray(nType)%mol(nMol)%z(Atm3)

        v2%x = molArray(nType)%mol(nMol)%x(Atm2) - molArray(nType)%mol(nMol)%x(Atm3)
        v2%y = molArray(nType)%mol(nMol)%y(Atm2) - molArray(nType)%mol(nMol)%y(Atm3)
        v2%z = molArray(nType)%mol(nMol)%z(Atm2) - molArray(nType)%mol(nMol)%z(Atm3)
        overlap = .false.
        call Rosen_BoltzWeight_Atom_Old(nType, nMol, Atm4, isIncluded,  E_Trial(1))
        do iRosen = 2, nRosenTrials(nType)
          call GenerateBondLength(r, k_bond, r_eq, Prob)
!          call GenerateBendAngle(bend_angle, k_bend, ang_eq, Prob)
          call GenerateBendAngle(bend_angle, bendType, Prob)
          call GenerateTorsAngle(tors_angle, torsType, Prob)
          call Generate_UnitTorsion(v1, v2, r, bend_angle, tors_angle, v3)
          trialPos%x = v3%x + molArray(nType)%mol(nMol)%x(Atm3)
          trialPos%y = v3%y + molArray(nType)%mol(nMol)%y(Atm3)
          trialPos%z = v3%z + molArray(nType)%mol(nMol)%z(Atm3)
          call Rosen_BoltzWeight_Atom_New(nType, Atm4, trialPos, isIncluded,  E_Trial(iRosen), overlap)
        enddo
        E_Max = minval(E_Trial)
        ProbRosen = 0d0
        do iRosen = 1, nRosenTrials(nType)
          ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen)-E_Max))         
        enddo
        rosenNorm = sum(ProbRosen)
        rosenRatio = rosenRatio*ProbRosen(1)*dble(nRosenTrials(nType))/rosenNorm
        regrown(Atm4) = .true.
      enddo


      end subroutine
!=======================================================================
      end module
