!========================================================================
      module CoodinateFunctions
      contains
!========================================================================
      subroutine AllocateCoordinateArrays
      use ForceField
      use SimParameters
      use Coords
      use EnergyTables
      use PairStorage
      implicit none
      integer :: i,j,k,cnt,AllocationStatus
      integer :: iType

      allocate( MolArray(1:nMolTypes), stat=AllocationStatus)
      do i = 1, nMolTypes
        allocate( MolArray(i)%mol(1:NMAX(i)), stat=AllocationStatus)   
        do j = 1, NMAX(i)
          allocate( MolArray(i)%mol(j)%x(1:nAtoms(i)), stat=AllocationStatus)
          allocate( MolArray(i)%mol(j)%y(1:nAtoms(i)), stat=AllocationStatus)
          allocate( MolArray(i)%mol(j)%z(1:nAtoms(i)), stat=AllocationStatus)
          allocate( MolArray(i)%mol(j)%globalIndx(1:nAtoms(i)), stat=AllocationStatus)
        enddo
      enddo

      nTotalAtoms = 0
      do iType = 1, nMolTypes
        nTotalAtoms = nTotalAtoms + NMAX(iType)*nAtoms(iType)
      enddo
      write(35,*) "Total Atoms in the System:", nTotalAtoms

      allocate(atomIndicies(1:nTotalAtoms), stat=AllocationStatus)
      cnt = 0
      do i = 1, nMolTypes
        do j = 1, NMAX(i)
          do k = 1, nAtoms(i)
            cnt = cnt + 1
            MolArray(i)%mol(j)%globalIndx(k) = cnt
            atomIndicies(cnt)%nType = i
            atomIndicies(cnt)%nMol = j
            atomIndicies(cnt)%nAtom = k
          enddo
        enddo
      enddo

!      do i = 1, size(atomIndicies)
!        write(35,*) i, atomIndicies(i)%nType, atomIndicies(i)%nMol, atomIndicies(i)%nAtom
!      enddo

      allocate( gasConfig(1:nMolTypes), stat=AllocationStatus)                
      do i=1,nMolTypes
        allocate( gasConfig(i)%x(1:nAtoms(i)), stat=AllocationStatus)   
        allocate( gasConfig(i)%y(1:nAtoms(i)), stat=AllocationStatus)   
        allocate( gasConfig(i)%z(1:nAtoms(i)), stat=AllocationStatus)           
      enddo      

      
      maxMol = sum(NMAX)
      allocate(typeList(1:maxMol), stat=AllocationStatus)
      allocate(subIndxList(1:maxMol), stat=AllocationStatus)
      allocate(ETable(1:maxMol), stat=AllocationStatus)
      allocate(NeiETable(1:maxMol), stat=AllocationStatus)
      allocate(neiCount(1:maxMol), stat=AllocationStatus)
      
      cnt = 0
      do i=1, nMolTypes      
        do j = 1, NMAX(i)      
          cnt = cnt + 1
          typeList(cnt) = i
          subIndxList(cnt) = j
        enddo
      enddo
      
      allocate(isActive(1:maxMol), stat=AllocationStatus)
      allocate(NeighborList(1:maxMol,1:maxMol), stat=AllocationStatus)
      maxAtoms = maxval(nAtoms)
      
      allocate(newMol%x(1:maxAtoms), stat=AllocationStatus)
      allocate(newMol%y(1:maxAtoms), stat=AllocationStatus)
      allocate(newMol%z(1:maxAtoms), stat=AllocationStatus)

      allocate(newMol2%x(1:maxAtoms), stat=AllocationStatus)
      allocate(newMol2%y(1:maxAtoms), stat=AllocationStatus)
      allocate(newMol2%z(1:maxAtoms), stat=AllocationStatus)

      call CreateDistArrays
      
      end subroutine
!=============================================================================
      subroutine CreateJointArray
      use SimParameters
      use Coords
      implicit none
      integer :: i,j
      integer :: cnt
      
      
      cnt=0
!      write(35,*) "Molecular Index List"
!      write(35,*) "Type, Mol, Index"      
      do i=1,nMolTypes
        do j=1,NMAX(i)
          cnt = cnt+1
!          JointArray(cnt)%molType = i          
          MolArray(i)%mol(j)%indx = cnt
!          write(35,*) i, j, cnt
!          do jj=1,nAtoms(i)
!            JointArray(cnt)%x(jj) => MolArray(i)%mol(j)%x(jj)
!            JointArray(cnt)%y(jj) => MolArray(i)%mol(j)%y(jj)
!            JointArray(cnt)%z(jj) => MolArray(i)%mol(j)%z(jj)
!          enddo
        enddo
      enddo

      end subroutine
!=============================================================================            
         
      subroutine ReadInitialConfiguration
      use ParallelVar
      use SimParameters
      use ForceField
      use Coords
      implicit none
      integer :: iType,iMol,iAtom,iIndx
      real(dp) :: x,y,z
      character(len=2) :: atmSymbol
      character(len=100) :: format_string ,fl_name, out1

      if(multipleInput) then
        if (myid .lt. 10) then
          format_string = "(A,I1,A)"
        elseif(myid .lt. 100) then
          format_string = "(A,I2,A)"
        elseif(myid .lt. 1000) then
          format_string = "(A,I3,A)"
        elseif(myid .lt. 1000) then
          format_string = "(A,I4,A)"          
        else
          format_string = "(A,I5,A)"      
        endif      
        write(fl_name,format_string) "configuration", myid,".dat"      
        open( unit=10, file=trim(adjustl(fl_name)), status = "Old" ) 
      else
        open( unit=10, file="configuration.dat", status = "Old")
      endif

      read(10,*) (NPART(iType),iType=1,nMolTypes)
      read(10,*)            
!     This block ensures the initial configuration is within the min/max bounds set in the input parameter file      
      do iType = 1, nMolTypes
        if(NPART(iType) .gt. NMAX(iType)) then
          write(nout,*) "Error! Number of particles in the initial"
          write(nout,*) "configuration are greater than the maxmium bounds set in the"
          write(nout,*) "input parameters."          
          write(nout,*) "Molecule Type:", iType
          write(nout,*) "Initial Size:", NPART(iType)          
          write(nout,*) "Maximum Allowed:", NMAX(iType)
          stop "Error! Configuration above maximum particle bounds!"
        endif
      enddo
      do iType = 1, nMolTypes
        if(NPART(iType) .lt. NMIN(iType)) then
          write(nout,*) "Error! Number of particles in the initial"
          write(nout,*) "configuration are less than the minimum bounds set in the"
          write(nout,*) "input parameters."          
          write(nout,*) "Molecule Type:", iType
          write(nout,*) "Initial Size:", NPART(iType)          
          write(nout,*) "Minimum Allowed:", NMIN(iType)
          stop "Error! Configuration below minimum particle bounds!"
        endif
      enddo      

      do iType=1,nMolTypes
        do iMol=1,NPART(iType)
          do iAtom=1,nAtoms(iType)
            read(10,*) atmSymbol, x, y, z
            MolArray(iType)%mol(iMol)%x(iAtom) = x
            MolArray(iType)%mol(iMol)%y(iAtom) = y
            MolArray(iType)%mol(iMol)%z(iAtom) = z
          enddo
        enddo
      enddo
      close(10)

     
      NTotal = sum(NPART)
      isActive = .false.
      do iType = 1,nMolTypes
        do iMol = 1,NPART(iType)
          iIndx = MolArray(iType)%mol(iMol)%indx
          isActive(iIndx) = .true.
        enddo
      enddo
      
      
      end subroutine
!========================================================            
      subroutine ReadInitialGasPhase
      use SimParameters
      use ForceField
      use Units
      use Coords
      implicit none
      integer :: iType,iAtom
      character(len=2) :: atmSymbol
      real(dp) :: x1, y1, z1
      
      open(unit=10, file="gasPhase.dat", status = "Old")
      do iType = 1,nMolTypes
        read(10,*)      
        do iAtom = 1,nAtoms(iType)
          read(10,*) atmSymbol, gasConfig(iType)%x(iAtom), gasConfig(iType)%y(iAtom), gasConfig(iType)%z(iAtom)
        enddo
      enddo
      
      close(10)


!     Shift the coordinates such that the first atom is located at (0,0,0) for convinence.  
      do iType = 1,nMolTypes
        x1 = gasConfig(iType)%x(1)
        y1 = gasConfig(iType)%y(1)
        z1 = gasConfig(iType)%z(1)
        do iAtom = 1,nAtoms(iType)
          gasConfig(iType)%x(iAtom) = gasConfig(iType)%x(iAtom) - x1
          gasConfig(iType)%y(iAtom) = gasConfig(iType)%y(iAtom) - y1
          gasConfig(iType)%z(iAtom) = gasConfig(iType)%z(iAtom) - z1
        enddo
      enddo

      
      end subroutine         
!========================================================            
      subroutine RecenterCoordinates
      use SimParameters
      use ForceField
      use Coords
      implicit none
      integer :: iType,iMol,iAtom
      real(dp) :: xcm, ycm, zcm
      
      xcm = MolArray(1)%mol(1)%x(1)
      ycm = MolArray(1)%mol(1)%y(1)
      zcm = MolArray(1)%mol(1)%z(1)
      
      do iType=1,nMolTypes
        do iMol=1,NPART(iType)
          do iAtom=1,nAtoms(iType)
            MolArray(iType)%mol(iMol)%x(iAtom) =  MolArray(iType)%mol(iMol)%x(iAtom) - xcm
            MolArray(iType)%mol(iMol)%y(iAtom) =  MolArray(iType)%mol(iMol)%y(iAtom) - ycm
            MolArray(iType)%mol(iMol)%z(iAtom) =  MolArray(iType)%mol(iMol)%z(iAtom) - zcm
          enddo
        enddo
      enddo

      
      end subroutine
!========================================================            
      subroutine Generate_UnitSphere(x,y,z)
      use VarPrecision
      implicit none
      real(dp), intent(out) :: x,y,z
      real(dp) :: u_12_sq, u1, u2, grnd
      
      u_12_sq = 2E0
      do while(u_12_sq .ge. 1)
       u1 = 2E0 * grnd() - 1E0
       u2 = 2E0 * grnd() - 1E0
       u_12_sq = u1 * u1 + u2 * u2
      enddo
 
      x = 2E0 * u1 * sqrt(1E0 - u_12_sq)
      y = 2E0 * u2 * sqrt(1E0 - u_12_sq)
      z = (1E0 - 2E0 * u_12_sq)
      
      
      end subroutine
!==========================================================================     
!     The purpose of this function is to generate a random position for an atom
!     given a fixed bond angle (bond_ang) with a given vector (v1) and a given distance (r2).
!     The torsional rotational angle is chosen randomly from 0 to 2pi.
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vector(v2) is calculated using a rotational matrix
 
      subroutine Generate_UnitCone(v1,r2,bond_ang,v2)
      use Constants      
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1
      type(SimpleAtomCoords), intent(out) :: v2
      real(dp), intent(in) :: bond_ang, r2
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj,grnd

      r1 = v1%x*v1%x + v1%y*v1%y + v1%z*v1%z
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()
        
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)      
      r_proj = sqrt(v1%x*v1%x + v1%y*v1%y)
        
      coeff1 = (r2/r1)*cos(bond_ang)
      coeff2 = (r2/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r1

      v2%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v2%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v2%z = coeff1*v1%z                      + coeff3*s_term*(r_proj*r_proj)

         
      end subroutine
	  
!==========================================================================     
!     The purpose of this function is to generate two random position for two atoms
!     given two fixed bond angles (bond_ang1 and bond_ang2), one dihedral angle (dihed) which 
!     is the angle between two planes made by two angles with a given vector (v1) and two 
!     given distances (r2 and r3).
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vectors(v2 and v3) is calculated using a rotational matrix

      subroutine Generate_UnitPyramid(v1, r2, r3, bond_ang1, bond_ang2, dihed, v2, v3)
      use Constants      
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1
      type(SimpleAtomCoords), intent(out) :: v2, v3
      real(dp), intent(in) :: bond_ang1, bond_ang2, dihed, r2, r3
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj,grnd

      r1 = v1%x*v1%x + v1%y*v1%y + v1%z*v1%z
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()
        
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)      
      r_proj = sqrt(v1%x*v1%x + v1%y*v1%y)
        
      coeff1 = (r2/r1)*cos(bond_ang1)
      coeff2 = (r2/r_proj)*sin(bond_ang1)
      coeff3 = coeff2/r1

      v2%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v2%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v2%z = coeff1*v1%z                      + coeff3*s_term*(r_proj*r_proj)
	  
      tors_angle = tors_angle + dihed
      s_term = sin(tors_angle)
      c_term = cos(tors_angle) 

      coeff1 = (r3/r1)*cos(bond_ang2)
      coeff2 = (r3/r_proj)*sin(bond_ang2)
      coeff3 = coeff2/r1
	  
      v3%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v3%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v3%z = coeff1*v1%z                      + coeff3*s_term*(r_proj*r_proj)
         
      end subroutine
	  
!==========================================================================     
!     The purpose of this function is similar to the UnitCone function however
!     in this case the torsional angle is fixed in addition to the distance and
!     bond angle which means there is exactly a single point in space which
!     corresponds to these constraints.
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x2,y2,z2)   w2=(-y2,x2,0)  w3=(-x2*z2, -y2*z2, x2^2 + y2^2)
      subroutine Generate_UnitTorsion(v1,v2,r3,bond_ang,tors_angle,v3)
      use Constants      
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1,v2
      type(SimpleAtomCoords), intent(out) :: v3
      real(dp), intent(in) :: bond_ang, tors_angle,r3
      real(dp) :: x1_u, y1_u, z1_u
      real(dp) :: x1_s, y1_s, z1_s
      real(dp) :: r2, rot_angle
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj

      r2 = v2%x*v2%x + v2%y*v2%y + v2%z*v2%z
      r2 = sqrt(r2)
      r_proj = sqrt(v2%x*v2%x + v2%y*v2%y)

!            
      x1_u =  v1%x - v2%x
      y1_u =  v1%y - v2%y
      z1_u =  v1%z - v2%z

!     Calculate the v1 vector's w2 and w3 components in the new orthonormal framework.
!      x1_s =  ( v2%x*x1_u + v2%y*y1_u + v2%z * z1_u)/(r2)
      y1_s =  (-v2%y*x1_u + v2%x*y1_u) / r_proj
      z1_s =  (-v2%x*v2%z*x1_u - v2%y*v2%z*y1_u + r_proj*r_proj*z1_u) / (r_proj*r2)


        
!     Calculate the torsional rotation angle for the new v3 vector from the v1 components
      rot_angle = atan2(z1_s, y1_s)
!      write(2,*) rot_angle*180E0/pi, tors_angle*180E0/pi, (rot_angle + tors_angle)*180E0/pi
      rot_angle = rot_angle + tors_angle

!     Rescale the angle      
!      do while(rot_angle .lt. 0E0) 
!        rot_angle = rot_angle + two_pi
!      enddo
!      do while(rot_angle .gt. two_pi) 
!        rot_angle = rot_angle - two_pi
!      enddo
        
      s_term = sin(rot_angle)
      c_term = cos(rot_angle)      

        
      coeff1 = (r3/r2)*cos(bond_ang)
      coeff2 = (r3/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r2

      v3%x = coeff1*v2%x - coeff2*c_Term*v2%y - coeff3*s_Term*v2%x*v2%z
      v3%y = coeff1*v2%y + coeff2*c_Term*v2%x - coeff3*s_Term*v2%y*v2%z
      v3%z = coeff1*v2%z                      + coeff3*s_Term*(r_proj*r_proj)

         
      end subroutine
!==========================================================================           
      subroutine GenerateBondLength(dist, k_bond, r_eq, ProbGen)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use AcceptRates, only: distGen_accpt, distGen_atmp
      implicit none
      logical acpt
      real(dp), intent(in) :: k_bond,r_eq
      real(dp), intent(out) :: dist  
      real(dp), intent(out), optional :: ProbGen   
      real(dp) :: normalSigma, grnd,eng
!      real(dp) ::

      if(k_bond .eq. 0E0) then
        dist = r_eq
        ProbGen = 1E0
        return
      endif
      
      normalSigma = 1E0/(0.5E0*k_bond*beta)
      acpt=.false.         
      do while(.not. acpt)
         distGen_atmp = distGen_atmp + 1E0
        dist = -log(1E0-grnd())
!        dist = sqrt(-2E0*log(grnd())) * cos(two_pi*grnd())
!        dist = dist*normalSigma + r_eq
        eng = Harmonic(dist, k_bond, r_eq)
        ProbGen = dist * dist * exp(-beta*eng+dist)
!        ProbGen = dist * dist
        if(ProbGen/3E0 .gt. grnd()) then        
          acpt=.true.
        endif
      enddo  
      distGen_accpt = distGen_accpt + 1E0   

      end subroutine
!==========================================================================           
      subroutine GenerateBendAngle(angle, bendType, ProbGen)
      use Constants
      use SimParameters
      use ForceField
      use ForceFieldFunctions
      use AcceptRates, only: angGen_accpt, angGen_atmp
      implicit none

      interface
        subroutine SelectBin_IntTable(integralTable, nBin)
          use VarPrecision
          implicit none
          real(dp), intent(in) :: integralTable(:)
          integer, intent(out) :: nBin
        end subroutine
      end interface

      logical :: acpt
      integer :: nSel, dummy
      integer, intent(in) :: bendType
      real(dp), intent(out) :: angle
      real(dp), intent(out), optional :: ProbGen
      real(dp) :: k_bend,theta_eq
      real(dp) :: grnd, eng, ranNum, sumInt
      real(dp) :: ProbSel

      k_bend = bendData(bendType)%k_eq
      theta_eq = bendData(bendType)%ang_eq  
      if(k_bend .eq. 0E0) then
         angle = theta_eq
         ProbGen = 1E0
         return
      endif
      
      
      acpt=.false.  
      do while(acpt .eqv. .false.)
         angGen_atmp = angGen_atmp + 1E0

!         angle = pi*grnd()
!         angle = acos(1E0-2E0*grnd())

         ranNum = grnd()
          !Since the majority of the probability density will be centered near the equilibrium angle
          !we can use this fact to find the correct bin with much fewer comparisons by simply starting
          !near the beginning of the gaussian curve unless the random number chosen is sufficiently small.
         if(ranNum .lt. startProb) then
           nSel = 1
         else 
           nSel = bendData(bendType)%startBin
         endif
         do while(bendData(bendType)%Prob(nSel) .lt. ranNum)
           nSel = nSel + 1
         enddo
 
          !Now that the bin has been chosen, select an angle uniformly from the bin and calculate the
          !acceptance probability.
         angle = ( dble(nSel)-grnd() ) * bendBinWidth
         if(nSel .ne. 1) then
           ProbSel = bendData(bendType)%Prob(nSel) - bendData(bendType)%Prob(nSel-1)
         else
           ProbSel = bendData(bendType)%Prob(nSel)
         endif

         eng = Harmonic(angle, k_bend, theta_eq)
         ProbGen = sin(angle)*exp(-beta*eng)/ProbSel
         ProbGen = ProbGen*bendData(bendType)%accptConstant
         if(ProbGen .gt. grnd()) then
           acpt=.true.
         endif
      enddo    
      angGen_accpt = angGen_accpt + 1E0   

      end subroutine
!==========================================================================  
      subroutine GenerateTwoBranches(ang1, ang2, dihedral, dihedType, bendType1, bendType2, bendType3, ProbGen)  
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use ForceField
      use AcceptRates, only: dihedGen_accpt, dihedGen_atmp
      use CBMC_Variables, only: diBinSize, dihedData
      use RandomTools
      implicit none

      logical acpt
      integer, intent(in) :: bendType1, bendType2, bendType3, dihedType
      real(dp), intent(out) :: ang1, ang2, dihedral
      real(dp), intent(out), optional :: ProbGen
      integer :: bin, nSel, cnt
      real(dp) :: ProbTemp, ProbDihed
      real(dp) :: k_bend1,theta_eq1,k_bend2,theta_eq2,k_bend3,theta_eq3
      real(dp) :: grnd, eng, ang3
      real(dp) :: ranNum, sumInt, weight
	  
      k_bend1 = bendData(bendType1)%k_eq
      theta_eq1 = bendData(bendType1)%ang_eq
      k_bend2 = bendData(bendType2)%k_eq
      theta_eq2 = bendData(bendType2)%ang_eq
      k_bend3 = bendData(bendType3)%k_eq
      theta_eq3 = bendData(bendType3)%ang_eq
	
      if((k_bend1 .eq. 0E0) .or. (k_bend2 .eq. 0E0) .or. (k_bend3 .eq. 0E0)) then
        ang1 = theta_eq1
        ang2 = theta_eq2
        ang3 = theta_eq3
        dihedral = acos((cos(ang3) - cos(ang1)*cos(ang2))/(sin(ang1)*sin(ang2)))
        ProbGen = 1E0
        return
      endif
      
      cnt = 0
      acpt=.false.         
      do while(acpt .eqv. .false.)
         call GenerateBendAngle(ang1, bendType1, ProbTemp)
         call GenerateBendAngle(ang2, bendType2, ProbTemp)
!         ang1 = pi*grnd()
!         ang2 = pi*grnd()
!         ang1 = acos(1E0-2E0*grnd())
!         ang2 = acos(1E0-2E0*grnd())

         dihedral = two_pi*grnd()


         ang3 = cos(ang1)*cos(ang2) + sin(ang1)*sin(ang2)*cos(dihedral)
         if (ang3 .ge. 1E0) then
           ang3 = 1E0
         elseif (ang3 .le. -1E0) then
           ang3 = -1E0
         endif
         ang3 = acos(ang3)
!         eng = Harmonic(ang1, k_bend1, theta_eq1)
!         eng = eng + Harmonic(ang2, k_bend2, theta_eq2)
!         eng = eng + Harmonic(ang3, k_bend3, theta_eq3)
         eng = Harmonic(ang3, k_bend3, theta_eq3)
         dihedGen_atmp = dihedGen_atmp + 1E0
!         ProbGen = exp(-beta*eng)
         ProbGen = exp(-beta*eng)
         if( ProbGen .gt. grnd() ) then
           acpt = .true.
         endif
      enddo    
      bin = floor(dihedral/diBinSize)
      dihedData(dihedType)%Hist(bin) = dihedData(dihedType)%Hist(bin) + 1E0
      dihedGen_accpt = dihedGen_accpt + 1E0   


      end subroutine
!========================================================================================
      subroutine GenerateTorsAngle(angle, torsType, ProbGen)
      use Constants
      use SimParameters
      use ForceField
      use ForceFieldFunctions
      implicit none
      logical :: acpt
      integer, intent(in) :: torsType
      real(dp), intent(out) :: angle
      real(dp), intent(out), optional :: ProbGen    
      real(dp) :: grnd, eng
      
      acpt=.false.         
      do while(acpt .eqv. .false.)
         angle = two_pi*grnd()
         eng = Trappe_CosNx(angle, torsData(torsType)%a)
         ProbGen = exp(-beta*eng)
         if(0.9999E0*ProbGen .gt. grnd()) acpt=.true.
      enddo    

      end subroutine      
!==========================================================================           
      subroutine FindBond(nType,mem1, mem2, bondType)
      use Constants
      use SimParameters
      use ForceField
      implicit none
      integer, intent(in) :: nType, mem1, mem2
      integer, intent(out) :: bondType
      integer :: iBond
      
      
      do iBond = 1, nBonds(nType)
        if(any(bondArray(nType,iBond)%bondMembr(1:2) .eq. mem1)) then
          if(any(bondArray(nType,iBond)%bondMembr(1:2) .eq. mem2)) then        
             bondType = bondArray(nType,iBond)%bondType
             return
          endif
        endif
      enddo

      write(6,*) "Error! FindBond function unable to find a bond"
      write(6,*) "containing memebers:", mem1, mem2
      stop
      
      end subroutine
!==========================================================================           
      subroutine FindAngle(nType, mem1, mem2, mem3, bendType)
      use Constants
      use SimParameters
      use ForceField
      implicit none
      integer, intent(in) :: nType, mem1, mem2, mem3
      integer, intent(out) :: bendType
      integer :: iBend
      
      
      do iBend = 1, nAngles(nType)
        if(any(bendArray(nType,iBend)%bendMembr(1:3) .eq. mem1)) then
          if(any(bendArray(nType,iBend)%bendMembr(1:3) .eq. mem2)) then        
            if(any(bendArray(nType,iBend)%bendMembr(1:3) .eq. mem3)) then  
              bendType = bendArray(nType,iBend)%bendType
              return
            endif
          endif
        endif
      enddo

      write(*,*) "Error! FindAngle function unable to find an angle"
      write(*,*) "containing memebers:", mem1, mem2, mem3
      stop
      
      end subroutine  
!==========================================================================           
      subroutine FindTorsion(nType, mem1, mem2, mem3, mem4, torsType)
      use Constants
      use SimParameters
      use ForceField
      implicit none
      integer, intent(in) :: nType, mem1, mem2, mem3, mem4
      integer, intent(out) :: torsType
      integer :: iTors
      
      
      do iTors = 1, nTorsional(nType)
        if(any(torsArray(nType,iTors)%torsMembr(1:4) .eq. mem1)) then
          if(any(torsArray(nType,iTors)%torsMembr(1:4) .eq. mem2)) then        
            if(any(torsArray(nType,iTors)%torsMembr(1:4) .eq. mem3)) then  
              if(any(torsArray(nType,iTors)%torsMembr(1:4) .eq. mem4)) then  
                torsType = torsArray(nType,iTors)%torsType
                return
              endif
            endif
          endif
        endif
      enddo

      write(*,*) "Error! FindTorsional function unable to find a torsional angle"
      write(*,*) "containing memebers:", mem1, mem2, mem3, mem4
      stop
      
      end subroutine 
!==========================================================================           
      subroutine FindSingleDihedral(nType, hubIndx, dihedType)
      use CBMC_Variables
      implicit none
      integer, intent(in) :: nType, hubIndx
      integer, intent(out) :: dihedType
      integer :: iDihed
      
      
      do iDihed = 1, totalDihed
        if(dihedData(iDihed)%molType .eq. nType) then
          if(dihedData(iDihed)%hubIndx .eq. hubIndx) then
            dihedType = iDihed
            return
          endif         
        endif
      enddo

      write(*,*) "Error! FindSingleDihedral function unable to find a dihedral angle"
      write(*,*) nType, hubIndx
      stop
      
      end subroutine 
!========================================================================
      end module
