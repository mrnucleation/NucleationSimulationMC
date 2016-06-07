!========================================================================
      subroutine AllocateCoordinateArrays
      use ForceField
      use SimParameters
      use Coords
      use EnergyTables
      implicit none
      integer :: i,j,cnt,AllocationStatus
      
      allocate( MolArray(1:nMolTypes), stat=AllocationStatus)
      do i=1,nMolTypes
        allocate( MolArray(i)%mol(1:NMAX(i)), stat=AllocationStatus)   
        do j=1,NMAX(i)
          allocate( MolArray(i)%mol(j)%x(1:nAtoms(i)), stat=AllocationStatus)
          allocate( MolArray(i)%mol(j)%y(1:nAtoms(i)), stat=AllocationStatus)
          allocate( MolArray(i)%mol(j)%z(1:nAtoms(i)), stat=AllocationStatus)
        enddo
      enddo

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

!       allocate( JointArray(1:NMaxMol), stat=AllocationStatus)
!       do i = 1,NMaxMol
!         allocate( JointArray(i)%x(1:maxAtoms),
!      &           stat=AllocationStatus)
!         allocate( JointArray(i)%y(1:maxAtoms),
!      &           stat=AllocationStatus)
!         allocate( JointArray(i)%z(1:maxAtoms),
!      &           stat=AllocationStatus)
!       enddo
      
      
      end subroutine
!=============================================================================
      subroutine CreateJointArray
      use SimParameters
      use Coords
      implicit none
      integer :: i,j
      integer :: cnt
      
      
      cnt=0
      write(35,*) "Molecular Index List"
      write(35,*) "Type, Mol, Index"      
      do i=1,nMolTypes
        do j=1,NMAX(i)
          cnt = cnt+1
!          JointArray(cnt)%molType = i          
          MolArray(i)%mol(j)%indx = cnt
          write(35,*) i, j, cnt
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
      real(kind(0.0d0)) :: x,y,z
      character(len=2) :: atmSymbol
      
      open(unit=10, file="configuration.dat")
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
      real(kind(0.0d0)) :: x1, y1, z1
      
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
      real(kind(0.0d0)) :: xcm, ycm, zcm
      
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
      implicit none
      real(kind(0.0d0)), intent(out) :: x,y,z
      real(kind(0.0d0)) :: u_12_sq, u1, u2, grnd
      
      u_12_sq = 2d0
      do while(u_12_sq .ge. 1)
       u1 = 2d0 * grnd() - 1d0
       u2 = 2d0 * grnd() - 1d0
       u_12_sq = u1 * u1 + u2 * u2
      enddo
 
      x = 2d0 * u1 * dsqrt(1d0 - u_12_sq)
      y = 2d0 * u2 * dsqrt(1d0 - u_12_sq)
      z = (1d0 - 2d0 * u_12_sq)
      
      
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
      real(kind(0.0d0)), intent(in) :: bond_ang, r2
      real(kind(0.0d0)) :: tors_angle
      real(kind(0.0d0)) :: r1
      real(kind(0.0d0)) :: s_term,c_term
      real(kind(0.0d0)) :: coeff1,coeff2,coeff3  
      real(kind(0.0d0)) :: r_proj,grnd

      r1 = v1%x*v1%x + v1%y*v1%y + v1%z*v1%z
      r1 = dsqrt(r1)

      tors_angle = two_pi*grnd()
        
      s_term = dsin(tors_angle)
      c_term = dcos(tors_angle)      
      r_proj = dsqrt(v1%x*v1%x + v1%y*v1%y)
        
      coeff1 = (r2/r1)*dcos(bond_ang)
      coeff2 = (r2/r_proj)*dsin(bond_ang)
      coeff3 = coeff2/r1

      v2%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v2%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v2%z = coeff1*v1%z                      + coeff3*s_term*(r_proj*r_proj)

         
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
      real(kind(0.0d0)), intent(in) :: bond_ang, tors_angle,r3
      real(kind(0.0d0)) :: x1_u, y1_u, z1_u
      real(kind(0.0d0)) :: x1_s, y1_s, z1_s
      real(kind(0.0d0)) :: r2, rot_angle
      real(kind(0.0d0)) :: s_term,c_term
      real(kind(0.0d0)) :: coeff1,coeff2,coeff3  
      real(kind(0.0d0)) :: r_proj

      r2 = v2%x*v2%x + v2%y*v2%y + v2%z*v2%z
      r2 = dsqrt(r2)
      r_proj = dsqrt(v2%x*v2%x + v2%y*v2%y)

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
!      write(2,*) rot_angle*180d0/pi, tors_angle*180d0/pi, (rot_angle + tors_angle)*180d0/pi
      rot_angle = rot_angle + tors_angle

!     Rescale the angle      
!      do while(rot_angle .lt. 0d0) 
!        rot_angle = rot_angle + two_pi
!      enddo
!      do while(rot_angle .gt. two_pi) 
!        rot_angle = rot_angle - two_pi
!      enddo
        
      s_term = dsin(rot_angle)
      c_term = dcos(rot_angle)      

        
      coeff1 = (r3/r2)*dcos(bond_ang)
      coeff2 = (r3/r_proj)*dsin(bond_ang)
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
      implicit none
      logical acpt
      real(kind(0.0d0)), intent(in) :: k_bond,r_eq
      real(kind(0.0d0)), intent(out) :: dist  
      real(kind(0.0d0)), intent(out), optional :: ProbGen   
      real(kind(0.0d0)) :: grnd,eng

      if(k_bond .eq. 0d0) then
         dist = r_eq
         ProbGen = 1d0
         return
      endif
        
      acpt=.false.         
      do while(.not. acpt)
        dist = -log(1d0-grnd())
        eng = Harmonic(dist, k_bond, r_eq)
        ProbGen = dist * dist * exp(-beta*eng+dist)
!        if(ProbGen/(3d0*exp(-dist)) .gt. grnd()) then
        if(ProbGen/3d0 .gt. grnd()) then        
          acpt=.true.
        endif
      enddo    
         
      end subroutine
!==========================================================================           
      subroutine GenerateBendAngle(angle,k_bend,theta_eq,ProbGen)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use AcceptRates, only: angGen_accpt, angGen_atmp
      implicit none
      logical acpt
      real(kind(0.0d0)), intent(out) :: angle
      real(kind(0.0d0)), intent(out), optional :: ProbGen
      real(kind(0.0d0)), intent(in) :: k_bend,theta_eq
      real(kind(0.0d0)) :: grnd, eng
      
      if(k_bend .eq. 0d0) then
         angle = theta_eq
         ProbGen = 1d0
         return
      endif
      
      
      acpt=.false.         
      do while(acpt .eqv. .false.)
         angGen_atmp = angGen_atmp + 1d0
!         angle = pi*grnd()
         angle = acos(1d0-2d0*grnd())
         eng = Harmonic(angle, k_bend, theta_eq)
!         ProbGen = dsin(angle)*exp(-beta*eng)
         ProbGen = exp(-beta*eng)
         if(0.9999d0*ProbGen .gt. grnd()) acpt=.true.
      enddo    
      angGen_accpt = angGen_accpt + 1d0   
      end subroutine

!==========================================================================           
      subroutine GenerateTorsAngle(angle, torsType, ProbGen)
      use Constants
      use SimParameters
      use ForceField
      use ForceFieldFunctions
      implicit none
      logical :: acpt
      integer, intent(in) :: torsType
      real(kind(0.0d0)), intent(out) :: angle
      real(kind(0.0d0)), intent(out), optional :: ProbGen    
      real(kind(0.0d0)) :: grnd, eng
      
      acpt=.false.         
      do while(acpt .eqv. .false.)
         angle = two_pi*grnd()
         eng = Trappe_CosNx(angle, torsData(torsType)%a)
         ProbGen = exp(-beta*eng)
         if(0.9999d0*ProbGen .gt. grnd()) acpt=.true.
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
!========================================================================
