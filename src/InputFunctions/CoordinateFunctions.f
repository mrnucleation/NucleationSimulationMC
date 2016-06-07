!========================================================================
      subroutine AllocateCoordinateArrays
      use ForceField
      use SimParameters
      use Coords
      implicit none
      integer :: i,j,AllocationStatus
      
      allocate( MolArray(1:nMolTypes),
     &           stat=AllocationStatus)
      do i=1,nMolTypes
        allocate( MolArray(i)%mol(1:NMAX(i)) ,
     &            stat=AllocationStatus)     
        do j=1,NMAX(i)
         allocate( MolArray(i)%mol(j)%x(1:nAtoms(i)),
     &            stat=AllocationStatus)
         allocate( MolArray(i)%mol(j)%y(1:nAtoms(i)),
     &            stat=AllocationStatus)
         allocate( MolArray(i)%mol(j)%z(1:nAtoms(i)),
     &            stat=AllocationStatus)
        enddo
      enddo

      maxMol = 0
      do i=1, nMolTypes      
        maxMol = maxMol + NMAX(i)
      enddo
      
      maxAtoms = maxval(nAtoms)
      
      allocate(newMol%x(1:maxAtoms),
     &         stat=AllocationStatus)
      allocate(newMol%y(1:maxAtoms),
     &         stat=AllocationStatus)
      allocate(newMol%z(1:maxAtoms),
     &         stat=AllocationStatus)
     
C       allocate( JointArray(1:NMaxMol),
C      &           stat=AllocationStatus)
C       do i = 1,NMaxMol
C         allocate( JointArray(i)%x(1:maxAtoms),
C      &           stat=AllocationStatus)
C         allocate( JointArray(i)%y(1:maxAtoms),
C      &           stat=AllocationStatus)
C         allocate( JointArray(i)%z(1:maxAtoms),
C      &           stat=AllocationStatus)
C       enddo
      
      
      end subroutine
!========================================================            
      subroutine ReadInitialConfiguration
      use SimParameters
      use ForceField
      use Units
      implicit none
      integer :: i,j,jj
      character(len=2) :: atmSymbol
      
      open(unit=10, file="configuration.dat")
      read(10,*) (NPART(j),j=1,nMolTypes)
      read(10,*)
      do i=1,nMolTypes
        do j=1,NPART(i)
          do jj=1,nAtoms(i)
            read(10,*) atmSymbol, MolArray(i)%mol(j)%x(jj),
     &                            MolArray(i)%mol(j)%y(jj),
     &                            MolArray(i)%mol(j)%z(jj)
          enddo
        enddo
      enddo
      
      close(10)
      
      
      end subroutine
!========================================================================