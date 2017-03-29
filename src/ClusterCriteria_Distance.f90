!=================================================================================
      module DistanceCriteria

      logical, allocatable :: ClusterMember(:), flipped(:) 
      contains
!=================================================================================     
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation. 
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation. 
      subroutine Detailed_DistanceCriteria(rejMove)
      use Coords
      use IndexingFunctions
      use ParallelVar
      use PairStorage
      use SimParameters
      implicit none     
      logical, intent(out) :: rejMove
!      real(dp), intent(inout) :: PairList(:, :)
      
!      logical :: ClusterMember(1:maxMol)
      integer :: h,cnt
      integer :: iType,jType, iMol, jMol, iIndx, jIndx
      integer :: globIndx1, globIndx2

      if(.not. allocated(ClusterMember) ) then
        allocate(ClusterMember(1:maxMol))
      endif
      if(.not. allocated(flipped) ) then
        allocate(flipped(1:maxMol))
      endif


      rejMove = .false.
      NeighborList = .false.
      if(NTotal .eq. 1) then
         return      
      endif

      
!      do iType=1,nMolTypes
!      do jType=iType,nMolTypes
!        do iMol = 1 ,NPART(iType)      
!        iIndx = MolArray(iType)%mol(iMol)%indx
!        do jMol = 1 ,NPART(jType) 
!          jIndx = MolArray(jType)%mol(jMol)%indx        
!          if(PairList(iIndx,jIndx) .le. Dist_Critr_sq ) then
!            NeighborList(iIndx,jIndx)=.true.         
!            NeighborList(jIndx,iIndx)=.true.          
!          endif
!        enddo
!        enddo
!      enddo
!      enddo

      do iType = 1, nMolTypes
        do iMol = 1 ,NPART(iType)      
          globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
          iIndx = MolArray(iType)%mol(iMol)%indx
          do jType = 1, nMolTypes
            do jMol = 1 ,NPART(jType) 
              globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
              jIndx = MolArray(jType)%mol(jMol)%indx
              if(jIndx .eq. iIndx) then
                cycle
              endif        
              if(rPair(globIndx1, globIndx2) % p % r_sq .le. Dist_Critr_sq ) then
!                write(*,*) globIndx1, globIndx2, rPair(globIndx1, globIndx2) % p % r_sq, rPair(globIndx1, globIndx2) % p % r
                NeighborList(iIndx, jIndx) = .true.         
                NeighborList(jIndx, iIndx) = .true.   
              else
                NeighborList(iIndx, jIndx) = .false.         
                NeighborList(jIndx, iIndx) = .false.   
              endif
            enddo
          enddo
        enddo
      enddo

      if( all(NeighborList .eqv. .false.) ) then
        write(nout,*) "------- Cluster Criteria Not Met!!!! --------"
        rejMove = .true.
        return
      endif



          
      do iType = 1, nMolTypes
        if(NPART(iType) .gt. 0) then
          iIndx = MolArray(iType)%mol(1)%indx
          ClusterMember(iIndx) = .true.
          exit
        endif
      enddo

      do h = 1, maxMol  
        do iIndx = 1, maxMol
          do jIndx = 1, maxMol
            if( NeighborList(iIndx, jIndx) ) then
              if( ClusterMember(iIndx) ) then
                ClusterMember(jIndx) = .true.
!                cnt = cnt + 1
              endif
              if( ClusterMember(jIndx) ) then
                ClusterMember(iIndx) = .true.
!                cnt = cnt + 1                
              endif
            endif
          enddo           
        enddo
      enddo     

      cnt = 0
      do iType = 1, nMolTypes
        if(NPART(iType) .lt. NMAX(iType)) then
          do iMol = NPART(iType)+1, NMAX(iType) 
            iIndx = MolArray(iType)%mol(iMol)%indx
            ClusterMember(iIndx) = .true.
            cnt = cnt + 1
          enddo
        endif   
      enddo
      
      if( any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
        write(nout,*) "------- Cluster Criteria Not Met!!!! --------"
      endif
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      subroutine Shift_DistanceCriteria(nIndx, rejMove)
      use Coords
      use IndexingFunctions
      use PairStorage
      use SimParameters     
      implicit none     
      
      logical, intent(out) :: rejMove      
!      real(dp), intent(in) :: PairList(:)
      integer,intent(in) :: nIndx
      
      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer :: iIndx, jIndx, h
      integer :: nType, i, jType, jMol, globIndx2
      integer :: curNeigh(1:60), neiMax
      
      rejMove=.false.
      if(NTotal .eq. 1) return        
      ClusterMember=.false.
      flipped=.false.
      
      nType = Get_MolType(nIndx,NMAX)
     
!     This section dermines which molecules are neighbored with the new trial position.  In the event
!     that the molecule's new location has no neghibors all further calcualtions are skipped and the move is
!     rejected.
   
      memberAdded = .false.
      do jType = 1, nMolTypes
        do jMol = 1, NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx
          if(nIndx .eq. jIndx) then
            cycle
          endif
          globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
          if(rPairNew(globIndx2)%p%r_sq .le. Dist_Critr_sq) then
            ClusterMember(jIndx) = .true.        
            memberAdded = .true.
          endif
        enddo
      enddo

      if(.not. memberAdded) then      
        rejMove = .true.
        return     
      endif    

!      This part of the code tabulates all the neighbors located around the particle's old position. 
      neiMax = 0
      curNeigh = 0
      do iIndx = 1, maxMol
        if(NeighborList(iIndx, nIndx)) then
          if(iIndx .ne. nIndx) then
            if(isActive(iIndx)) then
              neiMax = neiMax + 1
              curNeigh(neiMax) = iIndx
            endif
          endif      
        endif
      enddo      
      
      
!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.      
      neiFlipped = .true.
      do iIndx = 1, neiMax
        if(.not. clusterMember(curNeigh(iIndx))) then
          neiFlipped = .false.
          exit
        endif
      enddo

      if(neiFlipped) then
        rejMove = .false.
        return
      endif
      
    
      do h = 1, NTotal
!        cnt = 0
        memberAdded = .false.
        do iIndx = 1, maxMol
          if(ClusterMember(iIndx) .neqv. flipped(iIndx)) then
            do jIndx = 1, maxMol
              if(NeighborList(iIndx,jIndx)) then
                if(jIndx .ne. nIndx) then
                  ClusterMember(jIndx)=.true.   
                  memberAdded = .true.
                endif
              endif
            enddo        
            flipped(iIndx)=.true.
          endif
        enddo
 
        neiFlipped = .true.
        do i = 1, neiMax
          if(.not. clusterMember(curNeigh(i))) then
            neiFlipped = .false.
            exit
          endif
        enddo        
        if( neiFlipped ) then
          exit
        else 
          if(.not. memberAdded) then
            exit
          endif           
        endif
      enddo
  
       if( .not. neiFlipped ) then
         rejMove=.true.
       endif
     
      end subroutine

!=================================================================================     
!     This function determines if removing a particle from the cluster will result in the destruction of the cluster criteria. 
      subroutine SwapOut_DistanceCriteria(nSwap, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap
      integer :: iIndx, jIndx, h, cnt



      rejMove = .false.
      ClusterMember = .false.
      flipped = .false.

      do iIndx = 1, maxMol
        if( .not. isActive(iIndx) ) then      
          ClusterMember(iIndx) = .true.
          flipped(iIndx) = .true.
        endif
      enddo

!      In order to initialize the cluster criteria search, a single particle must be chosen as the starting point.
      do iIndx = 1, maxMol
        if( isActive(iIndx) ) then      
          if( nSwap .ne. iIndx ) then
            ClusterMember(iIndx) = .true.
            exit
          endif
        endif
      enddo

      cnt = 0
      do iIndx = 1, maxMol
        if( isActive(iIndx) .eqv. .false. ) then
          ClusterMember(iIndx) = .true.      
          flipped(iIndx) = .true.         
          cnt = cnt + 1
        endif
      enddo
      
      do h = 1, maxMol
        do iIndx = 1, maxMol
          if( iIndx .ne. nSwap ) then     
            if( ClusterMember(iIndx) .neqv. flipped(iIndx) ) then
              do jIndx = 1, maxMol
                if( NeighborList(iIndx, jIndx) ) then
                  ClusterMember(jIndx) = .true. 
                  cnt = cnt + 1
                endif
              enddo
              flipped(iIndx) = .true.
            endif
          endif
        enddo
        if(cnt .eq. maxMol-1) then
          exit
        endif         
      enddo
  
  
      ClusterMember(nSwap) = .true.
       
      if( any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
      endif      
     
      end subroutine
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine NeighborUpdate_Distance(nIndx)
      use Coords      
      use IndexingFunctions      
      use PairStorage
      use SimParameters
      implicit none     
      integer, intent(in) :: nIndx
!      real(dp), intent(in) :: PairList(:)
      integer :: nType, nMol, jMol, jType, jIndx, globIndx1, globIndx2

      nType = typeList(nIndx)
      nMol = subIndxList(nIndx)
      globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1) 
!      write(35,*) nIndx
      do jIndx=1,maxMol
        if(.not. isActive(jIndx)) then
          cycle
        endif
        if(jIndx .ne. nIndx) then  
          jType = typeList(jIndx)
          jMol = subIndxList(jIndx)
          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1) 
!          write(35,*) globIndx1, globIndx2, rPairNew(globIndx2) % p % r_sq, Dist_Critr_sq
!          if( rPair(globIndx1,globIndx2) % p % r_sq  .le.  Dist_Critr_sq ) then
          if( rPairNew(globIndx2) % p % r_sq  .le.  Dist_Critr_sq ) then
!            write(35,*) globIndx1, globIndx2, rPair(globIndx1,globIndx2) % p % r_sq 
            NeighborList(nIndx, jIndx) = .true.
            NeighborList(jIndx, nIndx) = .true.  
          else             
            NeighborList(nIndx, jIndx) = .false.
            NeighborList(jIndx, nIndx) = .false.            
          endif   
        endif
      enddo
!      NeighborList(nIndx,nIndx) = .false.


!      do jType = 1, nMolTypes
!        do jMol = 1, NPART(jType)
!          jIndx = MolArray(jType)%mol(jMol)%indx   
!          if(jIndx .ne. nIndx) then  
!            globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1) 
!            if( rPair(globIndx1,globIndx2) % p % r_sq  .le.  Dist_Critr_sq ) then
!              write(2,*) globIndx1, globIndx2, rPair(globIndx1,globIndx2) % p % r_sq 
!              NeighborList(nIndx, jIndx) = .true.
!              NeighborList(jIndx, nIndx) = .true.  
!            else             
!              NeighborList(nIndx, jIndx) = .false.
!              NeighborList(jIndx, nIndx) = .false.            
!            endif   
!          endif
!        enddo
!      enddo
!      NeighborList(nIndx,nIndx) = .false.
!      write(2,*) 

      end subroutine
!=================================================================================           
      end module
      
