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
      
!      logical :: ClusterMember(1:maxMol)      
!      logical :: flipped(1:maxMol)
      logical :: earlyExit
      integer :: iType, iMol
      integer :: iIndx, jIndx, h, cnt
      integer :: nType, jType, jMol
      integer :: globIndx2
      
      rejMove = .false.
      if(NTotal .eq. 1) then
        return
      endif
      ClusterMember = .false.
      flipped = .false.
      
      nType = Get_MolType(nIndx,NMAX)
      
!     This block performs a quick check to see if any neighbors were lost
!     in the process of making this move.  If all neighbors were retained
!     then the cluster remains in tact and the detailed calculations are no
!     longer needed. 

!      earlyExit = .true.
!      do jType = 1, nMolTypes     
!        do jMol = 1, NPART(jType)
!          jIndx = molArray(jType)%mol(jMol)%indx
!          if( jIndx .ne. nIndx ) then
!            if( NeighborList(nIndx, jIndx) ) then
!              globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
!              if( rPairNew(globIndx2) % p % r_sq .gt. Dist_Critr_sq ) then
!                earlyExit = .false. 
!                exit            
!              endif
!            endif
!          endif 
!        enddo
!        if(.not. earlyExit) then
!          exit
!        endif
!      enddo

!      if(earlyExit) then
!        write(*,*) sqrt(rPairNew(globIndx2) % p % r_sq)
!        write(*,*) "Early Exit"
!        return
!      endif

      
!     This section dermines which molecules are neighbored with the new trial position. These molecules are then activated as cluster members.
      cnt = 0
      do jType = 1, nMolTypes     
        do jMol = 1, NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx
          if(jIndx .ne. nIndx) then  
            globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
            if( rPairNew(globIndx2) % p % r_sq .le. Dist_Critr_sq ) then
              ClusterMember(jIndx) = .true.        
              cnt = cnt + 1
            endif
          endif
        enddo   
      enddo
      
      
!     This section checked to see if there were any neighbors for the molecule's new position.
!     If cnt is equal to 0 the new position has no neighbors which implies the cluster is broken.
      if(cnt .eq. 0) then
        rejMove = .true.
        return     
      endif    

      do iType = 1, nMolTypes     
        do iMol = NPART(iType)+1, NMAX(iType)
          iIndx = molArray(iType)%mol(iMol)%indx
          ClusterMember(iIndx) = .true.
          flipped(iIndx) = .true.
        enddo
      enddo

      ClusterMember(nIndx) = .true.        
      flipped(nIndx) = .true.     
      cnt = cnt + 1
      
      do h = 1, NTotal
        do iIndx = 1, maxMol
          if(iIndx .eq. nIndx) then
            cycle
          endif
          if( isActive(iIndx) ) then
            if( ClusterMember(iIndx) .neqv. flipped(iIndx) ) then
              do jIndx = 1, maxMol
                if(.not. isActive(jIndx)) then
                  cycle
                endif
                if(jIndx .eq. iIndx) then
                  cycle
                endif
                if(jIndx .eq. nIndx) then
                  cycle
                endif
                if( NeighborList(iIndx, jIndx) ) then
                  ClusterMember(jIndx) = .true.            
                  cnt = cnt+1
                endif
              enddo        
              flipped(iIndx)=.true.
            endif
          endif
        enddo
        if(cnt .eq. maxMol) then
          exit
        endif      
      enddo

     
  
      if( any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
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
      do iIndx= 1, maxMol
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
      do jType = 1, nMolTypes
        do jMol = 1, NPART(jType)
          jIndx = MolArray(jType)%mol(jMol)%indx   
          if(jIndx .ne. nIndx) then  
            globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1) 
            if( rPair(globIndx1,globIndx2) % p % r_sq  .le.  Dist_Critr_sq ) then
!              write(2,*) globIndx1, globIndx2, rPair(globIndx1,globIndx2) % p % r_sq 
              NeighborList(nIndx, jIndx) = .true.
              NeighborList(jIndx, nIndx) = .true.  
            else             
              NeighborList(nIndx, jIndx) = .false.
              NeighborList(jIndx, nIndx) = .false.            
            endif   
          endif
        enddo
      enddo
      NeighborList(nIndx,nIndx) = .false.
!      write(2,*) 

      end subroutine
!=================================================================================           
      end module
      
