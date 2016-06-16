!=================================================================================
      module EnergyCriteria
      contains
!=================================================================================     
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation. 
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation. 
      subroutine Detailed_EnergyCriteria(PairList,rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      use ParallelVar
      implicit none     
      logical,intent(out) :: rejMove
      real(kind(0.0d0)),intent(inout) :: PairList(1:maxMol,1:maxMol)
      
      logical :: ClusterMember(1:maxMol)
      integer :: i,j,h,cnt
      integer :: iType,jType, iMol, jMol, iIndx, jIndx


      rejMove = .false.
      NeighborList=.false.
      if(NTotal .eq. 1) then
         return      
      endif

      
      do iType=1,nMolTypes
      do jType=iType,nMolTypes
        do iMol = 1 ,NPART(iType)      
        iIndx = MolArray(iType)%mol(iMol)%indx
        do jMol = 1 ,NPART(jType) 
          jIndx = MolArray(jType)%mol(jMol)%indx        
          if(PairList(iIndx,jIndx) .le. Eng_Critr(iType,jType) ) then
            NeighborList(iIndx,jIndx)=.true.         
            NeighborList(jIndx,iIndx)=.true.          
          endif
        enddo
        enddo
      enddo
      enddo

      cnt = 0
      do i=1,maxMol
        if(isActive(i) .eqv. .false.) then
          ClusterMember(i) = .true.
          cnt = cnt + 1
        endif         
      enddo
      
      ClusterMember(1)=.true.      
      do h=1,maxMol  
        do i=1,maxMol
          do j=1,maxMol
            if(NeighborList(i,j)) then
              if(ClusterMember(i)) then
                ClusterMember(j)=.true.
!                cnt = cnt + 1
              endif
              if(ClusterMember(j)) then
                ClusterMember(i)=.true.
!                cnt = cnt + 1                
              endif
            endif
          enddo           
        enddo
      enddo     

      do i = 1, maxMol
        NeighborList(i,i) = .false.
      enddo
      
      if(any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
        write(nout,*) "------- Cluster Criteria Not Met! --------"
      endif
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      subroutine Shift_EnergyCriteria(PairList, nIndx, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(kind(0.0d0)), intent(in) :: PairList(1:maxMol)
      integer,intent(in) :: nIndx
      
      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer :: i,j,h,cnt
      integer :: nType, jType
      integer :: jlowerIndx      
      integer :: curNeigh(1:60), neiMax
      
      rejMove=.false.
      if(NTotal .eq. 1) return        
      ClusterMember=.false.
      flipped=.false.
      
      nType = Get_MolType(nIndx,NMAX)
     
!     This section dermines which molecules are neighbored with the new trial position

!      ClusterMember(nIndx) = .true.
!      flipped(nIndx) = .true.           
!      cnt = 1      
      memberAdded = .false.
      do j=1,maxMol
!        if(.not. isActive(j)) then
!          cycle
!        endif
!        if(j .ne. nIndx) then
          jType = typeList(j)          
          if(PairList(j) .le. Eng_Critr(jType,nType)) then
            ClusterMember(j) = .true.        
!            cnt = cnt + 1
            memberAdded = .true.
          endif
!        endif   
      enddo

!      jlowerIndx = 0
!      do jType = 1, nMolTypes
!       do j = jlowerIndx+1,jlowerIndx+NPART(jType)
!         if(PairList(j) .le. Eng_Critr(jType,nType)) then
!           if(j .ne. nIndx) then
!            ClusterMember(j) = .true.        
!            cnt = cnt + 1
!           endif
!         endif  
!        enddo
!        jlowerIndx = jlowerIndx + NMAX(jType)
!      enddo
      

      
!     This section checked to see if there were any neighbors for the molecule's new position.
!     If cnt is equal to 0 the new position has no neighbors which implies the cluster is broken.
!      if(cnt .eq. 1) then
      if(.not. memberAdded) then      
        rejMove = .true.
        return     
      endif    
      
      neiMax = 0
      curNeigh = 0
      do i=1,maxMol
        if(NeighborList(i,nIndx)) then
         if(i .ne. nIndx) then
           if(isActive(i)) then
            neiMax = neiMax + 1
            curNeigh(neiMax) = i
          endif
         endif      
        endif
      enddo      
      
      
!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.      
!      neiFlipped = .true.
!      do i=1,maxMol
!        if(.not. isActive(i)) then
!          cycle
!        endif
!        if(NeighborList(i,nIndx)) then
!          if(i .ne. nIndx) then
!            if(.not. ClusterMember(i)) then
!              neiFlipped = .false.
!              exit
!            endif
!          endif
!        endif
!      enddo

!      jlowerIndx = 0
!      do jType = 1, nMolTypes
!       do j = jlowerIndx+1,jlowerIndx+NPART(jType)
!         if(NeighborList(j,nIndx)) then
!           if(.not. ClusterMember(j)) then
!             neiFlipped = .false.
!             exit
!           endif
!         endif  
!        enddo
!        if(neiFlipped .eqv. .false.) then
!          exit
!        endif
!        jlowerIndx = jlowerIndx + NMAX(jType)
!      enddo
      neiFlipped = .true.
      do i = 1, neiMax
        if(.not. clusterMember(curNeigh(i))) then
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
        do i = 1, maxMol
!          if(isActive(i)) then
            if(ClusterMember(i) .neqv. flipped(i)) then
              do j=1,maxMol
!                if(.not. isActive(j)) then
!                  cycle
!                endif
                if(NeighborList(i,j)) then
                  if(j .ne. nIndx) then
                    ClusterMember(j)=.true.   
!                    cnt = cnt + 1                         
                    memberAdded = .true.
                  endif
                endif
              enddo        
              flipped(i)=.true.
            endif
!          endif
        enddo
!        if(h .gt. 1) then
 
        neiFlipped = .true.
        do i = 1, neiMax
          if(.not. clusterMember(curNeigh(i))) then
            neiFlipped = .false.
            exit
          endif
        enddo        

          
!        do i=1, maxMol
!          if(.not. isActive(i)) then
!            cycle
!          endif        
!          if( NeighborList(i,nIndx) ) then
!            if( .not. ClusterMember(i) ) then
!              neiFlipped = .false.
!              exit
!            endif
!          endif
!        enddo
!        jlowerIndx = 0
!        do jType = 1, nMolTypes
!         do j = jlowerIndx+1,jlowerIndx+NPART(jType)
!          if(NeighborList(j,nIndx)) then
!             if(.not. ClusterMember(j)) then
!               neiFlipped = .false.
!               exit
!             endif
!            endif  
!          enddo
!          if(neiFlipped .eqv. .false.) then
!            exit
!          endif
!          jlowerIndx = jlowerIndx + NMAX(jType)
!        enddo
          if( neiFlipped ) then
            exit
          else 
!            if(cnt .eq. 0) then
            if(.not. memberAdded) then
              exit
            endif           
          endif
!        endif
      enddo
  
       if( .not. neiFlipped ) then
         rejMove=.true.
       endif
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      pure subroutine SwapIn_EnergyCriteria(nType,PairList,rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(kind(0.0d0)), intent(in) :: PairList(1:maxMol)
      integer, intent(in) :: nType
      integer :: j, jType
      
      do j = 1, maxMol
        if(isActive(j)) then      
          jType = Get_MolType(j,NMAX)                 
          if(PairList(j) .le. Eng_Critr(nType,jType)) then
            rejMove = .false.
            return
          endif
        endif
      enddo
     
      rejMove = .true.
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      subroutine SwapOut_EnergyCriteria(nSwap,rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap

      logical :: neiFlipped, memberAdded
      logical :: clusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer :: i,j,h,cnt
      integer :: curNeigh(1:60), neiMax
      
      rejMove=.false.
      if(NTotal-1 .eq. 1) return  

      ClusterMember=.false.
      flipped=.false.
      
      do i=1,maxMol
        if(.not. isActive(i)) then
          cycle
        endif
        if( NeighborList(i,nSwap) ) then
          if(i .ne. nSwap) then
            clusterMember(i)=.true. 
            exit
          endif
        endif
      enddo

      neiMax = 0
      curNeigh = 0
      do i=1,maxMol
        if(.not. isActive(i)) then
          cycle
        endif
        if( NeighborList(i,nSwap) ) then
          if(i .ne. nSwap) then
            neiMax = neiMax + 1
            curNeigh(neiMax) = i
          endif
        endif        
      enddo
      
!      cnt = 0
!      do i=1,maxMol
!        if(isActive(i) .eqv. .false.) then
!          ClusterMember(i) = .true.      
!          flipped(i) = .true.         
!          cnt = cnt + 1
!        endif
!      enddo

     do h=1,NTotal
!        cnt = 0
        memberAdded = .false.
        do i=1,maxMol
!         if(.not. isActive(i)) then
!           cycle
!         endif
           if(ClusterMember(i) .neqv. flipped(i)) then
             if(i .ne. nSwap) then 
              do j=1,maxMol
!1                if(.not. isActive(j)) then
!                  cycle
!                endif
                if(NeighborList(i,j)) then
                  if(j .ne. nSwap) then
                    clusterMember(j) = .true. 
!                    cnt = cnt + 1                    
                    memberAdded = .true.
                  endif
                endif
              enddo
              flipped(i) = .true.
            endif
          endif
        enddo

!        if(h .gt. 1) then
          neiFlipped = .true.
          do i = 1, neiMax
            if(.not. clusterMember(curNeigh(i))) then
              neiFlipped = .false.
              exit
            endif
          enddo     
        
!        do i=1, maxMol
!          if( NeighborList(i,nSwap) ) then
!            if(i .ne. nSwap) then
!              if( .not. ClusterMember(i) ) then
!                neiFlipped = .false.
!                exit
!              endif
!            endif
!          endif
!        enddo
     
          if( neiFlipped ) then
            exit
          else
!            if(cnt .eq. 0) then
            if(.not. memberAdded) then            
              exit
            endif            
          endif
!        endif
     enddo
  
  
     if( .not. neiFlipped ) then
       rejMove = .true.
!     else
!       write(2,*) "Accepted!"
     endif
     
     end subroutine
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine NeighborUpdate(PairList, nIndx)
      use SimParameters
      use IndexingFunctions      
      use Coords      
      implicit none     
      integer iType,j,jType,nIndx
      real(kind(0.0d0)) :: PairList(1:maxMol)


!      do j=1,maxMol
!        if(.not. isActive(j)) cycle      
!        if(j .ne. nIndx) then
!            NeighborList(nIndx,j)=.true.
!            NeighborList(j,nIndx)=.true.  
!          else             
!            NeighborList(nIndx,j)=.false.
!            NeighborList(j,nIndx)=.false.            
!          endif   
!        endif
!      enddo


      iType = Get_MolType(nIndx,NMAX)  
      do j=1,maxMol
        if(.not. isActive(j)) then
            cycle
        endif
        if(j .ne. nIndx) then
          jType = Get_MolType(j,NMAX)        
          if(PairList(j) .le. Eng_Critr(iType,jType) ) then
            NeighborList(nIndx,j)=.true.
            NeighborList(j,nIndx)=.true.  
          else             
            NeighborList(nIndx,j)=.false.
            NeighborList(j,nIndx)=.false.            
          endif   
        endif
      enddo
      NeighborList(nIndx,nIndx) = .false.
      

      end subroutine
!=================================================================================     
!     This function updates the neighborlist if a move is accepted.
      subroutine NeighborUpdate_Delete(nIndx)
      use SimParameters
      use IndexingFunctions
      use Coords
      implicit none
      integer, intent(in) :: nIndx
      integer :: nType, nSwapIndx
      integer :: i

      nType = typeList(nIndx)     
      nSwapIndx = molArray(nType)%mol(NPART(nType))%indx

!      NeighborList(nIndx,:) = NeighborList(nSwapIndx,:)
!      NeighborList(:,nIndx) = NeighborList(:,nSwapIndx)

      if(nIndx .eq. nSwapIndx) then
        NeighborList(nSwapIndx,:) = .false.
        NeighborList(:,nSwapIndx) = .false.       
        return
      endif
     
      do i = 1, maxMol
        if(NeighborList(i,nSwapIndx)) then
          if(i .ne. nIndx) then
            NeighborList(i,nIndx) = .true.
            NeighborList(nIndx,i) = .true.            
          endif
        else
          NeighborList(i,nIndx) = .false.
          NeighborList(nIndx,i) = .false.
        endif        
      enddo
      
      NeighborList(nIndx,nIndx) = .false.

      NeighborList(nSwapIndx,:) = .false.
      NeighborList(:,nSwapIndx) = .false.      

      end subroutine      
!=================================================================================           
      end module
      
