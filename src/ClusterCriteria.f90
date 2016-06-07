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
      
      if(any(ClusterMember .eqv. .false.) ) then
        rejMove = .true.
        write(nout,*) "------- Cluster Criteria Not Met! --------"
      endif
     
      end subroutine
!=================================================================================     
!     This function determines if a given translational move will destroy a cluster. 
      pure subroutine Shift_EnergyCriteria(PairList, nIndx, rejMove)
      use SimParameters     
      use Coords
      use IndexingFunctions
      implicit none     
      
      logical, intent(out) :: rejMove      
      real(kind(0.0d0)), intent(in) :: PairList(1:maxMol)
      integer,intent(in) :: nIndx
      
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer :: i,j,h,cnt
      integer :: nType, jType
      
      rejMove=.false.
      if(NTotal .eq. 1) return        
      ClusterMember=.false.
      flipped=.false.
      
      nType = Get_MolType(nIndx,NMAX)
      
!     This block performs a quick check to see if any neighbors were lost
!     in the process of making this move.  If all neighbors were retained
!     then the cluster remains in tact and the detailed calculations are no
!     longer needed. 
      do j=1,maxMol+1
        if(j .eq. maxMol+1) then
          return      
        endif      
        if(j .ne. nIndx) then
          if(NeighborList(nIndx, j) ) then
!            jType = Get_MolType(j,NMAX)        
            jType = typeList(j)
            if(PairList(j) .gt. Eng_Critr(nType,jType)) then
              exit            
            endif
          endif
        endif 
      enddo
      

      
!     This section dermines which molecules are neighbored with the new trial position
      cnt = 0
      do j=1,maxMol
        if(j .ne. nIndx) then
!          jType = Get_MolType(j,NMAX)        
          jType = typeList(j)          
          if(PairList(j) .le. Eng_Critr(nType,jType)) then
            ClusterMember(j) = .true.        
            cnt = cnt + 1
          endif
        endif   
      enddo
      
      
!     This section checked to see if there were any neighbors for the molecule's new position.
!     If cnt is equal to 0 the new position has no neighbors which implies the cluster is broken.
      if(cnt .eq. 0) then
        rejMove = .true.
        return     
      endif    

!           
      do i=1,maxMol
        if(isActive(i) .eqv. .false.) then
         ClusterMember(i) = .true.
         flipped(i) = .true.
         cnt = cnt + 1
        endif
      enddo

      ClusterMember(nIndx) = .true.        
      flipped(nIndx) = .true.     
      cnt = cnt + 1
      
      do h=1,maxMol
        do i=1,maxMol
          if(isActive(i)) then
            if(ClusterMember(i) .neqv. flipped(i)) then
              do j=1,maxMol
                if(.not. isActive(j)) then
                  cycle
                endif
                if(NeighborList(i,j)) then
                  ClusterMember(j)=.true.            
                  cnt = cnt+1
                endif
              enddo        
              flipped(i)=.true.
            endif
          endif
        enddo
        if(cnt .eq. maxMol) exit      
      enddo

     
  
       if( any(ClusterMember .eqv. .false.) ) then
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
      
      logical :: ClusterMember(1:maxMol)      
      logical :: flipped(1:maxMol)
      integer i,j,h,cnt



      rejMove=.false.
      ClusterMember=.false.
      flipped=.false.


      
     do i=1,maxMol
       if(isActive(i) .eqv. .true.) then      
         if(nSwap .ne. i) then
           ClusterMember(i)=.true.
           exit
         endif
       endif
     enddo

     cnt = 0
     do i=1,maxMol
       if(isActive(i) .eqv. .false.) then
         ClusterMember(i) = .true.      
         flipped(i) = .true.         
         cnt = cnt + 1
       endif
     enddo
      
     do h=1,maxMol
        do i=1,maxMol
         if(i .ne. nSwap) then     
          if(ClusterMember(i) .neqv. flipped(i)) then
            do j=1,maxMol
              if(NeighborList(i,j)) then
                ClusterMember(j)=.true. 
                cnt = cnt + 1
              endif
            enddo
            flipped(i)=.true.
          endif
         endif
        enddo
        if(cnt .eq. maxMol-1) exit         
     enddo
  
  
     ClusterMember(nSwap) = .true.
       
     if( any(ClusterMember .eqv. .false.) ) then
       rejMove=.true.
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
        if(.not. isActive(j)) cycle
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
      
      NeighborList(nSwapIndx,:) = .false.
      NeighborList(:,nSwapIndx) = .false.      

      end subroutine      
!=================================================================================           
      end module
      
