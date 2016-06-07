      program BotAdjuster
      implicit none

      integer :: i,j,bin_Eq,bin_last,nbin, dummy
      real(kind(0.0d0)),Dimension(:),Allocatable:: Hist,ClusterBias,NewBias
      real(kind(0.0d0)) :: junk,d_the
      real(kind(0.0d0)) :: pi,ang_eq,Hist_eq,temp,min_new
      real(kind(0.0d0)) :: BiasLast, curHist, curbias

      write(6,*) "Wxhat is the smallest Clustersize?"
      read *,bin_Eq
      
      write(6,*) "What is the largest Clustersize?"
      read *,nbin
 
  
      Allocate(ClusterBias(1:nbin))
      Allocate(Hist(1:nbin))
      Allocate(NewBias(0:nBin))

      Clusterbias = 0d0
      Hist = 0d0
      NewBias = 0d0

      pi=4d0*datan(1d0)
      d_the=nBin/(pi)
      ang_eq=pi*70d0/180d0
      do i=1,nBin
        read(25,*,end=20) dummy, curBias
        ClusterBias(dummy) = curBias
      enddo

20    continue      
      Hist = 0d0
      
      open(unit=12,file="Histogram_ClusterSize.txt")
!      read(12,*)
      do i=1,nBin
        read(12,*, end=30) dummy, curHist
        if(dummy .ne. 0) then
          Hist(dummy) = Hist(dummy) + curHist
        endif
      enddo

30    continue

      write(6,*) "Histogram:"
      do i=1,nBin
        write(6,*) i,Hist(i)
      enddo
      
      do i=1,nBin
        NewBias(i) = ClusterBias(i)
      enddo
!      bin_Eq=1
      Hist_eq=Hist(bin_Eq)
!      if(Hist_eq .eq. 0) Hist_eq=5d-1
      bin_last=bin_Eq
!      write(6,*) bin_eq,Hist_Eq,nBin
      write(6,*) "Bias:"
      do i=bin_Eq,nBin
       if(Hist(i) .ne. 0d0) then
        NewBias(i)=ClusterBias(i)-log(Hist(i)/Hist_eq)
        if(i .gt. bin_last) then
         bin_last=i
         BiasLast=NewBias(i)
        endif
!       else
!        NewBias(i)=ClusterBias(i)
       endif
      enddo

      do i=bin_last,nBin
        NewBias(i)=BiasLast
      enddo  
  
      do i=1,nBin
!         Newbias(i)=Newbias(i)!-min_new
         write(6,*) i,ClusterBias(i),NewBias(i)
      enddo

!      NewBias(1)=0d0
      open(unit=35,file="NewBias.txt")
      do i=1,nBin
        write(35,*) i, NewBias(i)
      enddo

      end
