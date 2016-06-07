      program EngApprox
      implicit none
      integer, parameter :: nBin=1000
!      real(kind(0.0d0)), parameter :: cutOff = 3d0*(26d0/7d0)**(1d0/6d0)
      real(kind(0.0d0)), parameter :: cutOff = 3d0
      integer :: i,j,ii,jj
      integer :: atomType(1:4)
      real(kind(0.0d0)) :: cnt(1:nBin)
      real(kind(0.0d0)) :: x(1:nBin,1:4), y(1:nBin,1:4), z(1:nBin,1:4)
      real(kind(0.0d0)) :: Energy(1:nBin), AvgR(1:nBin,1:4,1:3)
      character(len=2) :: junk
      real(kind(0.0d0)) :: rx,ry,rz,r
      real(kind(0.0d0)) :: LJ
      real(kind(0.0d0)) :: ApproxEnergy(1:nBin)
      real(kind(0.0d0)) :: q(1:4), q_tab(1:4,1:4)
      
      write(6,*) "Cutoff:", cutOff
      write(6,*) "nBin:", nBin

      q(1) = 0d0
      q(2) = 0.52d0
      q(3) = 0.52d0      
      q(4) = 1.040d0
      
      atomType(1) = 1
      atomType(2) = 2
      atomType(3) = 2
      atomType(4) = 3
      
      do ii=1,4
        do jj=1,4
           q_tab(ii,jj) = q(ii)*q(jj)
        enddo
      enddo
      
      
      open(unit=12,file="NewConfig.xyz")
      read(12,*)
      read(12,*)
      do i=1,nBin
        do j=1,4      
          read(12,*) junk, x(i,j), y(i,j), z(i,j)
        enddo
      enddo

      Energy = 0d0
      AvgR = 0d0
      cnt = 0d0
      ApproxEnergy = 0d0
      do i=1,nBin-1
        do ii=2,4
         do j=i+1, nBin
          do jj=2,4
             rx = x(i,ii) - x(j,jj)
             ry = y(i,ii) - y(j,jj)
             rz = z(i,ii) - z(j,jj)
             r = rx*rx + ry*ry + rz*rz
             r = dsqrt(r)
             if(r .gt. cutOff) then
               AvgR(i,ii) = AvgR(i) + r
               AvgR(j) = AvgR(j) + r
!               LJ = 1d0/r**6
!               LJ = 4d0*50d0*LJ*(LJ-1d0)
               LJ = LJ-1.671d5/r
               Energy(i) = Energy(i) + LJ
               Energy(j) = Energy(j) + LJ
               cnt(i) = cnt(i) + 1d0
               cnt(j) = cnt(j) + 1d0
              endif
            enddo
          enddo          
        enddo
      enddo
     
      do i = 1, nBin
        AvgR(i) = AvgR(i)/cnt(i)
        LJ = 1d0/AvgR(i)**6
        LJ = 4d0*50d0*LJ*(LJ-1d0)
        LJ = LJ - 1.671d5/AvgR(i)
        ApproxEnergy(i) = LJ
      enddo 

      
      do i = 1,nBin
        write(6,*) i, cnt(i),AvgR(i),Energy(i)/300d0, 
     &                (cnt(i)*ApproxEnergy(i))/300d0
        write(2,*) Energy(i),ApproxEnergy(i)*cnt(i)
      enddo

 
      end program
