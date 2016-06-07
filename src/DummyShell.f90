       program DummyShell
	   implicit none
      logical :: screenEcho	  
      integer ::  i,j,seed
      integer(kind=8) :: ncycle,nmoves
      integer :: outFreq_Traj,outFreq_Screen,outFreq_GCD
      integer ::	  
      character(len=10) :: labelField 
      real(kind(0.0d0)) :: dummyCycle,max_dist,max_rot
	   
	   call ReadParameters(seed,ncycle,nmoves,max_dist,max_rot,
     &	                  outFreq_Traj,outFreq_Screen,outFreq_GCD,
     &	                  screenEcho)
	   
	   
	   end