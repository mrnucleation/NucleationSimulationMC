!========================================================================================
! ******************** Multiple Component Grand Canonical Monte Carlo Nucleation Code ********************
! Developed by: Troy Loeffler
! For use in the Bin Chen Research group at LSU
!
! The purpose of this module is to simulate the nucleation of a 
! many component chemical system using the grand canonical ensemble.  This simulation
! makes use of techniques such as the Aggregation-Volume-Bias Monte Carlo, 
! Configurational-Bias-Monte Carlo, etc. 
!========================================================================================
      program Swaper
      use AcceptRates
      use AVBMC_RejectionVar
      use CBMC_Variables
      use Coords
      use Constants
      use E_Interface
      use EnergyTables
      use Forcefield
      use Histogram
!      use MPI
      use ParallelVar
      use SimParameters
      use VarPrecision
      use WHAM_Module
      implicit none

      include 'mpif.h'
      
      logical :: errRtn
      logical :: screenEcho      
      integer(kind=8) :: ncycle,indx,ncycle2,indx2
      integer :: i,j,seed,AllocateStatus
      integer :: outFreq_Traj,outFreq_Screen,outFreq_GCD      
    
    
!      real(dp) :: max_dist, max_dist_single, max_rot
      real(dp) :: atmp_1,atmp_2,atmp_3,atmp_4
      real(dp) :: acc_1, acc_2, acc_3, acc_4
      real(dp) :: E_T, E_Final
      real(dp) :: grnd,ran_num
      real(dp) :: dist_limit,rot_limit
      real(dp) :: TimeStart,TimeFinish
      real(dp) :: check

      real(dp) :: E_Inter_Final, E_Bend_Final, E_Torsion_Final
      real(dp) :: E_Stretch_Final, E_NBond_Final
      
      character(len=30) :: format_string,fl_name, out1
      character(len=100) :: outFormat1, outFormat2
      
      logical,allocatable,dimension(:,:) :: FinalNeighborList
      real(dp),allocatable,dimension(:) :: FinalETable      
      real(dp),allocatable,dimension(:) :: FinalNeiETable
      
!      MPI CONTROL STATEMENTS.  This block initializes the MPI threads and assigns each process
!      a thread ID (myid) and collects the number of total processes (p_size)
      integer status(MPI_STATUS_SIZE)
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)      
!      This section of code generates the file number for each MPI thread's files
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
!     Assign screen output for each thread to fort.(100+myid)      
      nout = 100+myid

      write(fl_name,format_string) "Final_Report_", myid,".txt"      
      open( unit=35, file=trim(adjustl(fl_name)) ) 
      
!      This block calls the input functions required to set up the simulation. For more information
!      on the specific role of each function see the comments in the corresponding subfunction.
      call ReadParameters(seed,ncycle,ncycle2,outFreq_Traj, outFreq_Screen,outFreq_GCD,screenEcho)

!       This block assigns the root thread (myid=0) to output to the screen if the screen
!       echo input parameter is true.  Otherwise the screen data is exported to 100+myid      
      if(screenEcho) then      
        if(myid .eq. 0) then
          nout=6        
        endif
      endif         
      call ReadForcefield
      call AllocateCoordinateArrays
      call CreateJointArray      
      call ReadInitialConfiguration
      call RecenterCoordinates
      call ReadInitialGasPhase      
      
      
!      useWHAM = .true.
!      write(nout,*) "WHAM:", useWHAM
!      equilInterval = nint(0.02d0*intervalWHAM)
      

!   Counter for the number of Accepted Monte Carlo Moves   
      acc_1=0d0
      acc_2=0d0  
      acc_3=0d0
      acc_4=0d0 
      angGen_accpt = 0d0
!  Counter for the number of Attempted Monte Carlo Moves  
      atmp_1=0d0
      atmp_2=0d0
      atmp_3=0d0  
      atmp_4=0d0  
      angGen_atmp = 0d0

!  Counter for the number of moves rejected due to the cluster criteria
      NeighRej=0d0
      NHist=0d0
      E_Avg=0d0      
!  Detailed Counters for the swap move
      acptSwapIn = 0d0
      acptSwapOut = 0d0
      atmpSwapIn = 1d-30
      atmpSwapOut = 1d-30
      clusterCritRej = 0d0
!      Rejection Counters for the AVBMC Insertion move.
      totalRej = 0d0
      ovrlapRej = 0d0
      dbalRej = 0d0
      critriaRej = 0d0
      
!      Maximum Displacement used by the translational move
      max_dist = 0.1d0
!      Maximum Displacement used by the rotation move      
      max_rot = 0.01d0*pi
      max_dist_single = 0.01d0
!      Maximum Displacements allowed by the auto-tuning function. Failure to use this can result in
!      critical simulation errors in the auto-tuning function.
      dist_limit = 2d0
      rot_limit = pi
      
!      Histogram Initialization
      d_ang = dble(nBin_Hist)/two_pi
      d_r = dble(nBin_Hist)/(2d0)
      HistAngle = 0d0
      HistDist = 0d0

!      AVBMC related volume variables
      Dist_Critr_sq = Dist_Critr*Dist_Critr
      avbmc_vol = (4d0/3d0)*pi*Dist_Critr**3
      
!100   format(2x,I9,2x,I5,2x,E17.6,2x,F9.2,2x,F9.2,2x,F9.2)
!101   format(2x,I9,2x,I5,2x,F17.4,2x,F9.4,2x,F9.2,2x,F9.2)

      write(out1,"(A,I2,A)") "(",(8+nMolTypes),"(A))"
      write(outFormat1, out1) "(", "2x,I9", (",2x,I5",i=1,nMolTypes),",E17.6,2x",(",2x,F9.2",i=1,4), ")"
      write(outFormat2, out1) "(", "2x,I9", (",2x,I5",i=1,nMolTypes),",F17.4,2x",(",2x,F9.2",i=1,4), ")"  
!"
!      E_T is the total energy of the system
      E_T=0d0
      
!      Initialize random number generator      
      seed = p_size*seed + myid
      call sgrnd(seed)

      call CBMC_CreateTopology     

!      Perform the Intial Energy Calculations and perform the intial Cluster Criteria Check to ensure
!      the starting configuration is valid. 
      call Detailed_EnergyCalc(E_T,errRtn)
      if(errRtn) then
        stop      
      endif

      write(35,*) "Initial Energy Table:"
      do i=1,maxMol      
        if(isActive(i)) then
          if(ETable(i) .ne. 0d0) then
             write(35,*) i, ETable(i)
          endif            
        endif
      enddo      
      write(35,*)
!      Calculate the inital energy table used in the energy biased AVBMC algorithim
!      call Create_NeiETable
 
      write(fl_name,format_string) "Traj", myid,".xyz"    
      open(unit=30,file=trim(fl_name),FORM='FORMATTED')      

!      Print Dummy frame to VMD so VMD will correctly display varying cluster sizes.       
      call InitialTrajOutput
!      This function outputs current config to the trajectory visualization file.       
      call TrajOutput(indx)

!      Initialize the variables which keep track of the different energy types during the simulation.  
      E_Inter_T = 0d0
      E_NBond_T = 0d0
      E_Stretch_T = 0d0
      E_Bend_T = 0d0
      E_Torsion_T = 0d0
 

!      Output the initial parameters to the screen
      write(nout,*) "MPI threads:", p_size
      write(nout,*) "Thread ID:", myid
!      write(nout,*) "OpenMP threads:", omp_get_max_threads()
      write(nout,*) "Random Seed:", p_size*seed + myid
      write(nout,*) "Number of Molecule Types:", nMolTypes
      write(nout,*) "Number of Initial Particles:", NPART
      write(nout,*) "Number of Atoms per Molecule:", nAtoms
      write(nout,*) "Number of Cycles:", ncycle
      write(nout,*) "Number of Moves per Cycle:", ncycle2
      write(nout,*) "CBMC Regrow Type:", regrowType
      if(useWham) then
        write(nout,*) "   ----- WHAM ------- "
        write(nout,*) "    Expected number of WHAM itterations:", nWhamItter
        write(nout,*) "    Reference Bin:", refSizeNumbers
        write(nout,*) "    Number of Cycle before WHAM adjustment:", intervalWHAM
        write(nout,*) "    WHAM Convergence Tolerance:", tolLimit
        write(nout,*) "    Maximum WHAM Itterations:", maxSelfConsist
        write(nout,*) "    WHAM Number of Eqilibriation MC Cycles:", equilInterval
        write(nout,*) "   ---------  -------- "
      endif
      if(distCriteria) then
        write(nout,*) "Criteria Type: Distance" 
      else
        write(nout,*) "Criteria Type: Energy"       
      endif
      write(nout,*) "Distance Criteria:", Dist_Critr      
      write(nout,*) "Distance Criteria SQ:", Dist_Critr_sq
      write(nout,*) "AVBMC Volume:", avbmc_vol
      write(nout,*) "Energy Criteria:"      
      do i = 1, nMolTypes
        write(nout,*) (Eng_Critr(i,j), j=1, nMolTypes)
      enddo
      write(35,*) "Rejection Distance:", r_min
      write(nout,*) "AVBMC E_Bias Alpha:"      
      do i = 1, nMolTypes      
        write(nout,*) (biasAlpha(i,j),j=1,nMolTypes)       
      enddo
      write(nout,*) "Temperature:", temperature
      write(nout,*) "Gas Phase Density:", gas_dens
      write(nout,*) "Initial Energy:", E_T/outputEConv, outputEngUnits

      write(nout,*) "------------------------------------------------"
      write(nout,*) "         Simulation Start!"      
      write(nout,*) "------------------------------------------------"
      write(nout,*) "Cycle # ", "Particles ",  "Energy ", "Acceptance Rates"
!"
      atmp_1=0.00001d0
      atmp_2=0.00001d0      
      atmp_3=0.00001d0
      atmp_4=0.00001d0
      flush(35)
      call CPU_TIME(TimeStart)      
!--------------------------------------------------------------------------------------------------      
!      !**Begin simulation**!      
       do indx=1,nCycle
         do indx2=1,nCycle2
           ran_num=grnd()
           if(ran_num .lt. 1d0/2d0) then
             call AVBMC(E_T, acc_3, atmp_3)
!           elseif(ran_num .lt. 2d0/4d0) then
!             call CBMC(E_T, acc_2, atmp_2)
!             call SingleAtom_Translation(E_T,acc_2,atmp_2,max_dist_single)
!           elseif(ran_num .lt. 3d0/4d0) then
!             call Rotation(E_T,acc_4,atmp_4)  
           else
             call Translation(E_T,acc_1,atmp_1)
           endif

!            call SingleAtom_Translation(E_T,acc_2,atmp_2,0.05d0)           
!           call Simple_Intra(E_T, acc_2, atmp_2)  
!           call Calc_BendAngle(1, 1)
!           call Calc_BondStretch(1, 1)
!           call Calc_BondStretch(1, 2)           
           if(useWham) then
             if(mod(indx,intervalWham) .gt. equilInterval) then
               call NHistAdd(E_T) 
             endif
           else
             call NHistAdd(E_T) 
           endif
         enddo


            
         if(mod(indx,100) .eq. 0 ) then
           do i = 1, nMolTypes
             call AdjustMax(acptTrans(i), atmpTrans(i), max_dist(i), dist_limit)
!             call AdjustMax(acc_2,atmp_2,max_dist_single, 0.1d0)
             call AdjustMax(acptRot(i), atmpRot(i), max_rot(i), rot_limit)   
           enddo
         endif      

!         call TrajOutput(indx)
!        Mid Simulation Output Block
        if(mod(indx,outFreq_GCD) .eq. 0) then
          if(mod(indx,outFreq_Screen) .eq. 0) then
            if(abs(E_T) .lt. 1d6) then
             write(nout,outFormat2) indx,NPART, E_T,1d2*acc_1/atmp_1, 1d2*acc_3/atmp_3
            else
             write(nout,outFormat1) indx,NPART, E_T,1d2*acc_1/atmp_1, 1d2*acc_3/atmp_3
            endif
          endif
          if(mod(indx,outFreq_Traj) .eq. 0) call TrajOutput(indx)
        endif

        if(useWham) then
         if(mod(indx,intervalWHAM) .eq. 0) then
           call WHAM_AdjustHist
         endif
        endif

      enddo      
!--------------------------------------------------------------------------------------------------
      call CPU_TIME(TimeFinish)      
      write(nout,*) "--------------------------------------------"

!      This block initializes variables used in calculating averages as well as variables used
!      in checking for errors that may occur during the simulation. 
      ALLOCATE (FinalNeighborList(1:maxMol,1:maxMol),  STAT = AllocateStatus)
      ALLOCATE (FinalETable(1:maxMol), STAT = AllocateStatus)  
      ALLOCATE (FinalNeiETable(1:maxMol),STAT = AllocateStatus) 
      
!      Here the final neighbor list is placed into a temporary variable for error checking purposes.
!      The "NeighborList" variable will be recalculated by the detailed energy functions.  The two
!      lists will then be compared to see if there are any decrepancies.  If decrpancies exist
!      the cluster criteria was not properly maintained during the course of the simulation.
      FinalNeighborList = NeighborList
      FinalETable = ETable
      FinalNeiETable = NeiETable
      
      E_Inter_Final = E_Inter_T
      E_NBond_Final = E_NBond_T      
      E_Stretch_Final = E_Stretch_T
      E_Bend_Final = E_Bend_T
      E_Torsion_Final = E_Torsion_T

     
!      Calculate the final energy using the detailed energy function
!      and compare it to the culmative energy to check for possible errors during the simulation. 
!      If these two values do not match there is an error in the energy calculation routines. 
      call Detailed_EnergyCalc(E_Final,errRtn)      
      call Create_NeiETable      
!      Output final trajectory      
      call TrajOutput(indx)
      close(30)    
  
      write(35,*) "---------------------------------"
      write(35,*) "Culmative Energy,   Detailed Calc Final Energy"
      write(35,*) "Inter:", E_Inter_Final, E_Inter_T
      write(35,*) "Intra Nonbonded:", E_NBond_Final, E_NBond_T      
      write(35,*) "Stretch:", E_Stretch_Final, E_Stretch_T
      write(35,*) "Bend:", E_Bend_Final, E_Bend_T
      write(35,*) "Torsional:", E_Torsion_Final, E_Torsion_T
      
!      This block writes the final simulation report to the screen.  If any errors were detected 
!      those will written as well.

      write(nout,*) "Simulation Time:", TimeFinish - TimeStart
      write(nout,*) "Final Energy:", E_Final/outputEConv, outputEngUnits
      write(nout,*) "Final Cluster Size:", NPART
!     Check for errors in the energy calculation. 
      if(E_Final .eq. 0d0) then
        check = abs(E_Final - E_T)
      else
        check = abs((E_Final - E_T)/E_Final)      
      endif
      if(check .gt. 1d-6) then
        write(nout,*) "=========================================="
        write(nout,*) "Energy Disagreement"
        write(nout,*) "Culmative Energy:",E_T/outputEConv, outputEngUnits
        write(nout,*) "=========================================="
        write(35,*) "=========================================="
        write(35,*) "Energy Disagreement Error:"
        write(35,*) "Culmative Energy:",E_T/outputEConv,outputEngUnits
        write(35,*) "Final Energy:",E_Final/outputEConv,outputEngUnits          
        write(35,*) "=========================================="
      endif        
      write(nout,*) "Final Max Displacement", (max_dist(j), j=1,nMolTypes)
      write(nout,*) "Final Max Rotation", (max_rot(j), j=1,nMolTypes)
      write(nout,*) "Acceptance Rate Translation:", 1d2*acc_1/atmp_1
      write(nout,*) "Acceptance Rate CBMC:", 1d2*acc_2/atmp_2
      write(nout,*) "Acceptance Rate Swap:", 1d2*acc_3/atmp_3
      write(nout,*) "Acceptance Rate Rotate:", 1d2*acc_4/atmp_4
      write(nout,*) "Acceptance Translate (Mol Type):", (1d2*acptTrans(j)/atmpTrans(j),j=1,nMolTypes) 
      write(nout,*) "Acceptance Rotate (Mol Type):", (1d2*acptRot(j)/atmpRot(j),j=1,nMolTypes) 
      write(nout,*) "Angle Generation Success Rate:", 1d2*angGen_accpt/angGen_atmp
      write(nout,*) "Cluster Criteria Rejections (Total):", 1d2*clusterCritRej/atmp_3 
      write(nout,*) "********** AVBMC Insertion Rejection Breakdown ********"
      write(nout,*) "Percent of Moves Rejected due to Overlap:", 1d2*ovrlapRej/totalRej
      write(nout,*) "Percent of Moves Rejected due to Detailed Balance:", 1d2*dbalRej/totalRej
      write(nout,*) "Percent of Moves Rejected due to Cluster Criteria:", 1d2*critriaRej/totalRej
      write(nout,*) "**********"
      if(atmp_3 .ge. 1d0) then
        write(nout,*) "Acceptance Swap In (Mol Type):", (1d2*acptSwapIn(j)/atmpSwapIn(j),j=1,nMolTypes)
        write(nout,*) "Acceptance Swap Out (Mol Type):", (1d2*acptSwapOut(j)/atmpSwapOut(j),j=1,nMolTypes) 
        write(nout,*) "Cluster Criteria Rejections:",1d2*NeighRej/atmp_3
      endif

!      This block is the data collection block.  Variables such as the cluster histogram are collected by the root
!      processor. MPI_BARRIER ensures all threads have finished their simulations before data transfers occur
      write(nout,*) "Waiting for all processes to finish..."
      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
      write(nout,*) "Writting output..."      
      call CollapseN  
      if(myid .eq. 0) then
        call N_HistOutput    
      endif
      write(nout,*) "Histogram Outputted...."
!     Output Final Configuration to a visualization file that can be opened by a program like VMD or Avagadro.    
      if(myid .eq. 0) then
        call Output_VMD_Final
      endif

      
      write(35,*) "Energy Table:"
      do i=1,maxMol      
        if(isActive(i)) then
          if(FinalETable(i) .ne. 0d0) then
            if(abs((FinalETable(i) - ETable(i))/FinalETable(i)) .gt. 1d-6) then
              write(35,*) i, ETable(i),FinalETable(i),"<---ETable Error"
            else 
              write(35,*) i, ETable(i)
            endif            
          else
            if(ETable(i) .ne. 0) then
              write(35,*) i, ETable(i),FinalETable(i),"<---ETable Error"            
            else
              write(35,*) i, ETable(i)                     
            endif
          endif
        endif
      enddo      
       
!       If there were any disagreements between the culmative neighborlist and final neighborlist export to Debugfile 
      write(35,*) "NeighborList:"
      do i = 1, maxMol      
        do j = 1, maxMol
          if(FinalNeighborList(i,j) .neqv. NeighborList(i,j))then
            if(i .ne. j) then       
              write(35,*) i, j, NeighborList(i,j), FinalNeighborList(i,j) 
            endif
          endif
        enddo       
      enddo      

      if(any(HistAngle .ne. 0d0)) then
        write(35,*)
        write(35,*) "Angle Histogram:"      
        do i=0, nBin_Hist      
           write(35,*) i/d_ang, HistAngle(i) 
        enddo
      endif
      if(any(HistDist .ne. 0d0)) then
        write(35,*) "Dist Histogram:"      
        do i = 0, nBin_Hist      
          write(35,*) i/d_r, HistDist(i) 
        enddo
      endif

      if(any(atmpInSize .ne. 0d0)) then
        write(35,*) "Acceptance Rate by Cluster Size:"           
        do i=1,maxMol      
          if(atmpInSize(i) .ne. 0d0) then
           write(35,*) i, 1d2*acptInSize(i)/atmpInSize(i)
          endif
        enddo
      endif
      close(35)

      if(useWham) then
        call WHAM_Finalize
      endif

      
      write(nout,*) "Finished!"
      close(nout)
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
      call MPI_FINALIZE(ierror)      
      end program
!========================================================    
!     This function adjusts the maximum displacement for both the rotation and translational moves such that
!     a 50% acceptance rate is maitained through out the simulation.        
      subroutine AdjustMax(acc_x, atmp_x, max_x, limit)
      use VarPrecision
      implicit none
      real(dp), intent(in) :: acc_x,atmp_x,limit
      real(dp), intent(inout):: max_x
      
      if(atmp_x .eq. 0d0) then
        return
      endif

      if(acc_x/atmp_x .gt. 0.5d0) then
        if(max_x*1.01d0 .lt. limit) then
          max_x=max_x*1.01d0
        else 
          max_x=limit       
        endif
      else
        max_x=max_x*0.99d0      
      endif

 
      end subroutine
!===========================================================

      subroutine TrajOutput(indx)
      use VarPrecision
      use SimParameters
      use Coords
      use Forcefield
      implicit none
      integer(kind=8), intent(in) :: indx
      integer :: iType, iMol, iAtom
      integer ::  atmType
      integer :: cnt      
      real(dp) :: xcm,ycm,zcm


      xcm = 0d0
      ycm = 0d0
      zcm = 0d0
      cnt = 0
      do iType = 1,nMolTypes
        do iMol = 1,NPART(iType)
          xcm = xcm + MolArray(iType)%mol(iMol)%x(1)
          ycm = ycm + MolArray(iType)%mol(iMol)%y(1)
          zcm = zcm + MolArray(iType)%mol(iMol)%z(1)
          cnt = cnt + 1
        enddo
      enddo
      
      xcm = xcm/dble(cnt)
      ycm = ycm/dble(cnt)
      zcm = zcm/dble(cnt)

      write(30,*) vmdAtoms
      write(30,*) indx
      do iType = 1,nMolTypes
        do iMol = 1, NPART(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(30,*) atomData(atmType)%Symb, &
                        MolArray(iType)%mol(iMol)%x(iAtom)-xcm, &
                        MolArray(iType)%mol(iMol)%y(iAtom)-ycm, &
                        MolArray(iType)%mol(iMol)%z(iAtom)-zcm
          enddo
        enddo
        do iMol = NPART(iType)+1, NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(30,*) atomData(atmType)%Symb, 1d7, 1d7, 1d7
          enddo
        enddo
      enddo


      end subroutine
!===========================================================

      subroutine InitialTrajOutput
      use VarPrecision
      use SimParameters      
      use Coords
      use Forcefield
      implicit none
      real(dp), parameter :: xOffset = 50d0
      real(dp), parameter :: yOffset = 20d0     
      integer :: iType, iMol, iAtom
      integer :: atmType
      integer :: cnt      
      real(dp) :: xcm,ycm,zcm


      xcm = 0d0
      ycm = 0d0
      zcm = 0d0
      cnt = 0
      do iType = 1,nMolTypes
        do iMol = 1,NMAX(iType)
          do iAtom = 1, nAtoms(iType)
!            xcm = xcm + MolArray(iType)%mol(1)%x(iAtom) + 
!     &                xOffset*(dble(NMAX(iType)/2d0 - dble(iMol)))
!            ycm = ycm + MolArray(iType)%mol(1)%y(iAtom) + 
!     &                yOffset*dble(iType)
!            zcm = zcm + MolArray(iType)%mol(1)%z(iAtom)
            xcm = xcm + gasConfig(iType)%x(iAtom) + xOffset*(dble(NMAX(iType)/2d0 - dble(iMol)))
            ycm = ycm + gasConfig(iType)%y(iAtom) + yOffset*dble(iType)
            zcm = zcm + gasConfig(iType)%z(iAtom)
            cnt = cnt + 1
          enddo          
        enddo
      enddo

      xcm = xcm/dble(cnt)
      ycm = ycm/dble(cnt)
      zcm = zcm/dble(cnt)

      write(30,*) vmdAtoms
      write(30,*) 0
      do iType = 1,nMolTypes
        do iMol = 1, NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
!            write(30,*) atomData(atmType)%Symb
!     &                  ,MolArray(iType)%mol(1)%x(iAtom)
!     &                +xOffset*(dble(NMAX(iType)/2d0 - dble(iMol)))-xcm
!     &                  ,MolArray(iType)%mol(1)%y(iAtom)
!     &                +yOffset*dble(iType)-ycm
!     &                  ,MolArray(iType)%mol(1)%z(iAtom)-zcm            
            write(30,*) atomData(atmType)%Symb, &
                      gasConfig(iType)%x(iAtom) + xOffset*(dble(NMAX(iType)/2d0 - dble(iMol)))-xcm, &
                      gasConfig(iType)%y(iAtom) + yOffset*dble(iType) - ycm, &
                      gasConfig(iType)%z(iAtom) - zcm
          enddo
        enddo
      enddo


      end subroutine
      

