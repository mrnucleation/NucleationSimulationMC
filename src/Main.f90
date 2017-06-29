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
      use AnalysisMain
      use AVBMC_RejectionVar
      use CBMC_Variables
      use Coords
      use CoordinateFunctions
      use Constants
      use DistanceCriteria, only: Detailed_DistanceCriteria
      use PairStorage, only: CalcAllDistPairs
      use EnergyPointers, only: Detailed_ECalc
      use EnergyTables
      use Forcefield
      use MPI
      use MiscelaniousVars, only: CollectHistograms
      use MoveTypeModule
      use ParallelVar
      use Pressure_LJ_Electro, only: Detailed_PressCalc_Inter
      use ScriptInput, only: Script_ReadParameters
      use SelfAdaptive
      use SimParameters
      use VarPrecision
      use WHAM_Functions
      use UmbrellaSamplingNew, only: useUmbrella, curUIndx, UmbrellaHistAdd, OutputUmbrellaHist, &
                                     ScreenOutputUmbrella, CheckInitialValues, energyAnalytics, OutputUmbrellaAnalytics
      implicit none

!      include 'mpif.h'
      
      logical :: errRtn
      logical :: screenEcho      

      integer(kind=8) :: iCycle, iMove
      integer :: i,j,seed,AllocateStatus
      integer :: nSel, getBiasIndex
    
!      real(dp) :: max_dist, max_dist_single, max_rot
!      real(dp) :: atmp_1,atmp_2,atmp_3,atmp_4
!      real(dp) :: acc_1, acc_2, acc_3, acc_4
      real(dp) :: E_T, E_Final, E_Debug
      real(dp) :: grnd,ran_num
      real(dp) :: dist_limit,rot_limit
      real(dp) :: TimeStart,TimeFinish
      real(dp) :: check, norm

      real(dp) :: E_Inter_Final, E_Bend_Final, E_Torsion_Final
      real(dp) :: E_Stretch_Final, E_NBond_Final, P_Final
      integer :: pressLimit, curIndx, maxIndx
      real(dp), allocatable :: P_avg_Sum(:), NHist_Sum(:)
      
      character(len=100) :: format_string,fl_name, out1
      character(len=1500) :: outFormat1, outFormat2
      
      logical,allocatable,dimension(:,:) :: FinalNeighborList
      real(dp),allocatable,dimension(:) :: FinalETable      
      real(dp),allocatable,dimension(:) :: FinalNeiETable
      
!      MPI CONTROL STATEMENTS.  This block initializes the MPI threads and assigns each process
!      a thread ID (myid) and collects the number of total processes (p_size)
      integer status(MPI_STATUS_SIZE)
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  


      maxMol = 1
!      This section of code generates the file number for each MPI thread's files
      if (myid .lt. 10) then
        format_string = "(A,I1,A)"
      elseif(myid .lt. 100) then
        format_string = "(A,I2,A)"
      elseif(myid .lt. 1000) then
        format_string = "(A,I3,A)"
      elseif(myid .lt. 10000) then
        format_string = "(A,I4,A)"          
      else
        format_string = "(A,I5,A)"      
      endif      
!     Assign screen output for each thread to fort.(100+myid)      
      nout = 100 + myid
      if(p_size .eq. 1) then
        nout = 6
      endif

      write(fl_name,format_string) "Final_Report_", myid,".txt"      
      open( unit=35, file=trim(adjustl(fl_name)) ) 
      
!      This block calls the input functions required to set up the simulation. For more information
!      on the specific role of each function see the comments in the corresponding subfunction.
!      call ReadParameters(seed, screenEcho)
      
      call Script_ReadParameters(seed, screenEcho)
      if(p_size .eq. 1) then
        write(*,*) "Input Script Complete!"
      else
        write(nout,*) "Input Script Complete!"    
      endif

      
!       This block assigns the root thread (myid=0) to output to the screen if the screen
!       echo input parameter is true.  Otherwise the screen data is exported to 100+myid      
      if(screenEcho) then      
        if(myid .eq. 0) then
          nout = 6        
        endif
      endif         

      

!   Counter for the number of Accepted Monte Carlo Moves   
      movesAccepted = 0E0_dp
      distGen_accpt = 0E0_dp
      angGen_accpt = 0E0_dp
      dihedGen_accpt = 0E0_dp
      acptRot = 0E0_dp
!  Counter for the number of Attempted Monte Carlo Moves  
      movesAttempt = 1E-30_dp
      distGen_atmp = 0E0_dp
      angGen_atmp = 0E0_dp
      dihedGen_atmp = 0E0_dp
      acptTrans = 0E0_dp
      atmpTrans = 1E-40_dp
      acptRot = 0E0_dp
      atmpRot = 1E-40_dp
!  Counter for the number of moves rejected due to the cluster criteria
      NeighRej = 0E0_dp
!      NHist = 0E0
!      NHist = 0E0
!      E_Avg = 0E0      
!      Detailed Counters for the swap move
      acptSwapIn = 0E0_dp
      acptSwapOut = 0E0_dp
      atmpSwapIn = 1E-30_dp
      atmpSwapOut = 1E-30_dp
      clusterCritRej = 0E0_dp
!      Rejection Counters for the AVBMC Insertion move.
      totalRej = 0E0_dp
      ovrlapRej = 0E0_dp
      dbalRej = 0E0_dp
      critriaRej = 0E0_dp
      boundaryRej = 0E0_dp

      totalRej_out = 0E0_dp
      dbalRej_out = 0E0_dp
      critriaRej_out = 0E0_dp
      boundaryRej_out = 0E0_dp
!      Maximum Displacement used by the translational move
      max_dist = 0.05E0_dp
!      Maximum Displacement used by the rotation move      
      max_rot = 0.05E0_dp * pi
      max_dist_single = 0.01E0_dp
!      Maximum Displacements allowed by the auto-tuning function. Failure to use this can result in
!      critical simulation errors in the auto-tuning function.
      dist_limit = 2E0_dp
      rot_limit = pi

      prevMoveAccepted = .false.
      

!      AVBMC related volume variables
      Dist_Critr_sq = Dist_Critr*Dist_Critr
      avbmc_vol = (4E0_dp/3E0_dp)*pi*Dist_Critr**3
      
!100   format(2x,I9,2x,I5,2x,E17.6,2x,F9.2,2x,F9.2,2x,F9.2)
!101   format(2x,I9,2x,I5,2x,F17.4,2x,F9.4,2x,F9.2,2x,F9.2)

      write(out1,"(A,I2,A)") "(",(4+nMoveTypes+nMolTypes),"(A))"
      write(outFormat1, out1) "(", "2x,I9", (",2x,I5",i=1,nMolTypes),",E17.6,2x",(",2x,F6.2",i=1,nMoveTypes), ")"
      write(outFormat2, out1) "(", "2x,I9", (",2x,I5",i=1,nMolTypes),",F17.4,2x",(",2x,F6.2",i=1,nMoveTypes), ")"  
!"
!      E_T is the total energy of the system
      E_T = 0E0_dp
      
!      Initialize random number generator      
      seed = p_size*seed + myid
      call sgrnd(seed)
      call CBMC_CreateTopology     

!      Perform the Intial Energy Calculations and perform the intial Cluster Criteria Check to ensure
!      the starting configuration is valid. 
!      call Detailed_EnergyCalc(E_T,errRtn)
      call Detailed_ECalc(E_T, errRtn)
      if(calcPressure) then
        call Detailed_PressCalc_Inter(pressure)
        pressLimit = 1
        do i = 1, nMolTypes
          pressLimit = pressLimit * (NMAX(i) + 1)
        enddo
        allocate( NHist(1:pressLimit) )
        allocate( NHist_Sum(1:pressLimit) )
        allocate( P_Avg(1:pressLimit) )
        allocate( P_avg_Sum(1:pressLimit) )
        NHist = 0E0_dp
        NHist_Sum = 0E0_dp
        P_Avg = 0E0_dp
        P_Avg_Sum = 0E0_dp
      endif

      if(errRtn) then
        stop      
      endif

      write(35,*) "Initial Energy Table:"
      do i=1,maxMol      
        if(isActive(i)) then
          if(ETable(i) .ne. 0E0_dp) then
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
      call TrajOutput(iCycle, E_T)

!      Initialize the variables which keep track of the different energy types during the simulation.  
      E_Inter_T = 0E0_dp
      E_NBond_T = 0E0_dp
      E_Stretch_T = 0E0_dp
      E_Bend_T = 0E0_dp
      E_Torsion_T = 0E0_dp
 

!      Output the initial parameters to the screen
      write(nout,*) "MPI threads:", p_size
      write(nout,*) "Thread ID:", myid
!      write(nout,*) "OpenMP threads:", omp_get_max_threads()
      write(nout,*) "Random Seed:", seed
      write(nout,*) "Number of Molecule Types:", nMolTypes
      write(nout,*) "Number of Initial Particles:", NPART
      write(nout,*) "Minimum Number of Particles:", NMIN
      write(nout,*) "Maximum Number of Particles:", NMAX
      write(nout,*) "Number of Atoms per Molecule:", nAtoms
      write(nout,*) "Number of Cycles:", ncycle
      write(nout,*) "Number of Moves per Cycle:", ncycle2
      write(nout,*) "CBMC Regrow Type:", regrowType
      if(useWham) then
        write(nout,*) "   ----- WHAM ------- "
        write(nout,*) "    Expected number of WHAM itterations:", nWhamItter
        write(nout,*) "    Reference Bin:", refSizeNumbers
        write(nout,*) "    Reference Bin Indx:", refBin
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
        write(nout,*) (biasAlpha(i,j), j=1,nMolTypes)       
      enddo
      write(nout,*) "Temperature:", temperature
      write(nout,*) "Gas Phase Density:", gas_dens
      write(nout,*) "Initial Energy:", E_T/outputEConv, outputEngUnits
      write(nout,*) "Initial Energy (Per Molecule):", E_T/outputEConv/real(NTotal,dp), outputEngUnits

      write(nout,*) "------------------------------------------------"
      write(nout,*) "         Simulation Start!"      
      write(nout,*) "------------------------------------------------"
      write(nout,*) "Cycle # ", "Particles ",  "Energy ", "Acceptance Rates"
!"
!      Collect the initial values for the analysis variables and print out the . 
      if(useAnalysis) then
        call PostMoveAnalysis
      endif  

      if(abs(E_T) .lt. 1E6_dp) then
        write(nout,outFormat2) 0,NPART, E_T, (1E2_dp*movesAccepted(j)/movesAttempt(j), j=1, nMoveTypes)
      else
        write(nout,outFormat1) 0,NPART, E_T,  (1E2_dp*movesAccepted(j)/movesAttempt(j), j=1, nMoveTypes)
      endif
      if(useUmbrella) then
        call UmbrellaHistAdd(E_T) 
        call CheckInitialValues
!        call ScreenOutputUmbrella
      endif
      flush(nout)
      flush(35)
!      call Detailed_ECalc(E_Debug, errRtn)
      call CPU_TIME(TimeStart)      
!--------------------------------------------------------------------------------------------------      
!      !**Begin simulation**!      
       do iCycle = 1, nCycle
         do iMove = 1, nCycle2
            !This portion of the main loop is where the type of monte carlo move is selected.
           ran_num = grnd()
           nSel = 1
           do while(moveProbability(nSel) .lt. ran_num)
             nSel = nSel + 1
           enddo
           call mcMoveArray(nSel) % moveFunction(E_T, movesAccepted(nSel), movesAttempt(nSel))
!           write(*,*) nSel
!           if(prevMoveAccepted) then
!             write(35,*) "========================================================"
!             write(35,*) nSel
!            call DEBUG_Output_NeighborList
!            call DEBUG_NeighborQualityCheck(errRtn)
!             if(errRtn) then
!               write(*,*)  nSel
!               call TrajOutput(iCycle, E_T)
!               stop      
!             endif
!           endif


           if(calcPressure) then 
             curIndx = NPart(nMolTypes)
             maxIndx = NMAX(nMolTypes) + 1 
             do i = 1, nMolTypes-1
               curIndx = curIndx + maxIndx*NPart(nMolTypes-i)
               maxIndx = maxIndx * (NMAX(nMolTypes-i) + 1)
             enddo
             NHist(curIndx) = NHist(curIndx) + 1E0_dp
             P_Avg(curIndx) = P_Avg(curIndx) + pressure 
           endif

           if(useAnalysis) then !Property Calculation Variables
             call PostMoveAnalysis
           endif

           if(useWham) then !Free Energy Biasing Block
             if(mod(iCycle, intervalWham) .gt. equilInterval) then
               if(useUmbrella) then
                 call UmbrellaHistAdd(E_T)
               endif
             endif
           else
             if(useUmbrella) then
               call UmbrellaHistAdd(E_T) 
             endif
           endif

         enddo!  End Move Loop

         
         if(mod(iCycle, 100) .eq. 0 ) then
           do i = 1, nMolTypes
             if(NMAX(i) .le. 0) then
               cycle
             endif
             call AdjustMax(acptTrans(i), atmpTrans(i), max_dist(i), dist_limit)
!             call AdjustMax(acc_2,atmp_2,max_dist_single, 0.1E0)
             call AdjustMax(acptRot(i), atmpRot(i), max_rot(i), rot_limit)   
           enddo
         endif      

!        if(totalDihed .gt. 0) then
!          if(mod(iCycle, 5000) .eq. 0 ) then
!            call UpdateDihedralProbabilites
!          endif
!        endif

!         call TrajOutput(iCycle, E_T)

!        Mid Simulation Output Block
!        if(mod(iCycle, outFreq_GCD) .eq. 0) then
          if(mod(iCycle, outFreq_Screen) .eq. 0) then
            if(abs(E_T) .lt. 1E6_dp) then
             write(nout,outFormat2) iCycle,NPART, E_T, (1E2_dp*movesAccepted(j)/movesAttempt(j), j=1, nMoveTypes)
            else
             write(nout,outFormat1) iCycle,NPART, E_T,  (1E2_dp*movesAccepted(j)/movesAttempt(j), j=1, nMoveTypes)
            endif
            if(useUmbrella) then
              call ScreenOutputUmbrella
            endif
            flush(nout)
          endif
          if(mod(iCycle,outFreq_Traj) .eq. 0) then
            call TrajOutput(iCycle, E_T)
          endif
!        endif

        if(useWham) then
          if(mod(iCycle,intervalWHAM) .eq. 0) then
            call WHAM_AdjustHist
          endif
        endif

      enddo      
!--------------------------------------------------------------------------------------------------
      call CPU_TIME(TimeFinish)      
      write(nout,*) "--------------------------------------------"

!      Collect the initial values for the analysis variables and print out the . 
!      if(useAnalysis) then
!        call PostMoveAnalysis
!        call ScreenOutputUmbrella
!      endif 

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
!      call Detailed_EnergyCalc(E_Final,errRtn)      
      call Detailed_ECalc(E_Final,errRtn)
      if(calcPressure) then
        call Detailed_PressCalc_Inter(P_Final)
      endif 





!      Output final trajectory      
      call TrajOutput(iCycle, E_T)
      close(30)    
  
!      Begin epilogue
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
      write(nout,*) "Final Energy (Per Molecule):", E_Final/outputEConv/real(NTotal,dp), outputEngUnits
      write(nout,*) "Final Cluster Size:", NPART
!     Check for errors in the energy calculation. 
      if(E_Final .eq. 0E0) then
        check = abs(E_Final - E_T)
      else
        check = abs((E_Final - E_T)/E_Final)      
      endif
      if(check .gt. 1E-6) then
        write(nout,*) "=========================================="
        write(nout,*) "Energy Disagreement"
        write(nout,*) "Culmative Energy:",E_T/outputEConv, outputEngUnits
        write(nout,*) "Culmative Energy (Per Molecule)::",E_T/outputEConv/real(NTotal,dp), outputEngUnits
        write(nout,*) "=========================================="
        write(35,*) "=========================================="
        write(35,*) "Energy Disagreement Error:"
        write(35,*) "Culmative Energy:",E_T/outputEConv,outputEngUnits
        write(35,*) "Final Energy:",E_Final/outputEConv,outputEngUnits          
        write(35,*) "=========================================="
      endif

!     Check for errors in the energy calculation. 
      if(calcPressure) then 
        write(nout,*) "Final Pressure:", P_Final
        write(nout,*) "Final Pressure (Per Molecule):", P_Final/real(NTotal,dp)
        if(P_Final .eq. 0E0) then
          check = abs(P_Final - pressure)
        else
          check = abs((P_Final - pressure)/P_Final)      
        endif
        if(check .gt. 1E-6) then
          write(nout,*) "=========================================="
          write(nout,*) "Pressure Disagreement"
          write(nout,*) "Culmative Pressure:", pressure
          write(nout,*) "Culmative Pressure (Per Molecule)::", pressure/real(NTotal,dp)
          write(nout,*) "=========================================="
          write(35,*) "=========================================="
          write(35,*) "Pressure Disagreement Error:"
          write(35,*) "Culmative Pressure:",pressure/outputEConv,outputEngUnits
          write(35,*) "Final Pressure:",P_Final/outputEConv,outputEngUnits          
          write(35,*) "=========================================="
        endif
      endif  



      write(nout,*) "Final Max Displacement", (max_dist(j), j=1,nMolTypes)
      write(nout,*) "Final Max Rotation", (max_rot(j), j=1,nMolTypes)
      do i = 1, nMoveTypes
        if(movesAttempt(i) .ne. 0E0_dp ) then
          write(nout,"(1x,A,1x,A,A,F8.2)") "Acceptance Rate", trim(adjustl(moveName(i))), ": ",1E2*movesAccepted(i)/movesAttempt(i)
        endif
      enddo
      if(any(atmpTrans .gt. 1E0_dp)) then
        write(nout,*) "Acceptance Translate (Mol Type):", (1E2*acptTrans(j)/atmpTrans(j),j=1,nMolTypes) 
      endif
      if(any(atmpRot .gt. 1E0_dp)) then
        write(nout,*) "Acceptance Rotate (Mol Type):", (1E2*acptRot(j)/atmpRot(j),j=1,nMolTypes) 
      endif

      if(distGen_atmp .ne. 0E0_dp) then
        write(nout,*) "Distance Generation Success Rate:", 1E2*distGen_accpt/distGen_atmp
      endif
      if(angGen_atmp .ne. 0E0_dp) then
        write(nout,*) "Angle Generation Success Rate:", 1E2*angGen_accpt/angGen_atmp
      endif
      if(dihedGen_atmp .ne. 0E0_dp) then
        write(nout,*) "Dihedral Angle Generation Success Rate:", 1E2*dihedGen_accpt/dihedGen_atmp
      endif
!      write(nout,*) "Cluster Criteria Rejections (Total):", 1E2*clusterCritRej/atmp_3 
      if(avbmcUsed) then
        write(nout,*) "********** AVBMC Insertion Rejection Breakdown ********"
        write(nout,*) "Percent of Moves Rejected due to Overlap:", 1E2_dp*ovrlapRej/totalRej
        write(nout,*) "Percent of Moves Rejected due to Detailed Balance:", 1E2_dp*dbalRej/totalRej
        write(nout,*) "Percent of Moves Rejected due to Cluster Criteria:", 1E2_dp*critriaRej/totalRej
        write(nout,*) "Percent of Moves Rejected due to Boundary Condition:", 1E2_dp*boundaryRej/totalRej
        write(nout,*) "********** AVBMC Deletion Rejection Breakdown ********"
        write(nout,*) "Percent of Moves Rejected due to Detailed Balance:", 1E2_dp*dbalRej_out/totalRej_out
        write(nout,*) "Percent of Moves Rejected due to Cluster Criteria:", 1E2_dp*critriaRej_out/totalRej_out
        write(nout,*) "Percent of Moves Rejected due to Boundary Condition:", 1E2_dp*boundaryRej_out/totalRej_out
        write(nout,*) "**********"
        write(nout,*) "Acceptance Swap In (Mol Type):", (1E2_dp*acptSwapIn(j)/atmpSwapIn(j),j=1,nMolTypes)
        write(nout,*) "Acceptance Swap Out (Mol Type):", (1E2_dp*acptSwapOut(j)/atmpSwapOut(j),j=1,nMolTypes) 

!        write(nout,*) "Cluster Criteria Rejections:",1E2*NeighRej/atmp_3
      endif

!      This block is the data collection block.  Variables such as the cluster histogram are collected by the root
!      processor. MPI_BARRIER ensures all threads have finished their simulations before data transfers occur
      write(nout,*) "Waiting for all processes to finish..."
      call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
      write(nout,*) "Writting output..."      


      if(calcPressure) then 
        call MPI_REDUCE(NHist, NHist_Sum, pressLimit, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 
        call MPI_REDUCE(P_Avg, P_Avg_Sum, pressLimit, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 
        if(myid .eq. 0) then
          open(unit = 50, file = "ClusterPressure.dat")
          do i = 1, pressLimit
!            write(*,*) i, P_Avg(i), NHist(i)
            if(NHist_Sum(i) .ne. 0E0_dp) then
              write(50,*) i, P_Avg_Sum(i)/(3d0*NHist_Sum(i))
            endif
          enddo 
          close(50)
        endif
      endif


!     Output Final Configuration to a visualization file that can be opened by a program like VMD or Avagadro.    

      call CollectHistograms
      write(nout,*) "Histograms Condenced...."

      if(myid .eq. 0) then
        call Output_VMD_Final
        if(useAnalysis) then
          call OutputAnalysis
        endif
      endif

      if(useUmbrella) then
        call OutputUmbrellaHist
      endif
      write(nout,*) "Histograms Outputted...."
      
      write(35,*) "Energy Table:"
      do i=1,maxMol      
        if(isActive(i)) then
          if(FinalETable(i) .ne. 0E0_dp) then
            if(abs((FinalETable(i) - ETable(i))/FinalETable(i)) .gt. 1E-6_dp) then
              write(35,*) i, ETable(i)/outputEConv, FinalETable(i)/outputEConv,"<---ETable Error"
            else 
              write(35,*) i, ETable(i)/outputEConv
            endif            
          else
            if(ETable(i) .ne. 0) then
              write(35,*) i, ETable(i)/outputEConv, FinalETable(i)/outputEConv,"<---ETable Error"            
            else
              write(35,*) i, ETable(i)/outputEConv
            endif
          endif
        endif
      enddo      
       
!       If there were any disagreements between the culmative neighborlist and final neighborlist export to Debugfile 
      write(35,*) "NeighborList:"
      do i = 1, maxMol      
        if( .not. isActive(i) ) then
          cycle
        endif
        do j = 1, maxMol
          if( .not. isActive(j) ) then
            cycle
          endif
          if(FinalNeighborList(i,j) .neqv. NeighborList(i,j))then
            if(i .ne. j) then       
              write(35,*) "Error", i, j, NeighborList(i,j), FinalNeighborList(i,j) 
            endif
          else
            if(i .ne. j) then       
              write(35,*) i, j, NeighborList(i,j)
            endif             
          endif
        enddo       
      enddo      

!      if(any(atmpInSize .ne. 0E0_dp)) then
!        write(35,*) "Acceptance Rate by Cluster Size:"           
!        do i=1,maxMol      
!          if(atmpInSize(i) .ne. 0E0_dp) then
!           write(35,*) i, 1E2*acptInSize(i)/atmpInSize(i)
!          endif
!        enddo
!      endif

!      call OutputDihedral
      close(35)




!      do i=1,maxMol      
!        if(atmpInSize(i) .ne. 0E0_dp) then
!          write(35,*) i, 1E2*acptInSize(i)/atmpInSize(i)
!        endif
!      enddo

      if(useWham) then
        call WHAM_Finalize
      endif
      if(energyAnalytics) then
        call OutputUmbrellaAnalytics
      endif


      call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
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
      
      if(atmp_x .lt. 0.5E0_dp) then
        return
      endif

      if(acc_x/atmp_x .gt. 0.5E0_dp) then
        if(max_x*1.01E0_dp .lt. limit) then
          max_x = max_x * 1.01E0_dp
        else 
          max_x = limit       
        endif
      else
        max_x = max_x * 0.99E0_dp
      endif

 
      end subroutine
!===========================================================
      subroutine TrajOutput(iCycle, E_T)
      use VarPrecision
      use SimParameters
      use Coords
      use Forcefield
      implicit none
      integer(kind=8), intent(in) :: iCycle
      real(dp), intent(in) :: E_T
      integer :: iType, iMol, iAtom
      integer ::  atmType
      integer :: cnt      
      real(dp) :: xcm,ycm,zcm


      xcm = 0E0_dp
      ycm = 0E0_dp
      zcm = 0E0_dp
      cnt = 0
      do iType = 1,nMolTypes
        do iMol = 1,NPART(iType)
          xcm = xcm + MolArray(iType)%mol(iMol)%x(1)
          ycm = ycm + MolArray(iType)%mol(iMol)%y(1)
          zcm = zcm + MolArray(iType)%mol(iMol)%z(1)
          cnt = cnt + 1
        enddo
      enddo
      
      xcm = xcm/real(cnt, dp)
      ycm = ycm/real(cnt, dp)
      zcm = zcm/real(cnt, dp)

      write(30,*) vmdAtoms
      write(30,*) NPART, E_T
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
            write(30,*) atomData(atmType)%Symb, 1E7_dp, 1E7_dp, 1E7_dp
          enddo
        enddo
      enddo


      end subroutine
!=============================================================================================
!     This subroutine prints a dummy frame to the Trajectory file.  This
!     is done so that when the trajectory is loaded into VMD it will properly
!     display the cluster configuration.
      subroutine InitialTrajOutput
      use VarPrecision
      use SimParameters      
      use Coords
      use Forcefield
      implicit none
      real(dp), parameter :: xOffset = 50E0_dp
      real(dp), parameter :: yOffset = 20E0_dp
      integer :: iType, iMol, iAtom
      integer :: atmType
      integer :: cnt      
      real(dp) :: xcm,ycm,zcm


      xcm = 0E0_dp
      ycm = 0E0_dp
      zcm = 0E0_dp
      cnt = 0
      do iType = 1,nMolTypes
        do iMol = 1,NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            xcm = xcm + gasConfig(iType)%x(iAtom) + xOffset*( real(NMAX(iType),dp)/2E0_dp - real(iMol,dp) )
            ycm = ycm + gasConfig(iType)%y(iAtom) + yOffset*real(iType, dp)
            zcm = zcm + gasConfig(iType)%z(iAtom)
            cnt = cnt + 1
          enddo          
        enddo
      enddo

      xcm = xcm/real(cnt, dp)
      ycm = ycm/real(cnt, dp)
      zcm = zcm/real(cnt, dp)

      write(30,*) vmdAtoms
      write(30,*) NMAX
      do iType = 1,nMolTypes
        do iMol = 1, NMAX(iType)
          do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType,iAtom)
            write(30,*) atomData(atmType)%Symb, &
                      gasConfig(iType)%x(iAtom) + xOffset*( real(NMAX(iType),dp)/2E0_dp - real(iMol, dp) )-xcm, &
                      gasConfig(iType)%y(iAtom) + yOffset*real(iType, dp) - ycm, &
                      gasConfig(iType)%z(iAtom) - zcm
          enddo
        enddo
      enddo


      end subroutine
!========================================================            
!     The purpose of this subroutine is to lower case a given character string. 
      subroutine LowerCaseLine(line)
      implicit none
      character(len=*),intent(inout) :: line
      integer, parameter :: offset = ichar("a") - ichar("A")
      integer :: i,sizeLine
      integer :: curVal, newVal

      sizeLine = len(line)

      do i = 1, sizeLine
        curVal = ichar(line(i:i))
        if(curVal .le. ichar("Z")) then
          if(curVal .ge. ichar("A")) then
            newVal = curVal + offSet
            line(i:i) = char(newVal)
          endif
        endif
      enddo

      end subroutine
!=============================================================================================      

