!========================================================            
      module ScriptInput
      contains
!========================================================            
      subroutine ReadParameters(seed,ncycle,nmoves,outFreq_Traj, outFreq_Screen, outFreq_GCD, screenEcho)
      use SimParameters
      use Constants
      use ForceField
      use Units
      use ParallelVar
      use CBMC_Variables
      use Coords
      use EnergyTables
      implicit none
      integer :: indx,lineStat
      real(dp) :: varValue
      character(len=200) :: line
      character(len=15) :: command, dummy
      character(len=50) :: fileName
      
      
      read(5,*) fileName
      
!      open(unit=54,file="input_Parameters.dat",status='OLD')    
      open(unit=54,file=trim(adjustl(fileName)),status='OLD')    
      
      do indx = 1, nint(1d7)
        read(54,"(A)",iostat = lineStat) line
        write(35,*) line        
        if(lineStat .lt. 0) then
          exit
        endif
       
        lineStat = 0        
        call getCommand(line, command, lineStat)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 1) then
          cycle
        endif 
       
!        call lowercase(command1)
        
        select case(trim(adjustl( command )))
        case("setvar")
          call setVariable(line)
        case("forcefield_type")
          
        case default
          write(*,"(A,2x,I10)") "ERROR! Unknown Command on line", indx
          write(*,*) line
          stop 
        end select
        
      enddo
      
      

      end subroutine
!========================================================            
      subroutine getCommand(line, command, lineStat)
      implicit none
      character(len=*), intent(in) :: line
      character(len=15), intent(out) :: command
      integer, intent(out) :: lineStat
      integer :: i, sizeLine, lowerLim, upperLim

      sizeLine = len(line)
      lineStat = 0
      i = 1
!      Find the first character in the string
      while(i .le. sizeLine)
        if(line(i:i) .ne. " ") then
          exit
        endif
!        If the first character is a # then the line is a comment. Return linestat=1 to tell the parent function
!        to skip this line. 
        if(line(i:i) .eq. "#") then
          lineStat = 1
          return
        endif
        i = i + 1
      enddo
!      If no characters are found the line is empty
      if(i .ge. sizeLine) then
        lineStat = 1
        return
      endif
      lowerLim = i

!      
      while(i .le. sizeLine)
        if(line(i:i) .eq. " ") then
          exit
        endif
        i = i + 1
      enddo
      upperLim = i

      command = line(lowerLim:upperLim)
     
      end subroutine
!========================================================            
      subroutine setVariable(line)
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: line      

      character(len=15) :: stringValue
      character(len=15) :: fileName      
      logical :: logicValue
      integer :: intValue
      real(dp) :: realValue
      
      read(line,*) dummy, command
      select case(trim(adjustl(command)))
        case("avbmc_distance")
          read(line,*) dummy, command, realValue
          Dist_Critr = realValue   
        case("cycles")
          read(line,*) dummy, command, realValue  
          nCycle = nint(realValue)
        case("gasdensity")        
          if(.not. allocated(gas_dens)) then
            write(*,*) "INPUT ERROR! GasDensity is called before the number of molecular types has been assigned"
            stop
          endif
          read(line,*) dummy, command, (gas_dens(j), j=1, nMolTypes)  
        case("moves")
          read(line,*) dummy, command, realValue        
          nMoves = nint(realValue)      
        
        case("moleculetypes")
          read(line,*) dummy, command, realValue        
          nMolTypes = nint(realValue)
          allocate( NPART(1:nMolTypes), STAT = AllocateStatus )    
          allocate( NMIN(1:nMolTypes), STAT = AllocateStatus )     
          allocate( NMAX(1:nMolTypes), STAT = AllocateStatus )     
          allocate( gas_dens(1:nMolTypes), STAT = AllocateStatus )      
          allocate( nRosenTrials(1:nMolTypes), STAT = AllocateStatus )     
          NMIN = 0
          NMAX = 0
          gas_dens = 0
          nRosenTrials = 1
          allocate( Eng_Critr(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus )
          allocate( biasAlpha(1:nMolTypes,1:nMolTypes), STAT = AllocateStatus )
        case("molmin")        
          if(.not. allocated(NMIN)) then
            write(*,*) "INPUT ERROR! molmin is called before the number of molecular types has been assigned"
            stop
          endif
          read(line,*) dummy, command, (NMIN(j), j=1, nMolTypes)    
        
        case("molmax")        
          if(.not. allocated(NMAX)) then
            write(*,*) "INPUT ERROR! molmax is called before the number of molecular types has been assigned"
            stop
          endif
          read(line,*) dummy, command, (NMAX(j), j=1, nMolTypes)
          maxMol = sum(NMAX)        

        case("rosentrials")        
          if(.not. allocated(nRosenTrials)) then
            write(*,*) "INPUT ERROR! molmax is called before the number of molecular types has been assigned"
            stop
          endif
          read(line,*) dummy, command, (nRosenTrials(j), j=1, nMolTypes)
          maxRosenTrial = maxval(nRosenTrials)
          allocate(rosenTrial(1:maxRosenTrial))
        case("temperature")        
          read(line,*) dummy, command, realValue        
          temperature = realValue
          beta = 1d0/temperature
        
        case("screenecho")
          read(line,*) dummy, command, logicValue
          screenEcho = logicValue
        case("seed")
          read(line,*) dummy, command, intValue
          seed = intValue
        case("softcutoff")
          read(line,*) dummy, command, realValue
          softCutoff = realValue
      
        case("screen_outfreq")
          read(line,*) dummy, command, intValue
          outFreq_Screen = intValue             
        case("trajectory_outfreq")     
          read(line,*) dummy, command, intValue
          outFreq_Traj = intValue  
        case("multipleinputconfig")
          read(line,*) dummy, command, logicValue
          multipleInput = logicValue  
        case("out_energyunits")
          read(line,*) labelField, outputEngUnits
          outputEConv = FindEngUnit(outputEngUnits)
        case("out_distunits")
          read(line,*) labelField, outputLenUnits   
          outputLenConv = FindLengthUnit(outputLenUnits) 
        case("umbrellasampling")
          read(line,*) dummy, command, logicValue
          useBias = logicValue
          if(useBias) then
            read(line,*) labelField, useBias, fileName
            call AllocateUmbrellaBias(fileName)
          else 
            call BlankUmbrellaBias
          endif  
        case default
          write(*,*) "Invalid variable"
          write(*,*) line
          stop      
      end select

     
      end subroutine   
!========================================================            
      subroutine variableSafetyCheck
      implicit none

   

     
      end subroutine

!========================================================            
      end module
