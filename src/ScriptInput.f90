!========================================================            
      subroutine ReadParameters(seed,ncycle,nmoves,outFreq_Traj, &
                                outFreq_Screen,outFreq_GCD,screenEcho)
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
      real(kind(0.0d0)) :: varValue
      character(len=50) :: line
      character(len=15) :: command1, command2, dummy
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
        call getCommand(line, command1, lineStat)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 1) then
          cycle
        endif 
        if(lineStat .eq. -1) then
          write(*,*) "Invalid command on line", indx
          write(*,*) line
          stop
        endif 
        
!        call lowercase(command1)
        
        select case(trim(adjustl( command1 )))
        case("setvar")
          call setVariable(line)
        case("movetypes")
       
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
        if(line(i:i) .eq. "#") then
          lineStat = 1
          return
        endif
        i = i+1
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
        i = i+1
      enddo
      upperLim = i

      command = line(lowerLim:upperLim)
     
      end subroutine
!========================================================            
      subroutine setVariable(line)
      implicit none
      character(len=*), intent(in) :: line      

      character(len=15) :: stringValue
      character(len=15) :: fileName      
      logical :: logicValue
      integer :: intValue
      real(kind(0.0d0)) :: realValue
      
      read(line,*) dummy, command
      select case(trim(adjustl(command)))
        case("cycles")
          read(line,*) dummy, command, realValue  
          nCycle = nint(realValue)
        
        case("moves")
          read(line,*) dummy, command, realValue        
          nMoves = nint(realValue)      
        
        case("moleculetypes")
          read(line,*) dummy, command, realValue        
          nMolTypes = nint(realValue)
          allocate( NPART(1:nMolTypes),STAT = AllocateStatus )    
          allocate( NMIN(1:nMolTypes),STAT = AllocateStatus )     
          allocate( NMAX(1:nMolTypes),STAT = AllocateStatus )     
          allocate( gas_dens(1:nMolTypes),STAT = AllocateStatus )      
          allocate( nRosenTrials(1:nMolTypes),STAT = AllocateStatus )     
          NMIN = 0
          NMAX = 0
          gas_dens = 0
          nRosenTrials = 1
          allocate( Eng_Critr(1:nMolTypes,1:nMolTypes), stat = AllocateStatus )
          allocate( biasAlpha(1:nMolTypes,1:nMolTypes), stat = AllocateStatus )
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
        case("gasdensity")        
          if(.not. allocated(gas_dens)) then
            write(*,*) "INPUT ERROR! GasDensity is called before the number of molecular types has been assigned"
            stop
          endif
          read(line,*) dummy, command, (gas_dens(j), j=1, nMolTypes)  
        
        case("temperature")        
          read(line,*) dummy, command, realValue        
          temperature = realValue
          beta = 1d0/temperature
        
        case("screenecho")
          read(line,*) dummy, command, logicValue
          screenEcho = logicValue
        
        case("umbrellasampling")
          read(line,*) dummy, command, logicValue
          useBias = logicValue
          if(useBias) then
            read(line,*) labelField, useBias, fileName
            call AllocateUmbrellaBias(fileName)
          else 
            call BlankUmbrellaBias
          endif     
        case("softcut")
          read(line,*) dummy, command, realValue
          softCutoff = realValue
        case("avbmc_distance")
          read(line,*) dummy, command, realValue
          softCutoff = realValue          
        case("screen_outputfrequency")
          read(line,*) dummy, command, intValue
          outFreq_Screen = intValue             
        case("trajectory_outputfrequency")     
          read(line,*) dummy, command, intValue
          outFreq_Traj = intValue  
        case default
          write(*,*) "Invalid variable"
          write(*,*) line
          stop      
      end select

     
      end subroutine   
!========================================================            
