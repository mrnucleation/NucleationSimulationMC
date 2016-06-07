!==============================================================
      module Constants
 
      real(kind(0.0d0)),parameter :: pi=4d0*datan(1d0) 
      real(kind(0.0d0)),parameter :: two_pi=8d0*datan(1d0)
        
      end module
!===================================================================      
      module Units
      contains
        
!     !----------------------------------------------------------
        double precision function FindEngUnit(unitName)
        implicit none 
        character(len=10), intent(in) :: unitName       
        
        
        select case(trim(adjustl(unitName)))
          case("j-mol")
                FindEngUnit = 1d0/8.3144621d0
          case("kj-mol")
                FindEngUnit = 1d0/8.3144621d-3
          case("cal-mol")
                FindEngUnit = 1d0/1.9872041d0
          case("kcal-mol")
                FindEngUnit = 1d0/1.9872041d-3
          case("kB")
                FindEngUnit = 1d0
          case("kb")
                FindEngUnit = 1d0
          case default
            write(6,*) "Error! Invalid Energy Unit Type!"
            stop
        end select
        
        end function
!    !----------------------------------------------------------
        double precision function FindLengthUnit(unitName)
        implicit none 
        character(len=10) :: unitName       
        
        select case(trim(adjustl(unitName)))
          case("nm")
                FindLengthUnit = 1d-1
          case("a")
                FindLengthUnit = 1d0
          case("ang")
                FindLengthUnit = 1d0
          case("sigma")
                FindLengthUnit = 1d0                  
          case default
            write(6,*) "Error! Invalid Length Unit Type!"
            stop
        end select
 
        
        end function
!     !----------------------------------------------------------
        double precision function FindAngularUnit(unitName)
        use Constants         
        implicit none 
        character(len=10), intent(in) :: unitName       
        
        select case(trim(adjustl(unitName)))
          case("deg")
                FindAngularUnit = pi/180d0
          case("degree")
                FindAngularUnit = pi/180d0
          case("degrees")
                FindAngularUnit = pi/180d0                  
          case("rad")
                FindAngularUnit = 1d0
          case("radian")
                FindAngularUnit = 1d0
          case("radians")
                FindAngularUnit = 1d0                 
          case default
            write(6,*) "Error! Invalid Angular Unit Type!"
            stop
        end select
 
        
        end function          
!      !----------------------------------------------------------
      end module    
!===================================================================