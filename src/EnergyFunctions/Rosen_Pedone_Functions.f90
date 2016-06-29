      module Rosenbluth_Functions_Pedone
      contains
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_BoltzWeight_Molecule_New(nRosen, nType, included,  E_Trial, overlap)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: q
      real(dp) :: LJ, Ele, Morse
      real(dp) :: E_Ele, E_LJ, E_Morse
      real(dp) :: rmin_ij

      E_Trial = 0d0
      overlap = .false.
      
      E_LJ = 0d0
      E_Ele = 0d0      

      atmType1 = atomArray(nType, iAtom)
      do jType = 1, nMolTypes
        atmType2 = atomArray(jType, 1)
        r_eq = rEq_tab(atmType1, atmType2)
        q = q_tab(atmType1, atmType2)
        alpha = alpha_Tab(atmType1, atmType2)
        delta = D_Tab(atmType1, atmType2)
        repul_C = repul_tab(atmType1, atmType2)          
        rmin_ij = r_min_tab(atmType1, atmType2)  
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif              
          rx = rosenTrial(nRosen)%x(1) - MolArray(jType)%mol(jMol)%x(1)
          ry = rosenTrial(nRosen)%y(1) - MolArray(jType)%mol(jMol)%y(1)
          rz = rosenTrial(nRosen)%z(1) - MolArray(jType)%mol(jMol)%z(1)
          r = rx*rx + ry*ry + rz*rz
          if(r .lt. rmin_ij) then
            overlap = .true.
          endif           

          if(repul_C .ne. 0d0) then
            LJ = ( 1d0/r )**6
            LJ = repul_C * LJ
            E_LJ = E_LJ + LJ
          endif

          r = dsqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele

          if(delta .ne. 0d0) then
            Morse = exp(-alpha*(r-r_eq))
            Morse = delta*(Morse*Morse - 1d0)
            E_Morse = E_Morse + Morse
          endif
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele + E_Morse
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_BoltzWeight_Molecule_Old(mol_x, mol_y, mol_z, nType, included,  E_Trial)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_Trial = 0d0
      E_LJ = 0d0
      E_Ele = 0d0      

      atmType1 = atomArray(nType, 1)
      do jType = 1, nMolTypes
        atmType2 = atomArray(jType, 1)
        r_eq = rEq_tab(atmType1, atmType2)
        q = q_tab(atmType1, atmType2)
        alpha = alpha_Tab(atmType1, atmType2)
        delta = D_Tab(atmType1, atmType2)
        repul_C = repul_tab(atmType1, atmType2)          
        rmin_ij = r_min_tab(atmType1, atmType2) 
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif         
          rx = mol_x(1) - MolArray(jType)%mol(jMol)%x(1)
          ry = mol_y(1) - MolArray(jType)%mol(jMol)%y(1)
          rz = mol_z(1) - MolArray(jType)%mol(jMol)%z(1)
          r = rx*rx + ry*ry + rz*rz
          if(repul_C .ne. 0d0) then
            LJ = ( 1d0/r )**6
            LJ = repul_C * LJ
            E_LJ = E_LJ + LJ
          endif

          r = dsqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele

          if(delta .ne. 0d0) then
            Morse = exp(-alpha*(r-r_eq))
            Morse = delta*(Morse*Morse - 1d0)
            E_Morse = E_Morse + Morse
          endif
        enddo
      enddo

     
      E_Trial = E_LJ + E_Ele + E_Morse
      
      end subroutine 
!=====================================================================================    
      end module 


