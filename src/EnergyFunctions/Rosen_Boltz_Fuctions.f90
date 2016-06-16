      module Rosenbluth_Functions
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
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_Trial = 0d0
      overlap = .false.
      
      E_LJ = 0d0
      E_Ele = 0d0      

      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jMol = 1,NPART(jType)
            jIndx = molArray(jType)%mol(jMol)%indx              
            if(included(jIndx) .eqv. .false.) then
              cycle
            endif              
            do jAtom = 1,nAtoms(jType)        
              atmType2 = atomArray(jType, jAtom)
              ep = ep_tab(atmType2, atmType1)
              q = q_tab(atmType2, atmType1)
              if(q .eq. 0.0d0) then              
                if(ep .eq. 0.0d0) then              
                  cycle
                endif
              endif
              sig_sq = sig_tab(atmType2, atmType1)
              rmin_ij = r_min_tab(atmType2, atmType1)
        
              rx = rosenTrial(nRosen)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = rosenTrial(nRosen)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = rosenTrial(nRosen)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              if(r .lt. rmin_ij) then
                overlap = .true.
              endif              
              LJ = 0d0
              Ele = 0d0
              if(ep .ne. 0d0) then
                LJ = (sig_sq / r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ - 1d0)                
                E_LJ = E_LJ + LJ
              endif
              if(q .ne. 0d0) then
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
              endif
            enddo
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
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


      do iAtom = 1,nAtoms(nType)
        atmType1 = atomArray(nType, iAtom)
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            if(q .eq. 0.0d0) then              
              if(ep .eq. 0.0d0) then              
                cycle
              endif
            endif
            sig_sq = sig_tab(atmType2, atmType1)
            rmin_ij = r_min_tab(atmType2, atmType1)
            do jMol = 1,NPART(jType)
              jIndx = molArray(jType)%mol(jMol)%indx              
              if(included(jIndx) .eqv. .false.) then
                cycle
              endif         
              rx = mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
              ry = mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
              rz = mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
              r = rx*rx + ry*ry + rz*rz
              LJ = 0d0
              Ele = 0d0
              if(ep .ne. 0d0) then
                LJ = (sig_sq / r)
                LJ = LJ * LJ * LJ              
                LJ = ep * LJ * (LJ - 1d0)                
                E_LJ = E_LJ + LJ
              endif
              if(q .ne. 0d0) then
                r = sqrt(r)
                Ele = q / r
                E_Ele = E_Ele + Ele
              endif
            enddo
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      subroutine Rosen_BoltzWeight_Atom_New(nType, nAtom, trialPos, included,  E_Trial, overlap)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      overlap = .false.

      atmType1 = atomArray(nType, nAtom)
      do jType = 1, nMolTypes
        do jAtom = 1,nAtoms(jType)        
          atmType2 = atomArray(jType, jAtom)
          ep = ep_tab(atmType2, atmType1)
          q = q_tab(atmType2, atmType1)
          if(q .eq. 0.0d0) then              
            if(ep .eq. 0.0d0) then              
              cycle
            endif
          endif
          sig_sq = sig_tab(atmType2, atmType1)
          rmin_ij = r_min_tab(atmType2, atmType1)
          do jMol = 1,NPART(jType)
            jIndx = molArray(jType)%mol(jMol)%indx              
            if(included(jIndx) .eqv. .false.) then
              cycle
            endif              

            rx = trialPos%x - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = trialPos%y - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = trialPos%z - MolArray(jType)%mol(jMol)%z(jAtom)
            r = rx*rx + ry*ry + rz*rz
            if(r .lt. rmin_ij) then
              overlap = .true.
            endif              
!            if(r .ge. 25d0) then
!              cycle
!            endif
            LJ = 0d0
            Ele = 0d0
            if(ep .ne. 0d0) then
              LJ = (sig_sq / r)
              LJ = LJ * LJ * LJ              
              LJ = ep * LJ * (LJ - 1d0)                
              E_LJ = E_LJ + LJ
            endif
            if(q .ne. 0d0) then
              r = sqrt(r)
              Ele = q / r
              E_Ele = E_Ele + Ele
            endif
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_BoltzWeight_Atom_Old(nType, nMol, nAtom, included,  E_Trial)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(dp), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0

      atmType1 = atomArray(nType, nAtom)
      do jType = 1, nMolTypes
        do jMol = 1,NPART(jType)
          jIndx = molArray(jType)%mol(jMol)%indx              
          if(included(jIndx) .eqv. .false.) then
            cycle
          endif              
          do jAtom = 1,nAtoms(jType)        
            atmType2 = atomArray(jType, jAtom)
            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            if(q .eq. 0.0d0) then              
              if(ep .eq. 0.0d0) then              
                cycle
              endif
            endif
            sig_sq = sig_tab(atmType2, atmType1)
            rmin_ij = r_min_tab(atmType2, atmType1)
             
            rx = MolArray(nType)%mol(nMol)%x(nAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
            ry = MolArray(nType)%mol(nMol)%y(nAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
            rz = MolArray(nType)%mol(nMol)%z(nAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
            r = rx*rx + ry*ry + rz*rz
!            if(r .ge. 25d0) then
!              cycle
!            endif
            LJ = 0d0
            Ele = 0d0
            if(ep .ne. 0d0) then
              LJ = (sig_sq / r)
              LJ = LJ * LJ * LJ              
              LJ = ep * LJ * (LJ - 1d0)                
              E_LJ = E_LJ + LJ
            endif
            if(q .ne. 0d0) then
              r = sqrt(r)
              Ele = q / r
              E_Ele = E_Ele + Ele
            endif
          enddo
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      end subroutine 
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which each atom is regrown sequentially for the given trial.
      subroutine Rosen_BoltzWeight_Atom_Intra_New(nType, nAtom, trialPos, regrown,  E_Trial, overlap)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(dp), intent(out) :: E_Trial
      
      integer :: temp, iPair, jIndx, iAtom, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType,iPair)%nonMembr(1)
        jAtom = nonBondArray(nType,iPair)%nonMembr(2)
         
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif

        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif

        atmType1 = atomArray(nType,iAtom)
        atmType2 = atomArray(nType,jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        
        rx = trialPos%x - newMol%x(jAtom)
        ry = trialPos%y - newMol%y(jAtom)
        rz = trialPos%z - newMol%z(jAtom) 
        r = rx*rx + ry*ry + rz*rz
        if(ep .ne. 0d0) then
          LJ = (sig_sq/r)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1d0)
          E_LJ = E_LJ + LJ
        endif
        if(q .ne. 0d0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine 
!======================================================================================================
      pure subroutine Rosen_BoltzWeight_Atom_Intra_Old(nType, nMol, nAtom, regrown,  E_Trial)
      use ForceField
      use ForceFieldPara_LJ_Q
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: regrown(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(dp), intent(out) :: E_Trial
      
      integer :: temp, iPair, jIndx, iAtom, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      
      do iPair = 1,nIntraNonBond(nType)
        if(all(nonBondArray(nType,iPair)%nonMembr .ne. nAtom)) then
          cycle
        endif
        iAtom = nonBondArray(nType, iPair)%nonMembr(1)
        jAtom = nonBondArray(nType, iPair)%nonMembr(2)
         
        if(iAtom .ne. nAtom) then
          temp = iAtom
          iAtom = jAtom
          jAtom = temp
        endif

        if(regrown(jAtom) .eqv. .false.) then
          cycle
        endif

        atmType1 = atomArray(nType, iAtom)
        atmType2 = atomArray(nType, jAtom)
        ep = ep_tab(atmType1,atmType2)
        sig_sq = sig_tab(atmType1,atmType2)
        q = q_tab(atmType1,atmType2)
        
        rx = MolArray(nType)%mol(nMol)%x(iAtom) - MolArray(nType)%mol(nMol)%x(jAtom)
        ry = MolArray(nType)%mol(nMol)%y(iAtom) - MolArray(nType)%mol(nMol)%y(jAtom)
        rz = MolArray(nType)%mol(nMol)%z(iAtom) - MolArray(nType)%mol(nMol)%z(jAtom)
        r = rx*rx + ry*ry + rz*rz
        if(ep .ne. 0d0) then
          LJ = (sig_sq/r)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1d0)
          E_LJ = E_LJ + LJ
        endif
        if(q .ne. 0d0) then            
          r = sqrt(r)
          Ele = q/r
          E_Ele = E_Ele + Ele            
        endif
      enddo
      E_Trial = E_LJ + E_Ele
  
      end subroutine 
!=====================================================================================    
      end module 


