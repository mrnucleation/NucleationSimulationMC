      module Rosenbluth_Functions
      contains
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
      pure subroutine Rosen_BoltzWeight_Molecule_New(nRosen, nType, included,  E_Trial, overlap)
      use ForceField
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen      
      
      logical, intent(out) :: overlap
      real(kind(0.0d0)), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(kind(0.0d0)) :: rx,ry,rz,r
      real(kind(0.0d0)) :: ep,sig_sq,q
      real(kind(0.0d0)) :: LJ, Ele
      real(kind(0.0d0)) :: E_Ele,E_LJ
      real(kind(0.0d0)) :: rmin_ij

      
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      overlap = .false.

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
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      real(kind(0.0d0)), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      
      real(kind(0.0d0)), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(kind(0.0d0)) :: rx,ry,rz,r
      real(kind(0.0d0)) :: ep,sig_sq,q
      real(kind(0.0d0)) :: LJ, Ele
      real(kind(0.0d0)) :: E_Ele,E_LJ
      real(kind(0.0d0)) :: rmin_ij

      
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0

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
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      
      logical, intent(inout) :: overlap
      real(kind(0.0d0)), intent(out) :: E_Trial
      
      integer :: jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(kind(0.0d0)) :: rx,ry,rz,r
      real(kind(0.0d0)) :: ep,sig_sq,q
      real(kind(0.0d0)) :: LJ, Ele
      real(kind(0.0d0)) :: E_Ele,E_LJ
      real(kind(0.0d0)) :: rmin_ij

      
      E_LJ = 0d0
      E_Ele = 0d0      
      E_Trial = 0d0
      overlap = .false.

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
      use Coords
      use SimParameters
      implicit none
      
      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nAtom, nMol
      real(kind(0.0d0)), intent(out) :: E_Trial
      
      integer :: iAtom, jType, jIndx, jMol, jAtom
      integer(kind=2) :: atmType1,atmType2
      real(kind(0.0d0)) :: rx,ry,rz,r
      real(kind(0.0d0)) :: ep,sig_sq,q
      real(kind(0.0d0)) :: LJ, Ele
      real(kind(0.0d0)) :: E_Ele,E_LJ
      real(kind(0.0d0)) :: rmin_ij

      
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
!======================================================================================      
!      pure subroutine Rosen_BoltzWeight_IntraNonBonded_New(trialPos, regrown, E_Trial)
!      use ForceField
!      use Coords
!      use SimParameters
!      implicit none
!      real(kind(0.0d0)), intent(out) :: E_Trial
!      type(SimpleAtomCoords), intent(in) :: trialPos
!     
!      integer :: iType,iPair, iAtom, jAtom
!      integer(kind=2) :: atmType1,atmType2
!      real(kind(0.0d0)) :: rx,ry,rz,r
!      real(kind(0.0d0)) :: ep,sig_sq,q
!      real(kind(0.0d0)) :: LJ, Ele
!      real(kind(0.0d0)) :: E_Ele,E_LJ
!
!      E_LJ = 0d0
!      E_Ele = 0d0      
!      E_Trial = 0d0
!      iType = newMol%molType
!      do iPair = 1, nIntraNonBond(iType)
!        iAtom = nonBondArray(iType,iPair)%nonMembr(1)
!        jAtom = nonBondArray(iType,iPair)%nonMembr(2)
!        if(iAtom .ne. nAtom) then
!          if(jAtom .ne. nAtom) then
!            cycle
!          endif
!        endif
!
!        if(iAtom .eq. nAtom) then
!          if(regrown(jAtom) .eqv. .false.) then
!            cycle
!          endif
!        else
!          if(regrown(iAtom) .eqv. .false.) then
!            cycle
!          endif
!        endif
!      
!
!        atmType1 = atomArray(iType,iAtom)
!        atmType2 = atomArray(iType,jAtom)
!        ep = ep_tab(atmType1,atmType2)
!        sig_sq = sig_tab(atmType1,atmType2)
!        q = q_tab(atmType1,atmType2)
!        
!        rx = newMol%x(iAtom) - newMol%x(jAtom)
!        ry = newMol%y(iAtom) - newMol%y(jAtom)
!        rz = newMol%z(iAtom) - newMol%z(jAtom)
!        r = rx**2 + ry**2 + rz**2
!        if(ep .ne. 0d0) then
!          LJ = (sig_sq/r)
!          LJ = LJ * LJ * LJ
!          LJ = ep * LJ * (LJ-1d0)
!          E_LJ = E_LJ + LJ
!        endif
!        if(q .ne. 0d0) then            
!          r = dsqrt(r)
!          Ele = q/r
!          E_Ele = E_Ele + Ele            
!        endif
!      enddo
!
!      E_Trial = E_LJ + E_Ele
!      
!      
!      end subroutine    
!======================================================================================      
      end module 


