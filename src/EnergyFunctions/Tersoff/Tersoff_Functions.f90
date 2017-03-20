!*********************************************************************************************************************
!     This file contains the energy functions that work the Tersoff Model
!     The prefix naming scheme implies the following:
!           Detailed - Complete energy calculation inteded for use at the beginning and end
!                      of the simulation.  This function is not inteded for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves. 
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!           Exchange - Combines the Mol and New Mol routines for moves that simultaniously add and remove a particle at the same time.
!*********************************************************************************************************************
      module InterEnergy_Tersoff
      use VarPrecision


      contains
!======================================================================================      
      pure function Fc_Func(r, R_eq, D) result(val)
      use Constants, only: pi
      implicit none
      real(dp), intent(in) :: r, R_eq, D
      real(dp) :: val  
 
      if( r .lt. (R_eq-D) ) then
        val = 1E0_dp
      elseif( r .lt. (R_eq+D) ) then
        val = 0.5E0_dp * (1E0_dp - sin(pi*(r-R_eq)/(2d0*D)) )
      else
        val = 0E0_dp
      endif

      end function
!======================================================================================      
      pure function gik_Func(theta, c, d, h) result(val)
      implicit none
      real(dp), intent(in) :: theta, c, d, h
      real(dp) :: c_sq, d_sq
      real(dp) :: val  
 
      c_sq = c * c
      d_sq = d * d
      val = 1E0_dp + c_sq/d_sq - c_sq/(d_sq + (cos(theta) + h)**2) 

      end function

!======================================================================================
      pure function angleCalc(rx12, ry12, rz12, r12, rx23, ry23, rz23, r23) result(Angle)
      implicit none
      real(dp), intent(in) :: rx12, ry12, rz12, r12, rx23, ry23, rz23, r23
      real(dp) :: Angle  

      rx12 = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
      ry12 = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
      rz12 = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)

      rx23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
      ry23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
      rz23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)          
            
      Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
      Angle = Angle/(r12*r23)
      if(abs(Angle) .gt. 1E0_dp) then
        Angle = sign(1E0_dp, Angle)
      endif
      Angle = acos(Angle)

      end function
!======================================================================================      
      subroutine Detailed_ECalc_Inter(E_T, PairList)
      use ParallelVar
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use EnergyTables
      use PairStorage, only: rPair
      implicit none
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:,:)
      integer :: iAtom, jAtom
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, globIndx3, jMolMin
      real(dp) :: r_sq, r, rMax, rMax_sq
      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxjk, ryjk, rzjk, rjk
      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_Short
      real(dp) :: Zeta1, Zeta2
      real(dp) :: Beta, n
      real(dp) :: b1, b2, V1, V2
      real(dp) :: angijk, angjik

      E_T = 0E0_dp
      E_Short = 0E0_dp
      PairList = 0E0_dp
      ETable = 0E0_dp
      iType = 1
      jType = 1
      kType = 1
      atmType1 = atomIndicies(iType, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      beta = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n

      rMax = R_eq + D2
      rMax_sq = rMax * rMax

      do iMol = 1, nPart(iType)- 1
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx
        do jMol = iMol + 1, nPart(jType)
          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          r_sq = rPair(globIndx1, globIndx2)%p%r_sq
          if(r_sq .gt. rMax_sq) then
            cycle
          endif
          jIndx = MolArray(iType)%mol(iMol)%indx
          Zeta1 = 0E0_dp
          Zeta2 = 0E0_dp
          rij  = rPair(globIndx1, globIndx2)%p%r
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx1 .gt. globIndx2) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij
          endif
          do kMol = 1, nPart(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(jType)%mol(jMol)%globalIndx(1)
            rik  = rPair(globIndx1, globIndx3)%p%r
            rjk  = rPair(globIndx2, globIndx3)%p%r
            if(rik .gt. rMax) then
              if(rjk .gt. rMax) then
                cycle
              endif
            endif

            if(rjk .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx2 .gt. globIndx3) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, -rxjk, -ryjk, -rzjk, rjk)
              Zeta1 = Zeta1 + gik_Func(angijk, c, d, h) *  Fc_Func(rij, R_eq, D)
            endif

            if(rik .lt. rMax) then
              rxik  = rPair(globIndx1, globIndx3)%p%rx
              ryik  = rPair(globIndx1, globIndx3)%p%ry
              rzik  = rPair(globIndx1, globIndx3)%p%rz
              if(globIndx1 .gt. globIndx3) then
                rxik = -rxik
                ryik = -ryik
                rzik = -rzik
              endif
              angjik = angleCalc(-rxij, -ryij, -rzij, rij, -rxik, -ryik, -rzik, rik)
              Zeta1 = Zeta1 + gik_Func(angjik, c, d, h) *  Fc_Func(rik, R_eq, D)
            endif
          enddo
          if(Zeta1 .ne. 0E0_dp) then
            b1 = (1E0_dp + (beta*Zeta1)**n)**(-1d0/(2d0*n))
          endif
!          if(Zeta2 .ne. 0E0_dp) then
!            b2 = (1E0_dp + (beta*Zeta2)**n)**(-1d0/(2d0*n))
!          endif
          

          V1 = Fc_Func(rij, R_eq, D) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
          PairList(iIndx, jIndx) = V1
          PairList(jIndx, iIndx) = V1
          ETable(iIndx) = ETable(iIndx) + V1
          ETable(jIndx) = ETable(jIndx) + V1
          E_Short = E_Short + V1
        enddo
      enddo

      write(nout,*) "ShortRange Energy:", E_Short

      E_T = E_T + E_Short
      E_Inter_T = E_Short
      
      end subroutine
!======================================================================================      
      pure subroutine Shift_ECalc_Inter(E_Trial,disp,newDist, PairList,dETable,rejMove)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use PairStorage, only: distStorage, DistArrayNew, nNewDist, oldIndxArray
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove

      
      integer :: iType,jType,iMol,jMol,iAtom,jAtom, iPair
      integer(kind=atomIntType) :: atmType1,atmType2,iIndx,jIndx
      integer :: sizeDisp
!      integer, pointer :: oldIndx 
      integer :: gloIndx1, gloIndx2
      real(dp) :: r_new, r_new_sq   
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele, E_PairOld, E_Old
      real(dp) :: E_Ele,E_LJ


      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      E_Old = 0E0
      PairList = 0E0      
      rejMove = .false.
!      dETable = 0E0
!      if(NTotal .eq. 1) return
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      
!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary

      do iPair = 1, nNewDist
        gloIndx2 = newDist(iPair)%indx2
        jType = atomIndicies(gloIndx2)%nType
        jMol  = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(jIndx .ne. iIndx) then
          gloIndx1 = newDist(iPair)%indx1
          jAtom = atomIndicies(gloIndx2)%nAtom
          iAtom = atomIndicies(gloIndx1)%nAtom
       
          atmType1 = atomArray(iType,iAtom)
          atmType2 = atomArray(jType,jAtom)

          ep = ep_tab(atmType2, atmType1)
          q  = q_tab(atmType2, atmType1)
         
!          if(distCriteria) then
!            if(iAtom .eq. 1) then
!              if(jAtom .eq. 1) then
!                PairList(jIndx) = newDist(iPair)%r_sq
!              endif
!            endif
!          endif

          LJ = 0E0
          Ele = 0E0
          if(ep .ne. 0E0) then
            sig_sq = sig_tab(atmType2,atmType1)
            r_new_sq = newDist(iPair)%r_sq
            LJ = LJ_Func(r_new_sq, ep, sig_sq)             
            E_LJ = E_LJ + LJ
            if(.not. distCriteria) then
              PairList(jIndx) = PairList(jIndx) + LJ
            endif
!            dETable(iIndx) = dETable(iIndx) + LJ
!            dETable(jIndx) = dETable(jIndx) + LJ
!            newDist(iPair)%E_Pair = newDist(iPair)%E_Pair + LJ
          endif
          if(q .ne. 0E0) then
            r_new = newDist(iPair)%r
            Ele = q/r_new
!            Ele = Ele_Func(r_new, q)                
            E_Ele = E_Ele + Ele
            if(.not. distCriteria) then                
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
!            dETable(iIndx) = dETable(iIndx) + Ele
!            dETable(jIndx) = dETable(jIndx) + Ele
!            newDist(iPair)%E_Pair = newDist(iPair)%E_Pair + Ele
          endif
          newDist(iPair)%E_Pair = Ele + LJ
          E_PairOld = distStorage(oldIndxArray(iPair))%E_Pair
          dETable(iIndx) = dETable(iIndx) + LJ + Ele - E_PairOld
          dETable(jIndx) = dETable(jIndx) + LJ + Ele - E_PairOld  
          E_Old = E_Old + E_PairOld
        endif
      enddo




      sizeDisp = size(disp)
      if(.not. distCriteria) then      
        if(sizeDisp .lt. nAtoms(iType)) then
          call Shift_PairList_Correct(disp, PairList)
        endif
      endif
     
      E_Trial = E_LJ + E_Ele - E_Old
      
      
      end subroutine 
!======================================================================================      
      pure subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use PairStorage
      implicit none
      integer, intent(in) :: iType, iMol     
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:)
      
      integer :: iAtom,iIndx,jType,jIndx,jMol,jAtom
      integer  :: gloIndx1, gloIndx2
      real(dp) :: E_Pair

      E_Trial = 0E0
      dETable = 0E0
      iIndx = MolArray(iType)%mol(iMol)%indx

      do iAtom = 1,nAtoms(iType)
        gloIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom) 
        do jType = 1, nMolTypes
          do jAtom = 1,nAtoms(jType)        
            do jMol=1, NPART(jType)
              jIndx = MolArray(jType)%mol(jMol)%indx  
              if(iIndx .ne. jIndx) then
                gloIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom) 
                E_Pair = rPair(gloIndx1, gloIndx2)%p%E_Pair
                E_Trial = E_Trial + E_Pair
                dETable(iIndx) = dETable(iIndx) + E_Pair
                dETable(jIndx) = dETable(jIndx) + E_Pair
              endif
            enddo
          enddo
        enddo
      enddo

  
      
      end subroutine
!======================================================================================      
      subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use PairStorage
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iPair
      integer :: iType, iMol, iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1,atmType2
      integer :: gloIndx1, gloIndx2
      real(dp) :: r_sq, r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele

      real(dp) :: E_Ele,E_LJ

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      dETable = 0E0
      PairList = 0E0
      rejMove = .false.
      
      iType = newMol%molType
      iMol = NPART(iType)+1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iPair = 1, nNewDist
        gloIndx1 = newDist(iPair)%indx1
        gloIndx2 = newDist(iPair)%indx2

        jType = atomIndicies(gloIndx2)%nType
        jMol = atomIndicies(gloIndx2)%nMol
        jIndx = MolArray(jType)%mol(jMol)%indx
        if(iIndx .ne. jIndx) then
          iAtom = atomIndicies(gloIndx1)%nAtom
          jAtom = atomIndicies(gloIndx2)%nAtom

          atmType1 = atomArray(iType, iAtom)
          atmType2 = atomArray(jType, jAtom)

          ep = ep_tab(atmType2, atmType1)
          q = q_tab(atmType2, atmType1)

          LJ = 0d0
          Ele = 0d0
          if(ep .ne. 0E0) then
            r_sq = newDist(iPair)%r_sq
            sig_sq = sig_tab(atmType2,atmType1)
            LJ = LJ_Func(r_sq, ep, sig_sq)             
            E_LJ = E_LJ + LJ
            if(.not. distCriteria) then
              PairList(jIndx) = PairList(jIndx) + LJ
            endif
!            dETable(iIndx) = dETable(iIndx) + LJ
!            dETable(jIndx) = dETable(jIndx) + LJ
          endif
          if(q .ne. 0E0) then
            r = newDist(iPair)%r
            Ele = q/r
!            Ele = Ele_Func(r_sq, q)                
            E_Ele = E_Ele + Ele
            if(.not. distCriteria) then                
              PairList(jIndx) = PairList(jIndx) + Ele
            endif
!            dETable(iIndx) = dETable(iIndx) + Ele
!            dETable(jIndx) = dETable(jIndx) + Ele
          endif
          dETable(iIndx) = dETable(iIndx) + LJ + Ele
          dETable(jIndx) = dETable(jIndx) + LJ + Ele
          newDist(iPair)%E_Pair = LJ + Ele
        endif 
      enddo
     
      E_Trial = E_LJ + E_Ele
      
      
      end subroutine    
!======================================================================================      
      pure subroutine Exchange_ECalc_Inter(E_Trial, nType, nMol, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
      integer, intent(in) :: nType, nMol
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      
      integer :: iAtom, newIndx, jType, jIndx, jMol, jAtom
      integer :: iIndx2
      integer(kind=atomIntType) :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Ele,E_LJ
      real(dp) :: rmin_ij

      
      end subroutine    
!======================================================================================      
      subroutine QuickNei_ECalc_Inter_LJ_Q(jType, jMol, rejMove)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      implicit none
      integer, intent(in) :: jType, jMol     
      logical, intent(out) :: rejMove
      
      integer :: iAtom,jAtom
      integer(kind=atomIntType)  :: atmType1,atmType2
      real(dp) :: rx,ry,rz,r
      real(dp) :: ep,sig_sq,q
      real(dp) :: LJ, Ele
      real(dp) :: E_Trial,E_Ele,E_LJ
      real(dp) :: rmin_ij

      E_LJ = 0E0
      E_Ele = 0E0      
      E_Trial = 0E0
      rejMove = .false.
    
      do iAtom = 1,nAtoms(newMol%molType)
        atmType1 = atomArray(newMol%molType, iAtom)
        do jAtom = 1,nAtoms(jType)        
          atmType2 = atomArray(jType, jAtom)
          rmin_ij = r_min_tab(atmType2, atmType1)
          rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
          ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
          rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
          r = rx*rx + ry*ry + rz*rz

          if(r .lt. rmin_ij) then
            rejMove = .true.
            return
          endif          
          ep = ep_tab(atmType2, atmType1)
          q = q_tab(atmType2, atmType1)
          if(ep .ne. 0E0) then
            sig_sq = sig_tab(atmType2, atmType1)
            LJ = LJ_Func(r, ep, sig_sq)
            E_LJ = E_LJ + LJ
          endif
          if(q .ne. 0E0) then   
            r = sqrt(r)         
            Ele = Ele_Func(r, q)
            E_Ele = E_Ele + Ele
          endif
        enddo
      enddo
     
      E_Trial = E_LJ + E_Ele

      if( E_Trial .gt. Eng_Critr(newMol%molType,jType) ) then
        rejMove = .true.
      endif
!      write(35,*) newMol%molType, jType, E_Trial, Eng_Critr(newMol%molType,jType), rejMove

      
      end subroutine
!======================================================================================
      end module
      
       
