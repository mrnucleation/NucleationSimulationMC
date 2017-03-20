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
      integer :: iType, jType, kType
      integer :: iMol, jMol, kMol
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, globIndx3
      real(dp) :: r_sq, r, rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_Short
      real(dp) :: lam1, lam2
      real(dp) :: Zeta1, Zeta2
      real(dp) :: BetaPar, n, h
      real(dp) :: b1, b2, V1, V2
      real(dp) :: angijk, angjik

      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxjk, ryjk, rzjk, rjk
      real(dp) :: rxik, ryik, rzik, rik


      E_T = 0E0_dp
      E_Short = 0E0_dp
      PairList = 0E0_dp
      ETable = 0E0_dp
      iType = 1
      jType = 1
      kType = 1
      atmType1 = atomArray(iType, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax
!      write(*,*) rMax

      do iMol = 1, NPART(iType)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx
        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif
!          write(*,*) iMol, jMol
          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
!          write(*,*) iMol, jMol
          if(rij .gt. rMax) then
            cycle
          endif

          jIndx = MolArray(jType)%mol(jMol)%indx
          Zeta1 = 0E0_dp
          Zeta2 = 0E0_dp

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
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPair(globIndx1, globIndx3)%p%r
            if(rik .gt. rMax) then
              cycle
            endif

            if(rik .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx2 .gt. globIndx3) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, -rxjk, -ryjk, -rzjk, rjk)
              Zeta1 = Zeta1 + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif
          enddo
          if(Zeta1 .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta1)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif

          
          V1 = 0E0_dp
          V1 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
          PairList(iIndx, jIndx) = V1
          PairList(jIndx, iIndx) = V1
          ETable(iIndx) = ETable(iIndx) + 0.5d0*V1
          ETable(jIndx) = ETable(jIndx) + 0.5d0*V1
          E_Short = E_Short + V1
        enddo
      enddo
      E_Short = 0.5E0_dp*E_Short 
      write(nout,*) "ShortRange Energy:", E_Short

      E_T = E_T + E_Short
      E_Inter_T = E_Short
      
      end subroutine
!======================================================================================      
      subroutine Shift_ECalc_Inter(E_Trial, disp, newDist, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use PairStorage, only: distStorage, rPair, rPairNew, DistArrayNew, nNewDist, oldIndxArray
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove

      
      integer :: i, iType, jType, kType, iPair
      integer :: iMol, jMol, kMol
      integer :: iNei, jNei, kNei
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: iIndx, jIndx, globIndx1, globIndx2, globIndx3
      integer :: neiList(1:60), nNei
      integer :: pairIndxNew(1:6), nPair
      real(dp) :: r_sq, r, r_new,rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_Short, Short
      real(dp) :: lam1, lam2
      real(dp) :: Zeta
      real(dp) :: BetaPar, n, h
      real(dp) :: b1, b2, V1, V2
      real(dp) :: angijk, angjik

      real(dp) :: rxij, ryij, rzij, rij
      real(dp) :: rxjk, ryjk, rzjk, rjk
      real(dp) :: rxik, ryik, rzik, rik

      E_Trial = 0E0_dp
      E_Short = 0E0_dp
      PairList = 0E0_dp
      dETable = 0E0_dp
      iType = 1
      jType = 1
      kType = 1
      atmType1 = atomArray(iType, 1) 

      A = tersoffData(atmType1)%A
      B = tersoffData(atmType1)%B
      c = tersoffData(atmType1)%c
      d = tersoffData(atmType1)%d
      R_eq = tersoffData(atmType1)%R
      D2 = tersoffData(atmType1)%D2
      BetaPar = tersoffData(atmType1)%beta
      n = tersoffData(atmType1)%n
      h = tersoffData(atmType1)%h
      lam1 = tersoffData(atmType1)%lam1
      lam2 = tersoffData(atmType1)%lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax

      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx
      globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
      nNei = 0
      neiList = 0
      pairIndxNew = 0

!      nNei = 1
!      neiList(1) = globIndx1

      do iPair = 1, nNewDist
        r_new = newDist(iPair)%r
        if(r_new .lt. rMax) then
          nNei = nNei + 1
          neiList(nNei) = newDist(iPair)%indx2
          pairIndxNew(nNei) = iPair
        endif
      enddo

      do jMol = 1, NPART(jType)
        if(jMol .eq. iMol) then
          cycle
        endif
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        if(rPair(globIndx1, globIndx2)%p%r .lt. rMax ) then
          if( all(neiList(1:nNei) .ne. globIndx2) ) then
            nNei = nNei + 1
            neiList(nNei) = globIndx2
            do iPair = 1, nNewDist 
              if( (globIndx2 .eq. newDist(iPair)%indx1) .or. (globIndx2 .eq. newDist(iPair)%indx2) ) then
                pairIndxNew(nNei) = iPair
                exit
              endif
            enddo
          endif
        endif
      enddo

!      do i = 1, nNei
!        write(*,*) i, globIndx1, neiList(i), pairIndxNew(i)
!      enddo


          rij  = rPair(globIndx1, globIndx2)%p%r
          if(rij .gt. rMax) then
!            cycle
          endif

!      This portion of the code calculates the interactions that use
      do jNei = 1, nNei
        jMol = neiList(jNei)
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        nPair = pairIndxNew(jNei)
        jIndx = MolArray(jType)%mol(jMol)%indx
        Short = 0E0_dp
!        Compute new position
        rij  = newDist(nPair)%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          rxij = newDist(nPair)%rx
          ryij = newDist(nPair)%ry
          rzij = newDist(nPair)%rz
          rxij = -rxij
          ryij = -ryij
          rzij = -rzij
          do kMol = 1, NPART(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPairNew(globIndx3)%p%r
            if(rik .gt. rMax) then
              cycle
            endif
            rxjk  = rPair(globIndx2, globIndx3)%p%rx
            ryjk  = rPair(globIndx2, globIndx3)%p%ry
            rzjk  = rPair(globIndx2, globIndx3)%p%rz
            if(globIndx2 .gt. globIndx3) then
              rxjk = -rxjk
              ryjk = -ryjk
              rzjk = -rzjk
            endif
            angijk = angleCalc(rxij, ryij, rzij, rij, -rxjk, -ryjk, -rzjk, rjk)
            Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif
          V1 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
          if(.not. distCriteria) then                
            PairList(jIndx) = PairList(jIndx) + V1
          endif
          Short = Short + V1
        endif

        rij  = rPair(globIndx1, globIndx2)%p%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          rxij = newDist(nPair)%rx
          ryij = newDist(nPair)%ry
          rzij = newDist(nPair)%rz
          rxij = -rxij
          ryij = -ryij
          rzij = -rzij
          do kMol = 1, NPART(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPair(globIndx2, globIndx3)%p%r
            if(rik .gt. rMax) then
              cycle
            endif
            rxjk  = rPair(globIndx2, globIndx3)%p%rx
            ryjk  = rPair(globIndx2, globIndx3)%p%ry
            rzjk  = rPair(globIndx2, globIndx3)%p%rz
            if(globIndx2 .gt. globIndx3) then
              rxjk = -rxjk
              ryjk = -ryjk
              rzjk = -rzjk
            endif
            angijk = angleCalc(rxij, ryij, rzij, rij, -rxjk, -ryjk, -rzjk, rjk)
            Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif
          V1 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
          Short = Short - V1
        endif

        dETable(iIndx) = dETable(iIndx) + Short
        dETable(jIndx) = dETable(jIndx) + Short
        E_Short = E_Short + Short
      enddo
!      write(*,*) E_Short



!      Left to do:  Add the calculations that starts from the neighbors of the mobile particle
!                   

      E_Trial = E_Short


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
      subroutine QuickNei_ECalc_Inter_Tersoff(jType, jMol, rejMove)
      implicit none
      integer, intent(in) :: jType, jMol     
      logical, intent(out) :: rejMove
      
      rejMove = .false.

      
      end subroutine
!======================================================================================
      end module
      
       
