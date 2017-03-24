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
      val = 1E0_dp + c_sq/d_sq - c_sq/(d_sq + (cos(theta) - h)**2) 

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
      use Constants, only: pi
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
      real(dp) :: Zeta
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

!      write(*,*) A, B, c, d, R_eq, D2, BetaPar, n, h, lam1, lam2

      rMax = R_eq + D2
      rMax_sq = rMax * rMax
!      write(*,*) rMax, R_eq - D2

      do iMol = 1, NPART(iType)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx
        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif
          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
!          write(*,*) rij, Fc_Func(rij, R_eq, D2)
          jIndx = MolArray(jType)%mol(jMol)%indx
          if(distCriteria) then
            PairList(iIndx, jIndx) = rPair(globIndx1, globIndx2)%p%r_sq
            PairList(jIndx, iIndx) = PairList(iIndx,jIndx)
          endif
          if(rij .gt. rMax) then
            cycle
          endif




          Zeta = 0E0_dp

          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
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
            if(rik .lt. rMax) then
              rxik  = rPair(globIndx1, globIndx3)%p%rx
              ryik  = rPair(globIndx1, globIndx3)%p%ry
              rzik  = rPair(globIndx1, globIndx3)%p%rz

              if(globIndx3 .gt. globIndx1) then
                rxik = -rxik
                ryik = -ryik
                rzik = -rzik
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif      

          V1 = 0E0_dp
!          write(*,*) "rij:", rij
          V1 = 0.5d0*Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
!          write(*,*) "V1:",globIndx1, globIndx2, V1/outputEConv
          PairList(iIndx, jIndx) = PairList(iIndx, jIndx) + V1
          PairList(jIndx, iIndx) = PairList(jIndx, iIndx) + V1
          ETable(iIndx) = ETable(iIndx) + V1
          ETable(jIndx) = ETable(jIndx) + V1
          E_Short = E_Short + V1
!          write(*,*) 
        enddo
      enddo
!      E_Short = 0.5E0_dp*E_Short 
      write(nout,*) "ShortRange Energy:", E_Short/outputEConv

      E_T = E_T + E_Short
      E_Inter_T = E_Short
      
      end subroutine
!======================================================================================      
      subroutine Shift_ECalc_Inter(E_Trial, disp, newDist, PairList, dETable, rejMove)
      use ForceField
      use ForceFieldPara_Tersoff
      use Coords
      use SimParameters
      use Constants, only: pi
      use PairStorage, only: distStorage, rPair, DistArrayNew, nNewDist, oldIndxArray, rPairNew, nullPair 
      implicit none
      
      type(Displacement), intent(in) :: disp(:)  
      type(DistArrayNew), intent(inout) :: newDist(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(out) :: rejMove

      
      integer :: i, iType, jType, kType, nType, iPair
      integer :: iMol, jMol, kMol, nMol
      integer :: iNei, jNei, kNei
      integer(kind=atomIntType) :: atmType1, atmType2      
      integer :: iIndx, jIndx, nIndx 
      integer :: globIndx1, globIndx2, globIndx3
      integer :: neiList(1:60), nNei
      integer :: pairIndxNew(1:6), nPair
      real(dp) :: r_sq, r, r_new,rMax, rMax_sq

      real(dp) :: A, B, c, d, R_eq, D2 
      real(dp) :: E_Short, Short
      real(dp) :: lam1, lam2
      real(dp) :: Zeta, Zeta2
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
      nType = 1
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

      nType = disp(1)%molType
      nMol = disp(1)%molIndx
      nIndx = MolArray(nType)%mol(nMol)%indx
      globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1)
      nNei = 0
      neiList = 0
      pairIndxNew = 0




!      With the Tersoff model, 
      do iPair = 1, nNewDist
        r_new = newDist(iPair)%r
        if(r_new .lt. rMax) then
          nNei = nNei + 1
          neiList(nNei) = newDist(iPair)%indx2
          pairIndxNew(nNei) = iPair
        endif
      enddo

      do jMol = 1, NPART(jType)
        if(jMol .eq. nMol) then
          cycle
        endif
        jIndx = MolArray(jType)%mol(jMol)%indx

        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        if(distCriteria) then
          PairList(jIndx) = rPairNew(globIndx2)%p%r_sq
        endif
        if(globIndx2 .eq. globIndx1) then
          cycle
        endif
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





!      This portion of the code calculates the interactions that originate from the particle that just moved.
      do jNei = 1, nNei
        jMol = neiList(jNei)
        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
        jIndx = MolArray(jType)%mol(jMol)%indx
        Short = 0E0_dp
!        Compute new position
        rij  = rPairNew(globIndx2)%p%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp
          rxij = -rPairNew(globIndx2)%p%rx
          ryij = -rPairNew(globIndx2)%p%ry
          rzij = -rPairNew(globIndx2)%p%rz 
          do kMol = 1, NPART(kType)
            if((kMol .eq. nMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
!            write(*,*) "New"
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik  = rPairNew(globIndx3)%p%r
            rjk  = rPair(globIndx2, globIndx3)%p%r
            if(rik .lt. rMax) then
              rxik  = -rPairNew(globIndx3)%p%rx
              ryik  = -rPairNew(globIndx3)%p%ry
              rzik  = -rPairNew(globIndx3)%p%rz
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!              write(*,*) jMol, nMol, kMol, angijk*(180d0/pi)
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif     

!            Angle 

            if(rjk .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx2) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
!              write(*,*) nMol, jMol, kMol, angijk*(180d0/pi)
              Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rjk, R_eq, D2)
            endif  
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1d0/(2d0*n))
          else
            b2 = 1E0_dp
          endif
!          V1 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
!          V2 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b2*B*exp(-lam2*rij))
          V1 = Fc_Func(rij, R_eq, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij))
          if(.not. distCriteria) then                
            PairList(jIndx) = PairList(jIndx) + 0.5d0*V1
          endif
          Short = Short + V1
        endif

!        write(*,*) "Old"
        rij  = rPair(globIndx1, globIndx2)%p%r
        if(rij .lt. rMax) then
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij
          endif
          do kMol = 1, NPART(kType)
            if((kMol .eq. nMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
!            write(*,*) "Old"
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            rik = rPair(globIndx1, globIndx3)%p%r
            if(rik .lt. rMax) then
              rxik  = rPair(globIndx1, globIndx3)%p%rx
              ryik  = rPair(globIndx1, globIndx3)%p%ry
              rzik  = rPair(globIndx1, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx1) then
                rxik = -rxik
                ryik = -ryik
                rzik = -rzik
              endif
              angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!              write(*,*) nMol, jMol, kMol, angijk*(180d0/pi)
              Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2)
            endif

            rjk = rPair(globIndx2, globIndx3)%p%r
            if(rjk .lt. rMax) then
              rxjk  = rPair(globIndx2, globIndx3)%p%rx
              ryjk  = rPair(globIndx2, globIndx3)%p%ry
              rzjk  = rPair(globIndx2, globIndx3)%p%rz
              if(globIndx3 .gt. globIndx2) then
                rxjk = -rxjk
                ryjk = -ryjk
                rzjk = -rzjk
              endif
              angijk = angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
!              write(*,*) jMol, nMol, kMol, angijk*(180d0/pi)
              Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rjk, R_eq, D2)
            endif  
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1d0/(2d0*n))
          else
            b2 = 1E0_dp
          endif
!          V1 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
!          V2 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b2*B*exp(-lam2*rij))
          V1 = Fc_Func(rij, R_eq, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij))
          Short = Short - V1
        endif

!        write(*,*) globIndx1, globIndx2, Short
        dETable(nIndx) = dETable(nIndx) + 0.5E0_dp*Short
        dETable(jIndx) = dETable(jIndx) + 0.5E0_dp*Short
        E_Short = E_Short + Short
      enddo

!      write(*,*) Short
      if(NTotal .eq. 2) then
        E_Trial = 0.5E0_dp * E_Short
        return
      endif



!      This portion of the code computes the energy of the atoms which neighbored
!      the particle that moved. 
      do iNei = 1, nNei
        iMol = neiList(iNei)
        globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
        iIndx = MolArray(iType)%mol(iMol)%indx

        do jMol = 1, NPART(jType)
          if(iMol .eq. jMol) then
            cycle
          endif
          if(nMol .eq. jMol) then
            cycle
          endif

          globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
          rij  = rPair(globIndx1, globIndx2)%p%r
          if(rij .gt. rMax) then
            cycle
          endif

          Zeta = 0E0_dp        !Zeta for the New Configuration
          Zeta2 = 0E0_dp       !Zeta for the Old Configuration
          jIndx = MolArray(jType)%mol(jMol)%indx
 
          rxij = rPair(globIndx1, globIndx2)%p%rx
          ryij = rPair(globIndx1, globIndx2)%p%ry
          rzij = rPair(globIndx1, globIndx2)%p%rz
          if(globIndx2 .gt. globIndx1) then
            rxij = -rxij
            ryij = -ryij
            rzij = -rzij 
          endif

          do kMol = 1, nPart(kType)
            if((kMol .eq. iMol) .or. (kMol .eq. jMol)) then
              cycle
            endif
            globIndx3 = MolArray(kType)%mol(kMol)%globalIndx(1)
            if(kMol .eq. nMol) then
!              write(*,*) "New 2"
              rik  = rPairNew(globIndx1)%p%r
              if(rik .lt. rMax) then
                rxik  = rPairNew(globIndx1)%p%rx
                ryik  = rPairNew(globIndx1)%p%ry
                rzik  = rPairNew(globIndx1)%p%rz
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!                write(*,*) jMol, iMol, kMol, angijk*(180d0/pi)
                Zeta = Zeta + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
              endif

!              write(*,*) "Old 2"
              rik  = rPair(globIndx1, globIndx3)%p%r
              if(rik .lt. rMax) then
                rxik  = rPair(globIndx1, globIndx3)%p%rx
                ryik  = rPair(globIndx1, globIndx3)%p%ry
                rzik  = rPair(globIndx1, globIndx3)%p%rz
                if(globIndx3 .gt. globIndx1) then
                  rxik = -rxik
                  ryik = -ryik
                  rzik = -rzik
                endif
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!                write(*,*) iMol, jMol, kMol, angijk*(180d0/pi)
                Zeta2 = Zeta2 + gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
              endif
              
            else
              rik  = rPair(globIndx1, globIndx3)%p%r
              if(rik .lt. rMax) then
                rxik  = rPair(globIndx1, globIndx3)%p%rx
                ryik  = rPair(globIndx1, globIndx3)%p%ry
                rzik  = rPair(globIndx1, globIndx3)%p%rz
                rik   = rPair(globIndx1, globIndx3)%p%r
                if(globIndx3 .gt. globIndx1) then
                  rxik = -rxik
                  ryik = -ryik
                  rzik = -rzik
                endif
                angijk = angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                V1 = gik_Func(angijk, c, d, h) *  Fc_Func(rik, R_eq, D2) 
                Zeta = Zeta + V1
                Zeta2 = Zeta2 + V1 
              endif
            endif
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1d0/(2d0*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1d0/(2d0*n))
          else
            b2 = 1E0_dp
          endif


          V1 = Fc_Func(rij, R_eq, D2) * (B*exp(-lam2*rij))*(b2 - b1)
!          write(*,*) "Diff Formula", V1
!          V1 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
!          V2 = Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b2*B*exp(-lam2*rij))
!          write(*,*) "V1", 0.5E0_dp*V1
!          write(*,*) "V2", 0.5E0_dp*V2
!          write(*,*) "Diff", 0.5E0_dp*(V1-V2)


          dETable(iIndx) = dETable(iIndx) + 0.5d0*V1
          dETable(jIndx) = dETable(jIndx) + 0.5d0*V1
!          write(*,*) "V1", V1
          E_Short = E_Short + V1
        enddo
      enddo



      E_Trial = 0.5E0_dp*E_Short
!      write(*,*) E_Trial

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
      
       
