!===================================================================
     module SelfAdaptive
     use VarPrecision
     contains
!===================================================================
     subroutine UpdateDihedralProbabilites
     use Constants
     use CBMC_Variables
     implicit none
     integer :: iDihed, iBin
     real(dp) :: histNorm, probNorm
     real(dp) :: avg, tempConst
     real(dp) :: maxProb, minProb

     do iDihed = 1, totalDihed
       histNorm = sum(dihedData(iDihed)%Hist)
!       Compute the new distribution by averaging the new and old probability estimates
       do iBin = 0, nDihBins
         avg = dihedData(iDihed)%Prob(iBin) + dihedData(iDihed)%Hist(iBin)/histNorm
         avg = avg * 0.5E0
         dihedData(iDihed)%Prob(iBin) = avg
       enddo
!       Normalize the new distribution
       probNorm = sum(dihedData(iDihed)%Prob)
       do iBin = 0, nDihBins
         dihedData(iDihed)%Prob(iBin) = dihedData(iDihed)%Prob(iBin)/probNorm
       enddo



!       Tabulate the integral table
       dihedData(iDihed)%Integral(0) = dihedData(iDihed)%Prob(0)
       do iBin = 1, nDihBins
         dihedData(iDihed)%Integral(iBin) = dihedData(iDihed)%Integral(iBin-1) + dihedData(iDihed)%Prob(iBin)
       enddo

!       Determine the new acceptance constant
       dihedData(iDihed)%accConst = 0d0
       do iBin = 0, nDihBins
         tempConst = dihedData(iDihed)%Hist(iBin)/(dihedData(iDihed)%Prob(iBin)*histNorm)
         if(tempConst .gt. dihedData(iDihed)%accConst) then
           dihedData(iDihed)%accConst = tempConst
         endif
       enddo

 

       write(2,*) iDihed, dihedData(iDihed)%accConst 
       do iBin = 0, nDihBins
         write(2,*) iBin*diBinSize, dihedData(iDihed)%Prob(iBin)
       enddo
       write(2,*)
     enddo

     end subroutine
!=================================================================== 
     end module
