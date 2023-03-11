!!! recoded in f90 by aiyyer
!!! also see: ftp://texmex.mit.edu/pub/emanuel/TCMAX/
!!!======================================================================
SUBROUTINE PCMIN3(SST_IN,PSL_IN,P,T_IN,R_IN,IMAX,JMAX,NA,PMIN,VMAX,IFL)
       implicit none
!
!   Revised on 9/24/2005 to fix convergence problems at high pressure
!
!   ***   This subroutine calculates the maximum wind speed        ***
!   ***             and mimimum central pressure                   ***
!   ***    achievable in tropical cyclones, given a sounding       ***
!   ***             and a sea surface temperature.                 ***
!
!  INPUT:   SST: Sea surface temperature in K
!
!           PSL: Sea level pressure (Pa)
!
!           P: One-dimensional array of dimension NA containing
!              the list of pressure levels of the 3D data
!
!           T, R: 3D arrays (NA,IMAX,JMAX) of Temperature (K),
!             and mixing ratio (kg/kg). The arrays MUST be
!             arranged so that the lowest index corresponds
!             to the lowest model level, with increasing index
!             corresponding to decreasing pressure. The temperature
!             sounding should extend to at least the tropopause and 
!             preferably to the lower stratosphere, however the
!             mixing ratios are not important above the boundary
!             layer. Missing mixing ratios can be replaced by zeros.
!
!           NA: The vertical dimension of P,T and R
!           IMAX,JMAX: nx,ny dimensions of arrays T and R
!           N:  The actual number of points in the sounding
!                (N is less than or equal to NA)
!
!  OUTPUT:  PMIN is the minimum central pressure, in mb
!
!           VMAX is the maximum surface wind speed, in m/s
!                  (reduced to reflect surface drag)
!
!
!           IFL is a flag: A value of 1 means OK; a value of 0
!              indicates no convergence (hypercane); a value of 2
!              means that the CAPE routine failed.
!
!-----------------------------------------------------------------------------
       integer, intent(in)  :: NA, IMAX, JMAX
       real(4), intent(in)  :: SST_IN(IMAX,JMAX),PSL_IN(IMAX,JMAX),P(NA),T_IN(NA,IMAX,JMAX),R_IN(NA,IMAX,JMAX)
       real(4), intent(out) :: PMIN(IMAX,JMAX),VMAX(IMAX,JMAX)
       integer, intent(out) :: IFL
       real(4) :: SST, PSL, T(NA), R(NA)
       integer :: i,j,k, N


       N=NA
       do i=1,IMAX
        do j=1,JMAX                        !
          PSL=PSL_IN(i,j)/100.0            !SLP now in mb
          SST=SST_IN(i,j)-273.15           !SST now in C
          do k=1,NA
             T(k)=T_IN(k,i,j)-273.15       !T now in C
             R(k)=R_IN(k,i,j)*1000.0       !R now in g/kg
          enddo
          call PCMIN(SST,PSL,P,T,R,NA,N,PMIN(i,j),VMAX(i,j),IFL)
        enddo
       enddo
       RETURN
       END












!!$
!!$ f90 version of Emanuel's pcmin code
!!$
SUBROUTINE PCMIN(SST,PSL,P,T,R,NA,N,PMIN,VMAX,IFL)
!!$
!!$   Revised on 9/24/2005 to fix convergence problems at high pressure
!!$
!!$   ***   This subroutine calculates the maximum wind speed        ***
!!$   ***             and mimimum central pressure                   ***
!!$   ***    achievable in tropical cyclones, given a sounding       ***
!!$   ***             and a sea surface temperature.                 ***
!!$
!!$  INPUT:   SST: Sea surface temperature in C
!!$
!!$           PSL: Sea level pressure (mb)
!!$
!!$           P,T,R: One-dimensional arrays of dimension NA
!!$             containing pressure (mb), temperature (C),
!!$             and mixing ratio (g/kg). The arrays MUST be
!!$             arranged so that the lowest index corresponds
!!$             to the lowest model level, with increasing index
!!$             corresponding to decreasing pressure. The temperature
!!$             sounding should extend to at least the tropopause and 
!!$             preferably to the lower stratosphere, however the
!!$             mixing ratios are not important above the boundary
!!$             layer. Missing mixing ratios can be replaced by zeros.
!!$
!!$           NA: The dimension of P,T and R
!!$
!!$           N:  The actual number of points in the sounding
!!$                (N is less than or equal to NA)
!!$
!!$  OUTPUT:  PMIN is the minimum central pressure, in mb
!!$
!!$           VMAX is the maximum surface wind speed, in m/s
!!$                  (reduced to reflect surface drag)
!!$
!!$           IFL is a flag: A value of 1 means OK; a value of 0
!!$              indicates no convergence (hypercane); a value of 2
!!$              means that the CAPE routine failed.
!!$
!!$-----------------------------------------------------------------------------

  implicit none
  
  real(4), intent(in) :: SST, PSL
  integer(4), intent(in) :: NA, N
  real(4),  dimension(NA)  :: T, P, R
  real(4), intent(out) :: PMIN,VMAX
  integer(4), intent(out) ::IFL
  
!!$   ***   Adjustable constant: Ratio of C_k to C_D    ***  
  real(4), parameter :: CKCD=0.9
  
!!$ Sig: Adjustable constant for buoyancy of displaced parcels
!!$ 0=Reversible ascent;  1=Pseudo-adiabatic ascent
  real(4), parameter :: SIG=0.0
  
!!$   ***  Adjustable switch: if IDISS = 0, no dissipative heating is   ***
!!$   ***     allowed; otherwise, it is                                 ***
!!$
!!$   *** Set level (NK) from which parcels lifted   ***
!!$
  integer(4), parameter :: IDISS=1, NK=1
  
!!$
!!$   ***  Exponent, b, in assumed profile of azimuthal velocity in eye,   ***
!!$   ***   V=V_m(r/r_m)^b. Used only in calculation of central pressure   ***
!!$
  real(4), parameter ::  b=2.0
  
!!$
!!$   *** Factor to reduce gradient wind to 10 m wind
!!$
  real(4), parameter ::  VREDUC=0.8
  
  
  real :: ES0, SSTK, PM, TP, RP, PP, TVAV, TV1, TOM, TOA, TOMS, RS0
  real :: PNEW, FAC, CATFAC, CAT, CAPEMS, CAPEM, CAPEA, RAT
  integer(4):: NP, iflag
  logical :: TestConvergence
  
!!$------------------------------------------------------------------------------------
  
!!$
!!$   ***   Normalize certain quantities   ***
!!$
  SSTK=SST+273.15
  ES0=6.112*EXP(17.67*SST/(243.5+SST))
  
  R = R*.001
  T = T + 273.15
  
!!$
!!$   ***   Default values   ***
!!$
  
  
  VMAX=0.0
  PMIN=PSL 
  IFL=1
  
  NP=0
  PM=950.0
!!$
!!$   ***   Find environmental CAPE *** 
!!$
  TP=T(NK)
  RP=R(NK)
  PP=P(NK)
  
  CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEA,TOA,IFLAG)
  
  IF(IFLAG.NE.1) IFL=2
  
!!$  
!!$   ***   Begin iteration to find mimimum pressure   ***
!!$


!!$ a flag to check if we need to test for convergence  
  TestConvergence = .True. 
  do while (TestConvergence)
     
!!$  
!!$     ***  Find CAPE at radius of maximum winds   ***
!!$  
     
     
     TP=T(NK)
     PP=MIN(PM,1000.0)
     RP=0.622*R(NK)*PSL/(PP*(0.622+R(NK))-R(NK)*PSL)
     CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEM,TOM,IFLAG) 
     IF(IFLAG /= 1)IFL=2
     
!!$
!!$  ***  Find saturation CAPE at radius of maximum winds   ***
!!$
     
     TP=SSTK
     PP=MIN(PM,1000.0)
     RP=0.622*ES0/(PP-ES0)
     CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEMS,TOMS,IFLAG)
     IF(IFLAG.NE.1)IFL=2
     RAT=SSTK/TOMS
     IF (IDISS == 0) RAT=1.0
     
!!$
!!$  ***  Initial estimate of minimum pressure   ***
!!$
     
     RS0=RP
     TV1=T(1)*(1.+R(1)/0.622)/(1.+R(1))
     TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
!!$       CAT=0.5*CKCD*RAT*(CAPEMS-CAPEM)
     CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
     CAT=MAX(CAT,0.0)
     PNEW=PSL*EXP(-CAT/(287.04*TVAV))
!!$
!!$   ***  Test for convergence   ***
!!$
     
     IF(ABS(PNEW-PM) > 0.2) THEN
        PM=PNEW
        NP=NP+1
        IF(NP > 1000 .OR. PM< 400.0)THEN
           PMIN=PSL
           IFL=0
           TestConvergence = .False. ! exit loop failed convergence
        END IF
     ELSE
        CATFAC=0.5*(1.+1./b)
!!$        CAT=CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
        CAT=CAPEM-CAPEA+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
        CAT=MAX(CAT,0.0)
        PMIN=PSL*EXP(-CAT/(287.04*TVAV))
        TestConvergence = .False.    ! exit loop passed convergence
     END IF
     
  end do
  

!----------------------------------------------------------
  
  FAC=MAX(0.0,(CAPEMS-CAPEM))
  VMAX=VREDUC*SQRT(CKCD*RAT*FAC)
!!$
!!$   ***  Renormalize sounding arrays   ***
!!$
  
  
  R = R*1000.
  T = T - 273.15
  
  RETURN
END SUBROUTINE PCMIN



!!$ ==============================================================================


SUBROUTINE CAPE(TP,RP,PP,T,R,P,ND,N,SIG,CAPED,TOB,IFLAG)
!!$
!!$     This subroutine calculates the CAPE of a parcel with pressure PP (mb), 
!!$       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
!!$       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
!!$       of pressure (P in mb). ND is the dimension of the arrays T,R and P,
!!$       while N is the actual number of points in the sounding. CAPED is
!!$       the calculated value of CAPE and TOB is the temperature at the
!!$       level of neutral buoyancy.  IFLAG is a flag
!!$       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
!!$       not run owing to improper sounding (e.g.no water vapor at parcel level).
!!$       IFLAG=2 indicates that routine did not converge.                 
!!$
  
  implicit none
  integer(4), intent(in) :: ND, N
  
  real(4), intent (in) :: TP, RP, PP
  real(4), intent(in), dimension (ND) :: T,R,P
  real(4), dimension (100) :: TVRDIF
  real(4), intent (out) :: CAPED, TOB
  real(4) :: CPVMCL, EPS
  integer(4) :: IFLAG, NCMAX
  
!!$    CL=4190.0
  real(4), parameter :: CPD=1005.7, CPV=1870.0, CL=2500.0, RV=461.5
  real(4), parameter :: RD=287.04, ALV0=2.501E6
  
  real(4) :: TPC, ESP, EVP, RH, ALV, S, CHI, PLCL, ES, ENEW, EM, AP
  real(4) :: TLVR, TJC, TGNEW, TG, TC, SL, SG, SIG
  real(4) :: RMEAN, RG, PMA, PFAC, PINB, PAT, PA, NA
  integer(4) :: JMIN, NC,  J, INB

  
  EPS=RD/RV            
  CPVMCL=CPV-CL
  
!!$
!!$   ***   Default values   ***
!!$      
  CAPED=0.0
  TOB=T(1)
  IFLAG=1
!!$
!!$   ***   Check that sounding is suitable    ***
!!$


  
  
  IF(RP <= 1.0E-6 .OR. TP < 200.0) THEN
     IFLAG=0
     RETURN
  END IF
  
  
!!$
!!$   ***  Define various parcel quantities, including reversible   ***
!!$   ***                       entropy, S.                         ***
!!$                           
  
  
  TPC=TP-273.15
  ESP=6.112*EXP(17.67*TPC/(243.5+TPC))
  EVP=RP*PP/(EPS+RP)
  RH=EVP/ESP
  RH=MIN(RH,1.0)
  ALV=ALV0+CPVMCL*TPC
  S=(CPD+RP*CL)*LOG(TP)-RD*LOG(PP-EVP)+ALV*RP/TP-RP*RV*LOG(RH)            
  
  
!!$
!!$   ***  Find lifted condensation pressure, PLCL   ***
!!$     
  CHI=TP/(1669.0-122.0*RH-TP)
  PLCL=PP*(RH**CHI)
!!$
!!$   ***  Begin updraft loop   ***
!!$
  
  NCMAX=0
  TVRDIF(1:N) = 0.0
  
  JMIN=1E6
  
  do J=1,N
!!$
!!$    ***   Don't bother lifting parcel above 60 mb and skip sections of sounding below parcel level  ***
!!$  
     if (P(J).LT.59.0.OR.P(J).GE.PP) CYCLE
!!$
     JMIN=MIN(JMIN,J)
!!$
!!$    ***  Parcel quantities below lifted condensation level   ***
!!$        
     IF(P(J).GE.PLCL)THEN
        TG=TP*(P(J)/PP)**(RD/CPD)
        RG=RP
!!$
!!$   ***   Calculate buoyancy   ***
!!$  
        TLVR=TG*(1.+RG/EPS)/(1.+RG)
        TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
     ELSE
!!$
!!$   ***  Parcel quantities above lifted condensation level  ***
!!$        
        TG=T(J)          
        TJC=T(J)-273.15 
        ES=6.112*EXP(17.67*TJC/(243.5+TJC)) 
        RG=EPS*ES/(P(J)-ES)
!!$
!!$   ***  Iteratively calculate lifted parcel temperature and mixing   ***
!!$   ***                ratio for reversible ascent                    ***
!!$
        NC=0
120     CONTINUE
        NC=NC+1
!!$
!!$   ***  Calculate estimates of the rates of change of the entropy    ***
!!$   ***           with temperature at constant pressure               ***
!!$  
        ALV=ALV0+CPVMCL*(TG-273.15)
        SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
        EM=RG*P(J)/(EPS+RG)
        SG=(CPD+RP*CL)*LOG(TG)-RD*LOG(P(J)-EM)+ ALV*RG/TG
        IF(NC < 3)THEN
           AP=0.3
        ELSE
           AP=1.0
        END IF
        TGNEW=TG+AP*(S-SG)/SL  
!!$
!!$   ***   Test for convergence   ***
!!$
        IF(ABS(TGNEW-TG)>0.001)THEN
           TG=TGNEW
           TC=TG-273.15
           ENEW=6.112*EXP(17.67*TC/(243.5+TC))
!!$
!!$   ***   Bail out if things get out of hand   ***
!!$
           IF(NC>500.OR.ENEW>(P(J)-1.0))THEN
              IFLAG=2
              RETURN
           END IF
           RG=EPS*ENEW/(P(J)-ENEW)           
           GOTO 120
        END IF
        NCMAX=MAX(NC,NCMAX)
!!$
!!$   *** Calculate buoyancy   ***
!!$
        RMEAN=SIG*RG+(1.-SIG)*RP
        TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
        TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
     END IF
  end do
  
!!$
!!$  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
!!$
  NA=0.0
  PA=0.0
!!$
!!$   ***  Find maximum level of positive buoyancy, INB    ***
!!$
  INB=1
  do J=N,JMIN,-1
     IF(TVRDIF(J)>0.0)INB=MAX(INB,J)
  end do
  
  IF(INB.EQ.1)  then
     RETURN
  end if
!!$
!!$   ***  Find positive and negative areas and CAPE  ***
!!$
  IF (INB>1) THEN
     do J=JMIN+1,INB
        PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1))
        PA=PA+MAX(PFAC,0.0)
        NA=NA-MIN(PFAC,0.0)
     end DO
!!$
!!$   ***   Find area between parcel pressure and first level above it ***
!!$
     PMA=(PP+P(JMIN)) 
     PFAC=RD*(PP-P(JMIN))/PMA
     PA=PA+PFAC*MAX(TVRDIF(JMIN),0.0)
     NA=NA-PFAC*MIN(TVRDIF(JMIN),0.0)
!!$
!!$   ***   Find residual positive area above INB and TO  ***
!!$
     PAT=0.0
     TOB=T(INB)
     IF(INB<N)THEN
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/ &
             (TVRDIF(INB)-TVRDIF(INB+1))
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
        TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/ &
             (P(INB)-P(INB+1))
     END IF
!!$
!!$   ***   Find CAPE  ***
!!$            
     CAPED=PA+PAT-NA
     CAPED=MAX(CAPED,0.0)
  END IF
  
  RETURN





  
end SUBROUTINE CAPE

