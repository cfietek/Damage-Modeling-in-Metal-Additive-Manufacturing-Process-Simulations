CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                  UMAT OF JOHNSON-COOK HARDENING MODEL
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UMAT for use with Abaqus CAE 2021r in Fortran 77 (fixed) format         C
C     Last Edited by Carter Fietek 11.10.2022                                 C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'aba_param.inc' !Sets all variables as implicit real*8 (double precision)
C
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                  THE VARIABLES USED IN THE UMAT
C       
C--------------------------------------------------------------------
C        EELAS----ELASTIC STRAIN
C        EPLAS----PLASTIC STRAIN      
C        FLOW----FLOW STESS
C--------------------------------------------------------------------
      DIMENSION EELAS(6),EPLAS(6),FLOW(6)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0, HALF =0.5d0)
      DATA NEWTON,TOLER/1500,1.D-6/ !Default 150 and 1.D-6, change if not converging, or not results desired.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C        THE FOLLOWING IS UMAT FOR JOHNSON-COOK MODEL                         C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C                  CONSTANT USED IN THE UMAT
C
C--------------------------------------------------------------------
C     PROPS(1) - YOUNG'S MODULUS
C     PROPS(2) - POISSON RATIO
C     PROPS(3) - INELASTIC HEAT FRACTION
C     PARAMETERS OF JOHNSON-COOK MODEL:
C        PROPS(4) - A
C        PROPS(5) - B
C        PROPS(6) - n
C        PROPS(7) - C
C        PROPS(8) - m
C        PROPS(9) - rate0
C        PROPS(10) - Tr (Transition Temperature)
C        PROPS(11) - Tm (Melting Temperarure)
C        PROPS(12) - D1
C        PROPS(13) - D2
C        PROPS(14) - D3
C        PROPS(15) - D4
C        PROPS(16) - D5
C
C  State Variables
C     STATEV(1) = e11
C     STATEV(2) = e22
C     STATEV(3) = e33
C     STATEV(4) = e12
C     STATEV(5) = e13
C     STATEV(6) = e23
C     STATEV(7) = ep11
C     STATEV(8) = ep22
C     STATEV(9) = ep33
C     STATEV(10) = ep12
C     STATEV(11) = ep13
C     STATEV(12) = ep23
C     STATEV(13) = EQPS
C     STATEV(14) = EQPS Rate
C     STATEV(15) = Johnson-Cook Damage
C
C
C -------------------------------------------------------------------
C Checks that 3D Elements (Hexehedral, Tet...) are being used.
      IF (NDI.NE.3) THEN
         WRITE(6,1)
1        FORMAT(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',
     1     'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
      ENDIF
C--------------------------------------------------------------------
C
C                  ELASTIC PROPERTIES & Temperature Data
C
C--------------------------------------------------------------------
      EMOD=PROPS(1) !Young's Modulus E
      ENU=PROPS(2) !Poisson's Ratio
      A=PROPS(4) !For Johnson-Cook Model
      B=PROPS(5)
      n=PROPS(6)
      IF(ENU.GT.0.4999.AND.ENU.LT.0.5001) ENU=0.499
      D = STATEV(3+2*NTENS) 
      !EMOD = EMOD*(1-D) !Damaged elastic modulus for material degradation.
      EBULK3=EMOD/(ONE-TWO*ENU) !3*K, Bulk Modulus * 3          
      EG2=EMOD/(ONE+ENU) !Lame Constant mu*2.0              
      EG=EG2/TWO !Lame Constant mu
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE !Lame Constant lambda
      rate0 = PROPS(9)
      Tr = PROPS(10) 
      Tm = PROPS(11) 
      D1 = PROPS(12)
      D2 = PROPS(13)
      D3 = PROPS(14)
      D4 = PROPS(15)
      D5 = PROPS(16)

C--------------------------------------------------------------------
C
C                  ELASTIC STIFFNESS
C             (Jacobian = Stiffness Matrix for Elastic case)
C--------------------------------------------------------------------
      DO 20 K1=1,NTENS
        DO 10 K2=1,NTENS
          DDSDDE(K2,K1)=0.0 !Constructs Empty Jacobian Matrix (9x9 zero matrix)
10      CONTINUE
20    CONTINUE
C
      DO 40 K1=1,NDI
        DO 30 K2=1,NDI
          DDSDDE(K2,K1)=ELAM !Jacobian Normal Direction Triangle Terms: lambda
30      CONTINUE
        DDSDDE(K1,K1)=EG2+ELAM !Jacobian Normal Direction Diagonal Terms: (2mew + lambda) 
40    CONTINUE
C
      DO 50 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG !Jacobian Shear Direction Diagonal Terms: mew
50    CONTINUE
C--------------------------------------------------------------------
C
C                  RECOVER ELASTIC AND PLASTIC STRAINS
C
C--------------------------------------------------------------------
!This is important because we need to take rigid body rotation into account.
      CALL ROTSIG(STATEV(1),DROT, EELAS, 2, NDI, NSHR) !Recovers elastic strain and rotates it to original coordinate frame
      CALL ROTSIG(STATEV(NTENS+1),DROT, EPLAS, 2, NDI, NSHR) !Recovers plastic strain and rotates to original coordinates.
      EQPLAS=STATEV(1+2*NTENS) !Recovers equivalent plastic strain.
      EQRATE=STATEV(2+2*NTENS) !Recovers equivalent plastic strain rate from last increment.
C
      J = DFGRD1(1,1)*(DFGRD1(2,2)*DFGRD1(3,3)-DFGRD1(3,2)*DFGRD1(2,3))
     +   -DFGRD1(2,1)*(DFGRD1(1,2)*DFGRD1(3,3)-DFGRD1(3,2)*DFGRD1(1,3))
     +   +DFGRD1(3,1)*(DFGRD1(1,2)*DFGRD1(2,3)-DFGRD1(2,2)*DFGRD1(1,3)) !Determinant of Deformation Gradient Matrix
      J0 = DFGRD0(1,1)*(DFGRD0(2,2)*DFGRD0(3,3)-DFGRD0(3,2)*DFGRD0(2,3))
     +   -DFGRD0(2,1)*(DFGRD0(1,2)*DFGRD0(3,3)-DFGRD0(3,2)*DFGRD0(1,3))
     +   +DFGRD0(3,1)*(DFGRD0(1,2)*DFGRD0(2,3)-DFGRD0(2,2)*DFGRD0(1,3)) !Determinant of Previous Deformation Gradient Matrix     
C--------------------------------------------------------------------
C
C                  CALCULATE STRESS FROM ELASTIC STRAINS
C
C--------------------------------------------------------------------
C Stress(1:3) are Normal Stresses, STRESS(4:6) are shear stresses
      DO 70 K1=1,NTENS
        DO 60 K2=1,NTENS
          STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*DSTRAN(K1) !Calculates Trial Stress based on Hooke's Law 
60      CONTINUE
        EELAS(K1)=EELAS(K1)+DSTRAN(K1) !Updates strain, assuming totally elastic increment (this will be corrected at end)
70    CONTINUE
      sm = sum(STRESS(1:NDI))/3 !Hydrostatic Mean Stress (same as SHYDRO)

C--------------------------------------------------------------------
C
C                  CALCULATE VON MISES STRESS
C                  (Effective Trial Stress)
C--------------------------------------------------------------------
      CALL SINV(STRESS,SHYDRO,SMISES,NDI,NSHR) !Abaqus function that calculates the hydrostatic stress and J2 stress.
      !Effective Trial Stress
      !SMISES =  sqrt(1/2*[(S_xx-S_yy)^2+(S_yy-S_zz)^2+(S_zz-S_xx)^2]+3*(S_xy^2+S_yz^2+S_zx^2))
      !SHYDRO = 1/3*(S_xx+S_yy+S_zz)
C--------------------------------------------------------------------
C
C                  CALL USERHARD SUBROUTINE, GET HARDENING RATE AND YIELD STRESS, DETERMINE IF ACTIVELY YIELDING
C
C--------------------------------------------------------------------
C
      CALL USERHARD(SYIEL0,HARD,EQPLAS,PROPS(4),PROPS(5),PROPS(6),
     +      PROPS(7),PROPS(8),rate0,Tr,Tm,TEMP,DEQPL,DTIME,EQRATE)
C- USERHARD(SYIELD0,HARD,EQPLAS,A,B,EN,C,EM,rate0,Tr,Tm,T,deqpl,dt,eqrate)
C This gives the "flow stress" to evaluate against yield.
C Check Yield Criterion f = 0
      IF (SMISES.GT.(1+TOLER)*SYIEL0) THEN !Checks if material is actively yielding by comparing yield elastic predictor to yield.
C--------------------------------------------------------------------
C
C                  MATERIAL RESPONSE IS PLASTIC, DETERMINE FLOW DIRECTION
C
C--------------------------------------------------------------------
      ONESY=ONE/SMISES !1/S_vm
      DO 110 K1=1,NDI !Normal Stresses STRESS(1:3)
         FLOW(K1)=ONESY*(STRESS(K1)-SHYDRO) !Deviatoric Stress (Normal)/Von Mises Stress 
110   CONTINUE
      DO 120 K1=NDI+1,NTENS !Shear Stresses STRESS(4:6)
         FLOW(K1)=STRESS(K1)*ONESY !Deviatoric Stress (Shear)/Von Mises Stress
120   CONTINUE
C--------------------------------------------------------------------
C
C                  READ PARAMETERS OF JOHNSON-COOK MODEL
C
C--------------------------------------------------------------------
      A=PROPS(4)
      B=PROPS(5)
      EN=PROPS(6)
      C=PROPS(7)
      EM=PROPS(8)
C--------------------------------------------------------------------
C
C                  NEWTON ITERATION
C
C--------------------------------------------------------------------
C Single Equation Return Mapping, This equation only affects the deviatoric stress components.
! Elastic Predictor Return-Mapping Algorithm (Closest Point Projection)
      SYIELD=SYIEL0 !Sets yield to that of the elastic predictor
      DEQPL= 0.0
      DO 130 KEWTON=1, NEWTON
        RHS=SMISES-EG3*DEQPL-SYIELD !Effective Plastic Strain Increment (Residual Plastic Corrector)
        DEQPL=DEQPL+RHS/(EG3+HARD) !Update plastic strain increment with residual
        CALL USERHARD(SYIELD,HARD,EQPLAS+DEQPL,PROPS(4), 
     +                PROPS(5),PROPS(6),PROPS(7),PROPS(8),PROPS(9),
     +                PROPS(10),PROPS(11),TEMP,DEQPL,DTIME,EQRATE) !Update yield and hardening function for new deqpl.
C- USERHARD(SYIELD,HARD,EQPLAS,A,B,EN,C,EM,rate0,Tr,Tm,T,deqpl,dt,rate_old)
      IF(ABS(RHS).LT.TOLER*SYIEL0) GOTO 140 !Compares residual to initial difference of elastic trial stress.
130   CONTINUE
C
C WRITE WARNING MESSAGE TO THE .MSG FILE
C
      WRITE(6,2) NEWTON
2     FORMAT(//,30X,'***WARNING - PLASTICITY ALGORITHM DID NOT ',
     1     'CONVERGE AFTER ',I3,' ITERATIONS')
      PNEWDT = 0.25 !If it didn't converge, then stops step and cuts Time Period in Half
140    CONTINUE
      IF (KEWTON.LT.5) THEN 
            PNEWDT = 1.5 !If it converged in the first 3 iterations, then step will finish and and grow time by 1.5x to save computational time.                 
      END IF
C--------------------------------------------------------------------
C
C                  CALCULATE STRESS AND UPDATE STRAINS
C
C--------------------------------------------------------------------
      DO 150 K1=1,NDI
        STRESS(K1)=FLOW(K1)*SYIELD+SHYDRO !(S^pr/S_vm)*S_y + S_hydro Corrected Stress
        EPLAS(K1)=EPLAS(K1)+THREE*FLOW(K1)*DEQPL/TWO !dep = 3/2*DEQPL*S^pr/S_vm
        EELAS(K1)=EELAS(K1)-THREE*FLOW(K1)*DEQPL/TWO !de - dep
150   CONTINUE
      DO 160 K1=NDI+1,NTENS
        STRESS(K1)=FLOW(K1)*SYIELD !Corrected Normal Stresses.
        EPLAS(K1)=EPLAS(K1)+THREE*FLOW(K1)*DEQPL !dep = 3/2*DEQPL*S^pr/S_vm
        EELAS(K1)=EELAS(K1)-THREE*FLOW(K1)*DEQPL !de - dep
160   CONTINUE
      EQPLAS=EQPLAS+DEQPL !Updates EQPS for end of step.
      SPD=DEQPL*(SYIEL0+SYIELD)/TWO !Plastic Dissipation Calculation.
      RPL = PROPS(3)*SPD/DTIME !Volumetric Heat Generation
C--------------------------------------------------------------------
C                  Determine Damage
C--------------------------------------------------------------------
      CALL USERFAIL(sm,SMISES,EQPLAS,D1,D2,D3,D4,D5,rate0,Tr,Tm,
     +                    TEMP,DEQPL,DTIME,D)
C--------------------------------------------------------------------
C
C                  JACOBIAN
C          (Using Plastic Corrector)
C--------------------------------------------------------------------
      EFFG=EG*SYIELD/SMISES
      EFFG2=TWO*EFFG
      EFFG3=THREE*EFFG
      EFFLAM=(EBULK3-EFFG2)/THREE
      EFFHRD=EG3*HARD/(EG3+HARD)-EFFG3
      DO 220 K1=1,NDI
        DO 210 K2=1,NDI
          DDSDDE(K2,K1)=EFFLAM
210     CONTINUE
        DDSDDE(K1,K1)=EFFG2+EFFLAM
220   CONTINUE
      DO 230 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EFFG
230   CONTINUE
      DO 250 K1=1,NTENS
        DO 240 K2=1,NTENS
           DDSDDE(K2,K1)=DDSDDE(K2,K1)+FLOW(K2)*FLOW(K1)*EFFHRD
240     CONTINUE
250   CONTINUE
C
      ENDIF !Yield Criterian If Statement
C--------------------------------------------------------------------
C
C                  STORE STRAINS IN STATE VARIABLE ARRAY
C
C--------------------------------------------------------------------
      DO 310 K1=1,NTENS
        STATEV(K1)=EELAS(K1) !Elastic Strains in STATEV(1:6)
        STATEV(K1+NTENS)=EPLAS(K1) !Plastic Strains in STATEV(7:12)
310   CONTINUE
      STATEV(1+2*NTENS)=EQPLAS !EQPS STATEV(13)
      STATEV(2+2*NTENS)= DEQPL/DTIME !EQRATE STATEV(14)
      STATEV(3+2*NTENS)=D !Damage STATEV(15)
C
      RETURN !Call to return to Abaqus Solver
      END !End of main program



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C             THE FOLLOWING IS SUBROUTINE USERHARD                            C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE USERHARD(SYIELD,HARD,EQPLAS,A,B,EN,C,EM,rate0,Tr,Tm,
     +                    T,DEQPL,DTIME,EQRATE)
C
      INCLUDE 'aba_param.inc' !Defines all variables as implicit real*8 (double precision)
C
C Johnson-Cook Hardening Model S=(A+Be^n)(1+C*ln((de/dt)/(de/dt|0)))(1-T^m)
C--------------------------------------------------------------------
C 
C                  GET PARAMETERS, SET HARDENING TO ZERO
C
C--------------------------------------------------------------------
      HARD=0.0 !Sets initial hardening modulus to zero.
      if (isnan(EQPLAS)) then
            EQPLAS = 0.0 !If EQPLAS is Undefined, It is set to zero
      end if
      if (isnan(DEQPL)) then
            DEQPL = 0.0 !If DEQPL is undefined, it is set to zero.
      end if
      rate = rate0 !DEQPL/DTIME !Estimation of (dep/dt)

C--------------------------------------------------------------------
C
C                  CALCULATE CURRENT YIELD STRESS AND HARDENING RATE
C
C--------------------------------------------------------------------

      if (isnan(T)) then
            T = Tr !If Temperature is Undefined, Temperature Effects are Neglected
      end if

      if (T < Tr) then
            theta = 0.d0 !No Temp effect
      else if (T > Tm) then
            theta = 1.d0 !Max Temp Effect (No Stress)
      else 
            theta = (T - Tr)/(Tm - Tr) !Temp Function
      end if 
	
      if (rate <= 0.0) then
            rate = rate0 !If negative rate (non physical)
      else if (isnan(rate)) then
            rate = rate0 !If no EQPS yet.
      end if

      if (isnan(EQRATE)) then
            EQRATE = rate
      end if
      
      IF(EQPLAS.LE.0.0) THEN
        SYIELD=A*(1-theta**EM) !No rate effect because no dep yet.
        HARD = 0.0 !No Hardening slope because it just began deforming.
      else if (rate.LE.0.0) then
        SYIELD=A*(1-theta**EM) !No rate effect because no dep yet
        HARD = EN*B*EQPLAS**(EN-1.0)*(1-theta**EM) !No Hardening slope because it just began deforming.
      ELSE
        SYIELD=(A+B*EQPLAS**EN)*(1-theta**EM)*(1+C*log(rate/rate0)) !For Johnson-Cook Model S=(A+Be^n)(1+C*ln((de/dt)/(de/dt|0)))(1-T^m)
          TERMC = 1+C*log(rate/rate0) !Rate Multiplicative Term
          TERMT = 1-theta**EM !Temperature Multiplicative Term
          SYIELD0 = A+B*EQPLAS**EN !Yield with no rate or temp effects.
          HARD0 = EN*B*EQPLAS**(EN-1.0) !Mechanical Hardening Multiplicative Term.
        HARD = TERMT*(HARD0*TERMC +
     +    (A+B*EQPLAS**EN)*(C*(1/rate)*((ABS(EQRATE-rate))/DTIME)/(rate))) !Hardening Slope
      END IF
!      EQRATE = rate
C dS/de ~ {(n*B*e^(n-1))*(1+C*ln((de/dt)/(de/dt|0)))} *(1-T^m)
      RETURN !Return to main program
      END !End of USERHard subroutine



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C             THE FOLLOWING IS SUBROUTINE USERFAIL                            C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE USERFAIL(sm,SMISES,EQPLAS,D1,D2,D3,D4,D5,rate0,Tr,Tm,
     +                    T,DEQPL,DTIME,D)
C
      INCLUDE 'aba_param.inc' !Sets all variables to type implicit real*8 (double precision)

C Johnson-Cook Failure Model w=S(de/ef) 
C where ef = [D1+D2*exp(D3*S_hydro/S_vm)]*[1+D4*ln(((de/dt)/(de/dt|0))][1+D5*T]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



C      write(*,*) 'T,Tr,Tm=',T,Tr,Tm

      if (isnan(T)) then
            T = Tr !If Temperature is Undefined, Temperature Effects are Neglected
      end if
      if (isnan(EQPLAS)) then
            EQPLAS = 0.0 !If EQPLAS is Undefined, It is set to zero
      end if
      if (isnan(DEQPL)) then
            DEQPL = 0.0 !If DEQPL is undefined, it is set to zero.
      end if

      if (T < Tr) then
            theta = 0.d0
      else if (T > Tm) then
            theta = 1.d0
      else 
            theta = (T - Tr)/(Tm - Tr)
      end if 

      efail =(D1+D2*exp(D3*sm/SMISES))*(1+D4*log(DEQPL/DTIME/rate0))
     +                       *(1.d0+D5*theta) 
      D = (EQPLAS+DEQPL)/efail !Calculated Damage Value

      RETURN !Return to main program
      END !End of subroutine USERFAIL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C             END                                                             C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
