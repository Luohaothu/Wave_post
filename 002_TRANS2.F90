!C-----------------------------------------------------------------C
!C      FOURIER TRANSFORM FROM S.S TO P.S in x direction           C
!C-----------------------------------------------------------------C
!C
      SUBROUTINE FOURPX(U,MX,MZ)
      INCLUDE 'dims.h'
      PARAMETER (ND=NDX2)
      PARAMETER (LEN=2)
      PARAMETER (LENH=LEN/2)
      REAL UD(LENH,0:ND,2),WC(LENH,0:ND,2)
      COMMON /FFTX/NFAX,IFAX(20),TRIGX(2,0:NDX)
      COMPLEX U(0:NDX2-1,0:NDZ2-1)

      ISIGN=1
!C--------- X-DIRECTION TRANFORM ------------
      N=MX
      DO JZ=0,MZ-1
         DO K=0,N-1
            WC(1,K,1)=REAL(U(K,JZ))
            WC(1,K,2)=AIMAG(U(K,JZ))
         ENDDO
         CALL FFT2(WC,UD,N,ND,NFAX,IFAX,ISIGN,TRIGX,LENH)
         DO K=0,N-1
            U(K,JZ)=UD(1,K,1)+(0.,1.)*UD(1,K,2)
         ENDDO
      ENDDO
      RETURN
      END
!C-----------------------------------------------------------------C
!C      FOURIER TRANSFORM FROM S.S TO P.S in z direction           C
!C-----------------------------------------------------------------C
!C
      SUBROUTINE FOURPZ(U,MX,MZ)
      INCLUDE 'dims.h'
      PARAMETER (ND=NDZ2)
      PARAMETER (LEN=2)
      PARAMETER (LENH=LEN/2)
      REAL UD(LENH,0:ND,2),WC(LENH,0:ND,2)
      COMMON /FFTZ/NFAZ,IFAZ(20),TRIGZ(2,0:NDZ)
      COMPLEX U(0:NDX2-1,0:NDZ2-1)
   
      ISIGN=1
!C--------- Z-DIRECTION TRANSFORM -------------
!C
      N=MZ
      DO JX=0,MX-1
         DO K=0,N-1
            WC(1,K,1)=REAL(U(JX,K))
            WC(1,K,2)=AIMAG(U(JX,K))
         ENDDO
         CALL FFT2(WC,UD,N,ND,NFAZ,IFAZ,ISIGN,TRIGZ,LENH)
         DO K=0,N-1
            U(JX,K)=UD(1,K,1)+(0.,1.)*UD(1,K,2)
         ENDDO
      ENDDO
      RETURN
      END
!C----------------------------------------------------------------C
!C      FOURIER TRANSFORM FROM P.S TO S.S in x direction          C
!C----------------------------------------------------------------C
!C
      SUBROUTINE FOURSX(FX,MX,MZ)
      INCLUDE 'dims.h'
      PARAMETER (ND=NDX2)
      PARAMETER (LEN=2)
      PARAMETER (LENH=LEN/2)
      REAL UD(LENH,0:ND,2),WC(LENH,0:ND,2)
      COMMON /FFTX/NFAX,IFAX(20),TRIGX(2,0:NDX)
      COMPLEX FX(0:NDX2-1,0:NDZ2-1)

      ISIGN=-1
!C--------- X-DIRECTION TRANFORM ------------
      N=MX
      DO JZ=0,MZ-1
         DO K=0,N-1
            WC(1,K,1)=REAL(FX(K,JZ))
            WC(1,K,2)=AIMAG(FX(K,JZ))
         ENDDO
         CALL FFT2(WC,UD,N,ND,NFAX,IFAX,ISIGN,TRIGX,LENH)
         DO K=0,N-1
            FX(K,JZ)=UD(1,K,1)+(0.,1.)*UD(1,K,2)
         ENDDO
      ENDDO
      RETURN
      END
!C---------------------------------------------------------------C
!C      FOURIER TRANSFORM FROM P.S TO S.S in z direction         C
!C---------------------------------------------------------------C
!C
      SUBROUTINE FOURSZ(FX,MX,MZ)
      INCLUDE 'dims.h'
      PARAMETER (ND=NDZ2)
      PARAMETER (LEN=2)
      PARAMETER (LENH=LEN/2)
      REAL UD(LENH,0:ND,2),WC(LENH,0:ND,2)
      COMMON /FFTZ/NFAZ,IFAZ(20),TRIGZ(2,0:NDZ)
      COMPLEX FX(0:NDX2-1,0:NDZ2-1)

      ISIGN=-1
!C--------- Z-DIRECTION TRANSFORM -------------
!C
      N=MZ
      DO JX=0,MX-1
         DO K=0,N-1
            WC(1,K,1)=REAL(FX(JX,K))
            WC(1,K,2)=AIMAG(FX(JX,K))
         ENDDO
         CALL FFT2(WC,UD,N,ND,NFAZ,IFAZ,ISIGN,TRIGZ,LENH)
         DO K=0,N-1
            FX(JX,K)=UD(1,K,1)+(0.,1.)*UD(1,K,2)
         ENDDO
      ENDDO

      RETURN
      END
!C------------------------------------------------C
!C      CHEBSHEV TRANSFORM FROM S.S TO P.S OF     C
!C      U,V,W,OMX,OMY,OMZ                         C
!C------------------------------------------------C
!C
       SUBROUTINE CHEBP(FX,MX,N,MZ)
       INCLUDE 'dims.h'
       PARAMETER (ND=NDY)
       PARAMETER (LEN=NDX)
       PARAMETER (LT=LEN)
       COMPLEX FX(0:NDXH-1,0:NDY,0:NDZ-1)
       COMMON /FFTY/NFAY,IFAY(20),TRIGY(2,0:ND)              
       REAL U01((LEN+1)/2),U02((LEN+1)/2) &
            ,UN1((LEN+1)/2),UN2((LEN+1)/2) & 
            ,UR(LT,0:ND),UC(LT,0:ND) &
            ,A((LEN+1)/2,0:ND,2),C((LEN+1)/2,0:ND,2)

       PI=4.*ATAN(1.)
       DO IZ=0,MZ-1  

       DO IX=0,MX-1
          IR=IX+1
          II=MX+IR
          DO J=0,N
             UC(IR,J)=REAL(FX(IX,J,IZ))
             UC(II,J)=REAL(-(0.,1.)*FX(IX,J,IZ))
          ENDDO
       ENDDO

       ISIGN=1
       LENH = (LEN+1) / 2
       LP = LENH + 1
       LPP = LP - 1
       NH = N / 2

       DO IJ = 1,LENH
          U01(IJ) = (UC(IJ,0) + UC(IJ,N)) + UC(IJ,N-1)
          UN1(IJ) = (UC(IJ,0) + UC(IJ,N)) - UC(IJ,N-1)          
          U02(IJ) = (UC(LPP+IJ,0) + UC(LPP+IJ,N))  &
                  + UC(LPP+IJ,N-1)
          UN2(IJ) = (UC(LPP+IJ,0) + UC(LPP+IJ,N)) &
                  - UC(LPP+IJ,N-1)  
          A(IJ,0,1) = 2. * UC(IJ,0)
          A(IJ,0,2) = 2. * UC(LPP+IJ,0)          
          A(IJ,NH,1) = 2. * UC(IJ,N)          
          A(IJ,NH,2) = 2. * UC(LPP+IJ,N)          
       ENDDO
       DO K=1,N-2,2
          DO IJ=1,LENH
             U01(IJ) = U01(IJ) + UC(IJ,K) + UC(IJ,K+1)
             U02(IJ) = U02(IJ) + UC(LPP+IJ,K) + UC(LPP+IJ,K+1)          
             UN1(IJ) = UN1(IJ) - UC(IJ,K) + UC(IJ,K+1)
             UN2(IJ) = UN2(IJ) - UC(LPP+IJ,K) + UC(LPP+IJ,K+1)     
          ENDDO
       ENDDO

       L2N = 0
       DO K = 1,NH-1
          L2N = L2N + 2
          L2NM1 = L2N - 1
          L2NP1 = L2N + 1          
          DO IJ = 1,LENH
             A(IJ,K,1)=UC(IJ,L2N)-UC(LPP+IJ,L2NP1)+UC(LPP+IJ,L2NM1)
             A(IJ,K,2)=UC(LPP+IJ,L2N)+UC(IJ,L2NP1)-UC(IJ,L2NM1)        
          ENDDO
       ENDDO

       L2N = N
       DO K = NH+1,N-1
          L2N = L2N - 2
          L2NM1 = L2N + 1
          L2NP1 = L2N - 1          
          DO IJ = 1,LENH
             A(IJ,K,1)=UC(IJ,L2N)-UC(LPP+IJ,L2NP1)+UC(LPP+IJ,L2NM1)
             A(IJ,K,2)=UC(LPP+IJ,L2N)+UC(IJ,L2NP1)-UC(IJ,L2NM1)      
          ENDDO
       ENDDO
     
       CALL FFT2(A,C,N,ND,NFAY,IFAY,ISIGN,TRIGY,LENH)
       PIN=PI/N
       DO K=1,NH-1
          SSI = 0.125 / SIN(K*PIN)
          ALJ = 0.25 + SSI
          GLJ = 0.25 - SSI          
          DO IJ = 1,LENH
             UR(IJ,K)=ALJ*C(IJ,K,1)+GLJ*C(IJ,N-K,1)
             UR(IJ,N-K)=ALJ*C(IJ,N-K,1)+GLJ*C(IJ,K,1)             
          ENDDO
          DO IJ = 1,LENH
             UR(LPP+IJ,K)=ALJ*C(IJ,K,2)+GLJ*C(IJ,N-K,2)             
             UR(LPP+IJ,N-K)=ALJ*C(IJ,N-K,2)+GLJ*C(IJ,K,2)             
          ENDDO
       ENDDO
       DO IJ = 1,LENH
          UR(IJ,NH)=0.5*C(IJ,NH,1)
          UR(LPP+IJ,NH)=0.5*C(IJ,NH,2)          
       ENDDO
       DO IJ = 1,LENH
          UR(IJ,0) = U01(IJ)
          UR(IJ,N) = UN1(IJ)                    
       ENDDO
       DO IJ = 1,LENH
          UR(LPP+IJ,0) = U02(IJ)          
          UR(LPP+IJ,N) = UN2(IJ)          
       ENDDO
       DO IX=0,MX-1
          IR=IX+1
          II=MX+IR
          DO J=0,N
             FX(IX,J,IZ)= UR(IR,J)+(0.,1.)*UR(II,J)
          ENDDO
       ENDDO

       ENDDO
       RETURN
       END
      

!C-----------------------------------------------C
!C      CHEBSHEV TRANSFORM FROM P.S TO S.S OF    C
!C      FX,FY,FZ                                 C
!C-----------------------------------------------C
!C
       SUBROUTINE CHEBS(FX,MX,N,MZ)
       INCLUDE 'dims.h'
       PARAMETER (ND=NDY)
       PARAMETER (LEN=NDX)
       PARAMETER (LT=LEN)
       COMPLEX FX(0:NDXH-1,0:NDY,0:NDZ-1)
       COMMON /FFTY/NFAY,IFAY(20),TRIGY(2,0:ND)       
       REAL UR(LT,0:ND),UC(LT,0:ND) &
            ,A((LEN+1)/2,0:ND,2),C((LEN+1)/2,0:ND,2) &
            ,C01((LEN+1)/2),C02((LEN+1)/2) &
            ,CN1((LEN+1)/2),CN2((LEN+1)/2) 

       PI=4.*ATAN(1.)
       DO IZ=0,MZ-1  

       DO IX=0,MX-1
          IR=IX+1
          II=MX+IR
          DO J=0,N
             UR(IR,J)=REAL(FX(IX,J,IZ))
             UR(II,J)=REAL(-(0.,1.)*FX(IX,J,IZ))
          ENDDO
       ENDDO
     
       ISIGN=1
       LENH = (LEN+1) / 2
       LP = LENH + 1
       LPP = LP - 1

       XNI = 1./N
       NH = N / 2

       DO IJ = 1,LENH
          C(IJ,0,1) = UR(IJ,0)
          C(IJ,NH,1) = UR(IJ,N)          
          C01(IJ) = (0.5 * (UR(IJ,0) + UR(IJ,N))) + UR(IJ,N-1)
          CN1(IJ) = (0.5 * (UR(IJ,0) + UR(IJ,N))) - UR(IJ,N-1)          
          C(IJ,0,2) = UR(LPP+IJ,0)
          C(IJ,NH,2) = UR(LPP+IJ,N)          
          C02(IJ) = (0.5 * (UR(LPP+IJ,0) + UR(LPP+IJ,N))) &
                  + UR(LPP+IJ,N-1)
          CN2(IJ) = (0.5 * (UR(LPP+IJ,0) + UR(LPP+IJ,N))) &
                  - UR(LPP+IJ,N-1)  
       ENDDO
       L2N = 0
       DO K = 1,NH-1
          L2N = L2N + 2
          L2NM1 = L2N - 1
          L2NP1 = L2N + 1          
          DO IJ = 1,LENH
             C(IJ,K,1)=UR(IJ,L2N)-UR(LPP+IJ,L2NP1)+UR(LPP+IJ,L2NM1)
             C(IJ,K,2)=UR(LPP+IJ,L2N)+UR(IJ,L2NP1)-UR(IJ,L2NM1)       
          ENDDO
       ENDDO

       L2N = N
       DO K = NH+1,N-1
          L2N = L2N - 2
          L2NM1 = L2N + 1
          L2NP1 = L2N - 1          
          DO IJ = 1,LENH
             C(IJ,K,1)=UR(IJ,L2N)-UR(LPP+IJ,L2NP1)+UR(LPP+IJ,L2NM1)
             C(IJ,K,2)=UR(LPP+IJ,L2N)+UR(IJ,L2NP1)-UR(IJ,L2NM1)     
          ENDDO
       ENDDO

       DO K=1,N-2,2       
          DO IJ=1,LENH
             C01(IJ) = C01(IJ) + UR(IJ,K) + UR(IJ,K+1)
             C02(IJ) = C02(IJ) + UR(LPP+IJ,K) + UR(LPP+IJ,K+1)      
             CN1(IJ) = CN1(IJ) - UR(IJ,K) + UR(IJ,K+1)
             CN2(IJ) = CN2(IJ) - UR(LPP+IJ,K) + UR(LPP+IJ,K+1)   
          ENDDO
       ENDDO

       CALL FFT2(C,A,N,ND,NFAY,IFAY,ISIGN,TRIGY,LENH)
       PIN=PI/N
       DO K=1,NH-1
          SSI = 0.25 / SIN(K*PIN)
          ALJ = XNI * (0.5 + SSI)
          GLJ = XNI * (0.5 - SSI)          
          DO IJ = 1,LENH
             UC(IJ,K)=ALJ*A(IJ,K,1)+GLJ*A(IJ,N-K,1)
             UC(IJ,N-K)=ALJ*A(IJ,N-K,1)+GLJ*A(IJ,K,1)             
             UC(LPP+IJ,K)=ALJ*A(IJ,K,2)+GLJ*A(IJ,N-K,2)             
             UC(LPP+IJ,N-K)=ALJ*A(IJ,N-K,2)+GLJ*A(IJ,K,2)             
          ENDDO
       ENDDO
       DO IJ = 1,LENH
          UC(IJ,NH)=XNI*A(IJ,NH,1)
          UC(LPP+IJ,NH)=XNI*A(IJ,NH,2)          
       ENDDO
       DO IJ = 1,LENH
          UC(IJ,0) = XNI * C01(IJ)
          UC(LPP+IJ,0) = XNI * C02(IJ)          
          UC(IJ,N) = XNI * CN1(IJ)          
          UC(LPP+IJ,N) = XNI * CN2(IJ)          
       ENDDO
       DO IX=0,MX-1
          IR=IX+1
          II=MX+IR
          DO J=0,N
             FX(IX,J,IZ)= UC(IR,J)+(0.,1.)*UC(II,J)
          ENDDO
       ENDDO

       ENDDO
       RETURN
       END
!C-----------------------------------C
!C      PREFFT                       C
!C-----------------------------------C
       SUBROUTINE PREFFT(NX,NY,NZ)
       INCLUDE 'dims.h'
       COMMON /FFTX/NFAX,IFAX(20),TRIGX(2,0:NDX2)
       COMMON /FFTY/NFAY,IFAY(20),TRIGY(2,0:NDY)
       COMMON /FFTZ/NFAZ,IFAZ(20),TRIGZ(2,0:NDZ2)
       PI=4.*ATAN(1.)

       CALL FACTOR(NX,NFAX,IFAX)
       CALL FACTOR(NY,NFAY,IFAY)       
       CALL FACTOR(NZ,NFAZ,IFAZ)       
       DO K=0,NX-1
          ARG=2.*PI*K/NX
          TRIGX(1,K)=COS(ARG)
          TRIGX(2,K)=SIN(ARG)
       ENDDO
       DO K=0,NY-1
          ARG=2.*PI*K/NY
          TRIGY(1,K)=COS(ARG)
          TRIGY(2,K)=SIN(ARG)
       ENDDO
       DO K=0,NZ-1
          ARG=2.*PI*K/NZ
          TRIGZ(1,K)=COS(ARG)
          TRIGZ(2,K)=SIN(ARG)
       ENDDO
       RETURN
       END
!C-------------------------------------------------------
       SUBROUTINE FACTOR(N,NFAX,IFAX)
       INTEGER IFAX(20)       
       NFAX=0
       NN=N
!c
!c  extract factors of 3
!c
       DO II=1,20
          IF ( NN .EQ. 3*(NN/3) ) THEN
             NFAX=NFAX+1
             IFAX(NFAX)=3
             NN=NN/3
          ELSE
             GOTO 20
          ENDIF
       ENDDO
20     CONTINUE       
!c
!c  extract factors of 2
!c
       DO II=NFAX+1,20
          IF (NN .EQ. 2*(NN/2)) THEN
             NFAX=NFAX+1
             IFAX(NFAX)=2
             NN=NN/2
          ELSE
             GOTO 40
          ENDIF
       ENDDO
40     CONTINUE       
       IF (NN .NE. 1) THEN
          STOP
       ENDIF
       RETURN
       END
!CC       
!C
!C----------------------------------------------------
       SUBROUTINE FFT2(A,C,N,ND,NFAX,IFAX,ISIGN,TRIG,LEN)
       REAL A(LEN,0:ND,2),C(LEN,0:ND,2)
       REAL TRIG(2,0:ND)
       INTEGER IFAX(20)
       LOGICAL ODD
       
       LA = 1
       ODD = .TRUE.
       DO I = 1,NFAX
          IFAC = IFAX(I)
          IF(ODD) THEN
             CALL PASS2(A,C,N,ND,ISIGN,IFAC,LA,TRIG,LEN)
          ELSE   
             CALL PASS2(C,A,N,ND,ISIGN,IFAC,LA,TRIG,LEN)
          ENDIF   
          ODD = .NOT. ODD
          LA = LA * IFAC
       ENDDO   
       IF(ODD) THEN
          DO I = 0,N-1
             DO IJ = 1,LEN
                C(IJ,I,1) = A(IJ,I,1)
                C(IJ,I,2) = A(IJ,I,2)
             ENDDO
          ENDDO
       ENDIF
       IF ( ISIGN .EQ. -1 ) THEN
          XNI = 1./N
          DO I = 0,N-1
             DO IJ = 1,LEN
                C(IJ,I,1) = XNI * C(IJ,I,1)
                C(IJ,I,2) = XNI * C(IJ,I,2)             
             ENDDO
          ENDDO
       ENDIF
       RETURN
       END
!C----------------------------------------------       
       SUBROUTINE PASS2(A,C,N,ND,ISIGN,IFAC,LA,TRIG,LEN)

       REAL A(LEN,0:ND,2),C(LEN,0:ND,2),TRIG(2,0:ND)
       INTEGER IND(0:20),JND(0:20)
       PI=4.*ATAN(1.)
       ASN60=SIN(PI/3.)
       SN60 = ISIGN * ASN60
       M = N / IFAC

       DO K= 0,IFAC-1
          IND(K) = K * M
          JND(K) = K * LA
       ENDDO
       LLA = LA * LEN

       I = 0
       J = 0
       JUMP = (IFAC-1) * LA
       DO K = 0,M-LA,LA
          IF (IFAC .EQ. 2) THEN
             I0 = IND(0) + I
             I1 = IND(1) + I          
             J0 = JND(0) + J          
             J1 = JND(1) + J          
             CC = TRIG(1,K)
             SS = ISIGN * TRIG(2,K)
             IF (K .EQ. 0) THEN
                DO IJ= 1,LLA
                  C(IJ,J0,1) = A(IJ,I0,1) + A(IJ,I1,1)
                  C(IJ,J0,2) = A(IJ,I0,2) + A(IJ,I1,2)                 
                  C(IJ,J1,1) = A(IJ,I0,1) - A(IJ,I1,1)                 
                  C(IJ,J1,2) = A(IJ,I0,2) - A(IJ,I1,2)                 
                ENDDO   
             ELSE   
                DO IJ = 1,LLA
                  C(IJ,J0,1) = A(IJ,I0,1) + A(IJ,I1,1)
                  C(IJ,J0,2) = A(IJ,I0,2) + A(IJ,I1,2)                 
                  AM1 = A(IJ,I0,1) - A(IJ,I1,1)
                  AM2 = A(IJ,I0,2) - A(IJ,I1,2)                   
                  C(IJ,J1,1) = CC * AM1 - SS * AM2                   
                  C(IJ,J1,2) = SS * AM1 + CC * AM2                   
                ENDDO
             ENDIF
          ELSEIF (IFAC .EQ. 3) THEN
             I0 = IND(0) + I
             I1 = IND(1) + I          
             I2 = IND(2) + I
             J0 = JND(0) + J          
             J1 = JND(1) + J          
             J2 = JND(2) + J
             IF (K .EQ. 0) THEN
                DO IJ = 1,LLA
                   AP1 = A(IJ,I1,1) + A(IJ,I2,1)
                   AP2 = A(IJ,I1,2) + A(IJ,I2,2)                   
                   C(IJ,J0,1) = A(IJ,I0,1) + AP1
                   C(IJ,J0,2) = A(IJ,I0,2) + AP2                   
                   TA1 = A(IJ,I0,1) - 0.5 * AP1
                   TA2 = A(IJ,I0,2) - 0.5 * AP2                   
                   AM1 = SN60 * (A(IJ,I1,1) - A(IJ,I2,1))
                   AM2 = SN60 * (A(IJ,I1,2) - A(IJ,I2,2))              
                   C(IJ,J1,1) = TA1 - AM2                   
                   C(IJ,J1,2) = TA2 + AM1                   
                   C(IJ,J2,1) = TA1 + AM2                              
                   C(IJ,J2,2) = TA2 - AM1                              
                ENDDO   
             ELSE
                C1 = TRIG(1,K)             
                C2 = TRIG(1,2*K)
                S1 = ISIGN * TRIG(2,K)
                S2 = ISIGN * TRIG(2,2*K)                
                DO IJ = 1,LLA
                   AP1 = A(IJ,I1,1) + A(IJ,I2,1)
                   AP2 = A(IJ,I1,2) + A(IJ,I2,2)                   
                   C(IJ,J0,1) = A(IJ,I0,1) + AP1
                   C(IJ,J0,2) = A(IJ,I0,2) + AP2                   
                   TA1 = A(IJ,I0,1) - 0.5 * AP1
                   TA2 = A(IJ,I0,2) - 0.5 * AP2                   
                   AM1 = SN60 * (A(IJ,I1,1) - A(IJ,I2,1))
                   AM2 = SN60 * (A(IJ,I1,2) - A(IJ,I2,2))              
                   T1 = TA1 - AM2
                   T2 = TA2 + AM1
                   C(IJ,J1,1) = C1*T1 - S1*T2                   
                   C(IJ,J1,2) = S1*T1 + C1*T2               
                   T1 = TA1 + AM2
                   T2 = TA2 - AM1
                   C(IJ,J2,1) = C2*T1 - S2*T2
                   C(IJ,J2,2) = S2*T1 + C2*T2                         
                ENDDO
             ENDIF
          ENDIF
          I = I + LA
          J = J + LA
          J = J + JUMP
       ENDDO
       RETURN
       END

