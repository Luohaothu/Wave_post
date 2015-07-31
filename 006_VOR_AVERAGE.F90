!-------------------------------------------------C
!    Postprocessing Program for Channel Flow      C
!-------------------------------------------------C
      SUBROUTINE VOR_AVERAGE(IT,ITB,ITE,MT)
      
      INCLUDE 'HEAD.FI'

      REAL OMGXA,OMGYA,OMGZA,QA  
      COMMON /SDA/OMGXA(0:NDX2/4,0:NDY/2),OMGYA(0:NDX2/4,0:NDY/2) &
                  ,OMGZA(0:NDX2/4,0:NDY/2),QA(0:NDX2/4,0:NDY/2) 
                  
    integer I,J,K,II,III,I0

	COMPLEX DUDX(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DUDY(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DUDZ(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DVDX(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DVDY(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DVDZ(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DWDX(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DWDY(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DWDZ(0:NDXH-1,0:NDY,0:NDZ-1)

	REAL OMGX,OMGY,OMGZ,S11,S22,S33,S12,S13,S23,Q
	
    complex temp5(0:NDX2-1,0:NDZ2-1)
    complex temp6(0:NDX2-1,0:NDZ2-1)
    complex temp7(0:NDX2-1,0:NDZ2-1)
    complex temp8(0:NDX2-1,0:NDZ2-1)
    complex temp9(0:NDX2-1,0:NDZ2-1)
    complex temp10(0:NDX2-1,0:NDZ2-1)
    complex temp11(0:NDX2-1,0:NDZ2-1)
    complex temp12(0:NDX2-1,0:NDZ2-1)
    complex temp13(0:NDX2-1,0:NDZ2-1)
    
    REAL Y2(0:NDX2/4,0:NDY)
!===================================================  	
    IF(IT==ITB) THEN
	  OMGXA=0.0
	  OMGYA=0.0
	  OMGZA=0.0
	  QA=0.0
    END IF

      CALL PHI

      CALL PARTIAL_X(U,DUDX)
      CALL PARTIAL_Y(U,DUDY)
      CALL PARTIAL_Z(U,DUDZ)

      CALL PARTIAL_X(V,DVDX)
      CALL PARTIAL_Y(V,DVDY)
      CALL PARTIAL_Z(V,DVDZ)

      CALL PARTIAL_X(W,DWDX)
      CALL PARTIAL_Y(W,DWDY)
      CALL PARTIAL_Z(W,DWDZ)

      CALL CHEBP(DUDX,MX/2,MY,MZ)
      CALL CHEBP(DUDY,MX/2,MY,MZ)
      CALL CHEBP(DUDZ,MX/2,MY,MZ)

      CALL CHEBP(DVDX,MX/2,MY,MZ)
      CALL CHEBP(DVDY,MX/2,MY,MZ)
      CALL CHEBP(DVDZ,MX/2,MY,MZ)

      CALL CHEBP(DWDX,MX/2,MY,MZ)
      CALL CHEBP(DWDY,MX/2,MY,MZ)
      CALL CHEBP(DWDZ,MX/2,MY,MZ)

     DO J=0,MY
!---- transform from s.s to p.s ----C
     DO K=0,MZ2-1
	 DO I=0,MX2-1
	   TEMP5(I,K)=(0.0,0.0)
	   TEMP6(I,K)=(0.0,0.0)
	   TEMP7(I,K)=(0.0,0.0)
	   TEMP8(I,K)=(0.0,0.0)
	   TEMP9(I,K)=(0.0,0.0)
	   TEMP10(I,K)=(0.0,0.0)
	   TEMP11(I,K)=(0.0,0.0)
	   TEMP12(I,K)=(0.0,0.0)
	   TEMP13(I,K)=(0.0,0.0)
	ENDDO
	ENDDO
!---- z-direction
      DO I=0,MXH-1
	   DO K=0,MZH-1
	      TEMP5(I,K)=DUDX(I,J,K)
	      TEMP6(I,K)=DUDY(I,J,K)
	      TEMP7(I,K)=DUDZ(I,J,K)

	      TEMP8(I,K)=DVDX(I,J,K)
	      TEMP9(I,K)=DVDY(I,J,K)
	      TEMP10(I,K)=DVDZ(I,J,K)

	      TEMP11(I,K)=DWDX(I,J,K)
	      TEMP12(I,K)=DWDY(I,J,K)
	      TEMP13(I,K)=DWDZ(I,J,K)
	   ENDDO
	   DO K=MZH+1,MZ-1
          TEMP5(I,MZH+K)=DUDX(I,J,K)
	      TEMP6(I,MZH+K)=DUDY(I,J,K)
	      TEMP7(I,MZH+K)=DUDZ(I,J,K)

	      TEMP8(I,MZH+K)=DVDX(I,J,K)
	      TEMP9(I,MZH+K)=DVDY(I,J,K)
	      TEMP10(I,MZH+K)=DVDZ(I,J,K)

	      TEMP11(I,MZH+K)=DWDX(I,J,K)
	      TEMP12(I,MZH+K)=DWDY(I,J,K)
	      TEMP13(I,MZH+K)=DWDZ(I,J,K)
	   ENDDO
	ENDDO
      CALL FOURPZ(TEMP5,MXH,MZ2)
      CALL FOURPZ(TEMP6,MXH,MZ2)
      CALL FOURPZ(TEMP7,MXH,MZ2)
      CALL FOURPZ(TEMP8,MXH,MZ2)
      CALL FOURPZ(TEMP9,MXH,MZ2)
      CALL FOURPZ(TEMP10,MXH,MZ2)
      CALL FOURPZ(TEMP11,MXH,MZ2)
      CALL FOURPZ(TEMP12,MXH,MZ2)
      CALL FOURPZ(TEMP13,MXH,MZ2)
!---- x-direction
      DO K=0,MZ2-1
      DO I=1,MXH-1
        TEMP5(MX2-I,K)=CONJG(TEMP5(I,K))
        TEMP6(MX2-I,K)=CONJG(TEMP6(I,K))
        TEMP7(MX2-I,K)=CONJG(TEMP7(I,K))
        TEMP8(MX2-I,K)=CONJG(TEMP8(I,K))
        TEMP9(MX2-I,K)=CONJG(TEMP9(I,K))
        TEMP10(MX2-I,K)=CONJG(TEMP10(I,K))
        TEMP11(MX2-I,K)=CONJG(TEMP11(I,K))
        TEMP12(MX2-I,K)=CONJG(TEMP12(I,K))
        TEMP13(MX2-I,K)=CONJG(TEMP13(I,K))
	  ENDDO
	  ENDDO
      CALL FOURPX(TEMP5,MX2,MZ2)
      CALL FOURPX(TEMP6,MX2,MZ2)
      CALL FOURPX(TEMP7,MX2,MZ2)
      CALL FOURPX(TEMP8,MX2,MZ2)
      CALL FOURPX(TEMP9,MX2,MZ2)
      CALL FOURPX(TEMP10,MX2,MZ2)
      CALL FOURPX(TEMP11,MX2,MZ2)
      CALL FOURPX(TEMP12,MX2,MZ2)
      CALL FOURPX(TEMP13,MX2,MZ2)   
!---- END TRANSFORM FROM SS TO PS
    DO I=0,NDX2-2
        IF((ETAD(I,0).LT.0.) &
      .AND.(ETAD(I+1,0).GT.0.))THEN
            I0=I
            EXIT
        END IF
    END DO
    
    DO II=0,3
    
    DO I=0,NDX2/4
    III=I0+I+II*NDX2/4
    IF(III.GE.NDX2) III=III-NDX2
    IF(III.LT.0)    III=III+NDX2
    DO K=0,NDZ2-1
        OMGX=REAL(TEMP12(III,K)-TEMP10(III,K))
        OMGY=REAL(TEMP7(III,K)-TEMP11(III,K))
        OMGZ=REAL(TEMP8(III,K)-TEMP6(III,K))
        S11=REAL(TEMP5(III,K))
        S22=REAL(TEMP9(III,K))
        S33=REAL(TEMP13(III,K))
        S12=0.5*REAL(TEMP8(III,K)+TEMP5(III,K))
        S13=0.5*REAL(TEMP11(III,K)+TEMP7(III,K))
        S23=0.5*REAL(TEMP12(III,K)+TEMP10(III,K))
        Q=0.25*(OMGX**2+OMGY**2+OMGZ**2) &
               -0.5*(S11**2+S22**2+S33**2 &
               +2.0*S12*S12+2.0*S13*S13+2.0*S23*S23)   
                
        OMGXA(I,NDY-J)=OMGXA(I,NDY-J)+OMGX
        OMGYA(I,NDY-J)=OMGYA(I,NDY-J)+OMGY
        OMGZA(I,NDY-J)=OMGZA(I,NDY-J)+OMGZ
        QA(I,NDY-J)=QA(I,NDY-J)+Q
    END DO
    END DO    !END I

    END DO !END II

    ENDDO  ! END J LOOP
!=============================================  
!=============================================
    IF(IT.EQ.ITE) THEN
        OMGXA=OMGXA/(4.*MT*NDZ2)
        OMGYA=OMGYA/(4.*MT*NDZ2)
        OMGZA=OMGZA/(4.*MT*NDZ2)
        QA=QA/(4.*MT*NDZ2)     

    DO I=0,NDX2/4
    DO J=0,NDY
       ETA=-0.5*A*SIN(KW*X2(I))
	   ETA0=0.5*A*SIN(KW*X2(I))
       Y2(I,J)=(1.0+ETA)*Y(J)+ETA0
    END DO
    END DO
!=========================================!
!         Êä³ö                            !
!=========================================!  
    OPEN(10,FILE='CAL_VORA.DAT')
	WRITE(10,*) 'VARIABLES="X","Y","OMGXA","OMGYA","OMGZA","QA"'
	WRITE(10,*) 'ZONE  I=',NDX2/4+1,'J=',NDY+1
	DO J=0,NDY
	DO I=0,NDX2/4
		WRITE(10,1) X2(I),Y2(I,NDY-J) &
		            ,OMGXA(I,J),OMGYA(I,J),OMGZA(I,J),QA(I,J)
	END DO
	END DO
	CLOSE(10) 

    END IF	

1   FORMAT(1X,10E15.6)

	END