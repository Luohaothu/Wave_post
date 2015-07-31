!-------------------------------------------------C
!    Postprocessing Program for Channel Flow      C
!-------------------------------------------------C
      SUBROUTINE DRAG(IT,ITB,ITE,MT,DP)
      INCLUDE 'HEAD.FI'
      
!-------------------------------------------------------

    integer I,J,K,II,III,I0

    COMPLEX DUDX(0:NDXH-1,0:NDY,0:NDZ-1)
    COMPLEX DUDY(0:NDXH-1,0:NDY,0:NDZ-1)
    COMPLEX DUDZ(0:NDXH-1,0:NDY,0:NDZ-1)
    COMPLEX DVDX(0:NDXH-1,0:NDY,0:NDZ-1)
    COMPLEX DWDX(0:NDXH-1,0:NDY,0:NDZ-1)   
    
	COMPLEX ETADX(0:NDX2-1,0:NDZ2-1)
	COMPLEX ETADZ(0:NDX2-1,0:NDZ2-1)
	
    complex temp1(0:NDX2-1,0:NDZ2-1)
    complex temp2(0:NDX2-1,0:NDZ2-1)
    complex temp3(0:NDX2-1,0:NDZ2-1)
    complex temp4(0:NDX2-1,0:NDZ2-1)
    complex temp5(0:NDX2-1,0:NDZ2-1)
    complex temp6(0:NDX2-1,0:NDZ2-1)
    complex temp7(0:NDX2-1,0:NDZ2-1)
    complex temp8(0:NDX2-1,0:NDZ2-1)
    complex temp9(0:NDX2-1,0:NDZ2-1)
    
    REAL DP,EP,F1,F2,F3,F4,FN,P0
    REAL Y2(0:NDX2/4,0:NDY)
    
!===================================================  	
      CALL PHI

    IF(IT==ITB) THEN
	  DPA=0.0
	  DRPA=0.0
	  DRFA=0.0
    END IF
    
	DO I=0,NDX2-1
	DO K=0,NDZ2-1
	TEMP2(I,K)=ETAD(I,K)
	ENDDO
	ENDDO

      CALL FOURSX(TEMP2,NDX2,NDZ2)
      CALL FOURSZ(TEMP2,NDX2,NDZ2)

    DO K=0,NDZ2-1
    KK=K
    IF(K.GE.NDZ2/2) KK=K-NDZ2
    DO I=0,NDX2-1
    II=I
    IF(I.GE.NDX2/2) II=I-NDX2
        ETADX(I,K)=(0.0,1.0)*ARF*II*TEMP2(I,K)
        ETADZ(I,K)=(0.0,1.0)*BAT*KK*TEMP2(I,K)
    ENDDO
    ENDDO

      CALL FOURPZ(ETADX, NDX2,NDZ2)
      CALL FOURPZ(ETADZ, NDX2,NDZ2)
      
      CALL FOURPX(ETADX, NDX2,NDZ2)
      CALL FOURPX(ETADZ, NDX2,NDZ2)
!------------------------------------------------------

      CALL PARTIAL_X(U,DUDX)
      CALL PARTIAL_Y(U,DUDY)
      CALL PARTIAL_Z(U,DUDZ)
      CALL PARTIAL_X(V,DVDX)
      CALL PARTIAL_X(W,DWDX)

      CALL CHEBP(U,MX/2,MY,MZ)
      CALL CHEBP(V,MX/2,MY,MZ)
      CALL CHEBP(W,MX/2,MY,MZ)
      CALL CHEBP(P,MX/2,MY,MZ)
      CALL CHEBP(DUDX,MX/2,MY,MZ)
      CALL CHEBP(DUDY,MX/2,MY,MZ)
      CALL CHEBP(DUDZ,MX/2,MY,MZ)
      CALL CHEBP(DVDX,MX/2,MY,MZ)
      CALL CHEBP(DWDX,MX/2,MY,MZ)

     J=NDY
!---- transform from s.s to p.s ----C
     DO K=0,MZ2-1
	 DO I=0,MX2-1
	   TEMP1(I,K)=(0.0,0.0)
	   TEMP2(I,K)=(0.0,0.0)
	   TEMP3(I,K)=(0.0,0.0)
	   TEMP4(I,K)=(0.0,0.0)
	   TEMP5(I,K)=(0.0,0.0)
	   TEMP6(I,K)=(0.0,0.0)
	   TEMP7(I,K)=(0.0,0.0)
	   TEMP8(I,K)=(0.0,0.0)
	   TEMP9(I,K)=(0.0,0.0)
	ENDDO
	ENDDO
!---- z-direction
    DO I=0,MXH-1
    DO K=0,MZH-1
        TEMP1(I,K)=U(I,J,K)
        TEMP2(I,K)=V(I,J,K)
        TEMP3(I,K)=W(I,J,K)
        TEMP4(I,K)=P(I,J,K)

        TEMP5(I,K)=DUDX(I,J,K)
        TEMP6(I,K)=DUDY(I,J,K)
        TEMP7(I,K)=DUDZ(I,J,K) 
        TEMP8(I,K)=DVDX(I,J,K)
        TEMP9(I,K)=DWDX(I,J,K)
    ENDDO
    DO K=MZH+1,MZ-1
        TEMP1(I,MZH+K)=U(I,J,K)
        TEMP2(I,MZH+K)=V(I,J,K)
        TEMP3(I,MZH+K)=W(I,J,K)
        TEMP4(I,MZH+K)=P(I,J,K)

        TEMP5(I,MZH+K)=DUDX(I,J,K)
        TEMP6(I,MZH+K)=DUDY(I,J,K)
        TEMP7(I,MZH+K)=DUDZ(I,J,K) 
        TEMP8(I,MZH+K)=DVDX(I,J,K)
        TEMP9(I,MZH+K)=DWDX(I,J,K)
    ENDDO
    ENDDO
   
    CALL FOURPZ(TEMP1,MXH,MZ2)
    CALL FOURPZ(TEMP2,MXH,MZ2)
    CALL FOURPZ(TEMP3,MXH,MZ2)
    CALL FOURPZ(TEMP4,MXH,MZ2)
    CALL FOURPZ(TEMP5,MXH,MZ2)
    CALL FOURPZ(TEMP6,MXH,MZ2)
    CALL FOURPZ(TEMP7,MXH,MZ2)
    CALL FOURPZ(TEMP8,MXH,MZ2)
    CALL FOURPZ(TEMP9,MXH,MZ2)
!---- x-direction
      DO K=0,MZ2-1
      DO I=1,MXH-1
        TEMP5(MX2-I,K)=CONJG(TEMP5(I,K))
        TEMP4(MX2-I,K)=CONJG(TEMP4(I,K))
        TEMP1(MX2-I,K)=CONJG(TEMP1(I,K))
        TEMP2(MX2-I,K)=CONJG(TEMP2(I,K))
        TEMP3(MX2-I,K)=CONJG(TEMP3(I,K))
        TEMP9(MX2-I,K)=CONJG(TEMP9(I,K))        
        TEMP6(MX2-I,K)=CONJG(TEMP6(I,K))
        TEMP7(MX2-I,K)=CONJG(TEMP7(I,K))
        TEMP8(MX2-I,K)=CONJG(TEMP8(I,K))
	  ENDDO
	  ENDDO
        CALL FOURPX(TEMP5,MX2,MZ2)
        CALL FOURPX(TEMP4,MX2,MZ2)
        CALL FOURPX(TEMP1,MX2,MZ2)
        CALL FOURPX(TEMP2,MX2,MZ2)
        CALL FOURPX(TEMP3,MX2,MZ2)   
        CALL FOURPX(TEMP6,MX2,MZ2)
        CALL FOURPX(TEMP7,MX2,MZ2)
        CALL FOURPX(TEMP8,MX2,MZ2)
        CALL FOURPX(TEMP9,MX2,MZ2)
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
        EP=0.5*(REAL(TEMP1(III,K))*REAL(TEMP1(III,K)) &
               +REAL(TEMP2(III,K))*REAL(TEMP2(III,K)) &
               +REAL(TEMP3(III,K))*REAL(TEMP3(III,K)))
	    
        F1=-2.*REAL(TEMP5(III,K))*REAL(ETADX(III,K))/RE
        F2=REAL(TEMP8(III,K))/RE+REAL(TEMP6(III,K))/RE  
        F3=-(REAL(TEMP7(III,K))+REAL(TEMP9(III,K))) &
            *REAL(ETADZ(III,K))/RE
        F4=(REAL(TEMP4(III,K))-EP-PW0 & !+DP*X2(I)
          )*REAL(ETADX(III,K))
        
        FN=SQRT(1.+REAL(ETADX(III,K))**2+REAL(ETADZ(III,K))**2)
        
        DPA(I)=DPA(I)+REAL(TEMP4(III,K))-EP-PW0
        DRPA(I)=DRPA(I)+F4/FN
        DRFA(I)=DRFA(I)+(F1+F2+F3)/FN
    END DO
    END DO    !END I

    END DO !END II
!=============================================  
    IF(IT.EQ.ITE) THEN
        DPA=DPA/(4.*MT*NDZ2)
        DRPA=DRPA/(4.*MT*NDZ2)
        DRFA=DRFA/(4.*MT*NDZ2)  
        DRA=DRFA+DRPA

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
    OPEN(10,FILE='XDRAG.DAT')
	WRITE(10,*) 'VARIABLES="X","PA","DRA","DPA","DFA"'
	WRITE(10,*) 'ZONE  I=',NDX2/4+1
	DO I=0,NDX2/4
		WRITE(10,1) X2(I),DPA(I),DRA(I),DRPA(I),DRFA(I)
	END DO
	CLOSE(10) 
	
    OPEN(10,FILE='PROFILE_UA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","U1","Y2","U2","Y3","U3","Y4","U4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) (Y2(0,NDY-J)-Y2(0,NDY))*RE*SQRT(ABS(DRFA(0))) &
	            ,UA(0,J)/SQRT(ABS(DRFA(0))) &
	            ,(Y2(NDX2/16,NDY-J)-Y2(NDX2/16,NDY))*RE*SQRT(ABS(DRFA(NDX2/16))) &
	            ,UA(NDX2/16,J)/SQRT(ABS(DRFA(NDX2/16))) &
	            ,(Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            ,UA(NDX2/8,J)/SQRT(ABS(DRFA(NDX2/8))) &
	            ,(Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            ,UA(3*NDX2/16,J)/SQRT(ABS(DRFA(3*NDX2/16))) 
	END DO
	CLOSE(10) 
	
    OPEN(10,FILE='PROFILE_VA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","V1","Y2","V2","Y3","V3","Y4","V4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) (Y2(0,NDY-J)-Y2(0,NDY))*RE*SQRT(ABS(DRFA(0))) &
	            ,VA(0,J)/SQRT(ABS(DRFA(0))) &
	            ,(Y2(NDX2/16,NDY-J)-Y2(NDX2/16,NDY))*RE*SQRT(ABS(DRFA(NDX2/16))) &
	            ,VA(NDX2/16,J)/SQRT(ABS(DRFA(NDX2/16))) &
	            ,(Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            ,VA(NDX2/8,J)/SQRT(ABS(DRFA(NDX2/8))) &
	            ,(Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            ,VA(3*NDX2/16,J)/SQRT(ABS(DRFA(3*NDX2/16))) 
	END DO
	CLOSE(10) 
	
    OPEN(10,FILE='PROFILE_WA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","W1","Y2","W2","Y3","W3","Y4","W4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) (Y2(0,NDY-J)-Y2(0,NDY))*RE*SQRT(ABS(DRFA(0))) &
	            ,WA(0,J)/SQRT(ABS(DRFA(0))) &
	            ,(Y2(NDX2/16,NDY-J)-Y2(NDX2/16,NDY))*RE*SQRT(ABS(DRFA(NDX2/16))) &
	            ,WA(NDX2/16,J)/SQRT(ABS(DRFA(NDX2/16))) &
	            ,(Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            ,WA(NDX2/8,J)/SQRT(ABS(DRFA(NDX2/8))) &
	            ,(Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            ,WA(3*NDX2/16,J)/SQRT(ABS(DRFA(3*NDX2/16))) 
	END DO
	CLOSE(10) 
	
    OPEN(10,FILE='PROFILE_PA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","P1","Y2","P2","Y3","P3","Y4","P4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) (Y2(0,NDY-J)-Y2(0,NDY))*RE*SQRT(ABS(DRFA(0))) &
	            ,PA(0,J) & !/ABS(DRFA(0)) &
	            ,(Y2(NDX2/16,NDY-J)-Y2(NDX2/16,NDY))*RE*SQRT(ABS(DRFA(NDX2/16))) &
	            ,PA(NDX2/16,J) & !/ABS(DRFA(NDX2/16)) &
	            ,(Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            ,PA(NDX2/8,J) & !/ABS(DRFA(NDX2/8)) &
	            ,(Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            ,PA(3*NDX2/16,J)  !/ABS(DRFA(3*NDX2/16))
	END DO
	CLOSE(10) 
	
    END IF	

1   FORMAT(1X,10E15.6)

	END