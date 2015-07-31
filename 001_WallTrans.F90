!=============================================================================C
!     Calculate PHIX,PHIY,PHIZ,PHIT                                           C
!=============================================================================C
!
      SUBROUTINE PHI

      INCLUDE 'HEAD.FI'

	COMPLEX CETA(0:NDX2-1,0:NDZ2-1),CETA0(0:NDX2-1,0:NDZ2-1) !the same as temp
	COMPLEX ETAX(0:NDX2-1,0:NDZ2-1),ETA0X(0:NDX2-1,0:NDZ2-1)
	COMPLEX ETAZ(0:NDX2-1,0:NDZ2-1),ETA0Z(0:NDX2-1,0:NDZ2-1)
	
!-------------------------------------------------------------------------------------
      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   ETA =0.5*(ETAU(I,K)-ETAD(I,K))
	   ETA0=0.5*(ETAU(I,K)+ETAD(I,K))

	   PHIY(I,K,1)=0.0
	   PHIY(I,K,2)=1.0/(1.0+ETA)-1.0

	   CETA (I,K)=ETA
	   CETA0(I,K)=ETA0
	ENDDO
	ENDDO
	
	CALL FOURSX(CETA ,NDX2,NDZ2)
	CALL FOURSX(CETA0,NDX2,NDZ2)

      CALL FOURSZ(CETA ,NDX2,NDZ2)
      CALL FOURSZ(CETA0,NDX2,NDZ2)

!---- DERIVATIVE ----C

      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   ETAX(I,K)=(0.0,0.0)
	   ETAZ(I,K)=(0.0,0.0)
	   ETA0X(I,K)=(0.0,0.0)
	   ETA0Z(I,K)=(0.0,0.0)
	ENDDO
	ENDDO

      DO K=0,NDZ2-1
	   KK=K
	   IF(K.GT.NDZ2/2) KK=K-NDZ2
         DO I=0,NDX2-1
	      II=I
	      IF(I.GT.NDX2/2) II=I-NDX2
            ETAX(I,K)=(0.0,1.0)*ARF*II*CETA(I,K)
            ETAZ(I,K)=(0.0,1.0)*BAT*KK*CETA(I,K)

            ETA0X(I,K)=(0.0,1.0)*ARF*II*CETA0(I,K)
            ETA0Z(I,K)=(0.0,1.0)*BAT*KK*CETA0(I,K)
         ENDDO
      ENDDO

      CALL FOURPZ(ETAX, NDX2,NDZ2)
      CALL FOURPZ(ETAZ, NDX2,NDZ2)
      CALL FOURPZ(ETA0X,NDX2,NDZ2)
      CALL FOURPZ(ETA0Z,NDX2,NDZ2)

      CALL FOURPX(ETAX, NDX2,NDZ2)
      CALL FOURPX(ETAZ, NDX2,NDZ2)
      CALL FOURPX(ETA0X,NDX2,NDZ2)
      CALL FOURPX(ETA0Z,NDX2,NDZ2)

      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   ETA =0.5*(ETAU(I,K)-ETAD(I,K))

	   PHIX(I,K,1)=-REAL(ETAX (I,K))/(1.0+ETA)
	   PHIX(I,K,2)=-REAL(ETA0X(I,K))/(1.0+ETA)

	   PHIZ(I,K,1)=-REAL(ETAZ (I,K))/(1.0+ETA)
	   PHIZ(I,K,2)=-REAL(ETA0Z(I,K))/(1.0+ETA)
	ENDDO
	ENDDO



      RETURN
	END
!=============================================================================C
!     Calculate d_/dx                                                         C
!=============================================================================C
!
      SUBROUTINE PARTIAL_X(FI,FO)

      INCLUDE 'HEAD.FI'

	COMPLEX FI(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX FO(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DFDY(0:NDXH-1,0:NDY,0:NDZ-1),DF(0:NDY+1)

	COMPLEX TEMP1(0:NDX2-1,0:NDZ2-1)
!-------------------------------------------------------------------------------------
	DO K=0,NDZ-1
      DO I=0,NDXH-1
         DF(NDY)=(0.,0.)
         DF(NDY+1)=(0.,0.)
         DO JJ=NDY,2,-1
            DF(JJ-1)=DF(JJ+1)+2.*JJ*FI(I,JJ,K)
         ENDDO
         DF(0)=0.5*DF(2)+FI(I,1,K)

	   DO J=0,NDY
	      DFDY(I,J,K)=DF(J)
	   ENDDO
	ENDDO
	ENDDO
      CALL CHEBP(DFDY,NDXH,NDY,NDZ)

      DO J=0,NDY

!---- transform from s.s to p.s ----C
      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   TEMP1(I,K)=(0.0,0.0)
	ENDDO
	ENDDO
!---- z-direction
      DO I=0,NDXH-1
	   DO K=0,NDZH-1
	      TEMP1(I,K)=DFDY(I,J,K)
	   ENDDO
	   DO K=NDZH+1,NDZ-1
	      TEMP1(I,NDZH+K)=DFDY(I,J,K)
	   ENDDO
	ENDDO

      CALL FOURPZ(TEMP1,NDXH,NDZ2)
!---- x-direction
      DO K=0,NDZ2-1
         DO I=1,NDXH-1
            TEMP1(NDX2-I,K)=CONJG(TEMP1(I,K))
	   ENDDO
	ENDDO
      CALL FOURPX(TEMP1,NDX2,NDZ2)

      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   TEMP1(I,K)=TEMP1(I,K)*(PHIX(I,K,1)*Y(J)+PHIX(I,K,2))
	ENDDO
	ENDDO

!---- transform from p.s to s.s ----C

      CALL FOURSX(TEMP1,NDX2,NDZ2)
      CALL FOURSZ(TEMP1,NDXH,NDZ2)

!---- spectral truncation ----C

      DO I=0,NDXH-1
         DO K=0,NDZH-1
            FO(I,J,K)=TEMP1(I,K)
         ENDDO
         DO K=NDZH,NDZ-1
            FO(I,J,K)=TEMP1(I,K+NDZH)
         ENDDO
      ENDDO

      ENDDO ! END J LOOP


      CALL CHEBS(FO,NDXH,NDY,NDZ)

      DO K=0,NDZ-1
	DO J=0,NDY
	DO I=0,NDXH-1
	   FO(I,J,K)=FO(I,J,K)+(0.0,1.0)*ARF*I*FI(I,J,K)
	ENDDO
	ENDDO
	ENDDO

      RETURN
	END
!=============================================================================C
!     Calculate d_/dZ                                                         C
!=============================================================================C
!
      SUBROUTINE PARTIAL_Z(FI,FO)
      
      INCLUDE 'HEAD.FI'

	COMPLEX FI(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX FO(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DFDY(0:NDXH-1,0:NDY,0:NDZ-1),DF(0:NDY+1)

	COMPLEX TEMP1(0:NDX2-1,0:NDZ2-1)
!-------------------------------------------------------------------------------------
	
	DO K=0,NDZ-1
      DO I=0,NDXH-1
         DF(NDY)=(0.,0.)
         DF(NDY+1)=(0.,0.)
         DO JJ=NDY,2,-1
            DF(JJ-1)=DF(JJ+1)+2.*JJ*FI(I,JJ,K)
         ENDDO
         DF(0)=0.5*DF(2)+FI(I,1,K)

	   DO J=0,NDY
	      DFDY(I,J,K)=DF(J)
	   ENDDO
	ENDDO
	ENDDO

      CALL CHEBP(DFDY,NDXH,NDY,NDZ)

	DO J=0,NDY

!---- transform from s.s to p.s ----C
      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   TEMP1(I,K)=(0.0,0.0)
	ENDDO
	ENDDO
!---- z-direction
      DO I=0,NDXH-1
	   DO K=0,NDZH-1
	      TEMP1(I,K)=DFDY(I,J,K)
	   ENDDO
	   DO K=NDZH+1,NDZ-1
	      TEMP1(I,NDZH+K)=DFDY(I,J,K)
	   ENDDO
	ENDDO
      CALL FOURPZ(TEMP1,NDXH,NDZ2)

!---- x-direction
      DO K=0,NDZ2-1
         DO I=1,NDXH-1
            TEMP1(NDX2-I,K)=CONJG(TEMP1(I,K))
	   ENDDO
	ENDDO
      CALL FOURPX(TEMP1,NDX2,NDZ2)



      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   TEMP1(I,K)=TEMP1(I,K)*(PHIZ(I,K,1)*Y(J)+PHIZ(I,K,2))
	ENDDO
	ENDDO
!---- transform from p.s to s.s ----C

      CALL FOURSX(TEMP1,NDX2,NDZ2)
      CALL FOURSZ(TEMP1,NDXH,NDZ2)
!---- spectral truncation ----C
!
      DO I=0,NDXH-1
         DO K=0,NDZH-1
            FO(I,J,K)=TEMP1(I,K)
         ENDDO
         DO K=NDZH,NDZ-1
            FO(I,J,K)=TEMP1(I,K+NDZH)
         ENDDO
      ENDDO

      ENDDO ! END J LOOP

      CALL CHEBS(FO,NDXH,NDY,NDZ)

      DO K=0,NDZ-1
	   KK=K
	   IF(K.GT.NDZH) KK=K-NDZ
	   DO J=0,NDY
	   DO I=0,NDXH-1
	      FO(I,J,K)=FO(I,J,K)+(0.0,1.0)*BAT*KK*FI(I,J,K)
	   ENDDO
	   ENDDO
	ENDDO

      RETURN
	END
!=============================================================================C
!     Calculate d_/dy                                                         C
!=============================================================================C

      SUBROUTINE PARTIAL_Y(FI,FO)  !F IN, F OUT
      
      INCLUDE 'HEAD.FI'

	COMPLEX FI(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX FO(0:NDXH-1,0:NDY,0:NDZ-1)
	COMPLEX DFDY(0:NDXH-1,0:NDY,0:NDZ-1),DF(0:NDY+1)

	COMPLEX TEMP1(0:NDX2-1,0:NDZ2-1)
!-------------------------------------------------------------------------------------
      DO K=0,NDZ-1
      DO I=0,NDXH-1
         DF(NDY)=(0.,0.)
         DF(NDY+1)=(0.,0.)
         DO JJ=NDY,2,-1
            DF(JJ-1)=DF(JJ+1)+2.*JJ*FI(I,JJ,K)
         ENDDO
         DF(0)=0.5*DF(2)+FI(I,1,K)

	   DO J=0,NDY
	      DFDY(I,J,K)=DF(J)
	   ENDDO
	ENDDO
	ENDDO

      CALL CHEBP(DFDY,NDXH,NDY,NDZ)

      DO J=0,NDY

!---- transform from s.s to p.s ----C
      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   TEMP1(I,K)=(0.0,0.0)
	ENDDO
	ENDDO
!---- z-direction
      DO I=0,NDXH-1
	   DO K=0,NDZH-1
	      TEMP1(I,K)=DFDY(I,J,K)
	   ENDDO
	   DO K=NDZH+1,NDZ-1
	      TEMP1(I,NDZH+K)=DFDY(I,J,K)
	   ENDDO
	ENDDO
      CALL FOURPZ(TEMP1,NDXH,NDZ2)
!---- x-direction
      DO K=0,NDZ2-1
         DO I=1,NDXH-1
            TEMP1(NDX2-I,K)=CONJG(TEMP1(I,K))
	   ENDDO
	ENDDO
      CALL FOURPX(TEMP1,NDX2,NDZ2)

      DO K=0,NDZ2-1
	DO I=0,NDX2-1
	   TEMP1(I,K)=TEMP1(I,K)*(PHIY(I,K,1)*Y(J)+PHIY(I,K,2)+1.0)
	ENDDO
	ENDDO

!---- transform from p.s to s.s ----C

      CALL FOURSX(TEMP1,NDX2,NDZ2)

      CALL FOURSZ(TEMP1,NDXH,NDZ2)

!---- spectral truncation ----C
!
      DO I=0,NDXH-1
         DO K=0,NDZH-1
            FO(I,J,K)=TEMP1(I,K)
         ENDDO
         DO K=NDZH,NDZ-1
            FO(I,J,K)=TEMP1(I,K+NDZH)
         ENDDO
      ENDDO
!
      ENDDO ! END J LOOP

      CALL CHEBS(FO,NDXH,NDY,NDZ)

      RETURN
	END


