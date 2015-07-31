!-------------------------------------------------C
!    Postprocessing Program for Channel Flow      C
!-------------------------------------------------C
      SUBROUTINE RE_AVERAGE(IT,ITB,ITE,MT)
      
      INCLUDE 'HEAD.FI'
     
    integer I,J,K,II,III,I0
	
	REAL UP(0:NDX2-1,0:NDY,0:NDZ2-1)
	REAL VP(0:NDX2-1,0:NDY,0:NDZ2-1)
	REAL WP(0:NDX2-1,0:NDY,0:NDZ2-1)
	REAL PP(0:NDX2-1,0:NDY,0:NDZ2-1)
	
    complex temp1(0:NDX2-1,0:NDZ2-1)
    complex temp2(0:NDX2-1,0:NDZ2-1)
    complex temp3(0:NDX2-1,0:NDZ2-1)
    complex temp4(0:NDX2-1,0:NDZ2-1)
      
    REAL Y2(0:NDX2/4,0:NDY) , EP
    
!===================================================
    IF(IT==ITB) THEN  
	  RE_UUA=0.0
	  RE_VVA=0.0
	  RE_WWA=0.0
	  RE_UVA=0.0
	  PFA = 0.
    END IF

      CALL PHI
      CALL CHEBP(U,MX/2,MY,MZ)
      CALL CHEBP(V,MX/2,MY,MZ)
      CALL CHEBP(W,MX/2,MY,MZ)
      CALL CHEBP(P,MX/2,MY,MZ)
      
     DO J=0,MY
!
!---- transform from s.s to p.s ----C
     DO K=0,MZ2-1
	 DO I=0,MX2-1
	   TEMP1(I,K)=0.0
	   TEMP2(I,K)=0.0
	   TEMP3(I,K)=0.0
	   TEMP4(I,K)=0.0
	ENDDO
	ENDDO
!---- z-direction
      DO I=0,MXH-1
	   DO K=0,MZH-1
	      TEMP1(I,K)=U(I,J,K)
	      TEMP2(I,K)=V(I,J,K)
	      TEMP3(I,K)=W(I,J,K)
	      TEMP4(I,K)=P(I,J,K)
	   ENDDO
	   DO K=MZH+1,MZ-1
	      TEMP1(I,MZH+K)=U(I,J,K)
	      TEMP2(I,MZH+K)=V(I,J,K)
	      TEMP3(I,MZH+K)=W(I,J,K)
	      TEMP4(I,MZH+K)=P(I,J,K)
	   ENDDO
	ENDDO
      CALL FOURPZ(TEMP1,MXH,MZ2)
      CALL FOURPZ(TEMP2,MXH,MZ2)
      CALL FOURPZ(TEMP3,MXH,MZ2)
      CALL FOURPZ(TEMP4,MXH,MZ2)
!---- x-direction
      DO K=0,MZ2-1
         DO I=1,MXH-1
            TEMP1(MX2-I,K)=CONJG(TEMP1(I,K))
            TEMP2(MX2-I,K)=CONJG(TEMP2(I,K))
            TEMP3(MX2-I,K)=CONJG(TEMP3(I,K))
            TEMP4(MX2-I,K)=CONJG(TEMP4(I,K))
	   ENDDO
	ENDDO
    CALL FOURPX(TEMP1,MX2,MZ2)
    CALL FOURPX(TEMP2,MX2,MZ2)
    CALL FOURPX(TEMP3,MX2,MZ2)
    CALL FOURPX(TEMP4,MX2,MZ2)
!---- END TRANSFORM FROM SS TO PS
	DO K=0,MZ2-1
	   DO I=0,MX2-1
	      UP(I,J,K)=REAL(TEMP1(I,K))
	      VP(I,J,K)=REAL(TEMP2(I,K))
	      WP(I,J,K)=REAL(TEMP3(I,K))
	      PP(I,J,K)=REAL(TEMP4(I,K))
	   ENDDO
	ENDDO
	
    ENDDO  ! END J LOOP
    
    PP = PP - 0.5 * ( UP ** 2. + VP ** 2. + WP ** 2. ) - PW0
    
!=============================================
!    I0=NINT(MOD(C*IT/100.,2.*PAI/0.5)/(2.*PAI/0.5)*NDX2)
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
    DO J=0,NDY
    DO K=0,NDZ2-1
        RE_UUA(I,J)=RE_UUA(I,J)+(UP(III,NDY-J,K)-UA(I,J))**2
        RE_VVA(I,J)=RE_VVA(I,J)+(VP(III,NDY-J,K)-VA(I,J))**2
        RE_WWA(I,J)=RE_WWA(I,J)+(WP(III,NDY-J,K)-WA(I,J))**2
        RE_UVA(I,J)=RE_UVA(I,J)+(UP(III,NDY-J,K)-UA(I,J)) &
                               *(VP(III,NDY-J,K)-VA(I,J))
        PFA( I , J ) = PFA( I , J ) + ( PP( III , NDY-J , K ) - PA( I , J ) ) ** 2.
    END DO
    END DO
    END DO    !END I

    END DO !END II
    
!=============================================

    IF(IT.EQ.ITE) THEN
            RE_UUA=RE_UUA/(4.*MZ2*MT)
            RE_VVA=RE_VVA/(4.*MZ2*MT)
            RE_WWA=RE_WWA/(4.*MZ2*MT)
            RE_UVA=RE_UVA/(4.*MZ2*MT)
            PFA = SQRT( PFA / ( 4. * MZ2 * MT ) )
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

    OPEN(10,FILE='RE_UUA.DAT')
	WRITE(10,*) 'VARIABLES="X","Y","RE_UUA","RE_VVA","RE_WWA","RE_UVA","PFA"'
	WRITE(10,*) 'ZONE  I=',NDX2/4+1,'J=',NDY+1
	DO J=0,NDY
	DO I=0,NDX2/4
		WRITE(10,1) X2(I),Y2(I,NDY-J) &
                    ,RE_UUA(I,J),RE_VVA(I,J) &
			        ,RE_WWA(I,J),RE_UVA(I,J) &
			        , PFA( I , J )
	END DO
	END DO
	CLOSE(10)
	
    OPEN(10,FILE='PROFILE_UUA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","UF1","Y2","UF2","Y3","UF3","Y4","UF4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) ( Y2( 0 , NDY-J ) - Y2( 0 , NDY ) ) * RE * SQRT( ABS( DRFA( 0 ) ) ) &
	            , SQRT( RE_UUA( 0 , J ) ) / SQRT( ABS( DRFA( 0 ) ) ) &
	            ,( Y2( NDX2/16 , NDY-J ) - Y2( NDX2/16 , NDY ) ) * RE * SQRT( ABS( DRFA( NDX2/16 ) ) ) &
	            , SQRT( RE_UUA( NDX2/16 , J ) ) /SQRT(ABS(DRFA(NDX2/16))) &
	            ,( Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            , SQRT( RE_UUA( NDX2/8 , J ) ) /SQRT(ABS(DRFA(NDX2/8))) &
	            ,( Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            , SQRT( RE_UUA( 3 * NDX2/16 , J ) ) /SQRT(ABS(DRFA(3*NDX2/16)))
	END DO
	CLOSE(10)
	            
    OPEN(10,FILE='PROFILE_VVA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","VF1","Y2","VF2","Y3","VF3","Y4","VF4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) ( Y2( 0 , NDY-J ) - Y2( 0 , NDY ) ) * RE * SQRT( ABS( DRFA( 0 ) ) ) &
	            , SQRT( RE_VVA( 0 , J ) ) / SQRT( ABS( DRFA( 0 ) ) ) &
	            ,( Y2( NDX2/16 , NDY-J ) - Y2( NDX2/16 , NDY ) ) * RE * SQRT( ABS( DRFA( NDX2/16 ) ) ) &
	            , SQRT( RE_VVA( NDX2/16 , J ) ) /SQRT(ABS(DRFA(NDX2/16))) &
	            ,( Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            , SQRT( RE_VVA( NDX2/8 , J ) ) /SQRT(ABS(DRFA(NDX2/8))) &
	            ,( Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            , SQRT( RE_VVA( 3 * NDX2/16 , J ) ) /SQRT(ABS(DRFA(3*NDX2/16)))
	END DO
	CLOSE(10)
	            
	
	OPEN(10,FILE='PROFILE_WWA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","WF1","Y2","WF2","Y3","WF3","Y4","WF4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) ( Y2( 0 , NDY-J ) - Y2( 0 , NDY ) ) * RE * SQRT( ABS( DRFA( 0 ) ) ) &
	            , SQRT( RE_WWA( 0 , J ) ) / SQRT( ABS( DRFA( 0 ) ) ) &
	            ,( Y2( NDX2/16 , NDY-J ) - Y2( NDX2/16 , NDY ) ) * RE * SQRT( ABS( DRFA( NDX2/16 ) ) ) &
	            , SQRT( RE_WWA( NDX2/16 , J ) ) /SQRT(ABS(DRFA(NDX2/16))) &
	            ,( Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            , SQRT( RE_WWA( NDX2/8 , J ) ) /SQRT(ABS(DRFA(NDX2/8))) &
	            ,( Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            , SQRT( RE_WWA( 3 * NDX2/16 , J ) ) /SQRT(ABS(DRFA(3*NDX2/16)))
	END DO
	CLOSE(10)
   
	OPEN(10,FILE='PROFILE_UVA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","UV1","Y2","UV2","Y3","UV3","Y4","UV4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) ( Y2( 0 , NDY-J ) - Y2( 0 , NDY ) ) * RE * SQRT( ABS( DRFA( 0 ) ) ) &
	            , - RE_UVA( 0 , J ) / ABS( DRFA( 0 ) ) &
	            ,( Y2( NDX2/16 , NDY-J ) - Y2( NDX2/16 , NDY ) ) * RE * SQRT( ABS( DRFA( NDX2/16 ) ) ) &
	            , - RE_UVA( NDX2/16 , J ) / ABS( DRFA( NDX2/16 ) ) &
	            ,( Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            , - RE_UVA( NDX2/8  , J ) / ABS( DRFA( NDX2/8  ) ) &
	            ,( Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            , - RE_UVA( 3 * NDX2/16 , J ) / ABS( DRFA( 3*NDX2/16 ) ) 
    END DO
	CLOSE(10)
	
    OPEN(10,FILE='PROFILE_PFA.DAT')
	WRITE(10,*) 'VARIABLES="Y1","PF1","Y2","PF2","Y3","PF3","Y4","PF4"'
	WRITE(10,*) 'ZONE  J=',NDY+1
	DO J=0,NDY
	WRITE(10,1) ( Y2( 0 , NDY-J ) - Y2( 0 , NDY ) ) * RE * SQRT( ABS( DRFA( 0 ) ) ) &
	            , PFA( 0 , J )  & !/ ABS( DRFA( 0 ) ) &
	            ,( Y2( NDX2/16 , NDY-J ) - Y2( NDX2/16 , NDY ) ) * RE * SQRT( ABS( DRFA( NDX2/16 ) ) ) &
	            , PFA( NDX2/16 , J ) & ! / ABS( DRFA( NDX2/16 ) ) &
	            ,( Y2(NDX2/8,NDY-J)-Y2(NDX2/8,NDY))*RE*SQRT(ABS(DRFA(NDX2/8))) &
	            , PFA( NDX2/8 , J )  & !/ ABS( DRFA( NDX2/8 ) ) &
	            ,( Y2(3*NDX2/16,NDY-J)-Y2(3*NDX2/16,NDY))*RE*SQRT(ABS(DRFA(3*NDX2/16))) &
	            , PFA( 3*NDX2/16 , J )   !/ ABS( DRFA( 3*NDX2/16 ) ) 
	END DO
	CLOSE(10)
	
    OPEN(10,FILE='XPPLUS.DAT')
	WRITE(10,*) 'VARIABLES="X","PA","PFA"'
	WRITE(10,*) 'ZONE  I=',NDX2/4+1
	DO I=0,NDX2/4
		WRITE(10,1) X2(I),DPA(I) & !/ABS(DRFA(I))
		           ,PFA(I,0)/ABS(DRFA(I))
	END DO
	CLOSE(10) 
	
    END IF	

1   FORMAT(1X,10E15.6)

!==================================================================================

	END