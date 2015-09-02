    !-------------------------------------------------C
    !    Postprocessing Program for Channel Flow      C
    !-------------------------------------------------C
    SUBROUTINE UA_AVERAGE(IT,ITB,ITE,MT)

    INCLUDE 'HEAD.FI'

    integer I,J,K,II,III,I0

    REAL UP(0:NDX2-1,0:NDY,0:NDZ2-1)
    REAL VP(0:NDX2-1,0:NDY,0:NDZ2-1)
    REAL WP(0:NDX2-1,0:NDY,0:NDZ2-1)
    REAL PP(0:NDX2-1,0:NDY,0:NDZ2-1)

    REAL EP

    complex temp1(0:NDX2-1,0:NDZ2-1)
    complex temp2(0:NDX2-1,0:NDZ2-1)
    complex temp3(0:NDX2-1,0:NDZ2-1)
    complex temp4(0:NDX2-1,0:NDZ2-1)

    REAL Y2(0:NDX2/4,0:NDY)
    !===================================================
    IF(IT==ITB) THEN
        UA=0.0
        VA=0.0
        WA=0.0
        PA=0.0
        PW0=0.0
    END IF

    CALL PHI
    !
    CALL CHEBP(U,MX/2,MY,MZ)
    CALL CHEBP(V,MX/2,MY,MZ)
    CALL CHEBP(W,MX/2,MY,MZ)
    CALL CHEBP(P,MX/2,MY,MZ)

    DO J=0,MY
        !---- transform from s.s to p.s ----C
        DO K=0,MZ2-1
            DO I=0,MX2-1
                TEMP1(I,K)=(0.0,0.0)
                TEMP2(I,K)=(0.0,0.0)
                TEMP3(I,K)=(0.0,0.0)
                TEMP4(I,K)=(0.0,0.0)
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

        DO K=0,MZ2-1
            DO I=0,MX2-1
                EP=0.5*(UP(I,J,K)*UP(I,J,K) &
                    +VP(I,J,K)*VP(I,J,K) &
                    +WP(I,J,K)*WP(I,J,K))
                PP(I,J,K)=PP(I,J,K)-EP
            ENDDO
        ENDDO

    ENDDO  ! END J LOOP
    !=============================================
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
                    UA(I,J)=UA(I,J)+UP(III,NDY-J,K)
                    VA(I,J)=VA(I,J)+VP(III,NDY-J,K)
                    WA(I,J)=WA(I,J)+WP(III,NDY-J,K)
                    PA(I,J)=PA(I,J)+PP(III,NDY-J,K)
                END DO
            END DO
        END DO    !END I

    END DO !END II
    !=============================================
    IF(IT.EQ.ITE) THEN

        UA=UA/(4.*MT*NDZ2)
        VA=VA/(4.*MT*NDZ2)
        WA=WA/(4.*MT*NDZ2)
        PA=PA/(4.*MT*NDZ2)

        DO I=0,NDX2/4
            DO J=0,NDY
                ETA=-0.5*A*SIN(KW*X2(I))
                ETA0=0.5*A*SIN(KW*X2(I))
                Y2(I,J)=(1.0+ETA)*Y(J)+ETA0
            END DO
            PW0=PW0+PA(I,0)
        END DO
        PW0=PW0/(1.+NDX2/4.)
        PA=PA-PW0
        !=========================================!
        !         Êä³ö                            !
        !=========================================!
        OPEN(10,FILE='SLICE_UA.DAT')
        WRITE(10,*) 'VARIABLES="X","Y","UA","VA","WA","PA","MUA"'
        WRITE(10,*) 'ZONE  I=',NDX2/4+1,'J=',NDY+1
        DO J=0,NDY
            DO I=0,NDX2/4
                WRITE(10,1) X2(I),Y2(I,NDY-J) &
                    ,UA(I,J),VA(I,J),WA(I,J),PA(I,J) &
                    ,SQRT(UA(I,J)**2+VA(I,J)**2)
            END DO
        END DO
        CLOSE(10)

    END IF

1   FORMAT(1X,10E15.6)

    END