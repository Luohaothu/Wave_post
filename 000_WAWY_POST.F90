!-------------------------------------------------C
!    Postprocessing Program for Channel Flow      C
!-------------------------------------------------C
      INCLUDE 'HEAD.FI'

     
      REAL DP !,DPP(1:30000)        
!----------------------------------------------------------   
      OPEN(20,FILE='XINDAT')
         READ(20,*)
         READ(20,*)RE,ARF,BAT,DT,MX,MY,MZ,FR0,TIME0
         READ(20,*)
         READ(20,*)ITB,ITE,ITSTEP
         READ(20,*)
         READ(20,*)IFC,IFF,IFG
         READ(20,*)
         READ(20,*)DP0,CC0
         READ(20,*)
         READ(20,*)C,A,LMD,KW
	close(20)

      if(mx.ne.ndx.or.my.ne.ndy.or.mz.ne.ndz) then
	   write(*,*)'array dimension is incorrect'
	   stop
	endif

!      write(*,*)'input itb,ite,itstep'
!	read (*,*)itb,ite,itstep
      MT=(ITE-ITB)/ITSTEP+1
     
      PAI=ATAN(1.)*4.
      LMD=LMD*PAI

	MX2=MX*3/2
	MZ2=MZ*3/2
	MXH=MX/2
	MZH=MZ/2

!------ form coordinates ------C
!
      DO J=0,NDY
         Y(J)=COS(J*PAI/MY)
      ENDDO
      XL=2.*PAI/ARF
      DO I=0,NDX
         X(I)=XL*I/NDX
      ENDDO
      DO I=0,NDX2
         X2(I)=XL*I/NDX2
      ENDDO
      ZL=2.*PAI/BAT
      DO J=0,NDZ
         Z(J)=ZL*J/NDZ
      ENDDO
      DO J=0,NDZ2
         Z2(J)=ZL*J/NDZ2
      ENDDO

!------ do preFFT ------C 
!
      CALL PREFFT(NDX2,NDY,NDZ2)


!=====================================================
	DO IT=ITB,ITE,ITSTEP
	write(*,*)it
!--- TIME LOOP
      WRITE(CNAME,'(I5)')IT
      OPEN(15,FILE='../XS'//CNAME//'VP.BIN',FORM='BINARY')
         READ(15)TIME
         READ(15)U
         READ(15)V
         READ(15)W
         READ(15)P
      CLOSE(15)
      OPEN(15,FILE='../XS'//CNAME//'YWALL.BIN',FORM='BINARY')
         READ(15)ETAU
         READ(15)ETAD
      CLOSE(15)
 
	DO K=0,NDZ2-1
	DO I=0,NDX2-1
	DO J=0,MY
	   ETA =0.5*(ETAU(I,K)-ETAD(I,K))
	   ETA0=0.5*(ETAU(I,K)+ETAD(I,K))
	   YP(I,J,K)=(1.0+ETA)*Y(J)+ETA0
	ENDDO
	ENDDO
	ENDDO
	
	CALL UA_AVERAGE(IT,ITB,ITE,MT)
	
	END DO
!===========================================

	DO IT=ITB,ITE,ITSTEP
	write(*,*)it
!--- TIME LOOP
      WRITE(CNAME,'(I5)')IT
      OPEN(15,FILE='../XS'//CNAME//'VP.BIN',FORM='BINARY')
         READ(15)TIME
         READ(15)U
         READ(15)V
         READ(15)W
         READ(15)P
      CLOSE(15)
      OPEN(15,FILE='../XS'//CNAME//'YWALL.BIN',FORM='BINARY')
         READ(15)ETAU
         READ(15)ETAD
      CLOSE(15)
      
	DO K=0,NDZ2-1
	DO I=0,NDX2-1
	DO J=0,MY
	   ETA =0.5*(ETAU(I,K)-ETAD(I,K))
	   ETA0=0.5*(ETAU(I,K)+ETAD(I,K))
	   YP(I,J,K)=(1.0+ETA)*Y(J)+ETA0
	ENDDO
	ENDDO
	ENDDO
	
	DP=0. !DPP(ANINT(1./DT)*(IT-10000)/100)
      CALL DRAG(IT,ITB,ITE,MT,DP)

	END DO
	
!===========================================

	DO IT=ITB,ITE,ITSTEP
	write(*,*)it
!--- TIME LOOP
      WRITE(CNAME,'(I5)')IT
      OPEN(15,FILE='../XS'//CNAME//'VP.BIN',FORM='BINARY')
         READ(15)U
         READ(15)V
         READ(15)W
         READ(15)P
      CLOSE(15)
      OPEN(15,FILE='../XS'//CNAME//'YWALL.BIN',FORM='BINARY')
         READ(15)ETAU
         READ(15)ETAD
      CLOSE(15)
!      ETAU=(0.0,0.0)
!      ETAD=(0.0,0.0)      
	DO K=0,NDZ2-1
	DO I=0,NDX2-1
	DO J=0,MY
	   ETA =0.5*(ETAU(I,K)-ETAD(I,K))
	   ETA0=0.5*(ETAU(I,K)+ETAD(I,K))
	   YP(I,J,K)=(1.0+ETA)*Y(J)+ETA0
	ENDDO
	ENDDO
	ENDDO
	
	CALL RE_AVERAGE(IT,ITB,ITE,MT)
	
	END DO
!=================================================

	DO IT=ITB,ITE,ITSTEP
	write(*,*)it
!--- TIME LOOP
      WRITE(CNAME,'(I5)')IT
      OPEN(15,FILE='../XS'//CNAME//'VP.BIN',FORM='BINARY')
         READ(15)U
         READ(15)V
         READ(15)W
         READ(15)P
      CLOSE(15)
      OPEN(15,FILE='../XS'//CNAME//'YWALL.BIN',FORM='BINARY')
         READ(15)ETAU
         READ(15)ETAD
      CLOSE(15)
!      ETAU=(0.0,0.0)
!      ETAD=(0.0,0.0)      
	DO K=0,NDZ2-1
	DO I=0,NDX2-1
	DO J=0,MY
	   ETA =0.5*(ETAU(I,K)-ETAD(I,K))
	   ETA0=0.5*(ETAU(I,K)+ETAD(I,K))
	   YP(I,J,K)=(1.0+ETA)*Y(J)+ETA0
	ENDDO
	ENDDO
	ENDDO
	
	CALL VOR_AVERAGE(IT,ITB,ITE,MT)
	
	END DO
!===============================================================	
	END