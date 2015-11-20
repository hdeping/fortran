MODULE EVOLUTION_MODULE
USE PARAMETER_MODULE




CONTAINS

!--************************************************--!

SUBROUTINE EVOLUTION()
INTEGER, DIMENSION(1)				::				IT
INTEGER								::				I, J, II, JJ
REAL*8, DIMENSION(2)				::				RX, IRX, VXPI, VXPJ, VXVI, VXVJ, VX
REAL*8								::				AIJ, BIJ, CIJ, DIJ, TIJ
REAL*8								::				AS, PACKPHI						

DO 
	
	MOTION_STEP=MOTION_STEP+1

	IF (MOD(MOTION_STEP, 1000)==0) THEN
		AS=0.0D0
		DO I=1, PN
			AS=AS+PI*RA(I)**2.
		ENDDO
		PACKPHI=AS/DBLE(L*L)
		PRINT*, PACKPHI, TIME
	!	IF (PACKPHI>=0.80) THEN
	!		OPEN(1, FILE='INI.TXT')
	!		WRITE(1, *) PACKPHI
	!		DO I=1, PN
	!			WRITE(1, '(I5, 3F25.16)') I, RA(I), PP(I, :)
	!		ENDDO
	!		CLOSE(1)
	!		STOP
	!	ENDIF
	ENDIF



	IT(:)=MINLOC(COTIME)
	TIJ=COTIME(IT(1))
	I=IT(1)
	J=PARTNER(I)



	DO K=1, PN
		COTIME(K)=COTIME(K)-TIJ
		PP(K, :)=PP(K, :)+PV(K, :)*TIJ
		WHERE(PP(K, :)>L)
			PP(K, :)=PP(K, :)-L
		ELSEWHERE(PP(K, :)<=0.0D0)
			PP(K, :)=PP(K, :)+L
		ENDWHERE
		RA(K)=RA(K)+RIN(K)*TIJ	
	ENDDO
	
	TIME=TIME+TIJ
	IF (TIME>=2.0D0) THEN
		AS=0.0D0
		DO I=1, PN
			AS=AS+PI*RA(I)**2.
		ENDDO
		PACKPHI=AS/DBLE(L*L)
		OPEN(1, FILE='INI.TXT')
		WRITE(1, *) PACKPHI
		DO I=1, PN
			WRITE(1, '(I5, 3F25.16)') I, RA(I), PP(I, :)
		ENDDO
		CLOSE(1)
		STOP
	ENDIF

	II=I
	JJ=J
	
	RX(:)=PP(II, :)-PP(JJ, :)
	WHERE(ABS(RX)>L/2.0D0)
		RX=PP(II, :)-(PP(JJ, :)-L*DSIGN(1.0D0, PP(JJ, :)-L/2.0D0))
	ENDWHERE
	IRX(:)=RX(:)/DSQRT(RX(1)**2.+RX(2)**2.)	
	TVI=PV(II, 1)*IRX(1)+PV(II, 2)*IRX(2)
	TVJ=PV(JJ, 1)*IRX(1)+PV(JJ, 2)*IRX(2)
	VXPI(:)=TVI*IRX(:); VXPJ(:)=TVJ*IRX(:)
	VXVI(:)=PV(II, :)-VXPI(:); VXVJ(:)=PV(JJ, :)-VXPJ(:)
	
	PV(II, :)=(VXPJ(:)+(RIN(II)+RIN(JJ))*IRX(:))+VXVI(:)
	PV(JJ, :)=(VXPI(:)-(RIN(II)+RIN(JJ))*IRX(:))+VXVJ(:)
	
	DO I=1, PN
		IF ((I==II).OR.(I==JJ).OR.(PARTNER(I)==II).OR.(PARTNER(I)==JJ)) THEN
			COTIME(I)=BIGMAX
			DO J=1, PN
				IF (J/=I) THEN
					RX(:)=PP(I, :)-PP(J, :)
					WHERE(ABS(RX)>L/2.0D0)
						RX=PP(I, :)-(PP(J, :)-L*DSIGN(1.0D0, PP(J, :)-L/2.0D0))
					ENDWHERE
					VX(:)=PV(I, :)-PV(J, :)
					
					AIJ=(VX(1)**2.+VX(2)**2.)-(RIN(I)+RIN(J))**2.
					BIJ=(RX(1)*VX(1)+RX(2)*VX(2))-(RA(I)+RA(J))*(RIN(I)+RIN(J))
					CIJ=(RX(1)**2.+RX(2)**2.)-(RA(I)+RA(J))**2.
					DIJ=BIJ**2.-AIJ*CIJ					
					
					IF ((BIJ<=0.0D0.OR.AIJ<0.0D0).AND.DIJ>=0.0D0) THEN
						TIJ=(-BIJ-SQRT(DIJ))/AIJ
						IF (TIJ<COTIME(I)) THEN
							COTIME(I)=TIJ
							PARTNER(I)=J
						ENDIF
						IF (TIJ<COTIME(J)) THEN
							COTIME(J)=TIJ
							PARTNER(J)=I
						ENDIF
					ENDIF

				ENDIF
			ENDDO
		ENDIF
	ENDDO

ENDDO

!	CALL DIS_IJ()

END SUBROUTINE EVOLUTION

!--************************************************--!

SUBROUTINE DIS_IJ()
REAL*8					::				DISF, DIS
REAL*8, DIMENSION(2)	::				LX
INTEGER					::				I, J

	DISF=100000.0D0
	DO I=1, PN-1
		DO J=I+1, PN
			LX(:)=ABS(PP(I, :)-PP(J, :))
			WHERE(LX>L/2.0D0)
				LX=L-LX
			ENDWHERE
			DIS=DSQRT(LX(1)**2.+LX(2)**2.)
			IF (DIS<DISF) THEN
				DISF=DIS
			ENDIF
		ENDDO
	ENDDO

	PRINT*, DISF, 2*RA(1)

END SUBROUTINE DIS_IJ

!--************************************************--!

END MODULE EVOLUTION_MODULE