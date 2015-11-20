MODULE INITIAL_MODULE
USE PARAMETER_MODULE




CONTAINS

!--*************************************--!

SUBROUTINE INITIAL_PARAMETER()
INTEGER				::				I
REAL*8				::				R

	RA(:)=0.0D0; RIN(1:PNA)=0.1D0; RIN(PNA+1:1000)=0.070D0
	DO I=1, PN
		DO J=1, 2
			CALL RANDOM_NUMBER(R)
			PP(I, J)=L*R
			CALL GAUSSIAN_NOISE(R)
			PV(I, J)=R
		ENDDO
	ENDDO
	
	CALL INITIAL_TIME_LIST()
	MOTION_STEP=0
	TIME=0.0D0

END SUBROUTINE INITIAL_PARAMETER

!--*************************************--!

SUBROUTINE INITIAL_TIME_LIST()
INTEGER						::				I, J
REAL*8, DIMENSION(2)		::				RX, VX
REAL*8						::				AIJ, BIJ, CIJ, DIJ, TIJ

	COTIME=BIGMAX
	
	DO I=1, PN-1
		DO J=I+1, PN
			
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

		ENDDO
	ENDDO


END SUBROUTINE INITIAL_TIME_LIST


!--***********************************************--!

SUBROUTINE GAUSSIAN_NOISE(GAUSS_RAN)
REAL*8						::				U, G
INTEGER						::				N
REAL*8						::				BO, R
REAL*8						::				GAUSS_RAN

N=12; BO=0.0
G=1.0; U=0.0

	DO IJ=1, N
		CALL RANDOM_NUMBER(R)
		BO=BO+R
	ENDDO	

	GAUSS_RAN=U+G*(BO-REAL(N)/2.0)							  

END SUBROUTINE GAUSSIAN_NOISE

!--***********************************************--!


END MODULE INITIAL_MODULE