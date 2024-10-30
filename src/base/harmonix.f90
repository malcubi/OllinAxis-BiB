

! Addapted from the cactus module

!#include "cctk.h"

! MODULE: HARMONIX
! AUTHOR: JOSE MANUEL TORRES
!functions intended for development of spin-weighted spherical harmonics stuff
! VERSION: Seems like final
! DATE: Jan 29 2008
! Already tested to calculate coefficients of the expansion of a signal in terms
! of spherical harmonics within an error of 1% order 


MODULE HARMONIX

  IMPLICIT NONE


CONTAINS

  !I NEED A FACTORIAL
  FUNCTION FACT(N)
    IMPLICIT NONE
    integer, INTENT(IN):: N
    real*8:: FACT

    integer I

    FACT=1
    DO I=2,N
       FACT=FACT*I
    END DO
    RETURN
  END FUNCTION FACT


  !
  ! NOTE ON LEGENDRE ASSOCIATED FUNCTIONS
  !
  ! IS CONVENIENT TO TREAT THEM AS
  !
  !   m     m    m
  !  P = sin t* P*
  !   l          l
  !
  ! BY USING RECURRENCES FOR M FIXED
  ! THEN THE SAME APPLIES FOR THE NEW P*
  !


  !FROM NUMERICAL RECIPIES (with some additions)
  FUNCTION APLGNDR(L,M,X)
    IMPLICIT NONE
    integer, INTENT(IN):: L,M
    real*8, INTENT(IN):: X
    real*8:: APLGNDR
    !
    ! Computes the (almost)legendre associated function
    ! 
    !      m
    !     P* (x)
    !      l
    !
    ! Using the sign convention of Abramowitz
    !
    ! RECIVES
    !
    ! l,m integers POSITIVE BOTH
    ! x   real in the interval [-1,1]

    integer:: I,LL,MFLAG,MAUX
    real*8:: PLL,PMM,PMMP1,SOMX2,FAC

    ! I'm dropping the arguments validation
    !
    ! REMEMBER  M POSITIVE

    IF(ABS(M)>L) THEN
       APLGNDR= 0.D0
       RETURN
    END IF

    MFLAG=0
    MAUX=M
    !    IF(M<0) THEN
    !       MAUX=-M
    !       MFLAG= 1
    !    END IF

    PMM=1.D0
    IF(MAUX.GT.0) THEN
       !       SOMX2=SQRT((1.D0-X)*(1.D0+X))

       FAC=1.D0 !APPLY CLOSED FORMULA FOR Pmm
       DO I=1,MAUX
          !          PMM=-PMM*FAC*SOMX2
          PMM= -PMM*FAC
          FAC=FAC+2.D0
       END DO
    ENDIF
    IF(L==MAUX) THEN
       APLGNDR=PMM
    ELSE
       PMMP1=X*(2.D0*MAUX+1.D0)*PMM
       IF(L<=MAUX+1) THEN
          APLGNDR=PMMP1
       ELSE
          DO LL=MAUX+2,L
             PLL=(X*(2.D0*LL-1.D0)*PMMP1-(LL+MAUX-1.D0)*PMM)/(LL-MAUX)
             PMM=PMMP1
             PMMP1=PLL
          ENDDO
          APLGNDR=PLL
       ENDIF
    ENDIF
    !    !HERE WE CHECK IF M WAS NEGATIVE
    !    IF(MFLAG.EQ.1) PLGNDR = (-1)**MAUX*FACT(L-MAUX)/FACT(L+MAUX)*PLGNDR
    RETURN
  END FUNCTION APLGNDR

  ! FOR COMPLETNESS I RECONSTRUCT THE LEGENDRE ASSOCIATED FUNCTIONS 
  FUNCTION PLGNDR(L,M,X)
    IMPLICIT NONE
    integer, INTENT(IN):: L,M
    real*8, INTENT(IN):: X
    real*8 :: PLGNDR

    !BY THE USE OF THE P* FUNCTION 
    !THE ASSOCIATED LEGENDRE FUNCTIONS ARE GIVEN BY
    !
    !   m         m/2   m
    !  P = (1-X*X)   * P*  m>l
    !   l               l

    !   -m       m   (l-m)!    m
    !  P   = (-1)  *------- * P
    !   l            (l+m)!    l

    !    integer :: I
    !    real*8  :: PLA

    IF (ABS(M)>L) THEN
       PLGNDR=0.D0
       RETURN
    END IF

    PLGNDR= APLGNDR(L,ABS(M),X)

    IF(M>0) PLGNDR= SQRT((1.D0-X)*(1.D0+X))**M * PLGNDR
    IF(M<0) PLGNDR= (-1)**(-M)*FACT(L+M)/FACT(L-M) *SQRT((1.D0-X)*(1.D0+X))**(-M) * PLGNDR

    RETURN
  END FUNCTION PLGNDR


  FUNCTION GAUX(L,X)
    IMPLICIT NONE
    integer, INTENT(IN) :: L
    real*8, INTENT(IN)  :: X
    real*8 :: GAUX

    ! L MUST BE POSITIVE integer

    integer :: I

    GAUX=1.D0

    DO I=2,L
       GAUX=PLGNDR(I-1,0,X)+DBLE(I-1)/DBLE(I)*X*GAUX
    END DO

    RETURN
  END FUNCTION GAUX

  FUNCTION HAUX(L,X)
    IMPLICIT NONE
    integer, INTENT(IN):: L
    real*8,  INTENT(IN):: X
    real*8 :: HAUX

    ! L INTEGER EQUAL OR GREATER THAN 2

    ! HERE, INSTEAD OF CALLING GAUX AT EACH ITERATION
    ! IT IS CONSTRUCTED ON THE GO

    integer :: I
    real*8 :: GAU

    HAUX= 0.D0
    GAU = 1.D0

    DO I=2,L
       HAUX= DBLE(I-2)/DBLE(I-1)*HAUX-DBLE(2*(2*I-1))*GAU
       GAU=PLGNDR(I-1,0,X)+DBLE(I-1)/DBLE(I)*X*GAU
    END DO

    RETURN
  END FUNCTION HAUX

  ! NOW I IMPLEMENT THE SPIN-WEIGHTED SPHERICAL HARMONICS
  !
  ! FIRST, THE NORMALIZATION CONSTANT SO WE CAN TAKE IT
  ! OUT FROM LINEAR OPERATORS (I.E. TO AVOID CALCULATING A
  ! SQUARE ROOT AND 2 DIVISIONS FOR EACH VALUE WHEN INTEGRATING
  ! NUMERICALLY

  FUNCTION ASWSH(S,L,M)
    IMPLICIT NONE
    integer, INTENT(IN) :: S,L,M

    real*8 :: ASWSH

    ! JUST WORKS FOR S= -2,-1,0,1,2

    integer :: S1,M1,MF
    real*8, PARAMETER :: PI=4.D0*ATAN(1.D0)

    IF(L<ABS(S)) THEN
       ASWSH=0.D0
       RETURN
    END IF

    !CHECK FOR NEGARIVE M
    MF=0 !FLAG
    IF(M<0)THEN
       !SET NEGATIVE PARMS, TURN FLAG ON
       MF=1
       M1=-M
       S1=-S
    ELSE
       M1=M
       S1=S
    END IF

    SELECT CASE (S1)
    CASE (-2,2)
       IF (M1==0) THEN
          ASWSH= - SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-2))/DBLE(FACT(L+2)))*DBLE((L-1))
       ELSE  IF (M1==1) THEN
          ASWSH= SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-2))/DBLE(FACT(L+2))/DBLE(L*(L+1)))*DBLE((L-1))
       ELSE
          ASWSH= SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-2))/DBLE(FACT(L+2))*DBLE(FACT(L-M1))/DBLE(FACT(L+M1)))
       END IF
    CASE (-1)
       IF (M1==0) THEN
          ASWSH= -SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-1))/DBLE(FACT(L+1)))*DBLE(L)
       ELSE
          ASWSH= SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-1))/DBLE(FACT(L+1))*DBLE(FACT(L-M1))/DBLE(FACT(L+M1)))
       END IF
    CASE (0)
       ASWSH=SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-M1))/DBLE(FACT(L+M1)))
    CASE (1)
       IF (M1==0) THEN
          ASWSH= SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-1))/DBLE(FACT(L+1)))*DBLE(L)
       ELSE
          ASWSH= -SQRT(DBLE(2*L+1)/4.D0/PI*DBLE(FACT(L-1))/DBLE(FACT(L+1))*DBLE(FACT(L-M1))/DBLE(FACT(L+M1)))
       END IF
    CASE DEFAULT !YOU SHOULD NOT REACH HERE IF YOU KNOW WHAT YOU'RE DOING

    END SELECT

    !CHECK THE FLAG, PUT THE PHASE

    IF(MF==1) ASWSH= (-1)**(ABS(M1+S))*ASWSH

    RETURN
  END FUNCTION ASWSH

  !AND THE COMPLEX FACTOR THAT CARRIES THE DEPENDENCE ON THE ANGLE
  FUNCTION FSWSHR(S,L,M,T,P)
    IMPLICIT NONE
    integer, INTENT(IN) :: S,L,M
    real*8, INTENT(IN) :: T,P
    real*8 :: FSWSHR

    ! JUST WORKS FOR S= -2,-1,0,1,2

    integer :: S1,M1,MF
    real*8, PARAMETER :: PI=4.D0*ATAN(1.D0)
    real*8 :: X

    IF(L<ABS(S)) THEN
       FSWSHR=0.D0
       RETURN
    END IF

    !CHECK FOR NEGARIVE M
    MF=0
    IF(M<0)THEN
       !SET NEGATIVE PARMS, TURN FLAG ON
       MF=1
       M1=-M
       S1=-S
    ELSE
       M1=M
       S1=S
    END IF

    X=COS(T)

    SELECT CASE (S1)
    CASE (-2)
       IF (M1==0) THEN
          FSWSHR= DBLE(L)*PLGNDR(L,0,X)-2.D0*GAUX(L-1,X)
       ELSE  IF (M1==1) THEN
          FSWSHR= (-HAUX(L,X)+DBLE(L**2)*GAUX(L,X))*SIN(T)* COS(P)
       ELSE
          FSWSHR= ( (2.D0*M1**2-DBLE(L*(L+1)) + DBLE(2*M1*(L-1))*X + DBLE(L*(L-1))*X**2)*APLGNDR(L,M1,X) &
               + DBLE(2*(L+M1))*(X - DBLE(M1))*APLGNDR(L-1,M1,X))*SIN(T)**(M1-2)*COS(DBLE(M1)*P)
       END IF
    CASE (-1)
       IF (M1==0) THEN
          FSWSHR= GAUX(L,X)*SIN(T)
       ELSE
          FSWSHR= ( (DBLE(L)*X + DBLE(M1))*APLGNDR(L,M1,X)-DBLE(L+M1)*APLGNDR(L-1,M1,X)) &
               *SIN(T)**(M1-1)*COS(DBLE(M1)*P)
       END IF
    CASE (0)
       FSWSHR= APLGNDR(L,M1,X)*SIN(T)**(M1)*COS(DBLE(M1)*P)
    CASE (1)
       IF (M1==0) THEN
          FSWSHR= GAUX(L,X)*SIN(T)
       ELSE
          FSWSHR= ( (DBLE(L)*X - DBLE(M1))*APLGNDR(L,M1,X)-DBLE(L+M1)*APLGNDR(L-1,M1,X)) &
               *SIN(T)**(M1-1)*COS(DBLE(M1)*P)
       END IF
    CASE (2)
       IF (M1==0) THEN
          FSWSHR= DBLE(L)*PLGNDR(L,0,X)-2.D0*GAUX(L-1,X)
       ELSE  IF (M1==1) THEN
          FSWSHR= (HAUX(L,X)+DBLE(L**2)*GAUX(L,X))*SIN(T)* COS(P)
       ELSE
          FSWSHR= ( (2.D0*M1**2-DBLE(L*(L+1)) - DBLE(2*M1*(L-1))*X + DBLE(L*(L-1))*X**2)*APLGNDR(L,M1,X) &
               + DBLE(2*(L+M1))*(X + DBLE(M1))*APLGNDR(L-1,M1,X))*SIN(T)**(M1-2)*COS(DBLE(M1)*P)
       END IF

    CASE DEFAULT !YOU SHOULD NOT REACH HERE IF YOU KNOW WHAT YOU'RE DOING

    END SELECT


    RETURN
  END FUNCTION FSWSHR

  FUNCTION FSWSHI(S,L,M,T,P)
    IMPLICIT NONE
    integer, INTENT(IN) :: S,L,M
    real*8, INTENT(IN) :: T,P
    real*8 :: FSWSHI

    ! JUST WORKS FOR S= -2,-1,0,1,2

    integer :: S1,M1,MF
    real*8, PARAMETER :: PI=4.D0*ATAN(1.D0)
    real*8 :: X

    IF(L<ABS(S)) THEN
       FSWSHI=0.D0
       RETURN
    END IF

    !CHECK FOR NEGARIVE M
    MF=0
    IF(M<0)THEN
       !SET NEGATIVE PARMS, TURN FLAG ON
       MF=1
       M1=-M
       S1=-S
    ELSE
       M1=M
       S1=S
    END IF

    X=COS(T)

    SELECT CASE (S1)
    CASE (-2)
       IF (M1==0) THEN
          FSWSHI= 0.D0
       ELSE  IF (M1==1) THEN
          FSWSHI= (-HAUX(L,X)+DBLE(L**2)*GAUX(L,X))*SIN(T)* SIN(P)
       ELSE
          FSWSHI= ( (2.D0*M1**2-DBLE(L*(L+1)) + DBLE(2*M1*(L-1))*X + DBLE(L*(L-1))*X**2)*APLGNDR(L,M1,X) &
               + DBLE(2*(L+M1))*(X - DBLE(M1))*APLGNDR(L-1,M1,X))*SIN(T)**(M1-2)*SIN(DBLE(M1)*P)
       END IF
    CASE (-1)
       IF (M1==0) THEN
          FSWSHI= 0.D0
       ELSE
          FSWSHI= ( (DBLE(L)*X + DBLE(M1))*APLGNDR(L,M1,X)-DBLE(L+M1)*APLGNDR(L-1,M1,X)) &
               *SIN(T)**(M1-1)*SIN(DBLE(M1)*P)
       END IF
    CASE (0)
       FSWSHI= APLGNDR(L,M1,X)*SIN(T)**(M1)*SIN(DBLE(M1)*P)
    CASE (1)
       IF (M1==0) THEN
          FSWSHI= 0.D0
       ELSE
          FSWSHI= ( (DBLE(L)*X - DBLE(M1))*APLGNDR(L,M1,X)-DBLE(L+M1)*APLGNDR(L-1,M1,X)) &
               *SIN(T)**(M1-1)*SIN(DBLE(M1)*P)
       END IF
    CASE (2)
       IF (M1==0) THEN
          FSWSHI= 0.D0
       ELSE  IF (M1==1) THEN
          FSWSHI= (HAUX(L,X)+DBLE(L**2)*GAUX(L,X))*SIN(T)* SIN(P)
       ELSE
          FSWSHI= ( (2.D0*M1**2-DBLE(L*(L+1)) - DBLE(2*M1*(L-1))*X + DBLE(L*(L-1))*X**2)*APLGNDR(L,M1,X) &
               + DBLE(2*(L+M1))*(X + DBLE(M1))*APLGNDR(L-1,M1,X))*SIN(T)**(M1-2)*SIN(DBLE(M1)*P)
       END IF

    CASE DEFAULT !YOU SHOULD NOT REACH HERE IF YOU KNOW WHAT YOU'RE DOING

    END SELECT

    IF(MF==1) FSWSHI=-FSWSHI

    RETURN
  END FUNCTION FSWSHI


  FUNCTION SWSHR(S,L,M,T,P)
    IMPLICIT NONE
    integer, INTENT(IN) :: S,L,M
    real*8,  INTENT(IN) :: T,P
    real*8 :: SWSHR

    SWSHR= ASWSH(S,L,M)*FSWSHR(S,L,M,T,P)

    RETURN

  END FUNCTION SWSHR

  FUNCTION SWSHI(S,L,M,T,P)
    IMPLICIT NONE
    integer, INTENT(IN) :: S,L,M
    real*8,  INTENT(IN) :: T,P
    real*8 :: SWSHI

    SWSHI= ASWSH(S,L,M)*FSWSHI(S,L,M,T,P)

    RETURN

  END FUNCTION SWSHI


END MODULE HARMONIX




