
!-------------
!*** MLS PART
!-------------
REAL(8) FUNCTION dotprod(n,x,y) 
  IMPLICIT REAL(8) (a-h,o-z) 
  REAL(8),DIMENSION(*)::x,y 
  IF(n.GT.0)THEN
     dotprod=0.0d00
     DO in=1,n
        dotprod=dotprod+x(in)*y(in)
     END DO
  ELSE
     dotprod=0.0d00
  END IF
END FUNCTION dotprod

REAL(8) FUNCTION rnorm2(n,x) 
  IMPLICIT REAL(8) (a-h,o-z) 
  REAL(8),DIMENSION(*)::x 
  IF(n.LE.0)THEN 
     rnorm2=0.0d00 
  ELSE
     rnorm2=SQRT(dotprod(n,x,x))
  END IF
END FUNCTION rnorm2

SUBROUTINE nrmali(n,x)  
  IMPLICIT REAL(8) (a-h,o-z) 
  REAL(8),DIMENSION(*)::x 
  r=rnorm2(n,x(1)) 
  IF(r.GT.0.0d00)THEN 
     r=1.0d00/r
     DO in=1,n
        x(in)=r*x(in)
     END DO
  ENDIF
END SUBROUTINE nrmali

!-------------
!*** MLS PART
!-------------  
!-----------------------------
!*** GENERAL TRIANGULAR SOLVE
!*** MULTIPLE RIGHT-HAND-SIDES
!*** CHK0
!-----------------------------
  SUBROUTINE MLS_GENTRIANGSOLVE(UPPER,NRHS,N,R,B,X)
    IMPLICIT REAL(8) (a-h,o-z)
    LOGICAL::UPPER
    INTEGER::N
    REAL(8),DIMENSION(N,N)::R
    REAL(8),DIMENSION(N,NRHS)::B,X
    IF(UPPER)THEN
!-------------------------------------------
!*** SOLVE R.X=B WHERE R IS UPPER TRIANGULAR
!-------------------------------------------
       DO IR=1,NRHS
          DO ID=N,1,-1
             X(ID,IR)=B(ID,IR)
             DO JD=ID+1,N
                X(ID,IR)=X(ID,IR)-R(ID,JD)*X(JD,IR)
             END DO
             X(ID,IR)=X(ID,IR)/R(ID,ID)
          END DO
       END DO
    ELSE
!---------------------------------------------
!*** SOLVE R^T.X=B WHERE R IS UPPER TRIANGULAR
!---------------------------------------------
       DO IR=1,NRHS
          DO ID=1,N
             X(ID,IR)=B(ID,IR)
             DO JD=1,ID-1
                X(ID,IR)=X(ID,IR)-R(JD,ID)*X(JD,IR)
             END DO
             X(ID,IR)=X(ID,IR)/R(ID,ID)
          END DO
       END DO
    END IF
  END SUBROUTINE MLS_GENTRIANGSOLVE

!--------------------------
!*** PERFORMS A NAIVE QR
!*** DECOMPOSITION
!*** CHK0
!--------------------------
  SUBROUTINE MLS_PROJQR(N,U,A,UA)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::N
    REAL(8),DIMENSION(N)::U,A,UA
    UA=U*DOTPROD(N,U,A)/DOTPROD(N,U,U)
  END SUBROUTINE MLS_PROJQR

!--------------------------
!*** PERFORMS A NAIVE QR
!*** DECOMPOSITION
!*** CHK0
!--------------------------
  SUBROUTINE MLS_NAIVEQR(N,M,AT,R)
    IMPLICIT REAL(8) (a-h,o-z)
    INTEGER::N,M
    REAL(8),DIMENSION(N)::TEMP
    REAL(8),DIMENSION(N,M)::E
    REAL(8),DIMENSION(N,M)::A
    REAL(8),DIMENSION(M,M)::R
    REAL(8),DIMENSION(N,M)::AT
!*** SETS R TO ZERO
    R=0.0D00    
    DO IM=1,M
       E(1:N,IM)=AT(1:N,IM)
       DO JM=1,IM-1
          CALL MLS_PROJQR(N,E(1:N,JM),AT(1:N,IM),TEMP(1:N))
          E(1:N,IM)=E(1:N,IM)-TEMP
       END DO
       CALL NRMALI(N,E(1:N,IM))
    END DO
    DO IM=1,M
       A(1:N,IM)=0.0D00
       DO JM=1,IM
          A(1:N,IM)=A(1:N,IM)+DOTPROD(N,E(1:N,JM),AT(1:N,IM))*E(1:N,JM)
       END DO
    END DO
    DO IM=1,M
       DO JM=IM,M
          R(IM,JM)=DOTPROD(N,E(1:N,IM),A(1:N,JM))
       END DO
    END DO
  END SUBROUTINE MLS_NAIVEQR

!-------------------
!*** FROM W, P AND B
!*** DETERMINES U2
!*** CHK0 (BOTH)
!-------------------  
  SUBROUTINE MLS_DETERMU2(M,N,W,P,U2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::M,N
    REAL(8),DIMENSION(N)::W
    REAL(8),DIMENSION(M,N)::B,U1,U2
    REAL(8),DIMENSION(M,N)::P
    REAL(8),DIMENSION(N,M)::AT
    REAL(8),DIMENSION(M,M)::R,ATEMP
    DO IM=1,M
       DO IN=1,N
          AT(IN,IM)=SQRT(W(IN))*P(IM,IN)
       END DO
    END DO
    DO IN=1,N
       DO IM=1,M
          B(IM,IN)=P(IM,IN)*W(IN)
       END DO
    END DO
!    GOTO 133
    CALL MLS_NAIVEQR(N,M,AT,R)
    CALL MLS_GENTRIANGSOLVE(.FALSE.,N,M,R,B,U1)
    CALL MLS_GENTRIANGSOLVE(.TRUE.,N,M,R,U1,U2)
    RETURN
133 CONTINUE
    DO IM=1,M
       DO JM=1,M
          ATEMP(IM,JM)=0.0D00
          DO IN=1,N
             ATEMP(IM,JM)=ATEMP(IM,JM)+P(IM,IN)*W(IN)*P(JM,IN)
          END DO
       END DO
    END DO
    CALL MATINV(M,DET,ATEMP,ATEMP,.FALSE.)
    CALL MATMAT(M,M,N,ATEMP,B,U2,0)
  END SUBROUTINE MLS_DETERMU2
  
!--------------------------------------
!*** DETERMINES GENERIC SHAPE FUNCTIONS
!*** AND DERIVATIVES
!*** CHK0
!--------------------------------------
  SUBROUTINE MLS_SFDER(NDI,N,M,POLYN,DPOLYN,U2,FF,DFF)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::M,N
    REAL(8),DIMENSION(M)::POLYN
    REAL(8),DIMENSION(NDI,M)::DPOLYN
    REAL(8),DIMENSION(M,N)::U2
    REAL(8),DIMENSION(N)::FF
    REAL(8),DIMENSION(NDI,N)::DFF
    DO IN=1,N
       FF(IN)=DOTPROD(M,POLYN(1:M),U2(1:M,IN))
       DO ID=1,NDI
          DFF(ID,IN)=DOTPROD(M,DPOLYN(ID,1:M),U2(1:M,IN))
       END DO
    END DO
  END SUBROUTINE MLS_SFDER

!--------------------
!*** POLYNOMIAL BASES
!*** CHK0 
!--------------------
  SUBROUTINE MLS_POLYNQUAD1D(D,X,XBAR,POLYN,DPOLYN)
    IMPLICIT NONE
    DOUBLE PRECISION V(24),D,X(1),XBAR(1),POLYN(3),DPOLYN(1,3)
    V(19)=1D0/D
    V(18)=X(1)-XBAR(1)
    V(17)=1D0/D**2
    POLYN(1)=1D0
    POLYN(2)=V(18)*V(19)
    POLYN(3)=V(17)*(V(18)*V(18))
    DPOLYN(1,1)=0D0
    DPOLYN(1,2)=V(19)
    DPOLYN(1,3)=2D0*V(17)*V(18)
  END SUBROUTINE MLS_POLYNQUAD1D

  SUBROUTINE MLS_POLYNQUAD2D(D,X,XBAR,POLYN,DPOLYN)
    IMPLICIT NONE
    DOUBLE PRECISION V(48),D,X(2),XBAR(2),POLYN(6),DPOLYN(2,6)
    V(41)=1D0/D
    V(40)=X(2)-XBAR(2)
    V(39)=X(1)-XBAR(1)
    V(38)=1D0/D**2
    V(43)=2D0*V(38)
    V(42)=V(38)*V(39)
    POLYN(1)=1D0
    POLYN(2)=V(39)*V(41)
    POLYN(3)=V(40)*V(41)
    POLYN(4)=V(38)*(V(39)*V(39))
    POLYN(5)=V(38)*(V(40)*V(40))
    POLYN(6)=V(40)*V(42)
    DPOLYN(1,1)=0D0
    DPOLYN(1,2)=V(41)
    DPOLYN(1,3)=0D0
    DPOLYN(1,4)=V(39)*V(43)
    DPOLYN(1,5)=0D0
    DPOLYN(1,6)=V(38)*V(40)
    DPOLYN(2,1)=0D0
    DPOLYN(2,2)=0D0
    DPOLYN(2,3)=V(41)
    DPOLYN(2,4)=0D0
    DPOLYN(2,5)=V(40)*V(43)
    DPOLYN(2,6)=V(42)
  END SUBROUTINE MLS_POLYNQUAD2D

  SUBROUTINE MLS_POLYNQUAD3D(D,X,XBAR,POLYN,DPOLYN)
    IMPLICIT NONE
    DOUBLE PRECISION V(78),D,X(3),XBAR(*),POLYN(*),DPOLYN(3,*)
    V(73)=1D0/D
    V(72)=X(3)-XBAR(3)
    V(71)=X(2)-XBAR(2)
    V(70)=X(1)-XBAR(1)
    V(69)=1D0/D**2
    V(62)=V(69)*V(71)
    V(63)=V(69)*V(72)
    V(65)=V(69)*V(70)
    POLYN(1)=1D0
    POLYN(2)=V(70)*V(73)
    POLYN(3)=V(71)*V(73)
    POLYN(4)=V(72)*V(73)
    POLYN(5)=V(69)*(V(70)*V(70))
    POLYN(6)=V(69)*(V(71)*V(71))
    POLYN(7)=V(69)*(V(72)*V(72))
    POLYN(8)=V(65)*V(71)
    POLYN(9)=V(65)*V(72)
    POLYN(10)=V(63)*V(71)
    DPOLYN(1,1)=0D0
    DPOLYN(1,2)=V(73)
    DPOLYN(1,3)=0D0
    DPOLYN(1,4)=0D0
    DPOLYN(1,5)=2D0*V(65)
    DPOLYN(1,6)=0D0
    DPOLYN(1,7)=0D0
    DPOLYN(1,8)=V(62)
    DPOLYN(1,9)=V(63)
    DPOLYN(1,10)=0D0
    DPOLYN(2,1)=0D0
    DPOLYN(2,2)=0D0
    DPOLYN(2,3)=V(73)
    DPOLYN(2,4)=0D0
    DPOLYN(2,5)=0D0
    DPOLYN(2,6)=2D0*V(62)
    DPOLYN(2,7)=0D0
    DPOLYN(2,8)=V(65)
    DPOLYN(2,9)=0D0
    DPOLYN(2,10)=V(63)
    DPOLYN(3,1)=0D0
    DPOLYN(3,2)=0D0
    DPOLYN(3,3)=0D0
    DPOLYN(3,4)=V(73)
    DPOLYN(3,5)=0D0
    DPOLYN(3,6)=0D0
    DPOLYN(3,7)=2D0*V(63)
    DPOLYN(3,8)=0D0
    DPOLYN(3,9)=V(65)
    DPOLYN(3,10)=V(62)
  END SUBROUTINE MLS_POLYNQUAD3D
  
!-------------------------------
!*** MLS_POLYNQUADBASE DETERMINES
!*** POLY. BASE + DERIVATIVES
!*** 1D TO 3D
!*** CHK0
!-------------------------------  
  SUBROUTINE MLS_POLYNQBASE(NDI,DTEMP,X,XBAR,POLYN,DPOLYN,M)
    IMPLICIT REAL(8) (a-h,o-z)
!*** CHECK THIS PARAMETER FOR UPDATES (MPOL)
    REAL(8)::D,DTEMP
    INTEGER,PARAMETER::MPOL=10
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(4)::RS
    D=DTEMP
    SELECT CASE(NDI)
    CASE(1)
       M=3
       CALL MLS_POLYNQUAD1D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M))
    CASE(2)
       M=6
       CALL MLS_POLYNQUAD2D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M))
    CASE(3)      
       M=4!10       
       CALL MLS_POLYNQUAD3D(D,X,XBAR,POLYN(1:M),DPOLYN(1:NDI,1:M))
    CASE DEFAULT
       STOP "CASE NOT FOUND IN MLS_POLYNQBASE"
    END SELECT
  END SUBROUTINE MLS_POLYNQBASE
  
!----------------------------
!*** WEIGHT FUNCTION FOR MLS
!*** NDI
!*** CHK0
!----------------------------
  SUBROUTINE MLS_WEIGHT(NDI,D,TOL,XI,X,W)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8),DIMENSION(NDI)::XI,X
    REAL(8),PARAMETER::PI=4.0D00*ATAN(1.0D00)
    REAL(8),PARAMETER::BOTTOM=1.0D-16
    S=RNORM2(NDI,XI-X)
    WMAX=1.0D00/(TOL**2.0D00+TOL**4.0D00)
    WBOT=BOTTOM*WMAX
!!    W=MAX((1.0D00/((S*S)/(D*D)+TOL*TOL))-(1.0D00/(1.0D00+TOL*TOL)),WBOT)
!!    W=(1.0D00/((S*S)/(D*D)+TOL*TOL))-(1.0D00/(1.0D00+TOL*TOL))
!!    W=MAX(EXP(-S/(1.0D-6*D)),1.0d-13)
!!!    W=SQRT(1.0d00/(TOL*PI))*EXP(-((S/D)**2)/TOL)
    W=1.0d00/((S/D)+TOL)-1.0d00/(1.0d00+TOL)
!    W=MAX(W,1.0d-40)
!    W=1.0d00
!    w=1.0d00
!    W=(1.0D00/((S*S)/(D*D)+TOL*TOL))-(1.0D00/(1.0D00+TOL*TOL))
!    w=1.0d00
  END SUBROUTINE MLS_WEIGHT

!-------------------------------------------------
!*** DETERMINES U2 EVALUATED IN A GIVEN COORDINATE
!    CALL U2ATACOORDINATE(D,TOL,NDI,N,XN,X,U2)
!-------------------------------------------------
  SUBROUTINE MLS_U2ATACOORDINATE(D,TOL,NDI,N,XBAR,XN,X,U2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8)::D,TOL
    INTEGER,PARAMETER::MPOL=10
    REAL(8),DIMENSION(N)::W
    REAL(8),DIMENSION(MPOL,N)::P
    REAL(8),DIMENSION(MPOL,N)::U2
    REAL(8),DIMENSION(NDI,N)::XN
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(3,MPOL)::TRASH
!-------------------------------
!*** NOW WE MUST DEFINE MATRIX P
!*** DERIVATIVES ARE DISCARDED
!-------------------------------
    DO IN=1,N
       CALL MLS_POLYNQBASE(NDI,D,XN(1:NDI,IN),XBAR(1:NDI),P(1:MPOL,IN),TRASH(1:NDI,1:MPOL),M)
    END DO
!------------------------
!*** NOW WE MUST DEFINE W
!------------------------
    DO IN=1,N
       CALL MLS_WEIGHT(NDI,D,TOL,XN(1:NDI,IN),X(1:NDI),W(IN))
    END DO
!------------------------------
!*** FINALIZE IT WITH SPECIFICS
!------------------------------
    CALL MLS_DETERMU2(M,N,W(1:N),P(1:M,1:N),U2(1:M,1:N))
  END SUBROUTINE MLS_U2ATACOORDINATE

!------------------------------------------
!*** OBTAIN SHAPE FUNCTIONS AND DERIVATIVES
!*** FROM A GIVEN LIST OF NODES
!*** AT ONE POINT WITH COORDINATES X
!*** CHK0         
!------------------------------------------

  SUBROUTINE MLS_SHAPEFUNCTIONSQUAD(D,TOL,NDI,N,XN,X,FF,DFF)
    IMPLICIT REAL(8) (a-h,o-z)
    REAL(8)::D,TOL
    INTEGER,PARAMETER::MPOL=10
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(N)::FF,W
    REAL(8),DIMENSION(NDI,N)::DW
    REAL(8),DIMENSION(MPOL,N)::P
    REAL(8),DIMENSION(MPOL,N)::U2
    REAL(8),DIMENSION(NDI,N)::XN
    REAL(8),DIMENSION(NDI)::X,XBAR
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(3,MPOL)::TRASH
    REAL(8),DIMENSION(3,MPOL)::DPOLYN
    DO ID=1,NDI
       XBAR(ID)=0.0D00
    END DO
    TOP=1.0D00/(1.0D00*N)
    DO IN=1,N
       DO ID=1,NDI
          XBAR(ID)=XBAR(ID)+XN(ID,IN)*TOP
       END DO
    END DO
!--------------------
!*** MEAN COORDINATES
!--------------------    
    CALL MLS_U2ATACOORDINATE(D,TOL,NDI,N,XBAR(1:NDI),XN(1:NDI,1:N),X(1:NDI),U2(1:MPOL,1:N))
!------------------------------    
!*** POLYNOMIAL AND DERIVATIVES
!------------------------------
    CALL MLS_POLYNQBASE(NDI,D,X(1:NDI),XBAR(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),M)
    CALL MLS_SFDER(NDI,N,M,POLYN(1:M),DPOLYN(1:NDI,1:M),U2(1:M,1:N),FF,DFF)
  END SUBROUTINE MLS_SHAPEFUNCTIONSQUAD

!------------------------
!*** GETS INFLUENCE NODES
!------------------------
  SUBROUTINE MLS_GETSNODES(D,TOL,X,NDI,NTOT,XTOT,N,LISTN,XN,FF,DFF)
    IMPLICIT REAL(8)(A-H,O-Z)
    REAL(8),PARAMETER::CUTOFF=1.0D-3
    REAL(8)::D,TOL
    REAL(8),DIMENSION(NDI,*)::XTOT
    REAL(8),DIMENSION(NDI)::X
    INTEGER::NDI
    INTEGER,DIMENSION(:),ALLOCATABLE::LISTN
    REAL(8),DIMENSION(:,:),ALLOCATABLE::XN
    REAL(8),DIMENSION(:,:),ALLOCATABLE::DFF
    REAL(8),DIMENSION(:),ALLOCATABLE::FF
    N=0
    DO ITOT=1,NTOT
       IF(RNORM2(NDI,XTOT(1:NDI,ITOT)-X(1:NDI)).LE.D)THEN
          N=N+1
       END IF
    END DO
    ALLOCATE(LISTN(N))
    ALLOCATE(XN(NDI,N),DFF(NDI,N),FF(N))
!*** NOW INSERTS THE NODES IN THE LIST
    N=0
    DO ITOT=1,NTOT
       IF(RNORM2(NDI,XTOT(1:NDI,ITOT)-X(1:NDI)).LE.D)THEN
          N=N+1
          LISTN(N)=ITOT
          XN(1:NDI,N)=XTOT(1:NDI,ITOT)
       END IF
    END DO
  END SUBROUTINE MLS_GETSNODES
  
!-------------------------------------------------
!*** DETERMINES THE DEFORMATION GRADIENT
!*** AND GREEN-LAGRANGE STRAIN IN ENGINEERING FORM
!*** MLS STUFF
!*** CHK0,1
!-------------------------------------------------
  SUBROUTINE MLS_STRAIN(NDI,N,XDEF,DFF,F,STRAIN)
    IMPLICIT REAL(8) (a-h,o-z)
    REAL(8),DIMENSION(NDI,N)::XDEF
    REAL(8),DIMENSION(NDI,NDI)::C,F
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(NDI*(NDI+1)/2)::STRAIN
!*** DEFORMATION GRADIENT
    F=0.0D00
    DO IN=1,N
       DO JD=1,NDI
          DO ID=1,NDI
             F(ID,JD)=F(ID,JD)+DFF(JD,IN)*XDEF(ID,IN)
          END DO
       END DO
    END DO
!*** RIGHT CAUCHY-GREEN TENSOR
    CALL MATMAT(NDI,NDI,NDI,F,F,C,2)
    DO ID=1,NDI*(NDI+1)/2
       CALL APOMAT(I1,I2,ID,NDI)
       STRAIN(ID)=0.5D00*(C(I1,I2)-DELTAK(I1,I2))
    END DO
    DO ID=NDI+1,NDI*(NDI+1)/2
       STRAIN(ID)=2.0D00*STRAIN(ID)
    END DO
  END SUBROUTINE MLS_STRAIN

  SUBROUTINE MLS_2DFORCE(NK,F,S,NUCLEUS)
    IMPLICIT NONE
    DOUBLE PRECISION V(27),NK(2),F(2,2),S(3),NUCLEUS(2)
    V(22)=NK(2)*S(2)
    V(21)=NK(1)*S(1)
    NUCLEUS(1)=(F(1,2)*NK(1)+F(1,1)*NK(2))*S(3)+F(1,1)*V(21)+F(1,2)*V(22)
    NUCLEUS(2)=(F(2,2)*NK(1)+F(2,1)*NK(2))*S(3)+F(2,1)*V(21)+F(2,2)*V(22)
  END SUBROUTINE MLS_2DFORCE

  SUBROUTINE MLS_2D(NK,NL,F,S,DS,KERNEL)
    IMPLICIT NONE
    DOUBLE PRECISION V(67),NK(2),NL(2),F(2,2),S(3),DS(3,3),KERNEL(2,2)
    V(53)=NK(2)*(NL(2)*S(2)+NL(1)*S(3))+NK(1)*(NL(1)*S(1)+NL(2)*S(3))
    V(52)=F(2,2)*NK(1)+F(2,1)*NK(2)
    V(51)=F(1,2)*NK(1)+F(1,1)*NK(2)
    V(50)=F(2,2)*NK(2)
    V(49)=F(1,2)*NK(2)
    V(48)=F(2,1)*NK(1)
    V(47)=F(1,1)*NK(1)
    V(46)=F(2,2)*NL(1)+F(2,1)*NL(2)
    V(45)=F(1,2)*NL(1)+F(1,1)*NL(2)
    V(44)=F(2,2)*NL(2)
    V(60)=DS(2,2)*V(44)
    V(59)=DS(1,2)*V(44)+DS(1,3)*V(46)
    V(43)=F(1,2)*NL(2)
    V(54)=DS(1,2)*V(43)+DS(1,3)*V(45)
    V(42)=F(2,1)*NL(1)
    V(62)=DS(3,1)*V(42)+DS(3,2)*V(44)+DS(3,3)*V(46)
    V(61)=DS(2,1)*V(42)+DS(2,3)*V(46)
    V(58)=DS(1,1)*V(42)
    V(41)=F(1,1)*NL(1)
    V(56)=DS(3,1)*V(41)+DS(3,2)*V(43)+DS(3,3)*V(45)
    V(55)=DS(2,1)*V(41)+DS(2,3)*V(45)
    V(38)=V(47)*V(58)
    V(39)=V(49)*V(60)
    V(57)=V(38)+V(39)
    KERNEL(1,1)=V(53)+V(47)*(DS(1,1)*V(41)+V(54))+V(49)*(DS(2,2)*V(43)+V(55))+V(51)*V(56)
    KERNEL(1,2)=V(57)+V(47)*V(59)+V(49)*V(61)+V(51)*V(62)
    KERNEL(2,1)=V(48)*V(54)+V(50)*V(55)+V(52)*V(56)+V(57)
    KERNEL(2,2)=V(53)+V(48)*(V(58)+V(59))+V(50)*(V(60)+V(61))+V(52)*V(62)
  END SUBROUTINE MLS_2D

  SUBROUTINE MLS_3DFORCE(NK,F,S,NUCLEUS)
    IMPLICIT NONE
    DOUBLE PRECISION V(51),NK(3),F(3,3),S(6),NUCLEUS(3)
    V(46)=NK(3)*S(3)
    V(45)=NK(2)*S(2)
    V(44)=NK(1)*S(1)
    NUCLEUS(1)=(F(1,2)*NK(1)+F(1,1)*NK(2))*S(4)+(F(1,3)*NK(1)+F(1,1)*NK(3))*S(5)+(F(1,3)*NK(2)+F(1,2)*NK(3))*S(6)+F(1,1)*V&
         &(44)+F(1,2)*V(45)+F(1,3)*V(46)
    NUCLEUS(2)=(F(2,2)*NK(1)+F(2,1)*NK(2))*S(4)+(F(2,3)*NK(1)+F(2,1)*NK(3))*S(5)+(F(2,3)*NK(2)+F(2,2)*NK(3))*S(6)+F(2,1)*V&
         &(44)+F(2,2)*V(45)+F(2,3)*V(46)
    NUCLEUS(3)=(F(3,2)*NK(1)+F(3,1)*NK(2))*S(4)+(F(3,3)*NK(1)+F(3,1)*NK(3))*S(5)+(F(3,3)*NK(2)+F(3,2)*NK(3))*S(6)+F(3,1)*V&
         &(44)+F(3,2)*V(45)+F(3,3)*V(46)
  END SUBROUTINE MLS_3DFORCE
  
  SUBROUTINE MLS_3D(NK,NL,F,S,DS,KERNEL)
    IMPLICIT NONE
    DOUBLE PRECISION V(173),NK(3),NL(3),F(3,3),S(6),DS(6,6),KERNEL(3,3)
    V(144)=NK(1)*(NL(1)*S(1)+NL(2)*S(4)+NL(3)*S(5))+NK(3)*(NL(3)*S(3)+NL(1)*S(5)+NL(2)*S(6))+NK(2)*(NL(2)*S(2)+NL(1)*S(4)&
         &+NL(3)*S(6))
    V(143)=F(3,3)*NK(2)+F(3,2)*NK(3)
    V(142)=F(2,3)*NK(2)+F(2,2)*NK(3)
    V(141)=F(1,3)*NK(2)+F(1,2)*NK(3)
    V(140)=F(3,3)*NK(1)+F(3,1)*NK(3)
    V(139)=F(2,3)*NK(1)+F(2,1)*NK(3)
    V(138)=F(1,3)*NK(1)+F(1,1)*NK(3)
    V(137)=F(3,2)*NK(1)+F(3,1)*NK(2)
    V(136)=F(2,2)*NK(1)+F(2,1)*NK(2)
    V(135)=F(1,2)*NK(1)+F(1,1)*NK(2)
    V(134)=F(3,3)*NK(3)
    V(133)=F(2,3)*NK(3)
    V(132)=F(1,3)*NK(3)
    V(131)=F(3,2)*NK(2)
    V(130)=F(2,2)*NK(2)
    V(129)=F(1,2)*NK(2)
    V(128)=F(3,1)*NK(1)
    V(127)=F(2,1)*NK(1)
    V(126)=F(1,1)*NK(1)
    V(125)=F(3,3)*NL(2)+F(3,2)*NL(3)
    V(124)=F(2,3)*NL(2)+F(2,2)*NL(3)
    V(123)=F(1,3)*NL(2)+F(1,2)*NL(3)
    V(122)=F(3,3)*NL(1)+F(3,1)*NL(3)
    V(121)=F(2,3)*NL(1)+F(2,1)*NL(3)
    V(120)=F(1,3)*NL(1)+F(1,1)*NL(3)
    V(119)=F(3,2)*NL(1)+F(3,1)*NL(2)
    V(118)=F(2,2)*NL(1)+F(2,1)*NL(2)
    V(117)=F(1,2)*NL(1)+F(1,1)*NL(2)
    V(116)=F(3,3)*NL(3)
    V(115)=F(2,3)*NL(3)
    V(168)=DS(3,3)*V(115)
    V(114)=F(1,3)*NL(3)
    V(147)=DS(3,3)*V(114)
    V(113)=F(3,2)*NL(2)
    V(160)=DS(1,2)*V(113)+DS(1,3)*V(116)+DS(1,4)*V(119)+DS(1,5)*V(122)+DS(1,6)*V(125)
    V(112)=F(2,2)*NL(2)
    V(167)=DS(2,2)*V(112)
    V(154)=DS(1,2)*V(112)+DS(1,3)*V(115)+DS(1,4)*V(118)+DS(1,5)*V(121)+DS(1,6)*V(124)
    V(111)=F(1,2)*NL(2)
    V(148)=DS(1,2)*V(111)+DS(1,3)*V(114)+DS(1,4)*V(117)+DS(1,5)*V(120)+DS(1,6)*V(123)
    V(146)=DS(2,2)*V(111)
    V(110)=F(3,1)*NL(1)
    V(165)=DS(6,1)*V(110)+DS(6,2)*V(113)+DS(6,3)*V(116)+DS(6,4)*V(119)+DS(6,5)*V(122)+DS(6,6)*V(125)
    V(164)=DS(5,1)*V(110)+DS(5,2)*V(113)+DS(5,3)*V(116)+DS(5,4)*V(119)+DS(5,5)*V(122)+DS(5,6)*V(125)
    V(163)=DS(4,1)*V(110)+DS(4,2)*V(113)+DS(4,3)*V(116)+DS(4,4)*V(119)+DS(4,5)*V(122)+DS(4,6)*V(125)
    V(162)=DS(3,1)*V(110)+DS(3,2)*V(113)+DS(3,4)*V(119)+DS(3,5)*V(122)+DS(3,6)*V(125)
    V(161)=DS(2,1)*V(110)+DS(2,3)*V(116)+DS(2,4)*V(119)+DS(2,5)*V(122)+DS(2,6)*V(125)
    V(109)=F(2,1)*NL(1)
    V(166)=DS(1,1)*V(109)
    V(159)=DS(6,1)*V(109)+DS(6,2)*V(112)+DS(6,3)*V(115)+DS(6,4)*V(118)+DS(6,5)*V(121)+DS(6,6)*V(124)
    V(158)=DS(5,1)*V(109)+DS(5,2)*V(112)+DS(5,3)*V(115)+DS(5,4)*V(118)+DS(5,5)*V(121)+DS(5,6)*V(124)
    V(157)=DS(4,1)*V(109)+DS(4,2)*V(112)+DS(4,3)*V(115)+DS(4,4)*V(118)+DS(4,5)*V(121)+DS(4,6)*V(124)
    V(156)=DS(3,1)*V(109)+DS(3,2)*V(112)+DS(3,4)*V(118)+DS(3,5)*V(121)+DS(3,6)*V(124)
    V(155)=DS(2,1)*V(109)+DS(2,3)*V(115)+DS(2,4)*V(118)+DS(2,5)*V(121)+DS(2,6)*V(124)
    V(108)=F(1,1)*NL(1)
    V(153)=DS(6,1)*V(108)+DS(6,2)*V(111)+DS(6,3)*V(114)+DS(6,4)*V(117)+DS(6,5)*V(120)+DS(6,6)*V(123)
    V(152)=DS(5,1)*V(108)+DS(5,2)*V(111)+DS(5,3)*V(114)+DS(5,4)*V(117)+DS(5,5)*V(120)+DS(5,6)*V(123)
    V(151)=DS(4,1)*V(108)+DS(4,2)*V(111)+DS(4,3)*V(114)+DS(4,4)*V(117)+DS(4,5)*V(120)+DS(4,6)*V(123)
    V(150)=DS(3,1)*V(108)+DS(3,2)*V(111)+DS(3,4)*V(117)+DS(3,5)*V(120)+DS(3,6)*V(123)
    V(149)=DS(2,1)*V(108)+DS(2,3)*V(114)+DS(2,4)*V(117)+DS(2,5)*V(120)+DS(2,6)*V(123)
    V(145)=DS(1,1)*V(108)
    V(104)=V(127)*V(145)+V(130)*V(146)+V(133)*V(147)
    V(105)=V(128)*V(145)+V(131)*V(146)+V(134)*V(147)
    V(106)=V(128)*V(166)+V(131)*V(167)+V(134)*V(168)
    KERNEL(1,1)=V(144)+V(126)*(V(145)+V(148))+V(129)*(V(146)+V(149))+V(132)*(V(147)+V(150))+V(135)*V(151)+V(138)*V(152)+V&
         &(141)*V(153)
    KERNEL(1,2)=V(104)+V(126)*V(154)+V(129)*V(155)+V(132)*V(156)+V(135)*V(157)+V(138)*V(158)+V(141)*V(159)
    KERNEL(1,3)=V(105)+V(126)*V(160)+V(129)*V(161)+V(132)*V(162)+V(135)*V(163)+V(138)*V(164)+V(141)*V(165)
    KERNEL(2,1)=V(104)+V(127)*V(148)+V(130)*V(149)+V(133)*V(150)+V(136)*V(151)+V(139)*V(152)+V(142)*V(153)
    KERNEL(2,2)=V(144)+V(136)*V(157)+V(139)*V(158)+V(142)*V(159)+V(127)*(V(154)+V(166))+V(130)*(V(155)+V(167))+V(133)*(V&
         &(156)+V(168))
    KERNEL(2,3)=V(106)+V(127)*V(160)+V(130)*V(161)+V(133)*V(162)+V(136)*V(163)+V(139)*V(164)+V(142)*V(165)
    KERNEL(3,1)=V(105)+V(128)*V(148)+V(131)*V(149)+V(134)*V(150)+V(137)*V(151)+V(140)*V(152)+V(143)*V(153)
    KERNEL(3,2)=V(106)+V(128)*V(154)+V(131)*V(155)+V(134)*V(156)+V(137)*V(157)+V(140)*V(158)+V(143)*V(159)
    KERNEL(3,3)=V(144)+V(128)*(DS(1,1)*V(110)+V(160))+V(131)*(DS(2,2)*V(113)+V(161))+V(134)*(DS(3,3)*V(116)+V(162))+V(137&
         &)*V(163)+V(140)*V(164)+V(143)*V(165)
  END SUBROUTINE MLS_3D
