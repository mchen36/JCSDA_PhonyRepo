MODULE Durden_Model
  USE CSEM_Type_Kinds,       ONLY: fp => CSEM_fp

  IMPLICIT NONE
  PRIVATE
  REAL(fp), PARAMETER :: pi = 3.141592653589793
  
  REAL(fp), PARAMETER :: a=0.225, b=1.25, g=9.81
  REAL(fp), PARAMETER :: a0=0.008, gamma=7.25E-5
  REAL(fp), PARAMETER :: s=1.5E-4, kj=2.0
  REAL(fp)            :: ustar, kc, a0r, c
  
  PUBLIC :: W_ocean_Durden_init
  PUBLIC :: W_ocean_Durden
  PUBLIC :: search_kd
  PUBLIC :: slope_var_durden
  PUBLIC :: Pdf_Slope_Durden
  
! Author:Ming Chen, Dec-15-2015
!       ming.chen@noaa.gov
  
CONTAINS

!================================================================================
! This subroutine is used to initlize the Durden ocean wave model
! Author:Ming Chen, Dec-15-2015
!       ming.chen@noaa.gov
!================================================================================
  SUBROUTINE W_ocean_Durden_init(U,z)
    REAL(fp), INTENT(IN) :: u, z
    REAL(fp) :: z12 = 12.5, z19=19.5, u12, u19
    REAL(fp) :: z05 = 5.0,  z10=10.0, u05, u10
    REAL(fp) :: D, R
    REAL(fp),SAVE :: ul = -999.0, zl = -999.0
    
    IF(ul == U .AND. zl == z) RETURN
    ul = U; zl = z
    
    CALL Cal_Ustar(u, z, ustar)
    
    ! The wind speed at 19.5/12.5 meter height is calculated by
    ! the following Pierson's formula
    u19 = Ustar2Uz(ustar, z19)
    u12 = Ustar2Uz(ustar, z12)
    u10 = Ustar2Uz(ustar, z10)
    u05 = Ustar2Uz(ustar, z05)
    print*, "U05, U10", u05, u10
    
    kc=g/u19**2
    a0r=exp(0.74*(kc/kj)**2)  ! a0r=1

    ! A typo in Durden's paper for R
    R=(0.003+1.92E-3*u12)/(3.16E-3*u12)

    ! The integrals for sumD2 does converge.
    ! check with Stephen Durden to find out the cutoff
    ! frequency for his calculation.
    CALL D_durden(D)
    
    ! There is a typo in Stephen Durden's paper
    ! for the formula for c, which should be increased by a factor of 2.
    C=2*(1-R)/(1+R)/(1-D)
    
  END SUBROUTINE W_ocean_Durden_init
  
!================================================================================
! This subroutine is used to evaluate wind speed from the ocean roughnes
! Author:Ming Chen, Dec-15-2015
!       ming.chen@noaa.gov
!================================================================================
  FUNCTION Ustar2Uz(ustar, z) RESULT(uz)
  
    REAL(fp), INTENT(IN) :: ustar
    REAL(fp), INTENT(IN) :: z
    REAL(fp) :: a,b,c, z0,uz

    a = 6.84e-5 ; b = 4.28e-3 ; c = -4.43e-4
    Z0=a/ustar+b*ustar*ustar+c
    uz = ustar/0.4*log(Z/Z0)

  END FUNCTION Ustar2Uz
  
!================================================================================
! This subroutine is used to evaluate the ocean roughnes given the wind speed
! Author:Ming Chen, Dec-15-2015
!       ming.chen@noaa.gov
!================================================================================
  SUBROUTINE Cal_Ustar(u, z, ustar)
   
    REAL(fp), INTENT(IN)  :: u
    REAL(fp), INTENT(IN)  :: z
    REAL(fp), INTENT(OUT) :: ustar

    REAL(fp) :: up,ul,uz
    
    up=5.0; ul=0.0 ; uz = 0.0
    DO WHILE (abs(u-uz)/u > 0.02 )
      ustar = (up+ul)/2.0
      uz = Ustar2Uz(ustar, z)
      IF (u-uz > 0.0 ) THEN
        ul = ustar
      ELSE
        up = ustar
      ENDIF
      
    ENDDO

  END SUBROUTINE Cal_Ustar
  
!================================================================================
! This subroutine is used to numerically calculate
! the value of D.
! Input:
!	kd=high frequency cutoff for the large scale wave
! Output:

! Author: S. Yueh, Nov-3-92
!================================================================================

  SUBROUTINE D_durden(D)
    REAL(fp), INTENT(OUT) ::  D
    REAL(fp) ::  kd
    REAL(fp) ::  rlog2,rlog10,rlogkd, rlogkj, rlogbug
    REAL(fp) ::  dlogk, sumf, sumS, sumf_low, SumS_low
    REAL(fp) ::  x,y,gs, z, z1

    INTEGER :: i, NK,cn,n

    kd=(b*ustar**2/gamma+2)*2
    Nk=INT(kd/10)
    IF ( Nk .LE. 1000) Nk=1000
    rlog2=log(2.0)
    rlog10=log(10.0)
    rlogkd=log(kd)
    rlogkj=log(kj)
    dlogk=(rlogkd-rlogkj)/Nk

    sumf=0
    sumS=0

    DO i=1,Nk
      x=(i-0.5)*dlogk+rlogkj
      gs=g+gamma*exp(2*x)
      rlogbug=log(b*ustar**2/gs)
      y=exp(a*(x-rlog2)/rlog10*(x+rlogbug))
      sumS=sumS+y
      sumf=sumf+y*(1-exp(-s*exp(2*x)))
    ENDDO
    sumS=sumS*dlogk
    sumf=sumf*dlogk

    y=0.74*kc**2*s
    x=0.74*kc**2/kj**2

    sumf_low=0
    z=1.0
    cn=0
    z1 = 1.0
    DO while ( abs(z1) >= 1.E-7)
      cn=cn+1
      z=z*y/cn
      n=INT(cn+0.5)
      z1=z*Ein(-x,n)
      sumf_low=sumf_low+z1
    ENDDO

    sumf=sumf+sumf_low/2*a0r

    SumS_low=-Ei(-x)
    SumS=SumS+SumS_low/2*a0r

    SumS=SumS*a0
    Sumf=Sumf*a0

    D=1-Sumf/SumS
    RETURN
  END SUBROUTINE D_durden


  FUNCTION Ein(x,n) RESULT(E)

    REAL(fp) :: x
    REAL(fp) :: E
    INTEGER  :: i, n,l
    REAL(fp) ::  sum, xl
  
    E=Ei(x)
    DO i=1,n
      E=E/i
    ENDDO
    sum=0.0
    xl=1.0/x
    DO l=0,n-1
      xl=xl*x/(n-l)
      sum=sum+xl
    ENDDO
    E=E-exp(x)/x**n*sum
    RETURN
  END FUNCTION Ein

  FUNCTION Ei(x) RESULT(E)
 
    REAL(fp) :: x,xl,y
    REAL(fp) :: E

    REAL(fp),PARAMETER :: c = .577215
    INTEGER :: cl
    ! x < 0 is required.
    IF ( x .LE. 0.0) THEN
      E=C+log(-x)
      cl = 1
      xl = x
      y = x
      DO while (abs(y) .GE. 1.E-7)
        E=E+y
        cl=cl+1
        xl=xl*x/cl
        y=xl/cl
      END DO

    ELSE
      E = 0.0
      WRITE(*,*)'The input argument out of range'
    ENDIF
    RETURN
  END FUNCTION Ei
  
  
!=============================================================
! This routine searches the cutoff wavenumber
! which gives the specified value of the
! product of k0 and small scale rms height
! Input:
!	f: frequency in GHz
!	ksigma: k0*sigma
! output:
!	kd: cutoff wave number in unit 1/m
! coded by Simon H. Yueh, Jan-2-93
! Copyright: all rights reserved
! temp:
!	Kd1: Initial lower bound
!	Kd2: Initial upper bound
!=============================================================
  SUBROUTINE search_kd(f,kd,ksi)
    REAL(fp) :: sdd,ak,bb,dk,f,kd
    INTEGER :: i
    REAL(fp),optional :: ksi
    sdd=0.
    dk=1.0
    DO i=10000,1,-1
      ak=dk*i
      sdd=sdd+Sk_Durden(ak)*dk
      bb=(ak/((2.*pi*f/0.3)**2))**2
      IF (sdd .GT. bb) go to 281
    END DO
    IF (i .LE. 1) THEN
      PRINT *,' no data found '
      go to 282
    END IF
281 continue
282 continue
    kd=ak
  
    IF(PRESENT(ksi)) ksi = SQRT(sdd)
    RETURN
  END SUBROUTINE search_kd


!===============================================================
! This subroutine is used to numerically calculate
! the variances of the large scale surface with
! a given cutoff frequency.
! Input:
!	kd=high frequency cutoff for the large scale wave
! Output:
! 	Su: Slope variance in up wind direction
!	Sc: Slope variance in cross wind direction
!	St: Total slope variance
! Author: S. Yueh, Nov-2-92
!===============================================================
  SUBROUTINE slope_var_durden(kd,Su,Sc,St)
    REAL(fp) :: Su,Sc,St
    REAL(fp) :: kd
    REAL(fp) :: dk,k
    REAL(fp) :: Sumf, SumS, D, c_cal
    INTEGER ::  i, NK
    
    Nk=1000

    dk=kd/Nk
    SumS=0
    Sumf=0
    DO i=1,Nk
      k=(i-0.5)*dk
      SumS=SumS+Sk_durden(k)*k**2
      Sumf=Sumf+Sk_durden(k)*k**2*(1-exp(-s*k**2))
    ENDDO
    SumS=SumS*dk
    Sumf=Sumf*dk
    St=SumS
    Su=ABS(0.5*(St+0.5*c*Sumf))
    Sc=ABS(0.5*(St-0.5*c*Sumf))
    D=1-Sumf/SumS
    c_cal=2.0*(1.0-Sc/Su)/(1.0+Sc/Su)/(1.0-D)
   
  END SUBROUTINE slope_var_durden


  FUNCTION Sk_Durden(k) RESULT(sk)
    ! used together with W_ocean
    REAL(fp) :: k,sk
    REAL(fp) :: gs
    
    IF ( k < kj) THEN
      IF ( k <= kc*1.E-5)then
        Sk=0
      ELSE
        Sk=a0*a0r*k**(-3)*exp(-0.74*(kc/k)**2)
      ENDIF
    ELSE
      gs=g+gamma*k**2
      !Sk=a0*k**(-3)*(b*k*us**2/gs)**(a*log10(k/2.0))
      Sk=a0*k**(-3)*(b*k*ustar**2/gs)**(a*log10(k/kj))

    ENDIF
  END FUNCTION Sk_Durden
  

!==============================================================================
! This function returns the pdf of the slope of the large
! scale surface which is assumed to be Gaussian.
! This PDF was used by S. Durden for his two scale theory
! for a large scale surface.
! Input:
!	Zx=Slope in x (up-wind)
!	Zy=Slope in y (cross-wind)
!	Su=Slope variance in x (up-wind)
!	Sc=Slope variance in y (cross-wind)
! Author: S. Yueh, Nov-3-92
!==============================================================================
  FUNCTION Pdf_Slope_Durden(Su,Sc,Zx,Zy)
    REAL(fp) :: Su, Sc, Zx,Zy 
    REAL(fp) :: sx, sy
    REAL(fp) :: arg
    REAL(fp) :: Pdf_slope_durden
    
    Sx=sqrt(Su)
    Sy=sqrt(Sc)
    arg=(Zx**2/Su+Zy**2/Sc)*0.5
    Pdf_slope_durden=1/(2*pi*Sx*Sy)*exp(-arg)
   
  END FUNCTION Pdf_Slope_Durden
 
 
!===============================================================================  
! This function implements the ocean surface spectrum
! given in Li, F. K., Gregg Neumann, Scott Schaffer, and S. L. Durden,
! ``Studies of the location of azimuth modulation minimum for ku
! Band Ocean Radar Backscatter,'' JGR, Vol. 93, No. C7, pp. 8229-8238,
! July, 1988.
! The values of the parameters in the common block have to
! be initialized by executing W_ocean_init_Durden before running this
! function.
!================================================================================
  FUNCTION W_ocean_Durden(kx,ky)
    REAL(fp) :: k,kx,ky
    REAL(fp) :: cos2phi, SK, r, Psi
    REAL(fp) :: W_ocean_Durden

    k=sqrt(kx**2+ky**2)
    cos2phi=(kx**2-ky**2)/k**2
    Sk= Sk_Durden(k)
    r=c*(1-exp(-s*k**2))
    Psi=(1+r*cos2phi)/(2*pi)
    W_ocean_Durden=Sk*Psi/k
    RETURN
  END FUNCTION W_ocean_Durden
  
END MODULE Durden_Model
