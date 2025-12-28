      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock),
C
      character*80 cmname
C	  
      real*8 BULK_M, SHEAR_M
      real*8 LAMBDA1, LAMBDA2
      real*8 T1, T2
      real*8 TOL, DET_F
	  real*8 SIGMA11, SIGMA22, SIGMA33
	  
      parameter (TOL = 1.0D-10)
	  
      SHEAR_M = props(1)
      BULK_M = props(2)

      DO 100 km = 1,nblock
        LAMBDA1 = defgradNew(km, 1)
        LAMBDA2 = defgradNew(km, 3)
		
        IF (LAMBDA1 .LT. TOL) LAMBDA1 = TOL
        IF (LAMBDA2 .LT. TOL) LAMBDA2 = TOL
	  
        T1 = K_BULK * (1.0D0 - LAMBDA1 * LAMBDA2) + 
     &     MU_SHEAR * (LAMBDA1**2 - LAMBDA2**2) / 
     &     (2.0D0 * LAMBDA1**2 * LAMBDA2**2)
	 
        T2 = K_BULK * (1.0D0 - LAMBDA1 * LAMBDA2) + 
     &     MU_SHEAR * (LAMBDA2**2 - LAMBDA1**2) / 
     &     (2.0D0 * LAMBDA1**2 * LAMBDA2**2)
	 
	    DET_F = LAMBDA1 * LAMBDA2 * DEFGRADNEW(KM, 2)
		
		IF (DET_F .GT. TOL) THEN
          SIGMA11 = LAMBDA1 * T1 / DET_F   ! Axial Cauchy stress
          SIGMA33 = LAMBDA2 * T2 / DET_F   ! Hoop Cauchy stress
        ELSE
          SIGMA11 = 0.0D0
          SIGMA33 = 0.0D0
        END IF
		
		SIGMA22 = 0.0D0
		
		STRESSNEW(KM, 1) = SIGMA11
        STRESSNEW(KM, 2) = SIGMA22
        STRESSNEW(KM, 3) = SIGMA33
		
		IF (NSTATEV .GE. 4) THEN
          STATENEW(KM, 1) = LAMBDA1     ! Axial stretch
          STATENEW(KM, 2) = LAMBDA2     ! Hoop stretch
          STATENEW(KM, 3) = T1          ! Axial tension
          STATENEW(KM, 4) = T2          ! Hoop tension
        END IF
		
  100 continue

      return
      end