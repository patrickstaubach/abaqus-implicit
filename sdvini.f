      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      implicit none
C
      REAL(8) STATEV(NSTATV),COORDS(NCRDS)
      CHARACTER NAMES(2)*80
      real*8,DIMENSION(23395) :: ELLAB
      integer NOEL,ntens,i,ixx,jxx,nstatev,NCRDS,NPT,LAYER,KSPT,NSTATV

      statev(:)   = 0.0d0
      statev(1)   = 0.8278d0
      statev(1+3) = -0.0001d0
      RETURN
      END
