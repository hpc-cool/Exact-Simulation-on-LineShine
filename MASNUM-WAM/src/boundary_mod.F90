#include "wavedeff.h"
Module boundary_mod
  use varcommon_mod
  IMPLICIT NONE
  real(8) ,pointer::pvws(:,:) !计算优化 深度相关
  type  boundary_gbdy_type
     real(8) cosths(CJNTHET),sinths(CJNTHET);
     real(8) wk(CKL);
     real(8) windfield;
     integer nwpc
  end type boundary_gbdy_type
  type(boundary_gbdy_type) gbdy

CONTAINS
  SUBROUTINE InitSetspeca
#ifdef C_CALCULATE
    gbdy%cosths=cosths(1:jnthet)
    gbdy%sinths=sinths(1:jnthet)
    gbdy%wk=wk(1:kl)
    gbdy%nwpc=nwpc
    gbdy%windfield=windfield
    call initc_boundary(gbdy,pvws,nsp)
#endif
  END SUBROUTINE InitSetspeca
end module boundary_mod
SUBROUTINE InitSetspec(iacb,iace)
  use varcommon_mod
  use boundary_mod
  IMPLICIT NONE
  integer iacb,iace
  integer iac,k,ierr
  real*8 wkk,dk,tanhdk,wfk
  DO  iac=iacb,iace !{
    DO  k=1,kl !{
      wkk=wk(k)
      dk=dep(iac)*wkk
      tanhdk=1.
      if (dk.lt.4.) tanhdk=tanh(dk)
      wfk=sqrt(g*wkk*tanhdk)/zpi
      pvws(k,iac)=wfk*zpi
    ENDDO
  ENDDO
END SUBROUTINE
#ifndef C_CALCULATE
!风浪谱，n=1 全场初始化，n=2 边界
SUBROUTINE setspec(mode,e,nwpb,nwpe,n) !{
  use varcommon_mod
  use boundary_mod
  REALD e(kl,jnthet,0:*)
  integer nwpb,nwpe,mode
  real(8) gama,sq3,vx,vy,wv,xj0,xj,arlfa,wsj,wkj,wk0,wsk,wl,sigma,alpha
  integer n,j,k,iac
  gama=3.3
  sq3=sqrt(3.)
  do iac=nwpb,nwpe
    IF(nsp(iac)<n) CYCLE
    vx=wxy(1,iac)
    vy=wxy(2,iac)
    wv=wxy(3,iac)

    IF (wv<=0.) wv=0.9
    xj0=windfield*1000.
    xj=g*xj0/(wv**2)
    arlfa=(0.076*2*(xj**(-0.4)))/zpi
    wsj=22.*(xj**(-0.33))*g/wv
    wkj=wsj**2/g
    DO  j=1,jnthet !{
      DO  k=1,kl !{
        wk0=wk(k)
        wsk=pvws(k,iac)
        wl=vx*cosths(j)+vy*sinths(j)
        IF (wl>0) THEN !{
          IF (wsk<=wsj) THEN !{
            sigma=0.07
          ELSE !} {
            sigma=0.09
          end if !}
          alpha=arlfa/wk0**4*exp(-1.25*(wkj/wk0)**2)*gama**(exp(-0.5*((1.-wsk/wsj)/sigma)**2))*(wl/wv)**2
        ELSE !}  {
          alpha=0.0
        end if !}
        e(k,j,iac)=alpha

      end do !}
    end do !}
  end do !}
end SUBROUTINE setspec !}
#endif
