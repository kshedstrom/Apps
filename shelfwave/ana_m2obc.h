      SUBROUTINE ana_m2obc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets 2D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m2obc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     knew(ng),                                    &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     OCEAN(ng) % zeta)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(12)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m2obc
!
!***********************************************************************
      SUBROUTINE ana_m2obc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
     &                           angler, h, pm, pn, on_u,               &
#ifdef MASKING
     &                           umask,                                 &
#endif
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: kk, ll, alpha, omega
      real(r8) :: cos_wt, cos_ky, sin_wt, sin_ky
      real(r8), parameter :: Lx = 100.e+3
      real(r8), parameter :: Ly = 50.e+3
      real(r8), parameter :: my_amp = 1.e+3
      integer,  parameter :: jj = 2    ! cross-shore mode number

#include "set_bounds.h"

      alpha = 1. / Ly
      ll = 2. * pi / Lx
      kk = jj * pi / el(ng)
      omega = 2 * alpha * GRID(ng)%f(Istr,Jstr) * ll /                  &
     &           (kk*kk + alpha*alpha + ll*ll)
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
      IF (LBC(iwest,isUbar,ng)%acquire.and.                             &
     &    LBC(iwest,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          sin_wt = sin(ll*GRID(ng)%xu(1,j) - omega*time(ng))
          cos_wt = cos(ll*GRID(ng)%xu(1,j) - omega*time(ng))
          sin_ky = sin(kk * GRID(ng)%yu(1,j))
          cos_ky = cos(kk * GRID(ng)%yu(1,j))
          BOUNDARY(ng)%ubar_west(j) = exp(- alpha * GRID(ng)%yu(1,j))   &
     &           * my_amp * cos_wt * (alpha * sin_ky + kk * cos_ky)
        END DO

        DO j=JstrP,JendT
!         sin_wt = sin(ll*GRID(ng)%xv(0,j) - omega*time(ng))
!         cos_wt = cos(ll*GRID(ng)%xv(0,j) - omega*time(ng))
!         sin_ky = sin(kk * GRID(ng)%yv(0,j))
!         cos_ky = cos(kk * GRID(ng)%yv(0,j))
!         BOUNDARY(ng)%vbar_west(j) = ll * exp(- alpha *                &
!    &              my_amp * GRID(ng)%yv(0,j)) * sin_wt * sin_ky
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (LBC(isouth,isUbar,ng)%acquire.and.                            &
     &    LBC(isouth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_south(i)=0.0_r8
        END DO
      END IF
      IF (LBC(inorth,isUbar,ng)%acquire.and.                            &
     &    LBC(inorth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_north(i)=0.0_r8
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_m2obc_tile
