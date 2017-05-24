      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_fsobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 6)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
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
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF
      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          sin_wt = sin(ll*GRID(ng)%xr(0,j) - omega*time(ng))
          cos_wt = cos(ll*GRID(ng)%xr(0,j) - omega*time(ng))
          sin_ky = sin(kk * GRID(ng)%yr(0,j))
          cos_ky = cos(kk * GRID(ng)%yr(0,j))
          BOUNDARY(ng)%zeta_west(j) = exp(- alpha * GRID(ng)%yr(0,j))   &
     &            * my_amp * cos_wt * alpha / g * sin_ky
        END DO
      END IF
      IF (LBC(isouth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF
      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_fsobc_tile
