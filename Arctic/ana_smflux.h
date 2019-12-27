      SUBROUTINE ana_smflux (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface momentum flux (wind stress)     !
!  "sustr" and "svstr" (m2/s2) using an analytical expression.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_smflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % angler,                          &
#ifdef SPHERICAL
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
#else
     &                      GRID(ng) % xr,                              &
     &                      GRID(ng) % yr,                              &
#endif
#ifdef MASKING
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# ifdef MASK_HACK
     &                      GRID(ng) % mask2,                           &
# endif
#endif
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(24)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_smflux
!
!***********************************************************************
      SUBROUTINE ana_smflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            angler,                               &
#ifdef SPHERICAL
     &                            lonr, latr,                           &
#else
     &                            xr, yr,                               &
#endif
#ifdef MASKING
     &                            umask, vmask,                         &
# ifdef MASK_HACK
     &                            mask2,                                &
# endif
#endif
     &                            sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  ifdef MASK_HACK
      real(r8), intent(in) :: mask2(LBi:,LBj:)
#  endif
# endif
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  ifdef MASK_HACK
      real(r8), intent(in) :: mask2(LBi:UBi,LBj:UBj)
#  endif
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j, i0, j0
      real(r8) :: wind_min, wind_max, wind_amp, wind_dir, my_tdays

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  XI-direction (m2/s2) at horizontal U-points.
!-----------------------------------------------------------------------
!
      wind_min = 0.0_r8
      wind_max = 1.4e-4_r8
!     wind_max = 1.4e-4_r8 / sqrt(2.0)
      my_tdays = tdays(ng) - 41706.5
      if (my_tdays < 5.0) then
        wind_amp = wind_min
      else if (my_tdays < 7.0) then
        wind_amp = wind_min + 0.5*(my_tdays - 5.0)*(wind_max-wind_min)
      else if (my_tdays < 10.0) then
        wind_amp = wind_max
      else if (my_tdays < 12.0) then
        wind_amp = wind_min + 0.5*(12.0 - my_tdays)*(wind_max-wind_min)
      else
        wind_amp = wind_min
      endif
! This is now in degrees clockwise from north
      wind_dir = -45*pi/180._r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          sustr(i,j) = wind_amp*sin(wind_dir + angler(i,j))
!         sustr(i,j) = -wind_amp
#ifdef MASKING
          sustr(i,j) = sustr(i,j)*umask(i,j)
# ifdef MASK_HACK
! HACK for idealized Arctics
          i0 = max(i-1,IstrR)
          sustr(i,j) = sustr(i,j) * max(mask2(i,j),mask2(i0,j))
# endif
#endif
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          svstr(i,j) = wind_amp*cos(wind_dir + angler(i,j))
!         svstr(i,j) = wind_amp
#ifdef MASKING
          svstr(i,j) = svstr(i,j)*vmask(i,j)
# ifdef MASK_HACK
! HACK for idealized Arctics
          j0 = max(j-1,JstrR)
          svstr(i,j) = svstr(i,j) * max(mask2(i,j),mask2(i,j0))
# endif
#endif
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
#endif

      RETURN
      END SUBROUTINE ana_smflux_tile
