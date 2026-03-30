#include "cppdefs.h"
      MODULE set_zeta_mod

#ifdef SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine sets free-surface to its fast-time averaged value.     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
!!!Modified PWA 09/08/2025
#if defined SOLVE3D && defined OFFLINE
      PUBLIC  :: set_zeta, set_avg_zeta
#else
      PUBLIC  :: set_zeta
#endif
!!!End modified PWA 09/08/2025
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_zeta (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 12, __LINE__, MyFile)
# endif
      CALL set_zeta_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    COUPLING(ng) % Zt_avg1,                       &
     &                    OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE set_zeta
!
!***********************************************************************
      SUBROUTINE set_zeta_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Zt_avg1, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)

      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Prepare to time-step 2D equations:  set initial free-surface
!  to its fast-time averaged values (which corresponds to the time
!  step "n").
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,1)=Zt_avg1(i,j)
          zeta(i,j,2)=Zt_avg1(i,j)
        END DO
      END DO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,1))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,2))
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta(:,:,1),                                  &
     &                    zeta(:,:,2))
# endif
!
      RETURN
      END SUBROUTINE set_zeta_tile
#endif

!!!Inserted PWA 09/08/2025
#if defined SOLVE3D && defined OFFLINE
!
!***********************************************************************
      SUBROUTINE set_avg_zeta (ng, tile, FIRSTCALL)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      logical, intent(inout) :: FIRSTCALL
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 12, __LINE__, __FILE__)
# endif
      CALL set_avg_zeta_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    FIRSTCALL,                                    &
     &                    COUPLING(ng) % Zt_avg1,                       &
     &                    OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE set_avg_zeta
!
!***********************************************************************
      SUBROUTINE set_avg_zeta_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          FIRSTCALL,                              &
     &                          Zt_avg1,                                &
     &                          zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
      logical, intent(inout) :: FIRSTCALL
!
!  Local variable declarations.
!
      real(r8) :: cff, cff1
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  KMT added
!  For offline case, need to update Zt_avg1 so that set_depth 
!  case be used to update depths.
!-----------------------------------------------------------------------
!
      IF (FIRSTCALL) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,1)  ! or index 2 or 3
          END DO
        END DO
        FIRSTCALL = .FALSE.
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,2)  ! or index 2 or 3
          END DO
        END DO
      END IF

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Zt_avg1(:,:))
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Zt_avg1(:,:))
# endif

      RETURN
      END SUBROUTINE set_avg_zeta_tile

#endif  /* SOLVE3D and OFFLINE */

!!!End insert PWA 09/08/2025
      END MODULE set_zeta_mod
