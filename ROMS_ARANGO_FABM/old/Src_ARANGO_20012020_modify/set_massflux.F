#include "cppdefs.h"
      MODULE set_massflux_mod

#ifdef SOLVE3D
!
!svn $Id: set_massflux.F 995 2020-01-10 04:01:28Z arango $
!=======================================================================
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_massflux
!!!Inserted PWA 09/08/2025
# ifdef OFFLINE
      PUBLIC  :: set_massflux_avg
# endif
!!!End insert PWA 09/08/2025
# ifdef ADJOINT
      PUBLIC  :: reset_massflux
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_massflux (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 12, __LINE__, __FILE__)
# endif
      CALL set_massflux_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng),                                 &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
# ifdef NEARSHORE_MELLOR
     &                        OCEAN(ng) % u_stokes,                     &
     &                        OCEAN(ng) % v_stokes,                     &
# endif
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % om_v,                          &
     &                        GRID(ng) % on_u,                          &
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, model, 12, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE set_massflux
!!!Inserted PWA 09/08/2025

# ifdef OFFLINE
!
!***********************************************************************
      SUBROUTINE set_massflux_avg (ng, tile, model)
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  KMT Calculate mid-time step values for Huon, Hvom, and Akt (both indices)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
# if defined AKTCLIMATOLOGY || defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
      USE mod_mixing
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 12, __LINE__, __FILE__)
# endif
      CALL set_massflux_avg_tile (ng, tile, model,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng), nnew(ng),                       &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
# ifdef WEC
     &                        OCEAN(ng) % u_stokes,                     &
     &                        OCEAN(ng) % v_stokes,                     &
# endif
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % om_v,                          &
     &                        GRID(ng) % on_u,                          &
# if defined AKTCLIMATOLOGY || defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
     &                        MIXING(ng) % Akt, MIXING(ng) % gls,       & 
# endif
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, model, 12, __LINE__, __FILE__)
# endif

      RETURN
      END SUBROUTINE set_massflux_avg


!
!***********************************************************************
      SUBROUTINE set_massflux_avg_tile (ng, tile, model,                &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs, nnew,                         &
     &                              u, v,                               &
# ifdef WEC
     &                              u_stokes, v_stokes,                 &
# endif
     &                              Hz,om_v, on_u,                      &
# if defined AKTCLIMATOLOGY || defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
     &                              Akt, gls,                           &
# endif
     &                              Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#  ifdef WEC
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: om_v(LBi:,LBj:)
      real(r8), intent(inout) :: on_u(LBi:,LBj:)

      real(r8), intent(inout) :: Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: Hvom(LBi:,LBj:,:)
#  if defined AKTCLIMATOLOGY || defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
      real(r8), intent(out) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
#  endif
# else
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef WEC
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  if defined AKTCLIMATOLOGY || defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
      real(r8), intent(out) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:N(ng),3)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: Huonnew(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: Hvomnew(LBi:UBi,LBj:UBj,N(ng))

# include "set_bounds.h"

!-----------------------------------------------------------------------
!  KMT Calculate new time step of mass flux value using u, v corresponding to
!  the updated values for next time step nnew.
!-----------------------------------------------------------------------
      CALL set_massflux_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nnew,                                     &
     &                        u, v,                                     &
# ifdef WEC
     &                        u_stokes, v_stokes,                       &
# endif
     &                        Hz, om_v, on_u,                           &
     &                        Huonnew, Hvomnew)


!
!-----------------------------------------------------------------------
!  Compute average of  horizontal mass fluxes, Hz*u/n and Hz*v/m.
!  from previous and subsequent time steps. This is only for offline 
!  case.
!  Calculate average from previous and next time step values of Huon/Hvom
!  Previous time step held in GRID(ng)%Huon, next time step held in 
!  Huonnew. Newly calculated average will be put in GRID(ng)%Huon to be
!  used in other functions.
!-----------------------------------------------------------------------
!
!  Compute horizontal mass fluxes.
!
      DO k=1,N(ng)
        DO j=JstrT-1,JendT+1
          DO i=IstrP-1,IendT+1
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+Huonnew(i,j,k))
          END DO
        END DO
        DO j=JstrP-1,JendT+1
          DO i=IstrT-1,IendT+1
            Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+Hvomnew(i,j,k))
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!  Update Akt (both indices) to next time step, which is stored in GLS 
!-----------------------------------------------------------------------
      DO k=0,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
# if defined AKTCLIMATOLOGY || defined AKXCLIMATOLOGY
!  Update Akt to next time step, which is stored in GLS index nnew
!  Only deal with index 1 of Akt array since that corresponds to 
!  Akt instead of Aks.
!  Akt matters because it can be used with floats in VWALK 
             Akt(i,j,k,1)=gls(i,j,k,nnew)
# endif
# if defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
!  Update Aks to next time step, which is stored in GLS index 3
!  This is used in tracer advection in step3d_t
            Akt(i,j,k,2)=gls(i,j,k,3)
# endif
          END DO
        END DO
      END DO

!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Huon)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hvom)
! KMT I don't think these parts are used since the 
! boundary conditions are closed for offline simulations.
# if defined AKTCLIMATOLOGY || defined AKXCLIMATOLOGY
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akt(:,:,:,1))
# endif
# if defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akt(:,:,:,2))
# endif
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Huon, Hvom)
#  if defined AKTCLIMATOLOGY || defined AKSCLIMATOLOGY || defined AKXCLIMATOLOGY
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
#  endif
# endif

      RETURN
      END SUBROUTINE set_massflux_avg_tile

# endif  /* OFFLINE */
!!!End insert PWA 09/08/2025
!
!***********************************************************************
      SUBROUTINE set_massflux_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs,                               &
     &                              u, v,                               &
# ifdef NEARSHORE_MELLOR
     &                              u_stokes, v_stokes,                 &
# endif
     &                              Hz, om_v, on_u,                     &
     &                              Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(out) :: Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: Hvom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
!  Compute horizontal mass fluxes.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            Huon(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   &
     &                  on_u(i,j)
# ifdef NEARSHORE_MELLOR
            Huon(i,j,k)=Huon(i,j,k)+                                    &
     &                  0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                  u_stokes(i,j,k)*on_u(i,j)
# endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            Hvom(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   &
     &                  om_v(i,j)
# ifdef NEARSHORE_MELLOR
            Hvom(i,j,k)=Hvom(i,j,k)+                                    &
     &                  0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                 &
     &                  v_stokes(i,j,k)*om_v(i,j)
# endif
          END DO
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Huon)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hvom)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Huon, Hvom)
# endif

      RETURN
      END SUBROUTINE set_massflux_tile

# ifdef ADJOINT
!
!***********************************************************************
      SUBROUTINE reset_massflux (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, model, 12, __LINE__, __FILE__)
#  endif
      CALL reset_massflux_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nnew(ng),                               &
     &                          COUPLING(ng) % DU_avg2,                 &
     &                          COUPLING(ng) % DV_avg2,                 &
     &                          OCEAN(ng) % u,                          &
     &                          OCEAN(ng) % v,                          &
#  ifdef NEARSHORE_MELLOR
     &                          OCEAN(ng) % u_stokes,                   &
     &                          OCEAN(ng) % v_stokes,                   &
#  endif
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % on_u,                        &
     &                          GRID(ng) % Huon,                        &
     &                          GRID(ng) % Hvom)
#  ifdef PROFILE
      CALL wclock_off (ng, model, 12, __LINE__, __FILE__)
#  endif

      RETURN
      END SUBROUTINE reset_massflux
!
!***********************************************************************
      SUBROUTINE reset_massflux_tile (ng, tile, model,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                nnew,                             &
     &                                DU_avg2, DV_avg2,                 &
     &                                u, v,                             &
#  ifdef NEARSHORE_MELLOR
     &                                u_stokes, v_stokes,               &
#  endif
     &                                Hz, om_v, on_u,                   &
     &                                Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      integer, intent(in) :: nnew
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(inout) :: Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: Hvom(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(Ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute intermediate values of mass fluxes Huon and Hvom used by the
!  adjoint model.  The original values can be reinstated by calling
!  "set_massflux" after "ad_omega".
!-----------------------------------------------------------------------
!
!  Compute mass flux, Hz*u/n.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          DC(i,0)=0.0_r8
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
          END DO
        END DO
        DO k=N(ng),1,-1
          DO i=IstrP,IendT
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
#  ifdef NEARSHORE_MELLOR
            Huon(i,j,k)=Huon(i,j,k)+0.5_r8*u_stokes(i,j,k)*DC(i,k)
#  endif
            FC(i,0)=FC(i,0)+Huon(i,j,k)
          END DO
        END DO
!
!  Replace with correct vertical mean, DU_avg2.
!
        DO i=IstrP,IendT
          DC(i,0)=1.0_r8/DC(i,0)
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Compute mass flux, Hz*v/m.
!
        IF (j.ge.JstrP) THEN
          DO i=IstrT,IendT
            DC(i,0)=0.0_r8
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
            END DO
          END DO
          DO k=N(ng),1,-1
            DO i=IstrT,IendT
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
#  ifdef NEARSHORE_MELLOR
              Hvom(i,j,k)=Hvom(i,j,k)+0.5_r8*v_stokes(i,j,k)*DC(i,k)
#  endif
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
            END DO
          END DO
!
!  Replace with correct vertical mean, DV_avg2.
!
          DO i=IstrT,IendT
            DC(i,0)=1.0_r8/DC(i,0)
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        ENDIF
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Huon)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hvom)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Huon, Hvom)
#  endif

      RETURN
      END SUBROUTINE reset_massflux_tile
# endif
#endif
      END MODULE set_massflux_mod
