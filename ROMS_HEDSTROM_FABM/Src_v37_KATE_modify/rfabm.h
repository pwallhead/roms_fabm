      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group         Phil Wallhead   !
!                                                   Andre Staalstrom   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!  This routine computes the  biological sources and sinks for the     !
!  any model in the FABM family. Then, it adds those terms to the      !
!  global biological fields.                                           !
!  The code below is adapted from npzd_Franks.h + diagnostics code     !
!  from fennel.h                                                       !
!***********************************************************************

      USE mod_param     !PWA: Provides access to BOUNDS and FMODEL
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean     !PWA: Provides access to OCEAN (hence "t")
      USE mod_stepping
!!!Inserted PWA 08/03/2017
      USE mod_forces
      USE fabm
      USE fabm_config
!!!End insertion PWA 08/03/2017
!
      implicit none !Inserted PWA 30/07/2019 (missing from npzd_Franks.h, but present in fennel.h)
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask_full,                         &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
!!!PWA Inserted 08/03/2017
     &                   OCEAN(ng) % state1,                            &
     &                   OCEAN(ng) % state_sf,                          &
     &                   OCEAN(ng) % state_bt,                          &
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
#endif
!!!End insertion PWA 08/03/2017
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask_full,                              &
#endif
     &                         Hz, z_r, z_w,                            &
     &                         state1, state_sf, state_bt,              &
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                      &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology !PWA: This provides access to idbio, idfabmD2, idfabmD3
      USE mod_ncparam
      USE mod_scalars
!!!PWA Inserted 08/03/2017
      USE fabm
      USE fabm_config
!!!End insertion PWA 08/03/2017
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
!!!PWA Inserted 08/03/2017
      real(r8), intent(inout) :: state1(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: state_sf(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: state_bt(LBi:,LBj:,:,:)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
!!!End insertion PWA 08/03/2017
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
!!!PWA Inserted 08/03/2017
      real(r8), intent(inout) :: state1(LBi:UBi,LBj:UBj,UBk,NBT)
      real(r8), intent(inout) :: state_sf(LBi:UBi,LBj:UBj,3,NSAT)
      real(r8), intent(inout) :: state_bt(LBi:UBi,LBj:UBj,3,NBAT)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
!!!End insertion PWA 08/03/2017
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!!!PWA: The subroutine is passed chunks defined by tile boundaries including ghost points (LBi, UBi, etc.)
!!!     although it only modifies ranges Istr:Iend etc. which exclude ghost points.

!
!  Local variable declarations.
!
!!!PWA Modified 08/03/2017
!      integer, parameter :: Nsink = 1
!      ...
!      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      integer :: ibio, itrc, ivar, i, j, k, ks, nx, ny
      integer :: i1, j1 !These indices range strictly over (1,nx), (1,ny), while (i,j) range over (Istr,Iend), (Jstr,Jend)

!!!PWA: Rather than prescribe here the number of sinking variables,
!!!easier to just loop over NBT and take action for i where max(w(:,k,i))>0

! local variables
      real(r8)  yday, hour, dBdt1 !, dBdt1max
      integer :: iday, month, year
      real(r8) :: cff, cffL, cffR, cu, dltL, dltR

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: w      !ROMS vertical sinking velocity [m/s] positive downward
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      real(r8), dimension(IminS:ImaxS,N(ng),NBT) :: state_old !Old state (accounting for clipping)
      real(r8), dimension(IminS:ImaxS,N(ng),NBT) :: dBdt      !FABM SMS term (units [C/s])
      real(r8), dimension(IminS:ImaxS,NBT) :: flux_sf         !FABM surface fluxes (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,NSAT) :: sms_sf         !FABM surface-attached SMS terms (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,NBT) :: flux_bt         !FABM bottom fluxes (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,NBAT) :: sms_bt         !FABM bottom-attached SMS terms (units [Matter/m2/s])
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      logical,parameter :: repair = .TRUE.
      logical :: valid_int,valid_sf,valid_bt
#ifdef DIAGNOSTICS_BIO
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DiaBio2d1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: DiaBio3d1
#endif
#ifdef FABM_NONNEG_S
      real(r8), target, dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Snn    !Non-negative salinity for input to FABM
#endif
!!!End modification PWA 08/03/2017

#include "set_bounds.h" 
!PWA: Sets Istr,Iend etc. using BOUNDS (note that LBi,UBi etc. are passed through subroutine args)

! Chunksize variables for calls to FABM APIs
      nx = Iend-Istr+1
      ny = Jend-Jstr+1


!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
#ifdef DIAGNOSTICS_BIO
      IF (nDIA(ng).gt.0.and.                                            &
     &   (((iic(ng).gt.ntsDIA(ng)).and.(MOD(iic(ng),nDIA(ng)).eq.1)).or.&
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng))))) THEN
        DO ivar=1,NDbio2d
          IF (Dout(iDbio2(ivar),ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio2d(i,j,ivar)=0.0_r8
              END DO
            END DO
          END IF
        END DO
        DO ivar=1,NDbio3d
          IF (Dout(iDbio3(ivar),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  DiaBio3d(i,j,k,ivar)=0.0_r8
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
#endif


!-----------------------------------------------------------------------
!  If wetting and drying, re-point FABM to full mask on first time step
!-----------------------------------------------------------------------
# ifdef WET_DRY
      IF (iic(ng).eq.1) THEN
        CALL FMODEL(ng)%f(tile)%model%set_mask(                         &
     &                                rmask_full(Istr:Iend,Jstr:Jend))
      END IF
# endif


!-----------------------------------------------------------------------
! Relink environmental data as required (due to variable nstp).
!-----------------------------------------------------------------------
!Link environmental data to FABM (for variables that must be relinked for each time step since nstp varies)
!Note: We cannot use variable_needs_values to parse these calls, since these will return
!      .false. after the first linking calls in roms_fabm.F.
      IF (FMODEL(ng)%needsT(tile)) THEN
        CALL FMODEL(ng)%f(tile)%model%link_interior_data(               &
     &                    FMODEL(ng)%id_temp(tile),                     &
     &                    t(Istr:Iend,Jstr:Jend,1:N(ng),nstp,1))
      END IF
      IF (FMODEL(ng)%needsS(tile)) THEN
#ifdef FABM_NONNEG_S
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Snn(i,j,k) = MAX(0.0_r8, t(i,j,k,nstp,2))
            END DO
          END DO
        END DO
        CALL FMODEL(ng)%f(tile)%model%link_interior_data(               &
     &                    FMODEL(ng)%id_salt(tile),                     &
     &                    Snn(Istr:Iend,Jstr:Jend,1:N(ng)))
#else      
        CALL FMODEL(ng)%f(tile)%model%link_interior_data(               &
     &                    FMODEL(ng)%id_salt(tile),                     &
     &                    t(Istr:Iend,Jstr:Jend,1:N(ng),nstp,2))
#endif
      END IF



!-----------------------------------------------------------------------
! Relink surface/bottom state data as required (due to variable nstp).
!-----------------------------------------------------------------------
      DO itrc=1,NSAT
        CALL FMODEL(ng)%f(tile)%model%link_surface_state_data(          &
     &   itrc,state_sf(Istr:Iend,Jstr:Jend,nstp,itrc))
      END DO

      DO itrc=1,NBAT
        CALL FMODEL(ng)%f(tile)%model%link_bottom_state_data(           &
     &   itrc,state_bt(Istr:Iend,Jstr:Jend,nstp,itrc))
      END DO


!-----------------------------------------------------------------------------------
! Copy interior data from "t" to "state", capping if required using FABM check_state
!-----------------------------------------------------------------------------------
!  
!  Note: The combined FABM sms + ROMS sinking increment will be calculated as (state - state_old).
!        state_old is assigned within the J loop to save memory.
!  Note: If capping is applied here using FABM (cpp FABM_CHECK_STATE) then this capping
!        will be applied to BOTH "state_old" and "state".  This is consistent with the
!        capping applied to Bio and Bio_old in larger ROMS models (e.g. fennel.h, ecosim.h). 
!        However, in smaller ROMS biology modules e.g. npzd_Franks.h, nonnegativity is imposed 
!        more loosely as a correction that conserves total mass, and is not applied to Bio_old.
!        To reproduce the latter behaviour, FABM_CHECK_STATE should be deactivated (or parameter
!        "repair" set to "false") and the mass-conserving correction is done within the FABM module 
!        (e.g. niva_roms_npzd_Franks.F90).
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                state1(i,j,k,itrc)=t(i,j,k,nstp,ibio) !Note: FABM is permanently linked to state1 (in roms_fabm.F)
              END DO
            END DO
          END DO
        END DO
#ifdef FABM_CHECK_STATE
        DO j1=1,ny
          DO k=1,N(ng)
            CALL FMODEL(ng)%f(tile)%model%check_interior_state(1,nx,j1, &
     &                          k,repair,valid_int)   !This will cap state1(Istr:Iend,Jstr+j1,k,1:NBT)
          END DO
        END DO
        IF (NSAT.gt.0) THEN
          DO j1=1,ny
            CALL FMODEL(ng)%f(tile)%model%check_surface_state(1,nx,j1,  &
     &                          repair,valid_sf)      !This will cap state_sf(Istr:Iend,Jstr+j1,nstp,1:NSAT)
          END DO
        END IF
        IF (NBAT.gt.0) THEN
          DO j1=1,ny
            CALL FMODEL(ng)%f(tile)%model%check_bottom_state(1,nx,j1,   &
     &                          repair,valid_bt)      !This will cap state_bt(Istr:Iend,Jstr+j1,nstp,1:NBAT)
          END DO
        END IF
#endif


!-----------------------------------------------------------------------
!  Set date variables
!-----------------------------------------------------------------------
      CALL caldate(r_date, tdays(ng), year, yday, month, iday, hour)
      ydayc(ng) = yday-1.0_r8
      IF (iic(ng).le.icheckmax(ng)) THEN
        write(*,*) "tdays = ", tdays(ng)
        write(*,*) "yday, ydayc = ", yday, ydayc(ng)
      END IF
!     dBdt1max = 1e3 !Maximum absolute rate of change from FABM model (def = 1e3)
                     !Set to zero to examine initial FABM output
!     icheckmax = 0  !Maximum time step counter value (iig) for which MIN/MAX dBdt is checked
                     !Set to zero to switch off all checks
                     !Set to large value to check all time steps


!-----------------------------------------------------------------------
!  Prepare all fields FABM needs to compute SMS terms (e.g., light)
!-----------------------------------------------------------------------
      CALL FMODEL(ng)%f(tile)%model%prepare_inputs()




!-----------------------------------------------------------------------
!  Enter J loop to compute SMS terms
!-----------------------------------------------------------------------
      J_LOOP : DO j=Jstr,Jend
        j1 = j-Jstr+1 !FABM chunks are always indexed 1:chunksize
                      !Therefore, all calls to FABM APIs should use index j1 instead of j

!-----------------------------------------------------------------------
!  Compute inverse thickness to avoid repeated divisions.
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO


!-----------------------------------------------------------------------
!  Store the old interior state.
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO k=1,N(ng)
            DO i=Istr,Iend
              state_old(i,k,itrc)=state1(i,j,k,itrc)
            END DO
          END DO
        END DO


!-----------------------------------------------------------------------
!  Output grid parameters for this tile if req'd
!-----------------------------------------------------------------------
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Istr,Iend,Jstr,Jend = ",Istr,Iend,Jstr,Jend
          write(*,*) "iic,icheckmax,j,j1 = ",iic(ng),icheckmax(ng),j,j1
#ifdef MASKING
          write(*,*) "MIN,MAX(rmask_full(Istr:Iend,j)) = ",             &
     &  MINVAL(rmask_full(Istr:Iend,j)), MAXVAL(rmask_full(Istr:Iend,j))
#endif
        END IF


!-----------------------------------------------------------------------
!  Generate surface fluxes and SMS for surface-attached variables
!-----------------------------------------------------------------------
        flux_sf = 0.0_r8
        sms_sf = 0.0_r8
        !Note: surface fluxes may still be needed even if NSAT = 0
        CALL FMODEL(ng)%f(tile)%model%get_surface_sources(1,nx,j1,      &
     &    flux_sf(Istr:Iend,1:NBT),sms_sf(Istr:Iend,1:NSAT))
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do_surface"
          write(*,*) "MIN,MAX(flux_sf(Istr:Iend,:)) = ",                &
     &      MINVAL(flux_sf(Istr:Iend,:)), MAXVAL(flux_sf(Istr:Iend,:))
          IF (NSAT.gt.0) THEN
            write(*,*) "MIN,MAX(sms_sf(Istr:Iend,:)) = ",               &
     &        MINVAL(sms_sf(Istr:Iend,:)), MAXVAL(sms_sf(Istr:Iend,:))
          END IF
        END IF


!-----------------------------------------------------------------------
!  Generate bottom fluxes and SMS terms for bottom-attached variables
!-----------------------------------------------------------------------
        flux_bt = 0.0_r8
        sms_bt = 0.0_r8
        !Note: bottom fluxes may still be needed even if NBAT = 0
        CALL FMODEL(ng)%f(tile)%model%get_bottom_sources(1,nx,j1,       &
     &    flux_bt(Istr:Iend,1:NBT),sms_bt(Istr:Iend,1:NBAT))
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do_bottom"
          write(*,*) "MIN,MAX(flux_bt(Istr:Iend,:)) = ",                &
     &      MINVAL(flux_bt(Istr:Iend,:)), MAXVAL(flux_bt(Istr:Iend,:))
          IF (NBAT.gt.0) THEN
            write(*,*) "MIN,MAX(sms_bt(Istr:Iend,:)) = ",               &
     &      MINVAL(sms_bt(Istr:Iend,:)), MAXVAL(sms_bt(Istr:Iend,:))
          END IF
        END IF


!-----------------------------------------------------------------------
!  Generate SMS terms for pelagic variables
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          dBdt(Istr:Iend,k,1:NBT)=0.0_r8
          CALL FMODEL(ng)%f(tile)%model%get_interior_sources(1,nx,j1,k, &
     &              dBdt(Istr:Iend,k,1:NBT))
!Note: the biological tracer indices idbio range over (NAT+NPT+NCS+NNS)+1:NT = NT-NBT+1:NT, see rfabm_mod.h
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do"
          write(*,*) "MIN,MAX(dBdt(Istr:Iend,1:N(ng),1:NBT)) = ",       &
     &      MINVAL(dBdt(Istr:Iend,1:N(ng),1:NBT)),                      &
     &      MAXVAL(dBdt(Istr:Iend,1:N(ng),1:NBT))
        END IF


!-----------------------------------------------------------------------
!  Add contributions from surface and bottom fluxes
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO i=Istr,Iend
            dBdt(i,N(ng),itrc) = dBdt(i,N(ng),itrc) +                   &
     &            flux_sf(i,itrc)*Hz_inv(i,N(ng))
            dBdt(i,1,itrc) = dBdt(i,1,itrc) +                           &
     &            flux_bt(i,itrc)*Hz_inv(i,1)
          !E.g surface/bottom attached variables may have [mass/m2] while pelagic variables may have [mass/m3]
          END DO
        END DO


!-----------------------------------------------------------------------
!  Check the model state and FABM computations in detail if req'd
!-----------------------------------------------------------------------
        IF (iic(ng).le.icheckmax(ng)) THEN

          dBdt1 = MAXVAL(ABS(dBdt(Istr:Iend,:,1:NBT)))

          IF (dBdt1.ge.dBdt1max(ng)) THEN
            write(*,*) "dBdt1 = MAX(ABS(dBdt)) = ", dBdt1
#ifdef DIAGNOSTICS_BIO
            DiaBio3d1(Istr:Iend,Jstr:Jend,1:N(ng)) = FMODEL(ng)%        &
     &       f(tile)%model%get_interior_diagnostic_data(53)
            write(*,*) "MAXVAL(light_PAR0(Istr:Iend,j,N(ng))) = ",      &
     &        MAXVAL(DiaBio3d1(Istr:Iend,j,N(ng)))
            stop
#endif
            write(*,*) "MIN,MAX(Hz) = ", MINVAL(Hz), MAXVAL(Hz)
          END IF
        END IF


!-----------------------------------------------------------------------
!  Update tracers with rates of change from FABM (Euler step)
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO k=1,N(ng)
            DO i=Istr,Iend
              state1(i,j,k,itrc) = state1(i,j,k,itrc) +                 &
     &         dBdt(i,k,itrc)*dt(ng)
            END DO
          END DO
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Updated water column state with FABM dBdt"
          write(*,*) "MIN,MAX(state1(Istr:Iend,j,1:N(ng),1:NBT))=",     &
            MINVAL(state1(Istr:Iend,j,1:N(ng),1:NBT)),                  &
     &      MAXVAL(state1(Istr:Iend,j,1:N(ng),1:NBT))
        END IF


!-----------------------------------------------------------------------
!  Update surface states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NSAT
          DO i=Istr,Iend
            state_sf(i,j,nnew,itrc) = MAX(state_sf(i,j,nstp,itrc) +     &
     &                                 sms_sf(i,itrc)*dt(ng), 0.0_r8)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
          END DO
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          IF (NSAT.gt.0) THEN
            write(*,*) "Updated surface states with FABM sms_sf"
            write(*,*) "MIN,MAX(state_sf(Istr:Iend,j,nnew,1:NSAT))=",   &
              MINVAL(state_sf(Istr:Iend,j,nnew,1:NSAT)),                &
     &        MAXVAL(state_sf(Istr:Iend,j,nnew,1:NSAT))
          END IF
        END IF


!-----------------------------------------------------------------------
!  Update bottom states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NBAT
          DO i=Istr,Iend
            state_bt(i,j,nnew,itrc) = MAX(state_bt(i,j,nstp,itrc) +     &
     &                                 sms_bt(i,itrc)*dt(ng), 0.0_r8)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
          END DO
        END DO
        IF (iic(ng).le.icheckmax(ng)) THEN
          IF (NBAT.gt.0) THEN
            write(*,*) "Updated bottom states with FABM sms_bt"
            write(*,*) "MIN,MAX(state_bf(Istr:Iend,j,nnew,1:NBAT))=",   &
              MINVAL(state_bt(Istr:Iend,j,nnew,1:NBAT)),                &
     &        MAXVAL(state_bt(Istr:Iend,j,nnew,1:NBAT))
          END IF
        END IF


!-----------------------------------------------------------------------
!  Get vertical sinking velocities from FABM
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          CALL FMODEL(ng)%f(tile)%model%get_vertical_movement(          &
     &              1,nx,j1,k,w(Istr:Iend,k,NT(ng)-NBT+1:NT(ng)))
        END DO
        w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)) = -1*                  &
     &      w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))
        !FABM outputs vertical velocities in [m/s] positive upward
        !ROMS code below expects w in [m/s] positive downward
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm get_vertical_movement"
          write(*,*) "MIN,MAX(w(Istr:Iend,1:N(ng),NT-NBT+1:NT)) = ",    &
            MINVAL(w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))),           &
     &      MAXVAL(w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)))
        END IF


!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
!!!Modified PWA 08/03/2017
!        SINK_LOOP: DO isink=1,Nsink !!!PWA: This is not convenient for arbitrary fabm model
!                                    !!!PWA: Instead we loop over all biol tracers and use IF statement (leaves indentation unchanged since no BioIter loop)
!        ibio=idsink(isink)
        SINK_LOOP: DO itrc=1,NBT
          ibio=idbio(itrc)
          IF (MAXVAL(ABS(w(Istr:Iend,1:N(ng),ibio))).gt.1.0E-12_r8) THEN
!!!Note: 1.0E-12 m/s => 0.3 mm/decade (safely negligible)
!!!End modification PWA 08/03/2017
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=state1(i,j,k,itrc) !PWA: Note this uses the bgc-updated state (hence operator splitting), consistent with ROMS biology modules
              END DO
            END DO
  !
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
!!!Modified PWA 08/03/2017
!            cff=dtdays*ABS(Wbio(isink))
!!!End modification PWA 08/03/2017
            DO k=1,N(ng)
              DO i=Istr,Iend
!!!Inserted PWA 08/03/2017
                cff=dt(ng)*w(i,k,ibio) !Modified 12/11/2020 to allow negative w => buoyant particles
!!!End insertion PWA 08/03/2017
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                state1(i,j,k,itrc)=qc(i,k)+                             &
     &           (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
!!!Modified PWA 08/03/2017
!          END DO SINK_LOOP
!        END DO ITER_LOOP
          END IF
        END DO SINK_LOOP
!!!End modification PWA 08/03/2017
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "state_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "state"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of nutrients
!  when advection causes tracer concentration to go negative.
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff = state1(i,j,k,itrc) - state_old(i,k,itrc)
#ifdef MASKING
              cff = cff*rmask_full(i,j)
!It seems safer to ensure that no cells are not updated when dry.
!(rmask_full = rmask_wet*rmask, where rmask is the imposed time-independent mask (0=dry, 1=wet))
#endif
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP


!-----------------------------------------------------------------------
!  Finalize outputs
!-----------------------------------------------------------------------
      CALL FMODEL(ng)%f(tile)%model%finalize_outputs()


#ifdef DIAGNOSTICS_BIO
      IF (nDIA(ng).gt.0) THEN
!-----------------------------------------------------------------------
!  Extract diagnostics from FABM
!-----------------------------------------------------------------------
        DO ivar=1,NDbio2d
          IF (Dout(iDbio2(ivar),ng)) THEN
            i=idfabmD2(ivar)
            DiaBio2d1(Istr:Iend,Jstr:Jend) = FMODEL(ng)%f(tile)%model%  &
     &       get_horizontal_diagnostic_data(i)
            DO j=Jstr,Jend
              DO i=Istr,Iend
# ifdef MASKING
                DiaBio2d(i,j,ivar) = DiaBio2d(i,j,ivar) +               &
     &           rmask_full(i,j) * DiaBio2d1(i,j)
# else
                DiaBio2d(i,j,ivar) = DiaBio2d(i,j,ivar) + DiaBio2d1(i,j)
# endif
              END DO
            END DO
          END IF
        END DO
        DO ivar=1,NDbio3d
          IF (Dout(iDbio3(ivar),ng)) THEN
            i=idfabmD3(ivar)
            DiaBio3d1(Istr:Iend,Jstr:Jend,1:N(ng)) =                    &
     &       FMODEL(ng)%f(tile)%model%get_interior_diagnostic_data(i)
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef MASKING
                  DiaBio3d(i,j,k,ivar) = DiaBio3d(i,j,k,ivar) +         &
     &             rmask_full(i,j) * DiaBio3d1(i,j,k)
# else
                  DiaBio3d(i,j,k,ivar) = DiaBio3d(i,j,k,ivar) +         &
     &                                    DiaBio3d1(i,j,k)
# endif
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE biology_tile
