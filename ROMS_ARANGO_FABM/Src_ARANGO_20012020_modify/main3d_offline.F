#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D && (defined OFFLINE)
      SUBROUTINE main3d_offline (RunInterval)
!
!svn $Id: main3d.F 995 2020-01-10 04:01:28Z arango $
!=======================================================================
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for all  nested  grids, if any,     !
!  for the specified time interval (seconds), RunInterval.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_coupler
# endif
      USE mod_iounits
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_stepping
!
# ifdef ANA_VMIX
      USE analytical_mod,       ONLY : ana_vmix
# endif
# ifdef BIOLOGY
      USE biology_mod,          ONLY : biology
# endif
# ifdef BBL_MODEL
      USE bbl_mod,              ONLY : bblm
# endif
# ifdef BULK_FLUXES
      USE bulk_flux_mod,        ONLY : bulk_flux
# endif
# ifdef BVF_MIXING
      USE bvf_mix_mod,          ONLY : bvf_mix
# endif
      USE dateclock_mod,        ONLY : time_string
      USE diag_mod,             ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod,       ONLY : nl_dotproduct
# endif
# if defined NLM_OUTER           || defined W4DPSAS                 ||\
     defined W4DPSAS_SENSITIVITY || defined W4DPSAS_FCT_SENSITIVITY
      USE forcing_mod,          ONLY : forcing
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE frc_adjust_mod,       ONLY : frc_adjust
# endif
# ifdef GLS_MIXING
      USE gls_corstep_mod,      ONLY : gls_corstep
      USE gls_prestep_mod,      ONLY : gls_prestep
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
      USE hmixing_mod,          ONLY : hmixing
# endif
      USE ini_fields_mod,       ONLY : ini_fields, ini_zeta
# ifdef LMD_MIXING
      USE lmd_vmix_mod,         ONLY : lmd_vmix
# endif
# ifdef MY25_MIXING
      USE my25_corstep_mod,     ONLY : my25_corstep
      USE my25_prestep_mod,     ONLY : my25_prestep
# endif
# ifdef NESTING
      USE nesting_mod,          ONLY : nesting
#  ifndef ONE_WAY
      USE nesting_mod,          ONLY : do_twoway
#  endif
# endif
# if defined ADJUST_BOUNDARY
      USE obc_adjust_mod,       ONLY : obc_adjust, load_obc
# endif
# if defined ATM_COUPLING && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : ocn2atm_coupling
# endif
# if defined WAV_COUPLING && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : ocn2wav_coupling
# endif
      USE omega_mod,            ONLY : omega
# ifdef NEARSHORE_MELLOR
      USE radiation_stress_mod, ONLY : radiation_stress
# endif
# ifndef TS_FIXED
      USE rho_eos_mod,          ONLY : rho_eos
# endif
      USE rhs3d_mod,            ONLY : rhs3d
# ifdef SEDIMENT
      USE sediment_mod,         ONLY : sediment
# endif
# ifdef AVERAGES
      USE set_avg_mod,          ONLY : set_avg
# endif
      USE set_depth_mod,        ONLY : set_depth
      USE set_massflux_mod,     ONLY : set_massflux, set_massflux_avg
# if defined SSH_TIDES || defined UV_TIDES
      USE set_tides_mod,        ONLY : set_tides
# endif
      USE set_vbc_mod,          ONLY : set_vbc
      USE set_zeta_mod,         ONLY : set_zeta, set_avg_zeta
      USE step2d_mod,           ONLY : step2d
# ifndef TS_FIXED
      USE step3d_t_mod,         ONLY : step3d_t
# endif
      USE step3d_uv_mod,        ONLY : step3d_uv
# ifdef FLOATS
      USE step_floats_mod,      ONLY : step_floats
# endif
      USE strings_mod,          ONLY : FoundError
      USE wvelocity_mod,        ONLY : wvelocity
!
      implicit none
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoNestLayer, Time_Step, FIRSTCALL

      integer :: Nsteps, Rsteps
      integer :: ig, il, istep, ng, nl, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations by the specified time.
!=======================================================================
!
!  Time-step the 3D kernel for the specified time interval (seconds),
!  RunInterval.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
!
      KERNEL_LOOP : DO WHILE (Time_Step)
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
        nl=0
#ifdef NESTING
        TwoWayInterval(1:Ngrids)=0.0_r8
#endif
!
        NEST_LAYER : DO WHILE (DoNestLayer)
!
!  Determine number of time steps to compute in each nested grid layer
!  based on the specified time interval (seconds), RunInterval. Non
!  nesting applications have NestLayers=1. Notice that RunInterval is
!  set in the calling driver. Its value may span the full period of the
!  simulation, a multi-model coupling interval (RunInterval > ifac*dt),
!  or just a single step (RunInterval=0).
!
          CALL ntimesteps (iNLM, RunInterval, nl, Nsteps, Rsteps)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
          IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
!
!  Time-step governing equations for Nsteps.
!
          STEP_LOOP : DO istep=1,Nsteps
!
!  Set time indices and time clock.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)+1
              nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
              nnew(ng)=3-nstp(ng)
              nrhs(ng)=nstp(ng)
              time(ng)=time(ng)+dt(ng)
              tdays(ng)=time(ng)*sec2day
              CALL time_string (time(ng), time_code(ng))
              IF (step_counter(ng).eq.Rsteps) Time_Step=.FALSE.
            END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
!$OMP MASTER
              CALL get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL set_data (ng, tile)
              END DO
!$OMP BARRIER
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntstart(ng)) THEN

                ! Set this here so that only occurs the first time
                ! set_avg_zeta is called, the first time step
                FIRSTCALL = .TRUE.


!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  ! commented by KMT since I am inputting zeta. Revisit this
                  ! if want to be able to calculate zeta instead of
                  ! inputting.
                  !CALL ini_zeta (ng, tile, iNLM)
                  CALL set_avg_zeta(ng, tile, FIRSTCALL)  ! to initialize Zt_avg1
                  CALL set_depth (ng, tile, iNLM)
                END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ini_fields (ng, tile, iNLM)
                END DO
!$OMP BARRIER

# ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL nesting (ng, iNLM, ngetD)
                END IF
# endif
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL set_massflux (ng, tile, iNLM)
# ifndef TS_FIXED
                CALL rho_eos (ng, tile, iNLM)
# endif
                CALL diag (ng, tile)
              END DO
!$OMP BARRIER
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN

# if defined OFFLINE_TPASSIVE || defined OFFLINE_BIOLOGY
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
#  ifdef BULK_FLUXES
                CALL bulk_flux (ng, tile)
#  endif
#  ifdef BBL_MODEL
                CALL bblm (ng, tile)
#  endif
                CALL set_vbc (ng, tile)
#  if defined SSH_TIDES || defined UV_TIDES
                CALL set_tides (ng, tile)
#  endif
              END DO
!$OMP BARRIER
            END DO
#  ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for bottom stress variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, nbstr)
              END IF
            END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
# if defined OFFLINE_TPASSIVE || defined OFFLINE_BIOLOGY
#  if defined ANA_VMIX
                CALL ana_vmix (ng, tile, iNLM)
#  elif defined LMD_MIXING
                CALL lmd_vmix (ng, tile)
#  elif defined BVF_MIXING
                CALL bvf_mix (ng, tile)
#  endif
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
                CALL hmixing (ng, tile)
# endif
# if !defined OCLIMATOLOGY
                CALL omega (ng, tile, iNLM)
# endif
                CALL wvelocity (ng, tile, nstp(ng))
              END DO
!$OMP BARRIER
            END DO
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1     ! irreversible
                ! KMT commented for all case since will be inputs
                !CALL set_zeta (ng, tile)
# ifdef DIAGNOSTICS
                CALL set_diags (ng, tile)
# endif
# ifdef AVERAGES
                CALL set_avg (ng, tile)
# endif
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for 3D kernel free-surface.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, nzeta)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
!$OMP MASTER
              CALL output (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF ((FoundError(exit_flag, NoError, __LINE__,             &
     &            __FILE__)).or.                                        &
     &            ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) THEN
                RETURN
              END IF
            END DO

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  If refinement grid, interpolate (space, time) state variables
!  contact points from donor grid extracted data.
#  ifdef NESTING_DEBUG
!
!  Also, fill BRY_CONTACT(:,:)%Mflux to check for mass conservation
!  between coarse and fine grids.  This is only done for diagnostic
!  purposes. Also, debugging is possible with very verbose output
!  to fort.300 is allowed by activating uppercase(nesting_debug).
#  endif
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                CALL nesting (ng, iNLM, nputD)
#  ifdef NESTING_DEBUG
                CALL nesting (ng, iNLM, nmflx)
#  endif
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
# if defined OFFLINE_TPASSIVE || defined OFFLINE_BIOLOGY
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rhs3d (ng, tile)
#  if ((!defined OFFLINE && \
        ( !defined MIXCLIMATOLOGY || !defined AKXCLIMATOLOGY )) || \
       ((defined CALC_TDR || defined CALC_TSR) && \
       !defined MIXCLIMATOLOGY))
#   ifdef MY25_MIXING
                CALL my25_prestep (ng, tile)
#   elif defined GLS_MIXING
                CALL gls_prestep (ng, tile)
#   endif
#  endif
              END DO
!$OMP BARRIER
            END DO

#  ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for right-hand-side terms
!  (tracers).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, nrhst)
              END IF
            END DO
#  endif
# endif

# ifdef NESTING
#  if defined MASKING && defined WET_DRY
!
!  If nesting and wetting and drying, scale horizontal interpolation
!  weights to account for land/sea masking in contact areas. This needs
!  to be done at very time-step since the Land/Sea masking is time
!  dependent.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL nesting (ng, iNLM, nmask)
            END DO
#  endif
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the time-averaged
!  momentum fluxes (DU_avg1, DV_avg1) and free-surface (Zt_avg).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, n2dfx)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.
!  KMT In this offline case, the fast time steps are not run. Instead,
!  zeta is updated here with the average of this and the subsequent time
!  step values. Then, depths are updated to use this new zeta value
!  (stored in Zt_avg1). 
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL set_avg_zeta (ng, tile, FIRSTCALL)
                CALL set_depth (ng, tile, iNLM)
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If nesting, determine vertical indices and vertical interpolation
!  weights in the contact zone using new depth arrays.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL nesting (ng, iNLM, nzwgt)
            END DO
# endif


! KMT this is where step3d_uv would occur in online case. Need to update
! Huon here using the newly calculated depths and the next time step's
! u field?
!-----------------------------------------------------------------------
!  KMT This calculates values for Huon, Hvom, Akt (indices 1 and 2)
!  using updated averaged Zt_avg1 (zeta) value so that offline tracer is
!  calculated correctly. 
!-----------------------------------------------------------------------

            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL set_massflux_avg (ng, tile, iNLM)
              END DO
            END DO

!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                ! KMT run this even if OCLIMATOLOGY is on because omega
                ! needs to be updated mid-time loop for tracer advection.
                CALL omega (ng, tile, iNLM)
# if ((!defined OFFLINE && \
       ( !defined MIXCLIMATOLOGY || !defined AKXCLIMATOLOGY )) || \
      (defined CALC_TDR || defined CALC_TSR))
! PWA: Note that we allow the call to my25/gls_corstep if CALC_TDR/TSR is defined
!      even if we have MIXCLIMATOLOGY defined. This is to enable calculation
!      of tdr/tsr within my25/gls_corstep -- any modifications to tke/gls will
!      be avoided within my25/gls_corstep if MIXCLIMATOLOGY is defined.
#  ifdef MY25_MIXING
                CALL my25_corstep (ng, tile)
#  elif defined GLS_MIXING
                CALL gls_corstep (ng, tile)
#  endif
# endif

# ifdef BIOLOGY
                CALL biology (ng, tile)
# endif
# ifdef SEDIMENT
                CALL sediment (ng, tile)
# endif
              END DO
!$OMP BARRIER
            END DO

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!

#  if defined OFFLINE_TPASSIVE || defined OFFLINE_BIOLOGY
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL step3d_t (ng, tile)
              END DO
!$OMP BARRIER
            END DO

#   ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for Tracer Variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, n3dTV)
              END IF
            END DO
#   endif
#  endif
# endif

# ifdef NESTING
#  ifndef ONE_WAY
!
!-----------------------------------------------------------------------
!  If refinement grids, perform two-way coupling between fine and
!  coarse grids. Correct coarse grid tracers values at the refinement
!  grid with refined accumulated fluxes.  Then, replace coarse grid
!  state variable with averaged refined grid values (two-way nesting).
!  Update coarse grid depth variables.
!
!  The two-way exchange of infomation between nested grids needs to be
!  done at the correct time-step and in the right sequence.
!-----------------------------------------------------------------------
!
            DO il=NestLayers,1,-1
              DO ig=1,GridsInLayer(il)
                ng=GridNumber(ig,il)
                IF (do_twoway(iNLM, nl, il, ng, istep)) THEN
                  CALL nesting (ng, iNLM, n2way)
                END IF
              END DO
            END DO
#  endif
!
!-----------------------------------------------------------------------
!  If donor to a finer grid, extract data for the external contact
!  points. This is the latest solution for the coarser grid.
!
!  It is stored in the REFINED structure so it can be used for the
!  space-time interpolation when "nputD" argument is used above.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (DonorToFiner(ng)) THEN
                CALL nesting (ng, iNLM, ngetD)
              END IF
            END DO
# endif
          END DO STEP_LOOP

        END DO NEST_LAYER

      END DO KERNEL_LOOP

      RETURN
      END SUBROUTINE main3d_offline
#else
      SUBROUTINE main3d_offline
      RETURN
      END SUBROUTINE main3d_offline
#endif
