#include "cppdefs.h"
      MODULE roms_fabm_mod
!
! This module houses all the subroutines needed to link ROMS to FABM.
! It should be independent of ROMS version.
!
!!Phil Wallhead 2020

#ifdef RFABM

      implicit none

      PRIVATE
      PUBLIC  :: roms_fabm_init


      CONTAINS


      SUBROUTINE roms_fabm_init

      USE fabm
      USE fabm_config

      USE mod_param   !Provides access to BOUNDS, FMODEL etc.
      USE mod_ncparam !Provides access to Dout etc.
      USE mod_parallel!Provides access to MyRank
      USE mod_forces  !Provides access to FORCES
      USE mod_grid    !Provides access to GRID
      USE mod_ocean   !Provides access to OCEAN (hence "t")
      USE mod_biology !Provides access to idbio, iDbio2 etc. via rfabm_mod.h
      USE mod_scalars !Provides access to ydayc = No. days since start of year
      USE mod_iounits !Provides access to stdout

      integer :: ng, tile, Istr, Iend, Jstr, Jend, nx, ny, i, j, k
      integer :: itrc, ibio, NDbio2dmax, NDbio3dmax 

# ifdef DEBUGFABM
  20  FORMAT (1x,a,t26,a)
# endif
# ifdef DISTRIBUTE
      tile = MyRank
# else
      tile = 0  !Zero to access the full domain (as for "t", see mod_arrays -> allocate ocean)
# endif
      DO ng=1,Ngrids
        !  Read FABM's runtime config from fabm.yaml in current working directory
!        CALL fabm_create_model_from_yaml_file(FMODEL(ng)%f(tile)%model)
        FMODEL(ng)%f(tile)%model => fabm_create_model()
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm_create_model"
# endif

# if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
        !  All diagnostics unsaved by default
        NDbio2dmax =                                                    &
     &  size(FMODEL(ng)%f(tile)%model%horizontal_diagnostic_variables)
        DO i=1,NDbio2dmax
          FMODEL(ng)%f(tile)%model%                                     &
     &    horizontal_diagnostic_variables(i)%save = .FALSE.
        END DO
        !  Switch on the specified diagnostics
        DO i=1,NDbio2d
          IF (Dout(iDbio2(i),ng)) THEN          
            ibio=idfabmD2(i)
            FMODEL(ng)%f(tile)%model%                                   &
     &      horizontal_diagnostic_variables(ibio)%save = .TRUE.
#  ifdef DEBUGFABM
            IF (Master) write(stdout,20) "Save horiz. diag. ",          &
     & TRIM(FMODEL(ng)%f(tile)%model%                                   &
     &  horizontal_diagnostic_variables(ibio)%name)
#  endif
          END IF
        END DO

        !  All diagnostics unsaved by default
        NDbio3dmax =                                                    &
     &  size(FMODEL(ng)%f(tile)%model%interior_diagnostic_variables)
        DO i=1,NDbio3dmax
          FMODEL(ng)%f(tile)%model%                                     &
     &    interior_diagnostic_variables(i)%save = .FALSE.
        END DO
        !  Switch on the specified diagnostics
        DO i=1,NDbio3d
          IF (Dout(iDbio3(i),ng)) THEN          
            ibio=idfabmD3(i)
            FMODEL(ng)%f(tile)%model%                                   &
     &      interior_diagnostic_variables(ibio)%save = .TRUE.
#  ifdef DEBUGFABM
            IF (Master) write(stdout,20) "Save interior diag. ",TRIM(   &
     &FMODEL(ng)%f(tile)%model%interior_diagnostic_variables(ibio)%name)
#  endif
          END IF
        END DO
# endif

        !  Provide FABM with the extents of the spatial domain
        Istr = BOUNDS(ng)%Istr(tile) !Note: this is the lower bound index NOT INCLUDING GHOST POINTS (cf. LBi)
        Iend = BOUNDS(ng)%Iend(tile)
        Jstr = BOUNDS(ng)%Jstr(tile)
        Jend = BOUNDS(ng)%Jend(tile)

        nx = Iend-Istr+1
        ny = Jend-Jstr+1

        CALL FMODEL(ng)%f(tile)%model%set_domain(nx,ny,N(ng))
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm set_domain"
# endif

# ifdef MASKING
        CALL FMODEL(ng)%f(tile)%model%set_mask(                         &
     &                  GRID(ng)%rmask(Istr:Iend,Jstr:Jend))
        ! Note: Only the permanent land mask rmask should be passed to FABM at this point
#  ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm_set_mask"
#  endif
# endif

        !! Link the FABM state variables
        ! First link the interior state variables to FABM (using temporary variable "state1")
        ! This is a permanent link.
        DO itrc=1,NBT
          CALL FMODEL(ng)%f(tile)%model%link_interior_state_data(       &
     &     itrc,OCEAN(ng)%state1(Istr:Iend,Jstr:Jend,1:N(ng),itrc))
        END DO
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm link_interior_state_data"
# endif

        ! Next link the surface-attached state variables (using "state_sf" for nstp=1, cf. ana_initial.h)
        ! This link will have to be adjusted every time step in rfabm.h as nstp varies
        DO itrc=1,NSAT
          CALL FMODEL(ng)%f(tile)%model%link_surface_state_data(        &
     &     itrc,OCEAN(ng)%state_sf(Istr:Iend,Jstr:Jend,1,itrc))
        END DO
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm link_surface_state_data"
# endif

        ! Next link the bottom-attached (benthic) state variables (using "state_bt" for nstp=1, cf. ana_initial.h)
        ! This link will have to be adjusted every time step in rfabm.h as nstp varies
        DO itrc=1,NBAT
          CALL FMODEL(ng)%f(tile)%model%link_bottom_state_data(         &
     &     itrc,OCEAN(ng)%state_bt(Istr:Iend,Jstr:Jend,1,itrc))
        END DO
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm link_bottom_state_data"
# endif

        ! Set the pointers for cell thickness and bgc forcings
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%cell_thickness)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !cell thickness (m)
     &     fabm_standard_variables%cell_thickness,                      &
     &     GRID(ng)%Hz(Istr:Iend,Jstr:Jend,:))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%depth)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !depth of layer midpoints (m)
     &     fabm_standard_variables%depth,                               &
     &     GRID(ng)%z_r(Istr:Iend,Jstr:Jend,:))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%temperature)) THEN
          FMODEL(ng)%needsT(tile) = .true. !This will be used to parse subsequent linking calls
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !temperature [degC]
     &     fabm_standard_variables%temperature,                         &
     &     OCEAN(ng)%t(Istr:Iend,Jstr:Jend,1:N(ng),1,1))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%practical_salinity)) THEN
          FMODEL(ng)%needsS(tile) = .true. !This will be used to parse subsequent linking calls
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !salinity [psu]
     &     fabm_standard_variables%practical_salinity,                  &
     &     OCEAN(ng)%t(Istr:Iend,Jstr:Jend,1:N(ng),1,2))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%density)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !density [kg/m3]
     &     fabm_standard_variables%density,                             &
     &     OCEAN(ng)%rhof(Istr:Iend,Jstr:Jend,1:N(ng)))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%pressure)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !pressure due to water [dbar]
     &     fabm_standard_variables%pressure,                            &
     &     OCEAN(ng)%pres(Istr:Iend,Jstr:Jend,1:N(ng)))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &fabm_standard_variables%surface_downwelling_shortwave_flux)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !surface downwelling shortwave flux [W/m2]
     &     fabm_standard_variables%surface_downwelling_shortwave_flux,  &
     &     FORCES(ng)%srflxWm2(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%wind_speed)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !wind speed [m/s]
     &     fabm_standard_variables%wind_speed,                          &
     &     FORCES(ng)%Swind(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%bottom_stress)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !bottom stress [Pa]
     &     fabm_standard_variables%bottom_stress,                       &
     &     FORCES(ng)%bstrPa(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%bottom_depth_below_geoid)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !bottom depth below geoid [m]
     &     fabm_standard_variables%bottom_depth_below_geoid,            &
     &     GRID(ng)%h(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%longitude)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !longitude [degrees]
     &     fabm_standard_variables%longitude,                           &
     &     GRID(ng)%lonr(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &   fabm_standard_variables%latitude)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !latitude [degrees]
     &     fabm_standard_variables%latitude,                            &
     &     GRID(ng)%latr(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &  fabm_standard_variables%number_of_days_since_start_of_the_year)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_scalar(                    & !year day (0-365) [days]
     &fabm_standard_variables%number_of_days_since_start_of_the_year,   &
     &     ydayc(ng))
        END IF

!!! Additional forcing variables that require cpp activation
# ifdef FABM_ASILT
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &  fabm_standard_variables%absorption_of_silt)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_interior_data(             & !absorption due to silt [m-1]
     &     fabm_standard_variables%absorption_of_silt,                  &
     &     FORCES(ng)%Asilt(Istr:Iend,Jstr:Jend,1:N(ng)))
        END IF
# endif
!! Note: absorption of silt (needed by light_iop) is a custom variable for ERSEM, not a standard variable
!! It is actually the PAR attenuation due to silt
!! Constant values can be implemented more simply by adding to the fabm.yaml:
! silt_abs:
!   model: bulk_constant
!   parameters:
!     standard_name: absorption_of_silt
!     value: 0.1
# ifdef FABM_PCO2ATM
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &  fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !partial pressure of atmospheric CO2 [uatm]
     &  fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air, & !note: the name in FABM is not quite correct (only equivalent if total atmospheric pressure = 1 atm)
     &     FORCES(ng)%pCO2atm(Istr:Iend,Jstr:Jend))
        END IF
# endif
# ifdef FABM_ADYTRACER
        IF (FMODEL(ng)%f(tile)%model%variable_needs_values(             &
     &  fabm_standard_variables%gelbstoff_absorption_satellite)) THEN
          CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(           & !absorption due to dissolved and detrital matter (adg) [m-1]
     &  fabm_standard_variables%gelbstoff_absorption_satellite,         &
     &     FORCES(ng)%ADY_0(Istr:Iend,Jstr:Jend))
        END IF
# endif
# ifdef FABM_N3ATMDEPO
        CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(             & !atmospheric deposition of oxidized nitrogen [mmolN/m2/s]
     &   'n3_deposition/flux',FORCES(ng)%N3atmd(Istr:Iend,Jstr:Jend))
# endif
# ifdef FABM_N4ATMDEPO
        CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(             & !atmospheric deposition of reduced nitrogen [mmolN/m2/s]
     &   'n4_deposition/flux',FORCES(ng)%N4atmd(Istr:Iend,Jstr:Jend))
# endif
# ifdef FABM_ISOURCES
        IF (FMODELS(ng)%fmodel(tile)%variable_needs_values(             &
     &  standard_variables%absorption_of_silt)) THEN
          CALL FMODELS(ng)%fmodel(tile)%link_interior_data(             & !internal sources via FABM (work in progress, not yet usable)
     &     standard_variables%absorption_of_silt,                       &
     &     FORCES(ng)%Asilt(Istr:Iend,Jstr:Jend,1:N(ng)))
        END IF
# endif
# ifdef FABM_AICE
        CALL FMODEL(ng)%f(tile)%model%link_horizontal_data(             & !fraction of cell covered by ice from internal sea-ice module (Hedstrom code only) []
     &   'aice',ICE(ng)%ai(Istr:Iend,Jstr:Jend,:))
# endif

# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm_link_data calls"
# endif

        ! Set variable ids for pointers that will change depending on time step (via nstp)
        !Note: We cannot use variable_needs_values to parse these calls, since these will return
        !      .false. after the first linking calls above.
        IF (FMODEL(ng)%needsT(tile)) THEN
          FMODEL(ng)%id_temp(tile) =                                    &
     &      FMODEL(ng)%f(tile)%model%get_interior_variable_id(          &
     &      fabm_standard_variables%temperature)
        END IF
        IF (FMODEL(ng)%needsS(tile)) THEN
          FMODEL(ng)%id_salt(tile) =                                    &
     &      FMODEL(ng)%f(tile)%model%get_interior_variable_id(          &
     &       fabm_standard_variables%practical_salinity)
        END IF
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done setting runtime variable ids"
# endif

        ! Check this fabm model is ready
!        CALL fabm_check_ready(FMODEL(ng)%f(tile)%model)
        CALL FMODEL(ng)%f(tile)%model%start()
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm start"
# endif

# ifdef FABM_INITIAL
        ! Call fabm_initialize_state for each layer and j coordinate
        ! NOTE: This will OVERRIDE any ROMS tracer initialization via .nc input
        !       --- call to ROMS_initialise occurs before call to roms_fabm_init in ocean.h
        DO k=1,N(ng)
          DO j=1,ny
            CALL FMODEL(ng)%f(tile)%model%initialize_interior_state(    &
     &                                 1,nx,j,k)
          END DO
        END DO
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                OCEAN(ng)%t(i,j,k,1,ibio) = OCEAN(ng)%state1(i,j,k,itrc)
              END DO
            END DO
          END DO
        END DO
#  ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done initialize_interior_state"
#  endif
# endif

# if defined FABM_INITIAL || defined FABM_INITIAL_SB
        IF (NSAT.gt.0) THEN
          DO j=1,ny
            CALL FMODEL(ng)%f(tile)%model%initialize_surface_state(     &
     &                                        1,nx,j)
          END DO
#  ifdef DEBUGFABM
          IF (Master) write(stdout,*) "Done initialize_surface_state"
#  endif
        END IF
        IF (NBAT.gt.0) THEN
          DO j=1,ny
            CALL FMODEL(ng)%f(tile)%model%initialize_bottom_state(      &
     &                                        1,nx,j)
          END DO
#  ifdef DEBUGFABM
          IF (Master) write(stdout,*) "Done initialize_bottom_state"
#  endif
        END IF
# endif

      END DO

      RETURN
      END SUBROUTINE roms_fabm_init

#endif

      END MODULE roms_fabm_mod
