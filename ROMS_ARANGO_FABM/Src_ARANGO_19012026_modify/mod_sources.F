#include "cppdefs.h"
      MODULE mod_sources
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Msrc       Maximum number of analytical point Sources/Sinks.        !
!  Nsrc       Number of point Sources/Sinks.                           !
!  Dsrc       Direction of point Sources/Sinks:                        !
!               Dsrc(:) = 0,  Along XI-direction.                      !
!               Dsrc(:) > 0,  Along ETA-direction.                     !
!  Fsrc       Point Source/Sinks identification flag:                  !
!               Fsrc(:) = 0,  All Tracer source/sink are off.          !
!               Fsrc(:) = 1,  Only temperature is on.                  !
!               Fsrc(:) = 2,  Only salinity is on.                     !
!               Fsrc(:) = 3,  Both temperature and salinity are on.    !
!               Fsrc(:) = 4,  Both nitrate and salinity are on.        !
!               Fsrc(:) = ... And other combinations.                  !
!                             (We need a more robust logic here)       !
!  Isrc       I-grid location of point Sources/Sinks,                  !
!               1 =< Isrc =< Lm(ng).                                   !
!  Jsrc       J-grid location of point Sources/Sinks,                  !
!               1 =< Jsrc =< Mm(ng).                                   !
!  Qbar       Vertically integrated mass transport (m3/s) of point     !
!               Sources/Sinks at U- or V-points:                       !
!               Qbar -> positive, if the mass transport is in the      !
!                       positive U- or V-direction.                    !
!               Qbar -> negative, if the mass transport is in the      !
!                       negative U- or V-direction.                    !
!  QbarG      Latest two-time snapshots of vertically integrated       !
!               mass transport (m3/s) of point Sources/Sinks.          !
!  Qshape     Nondimensional shape function to distribute mass         !
!               mass point Sources/Sinks vertically.                   !
!  Qsrc       Mass transport profile (m3/s) of point Sources/Sinks.    !
!  Tsrc       Tracer (tracer units) point Sources/Sinks.               !
!  TsrcG      Latest two-time snapshots of tracer (tracer units)       !
!               point Sources/Sinks.                                   !
!                                                                      !
!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
!  Nisrc      Number of internal tracer point Sources/Sinks.           !
!  Iisrc      I-grid location of internal point Sources/Sinks,         !
!               1 =< Iisrc =< Lm(ng).                                  !
!  Jisrc      J-grid location of internal point Sources/Sinks,         !
!               1 =< Jisrc =< Mm(ng).                                  !
!  Tisrc      Tracer (tracer units) internal point Sources/Sinks.      !
!  TisrcG     Latest two-time snapshots of tracer (tracer units)       !
!               internal point Sources/Sinks.                          !
#endif
!!!End insert PWA 18/10/2021
!=======================================================================
!
        USE mod_kinds
        USE mod_param
!
        implicit none
!
        PUBLIC :: allocate_sources
        PUBLIC :: deallocate_sources
!
!-----------------------------------------------------------------------
!  Define T_SOURCES structure.
!-----------------------------------------------------------------------
!
        TYPE T_SOURCES

          integer, pointer :: Isrc(:)
          integer, pointer :: Jsrc(:)
!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
          integer, pointer :: Iisrc(:)
          integer, pointer :: Jisrc(:)
          integer, pointer :: Nisrct(:)
          integer, pointer :: Iisrct(:,:)
          integer, pointer :: Jisrct(:,:)
          integer, pointer :: isrct(:,:)
#endif
!!!End insert PWA 18/10/2021

          real(r8), pointer :: Dsrc(:)
          real(r8), pointer :: Fsrc(:)
          real(r8), pointer :: Qbar(:)
          real(r8), pointer :: Qshape(:,:)
          real(r8), pointer :: Qsrc(:,:)
          real(r8), pointer :: Tsrc(:,:,:)
          real(r8), pointer :: Xsrc(:)
          real(r8), pointer :: Ysrc(:)

#ifndef ANA_PSOURCE
          real(r8), pointer :: QbarG(:,:)
          real(r8), pointer :: TsrcG(:,:,:,:)
#endif

!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
          real(r8), pointer :: Tisrc(:,:,:)
          real(r8), pointer :: Xisrc(:)
          real(r8), pointer :: Yisrc(:)
          real(r8), pointer :: TisrcG(:,:,:,:)
#endif
!!!End insert PWA 18/10/2021

#ifdef ADJOINT
          real(r8), pointer :: ad_Qbar(:)
          real(r8), pointer :: ad_Qsrc(:,:)
          real(r8), pointer :: ad_Tsrc(:,:,:)
#endif

#ifdef TANGENT
          real(r8), pointer :: tl_Qbar(:)
          real(r8), pointer :: tl_Qsrc(:,:)
          real(r8), pointer :: tl_Tsrc(:,:,:)
#endif

        END TYPE T_SOURCES
!
        TYPE (T_SOURCES), allocatable :: SOURCES(:)
!
!-----------------------------------------------------------------------
!  Define other variables in module.
!-----------------------------------------------------------------------
!
        integer, allocatable :: Msrc(:)
        integer, allocatable :: Nsrc(:)
!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
        integer, allocatable :: Nisrc(:)
#endif
!!!End insert PWA 18/10/2021
!
      CONTAINS
!
      SUBROUTINE allocate_sources (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in the module  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!!!Modified PWA 18/10/2021
#if !defined ANA_PSOURCE || defined TS_ISOURCE
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
      USE mod_scalars
#endif
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
!!!Modified PWA 20/11/2019 following Hedstrom 14092019 code
      integer :: i, ng 
!!!End modified PWA 20/11/2019
!
!  Local variable declarations.
!
!!!Modified PWA 18/10/2021
#if !defined ANA_PSOURCE || defined TS_ISOURCE
      logical :: foundit
!
      integer :: Vid, ifile, nvatt, nvdim
#endif
      integer :: is, itrc, k, mg
#ifdef TS_ISOURCE
      integer :: Ntiles 
#endif

      real(r8), parameter :: IniVal = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", allocate_sources"

#ifndef ANA_PSOURCE
# if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (Var_desc_t)  :: my_pioVar
# endif
#endif
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(Msrc)) THEN
        allocate ( Msrc(Ngrids) )
      END IF

      IF (.not.allocated(Nsrc)) THEN
        allocate ( Nsrc(Ngrids) )
      END IF

#ifndef ANA_PSOURCE
!
!  Inquire about the number of point Sources/Sinks.
!
      IF (ng.eq.1) THEN
        DO mg=1,Ngrids
          foundit=.FALSE.
          IF (LuvSrc(mg).or.LwSrc(mg).or.ANY(LtracerSrc(:,mg))) THEN
!!!Modified PWA 21/01/2026 to allow for multiple files (use only first)
            SELECT CASE (SSF(1,ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_inq_var (ng, iNLM, SSF(1,mg)%name,          &
     &                               MyVarName = Vname(1,idRxpo),       &
     &                               SearchVar = foundit,               &
     &                               VarID = Vid,                       &
     &                               nVardim = nvdim,                   &
     &                               nVarAtt = nvatt)

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_inq_var (ng, iNLM, SSF(1,mg)%name,      &
     &                                   MyVarName = Vname(1,idRxpo),   &
     &                                   SearchVar = foundit,           &
     &                                   pioVar = my_pioVar,            &
     &                                   nVardim = nvdim,               &
     &                                   nVarAtt = nvatt)
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (foundit) THEN
              Nsrc(mg)=var_Dsize(1)         ! first dimension
              Msrc(mg)=Nsrc(mg)
              CALL check_sources (mg, SSF(1,mg)%name, Nsrc(mg))
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                         MyFile)) RETURN

            END IF
!!!End modified PWA 21/01/2026
          END IF
        END DO
      END IF
#else
!
!  Set number of point sources to maximum number of analytical sources.
!  Notice that a maximum of 200 analytical sources are set-up here.
!
      Msrc(ng)=200
      Nsrc(ng)=Msrc(ng)
#endif

!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
      IF (ng.eq.1) THEN
        allocate ( Nisrc(Ngrids) )
      END IF
!
!  Inquire about the number of internal tracer point sources.
!
      IF (ng.eq.1) THEN
        DO mg=1,Ngrids
          foundit=.FALSE.
          IF (ANY(LtracerISrc(:,mg))) THEN
            SELECT CASE (ISSF(1,ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_inq_var (ng, iNLM, ISSF(1,mg)%name,         &
     &                               MyVarName = Vname(1,idIxpo),       &
     &                               SearchVar = foundit,               &
     &                               VarID = Vid,                       &
     &                               nVardim = nvdim,                   &
     &                               nVarAtt = nvatt)

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_inq_var (ng, iNLM, ISSF(1,mg)%name,     &
     &                                   MyVarName = Vname(1,idIxpo),   &
     &                                   SearchVar = foundit,           &
     &                                   pioVar = my_pioVar,            &
     &                                   nVardim = nvdim,               &
     &                                   nVarAtt = nvatt)
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (foundit) THEN
              Nisrc(mg)=var_Dsize(1)         ! first dimension
            END if
          END IF
        END DO
      END IF
#endif
!!!End insert PWA 18/10/2021

!
!  Allocate structure.
!
      IF (ng.eq.1) allocate ( SOURCES(Ngrids) )
!
!  Allocate point Sources/Sinks variables.
!
      allocate ( SOURCES(ng) % Isrc(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % Jsrc(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % Dsrc(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % Fsrc(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % Qbar(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % Qshape(Nsrc(ng),N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng),r8)

      allocate ( SOURCES(ng) % Qsrc(Nsrc(ng),N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng),r8)

      allocate ( SOURCES(ng) % Tsrc(Nsrc(ng),N(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*NT(ng),r8)

      allocate ( SOURCES(ng) % Xsrc(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % Ysrc(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

#ifndef ANA_PSOURCE
      allocate ( SOURCES(ng) % QbarG(Nsrc(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % TsrcG(Nsrc(ng),N(ng),2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nsrc(ng)*N(ng)*NT(ng),r8)
#endif

!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
      allocate ( SOURCES(ng) % Iisrc(Nisrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng),r8)

      allocate ( SOURCES(ng) % Jisrc(Nisrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng),r8)

      Ntiles=NtileI(ng)*NtileJ(ng)
      allocate ( SOURCES(ng) % Nisrct(Ntiles) )
      Dmem(ng)=Dmem(ng)+REAL(Ntiles,r8)

      allocate ( SOURCES(ng) % Iisrct(Nisrc(ng),Ntiles) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng)*Ntiles,r8)

      allocate ( SOURCES(ng) % Jisrct(Nisrc(ng),Ntiles) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng)*Ntiles,r8)

      allocate ( SOURCES(ng) % isrct(Nisrc(ng),Ntiles) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng)*Ntiles,r8)

      allocate ( SOURCES(ng) % Tisrc(Nisrc(ng),N(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng)*N(ng)*NT(ng),r8)

      allocate ( SOURCES(ng) % Xisrc(Nisrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng),r8)

      allocate ( SOURCES(ng) % Yisrc(Nisrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nisrc(ng),r8)

      allocate ( SOURCES(ng) % TisrcG(Nisrc(ng),N(ng),2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nisrc(ng)*N(ng)*NT(ng),r8)
#endif
!!!End insert PWA 18/10/2021

#ifdef ADJOINT
      allocate ( SOURCES(ng) % ad_Qbar(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % ad_Qsrc(Nsrc(ng),N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng),r8)

      allocate ( SOURCES(ng) % ad_Tsrc(Nsrc(ng),N(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*NT(ng),r8)
#endif

#ifdef TANGENT
      allocate ( SOURCES(ng) % tl_Qbar(Nsrc(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng),r8)

      allocate ( SOURCES(ng) % tl_Qsrc(Nsrc(ng),N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng),r8)

      allocate ( SOURCES(ng) % tl_Tsrc(Nsrc(ng),N(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nsrc(ng)*N(ng)*NT(ng),r8)
#endif

!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO is=1,Nsrc(ng)
        SOURCES(ng) % Isrc(is) = 0
        SOURCES(ng) % Jsrc(is) = 0
        SOURCES(ng) % Dsrc(is) = IniVal
        SOURCES(ng) % Fsrc(is) = IniVal
        SOURCES(ng) % Xsrc(is) = IniVal
        SOURCES(ng) % Ysrc(is) = IniVal
        SOURCES(ng) % Qbar(is) = IniVal
#ifndef ANA_PSOURCE
        SOURCES(ng) % QbarG(is,1) = IniVal
        SOURCES(ng) % QbarG(is,2) = IniVal
#endif
#ifdef ADJOINT
        SOURCES(ng) % ad_Qbar(is) = IniVal
#endif
#ifdef TANGENT
        SOURCES(ng) % tl_Qbar(is) = IniVal
#endif
      END DO
      DO k=1,N(ng)
        DO is=1,Nsrc(ng)
          SOURCES(ng) % Qshape(is,k) = IniVal
          SOURCES(ng) % Qsrc(is,k) = IniVal
#ifdef ADJOINT
          SOURCES(ng) % ad_Qsrc(is,k) = IniVal
#endif
#ifdef TANGENT
          SOURCES(ng) % tl_Qsrc(is,k) = IniVal
#endif
        END DO
      END DO
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO is=1,Nsrc(ng)
            SOURCES(ng) % Tsrc(is,k,itrc) = IniVal
#ifdef ADJOINT
            SOURCES(ng) % ad_Tsrc(is,k,itrc) = IniVal
#endif
#ifdef TANGENT
            SOURCES(ng) % tl_Tsrc(is,k,itrc) = IniVal
#endif
#ifndef ANA_PSOURCE
            SOURCES(ng) % TsrcG(is,k,1,itrc) = IniVal
            SOURCES(ng) % TsrcG(is,k,2,itrc) = IniVal
#endif
          END DO
        END DO
      END DO

!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
      DO is=1,Nisrc(ng)
        SOURCES(ng) % Iisrc(is) = 0
        SOURCES(ng) % Jisrc(is) = 0
        SOURCES(ng) % Xisrc(is) = IniVal
        SOURCES(ng) % Yisrc(is) = IniVal
      END DO
      DO k=1,Ntiles
        SOURCES(ng) % Nisrct(k) = 0
        DO is=1,Nisrc(ng)
          SOURCES(ng) % Iisrct(is,k) = 0
          SOURCES(ng) % Jisrct(is,k) = 0
          SOURCES(ng) % isrct(is,k) = 0
        END DO
      END DO
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO is=1,Nisrc(ng)
            SOURCES(ng) % Tisrc(is,k,itrc) = IniVal
            SOURCES(ng) % TisrcG(is,k,1,itrc) = IniVal
            SOURCES(ng) % TisrcG(is,k,2,itrc) = IniVal
          END DO
        END DO
      END DO
#endif
!!!End insert PWA 18/10/2021
!
      RETURN
      END SUBROUTINE allocate_sources
!
      SUBROUTINE deallocate_sources (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================

#ifdef SUBOBJECT_DEALLOCATION
!
      USE destroy_mod, ONLY : destroy
#endif
!
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_mixing"

#ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_SOURCES structure
!  separately.
!-----------------------------------------------------------------------
!
      IF (.not.destroy(ng, SOURCES(ng)%Isrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Isrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Jsrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Jsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Dsrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Dsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Fsrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Fsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Qbar, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Qbar')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Qshape, MyFile,                  &
     &                 __LINE__, 'SOURCES(ng)%Qshape')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Qsrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Qsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Tsrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Tsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Xsrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Xsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%Ysrc, MyFile,                    &
     &                 __LINE__, 'SOURCES(ng)%Ysrc')) RETURN

# ifndef ANA_PSOURCE
      IF (.not.destroy(ng, SOURCES(ng)%QbarG, MyFile,                   &
     &                 __LINE__, 'SOURCES(ng)%QbarG')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%TsrcG, MyFile,                   &
     &                 __LINE__, 'SOURCES(ng)%TsrcG')) RETURN
# endif

# ifdef ADJOINT
      IF (.not.destroy(ng, SOURCES(ng)%ad_Qbar, MyFile,                 &
     &                 __LINE__, 'SOURCES(ng)%ad_Qbar')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%ad_Qsrc, MyFile,                 &
     &                 __LINE__, 'SOURCES(ng)%ad_Qsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%ad_Tsrc, MyFile,                 &
     &                 __LINE__, 'SOURCES(ng)%ad_Tsrc')) RETURN
# endif

# ifdef TANGENT
      IF (.not.destroy(ng, SOURCES(ng)%tl_Qbar, MyFile,                 &
     &                 __LINE__, 'SOURCES(ng)%tl_Qbar')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%tl_Qsrc, MyFile,                 &
     &                 __LINE__, 'SOURCES(ng)%tl_Qsrc')) RETURN

      IF (.not.destroy(ng, SOURCES(ng)%tl_Tsrc, MyFile,                 &
     &                 __LINE__, 'SOURCES(ng)%tl_Tsrc')) RETURN
# endif
#endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type SOURCES structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(SOURCES)) deallocate ( SOURCES )
      END IF
!
!-----------------------------------------------------------------------
!  Deallocate other variables in module.
!-----------------------------------------------------------------------
!
      IF (allocated(Msrc)) THEN
        deallocate ( Msrc )
      END IF

      IF (allocated(Nsrc)) THEN
        deallocate ( Nsrc )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_sources
!
      SUBROUTINE check_sources (ng, ncname, Npsrc)
!
!=======================================================================
!                                                                      !
!  It checks input NetCDF file for correctness in the spefication of   !
!  the point Source/Sink grid-cell face flag values (0, 1, 2) and the  !
!  implementation methodology using "LuvSrc" and or "LwSrc".           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError

!  Imported variable declarations.
!
      integer, intent(in) :: ng, Npsrc
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: i, ic_u, ic_v, ic_w
!
      real(r8) :: Dsrc_min, Dsrc_max
!
      real(r8) :: Dsrc(Npsrc)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", check_sources"
!
!-----------------------------------------------------------------------
!  Read in Point Source/Sink cell-face flag and check for correct values
!  according to the implementation methodology:
!
!    If LuvSrc = T,  needs Dsrc = 0 (flow across grid cell u-face) or
!                          Dsrc = 1 (flow across grid cell v-face)
!
!    If LwSrc = T,   needs Dsrc = 2 (flow across grid cell w-face)
!-----------------------------------------------------------------------
!
!  Read in Source/Sink cell-face flag.
!
      SELECT CASE (SSF(1,ng)%IOtype) !!!Modified PWA 21/01/2026 to allow for multifiles
        CASE (io_nf90)
          CALL netcdf_get_fvar (ng, iNLM, ncname,                       &
     &                          Vname(1,idRdir),  Dsrc,                 &
     &                          min_val = Dsrc_min,                     &
     &                          max_val = Dsrc_max)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_fvar (ng, iNLM, ncname,                   &
     &                              Vname(1,idRdir),  Dsrc,             &
     &                              min_val = Dsrc_min,                 &
     &                              max_val = Dsrc_max)
# endif
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Count and report Source/Sink for each cell-face flag value.
!
      ic_u=0
      ic_v=0
      ic_w=0
!
      DO i=1,NpSrc
        IF (INT(Dsrc(i)).eq.0) ic_u=ic_u+1
        IF (INT(Dsrc(i)).eq.1) ic_v=ic_v+1
        IF (INT(Dsrc(i)).eq.2) ic_w=ic_w+1
      END DO
!
      IF (Master) THEN
        IF (ng.eq.1) WRITE (stdout,10)
        WRITE (stdout,20) ng, TRIM(ncname),                             &
                          LuvSrc(ng), ic_u,                             &
                          LuvSrc(ng), ic_v,                             &
                          LwSrc(ng),  ic_w
      END IF
!
!  Stop if illegal configuration.
!
      IF (.not.LwSrc(ng).and.LuvSrc(ng).and.                            &
          (ic_u.eq.0).and.(ic_v.eq.0)) THEN
        IF (Master) WRITE (stdout,30) 'LuvSrc'
        exit_flag=5
      END IF
      IF (.not.LuvSrc(ng).and.LwSrc(ng).and.(ic_w.eq.0)) THEN
        IF (Master) WRITE (stdout,30) 'LwSrc'
        exit_flag=5
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,1x,'Point Sources/Sinks grid-cell flag locations ',     &
     &        'counter:',/)
  20  FORMAT (4x,'Grid: ',i0,', file: ',a,/,                            &
     &        19x,'LuvSrc = ',l1,2x,'u-face = ',i0,/,                   &
     &        19x,'LuvSrc = ',l1,2x,'v-face = ',i0,/,                   &
     &        19x,'LwSrc  = ',l1,2x,'w-face = ',i0)
  30  FORMAT (/,' CHECK_SOURCES - Cannot find point Souces/Sinks ',     &
     &        "the '",a,"' method.")
!
      RETURN
      END SUBROUTINE check_sources
!
      END MODULE mod_sources
