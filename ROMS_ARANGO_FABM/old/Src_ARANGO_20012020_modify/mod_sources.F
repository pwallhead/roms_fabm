#include "cppdefs.h"
      MODULE mod_sources
!
!svn $Id: mod_sources.F 995 2020-01-10 04:01:28Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
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

        implicit none

        integer, allocatable :: Msrc(:)
        integer, allocatable :: Nsrc(:)
!!!Inserted PWA 18/10/2021
#ifdef TS_ISOURCE
        integer, allocatable :: Nisrc(:)
#endif
!!!End insert PWA 18/10/2021

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

        END TYPE T_SOURCES

        TYPE (T_SOURCES), allocatable :: SOURCES(:)

      CONTAINS

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

      integer :: Vid, ifile, nvatt, nvdim
#endif
      integer :: is, itrc, k, mg
#ifdef TS_ISOURCE
      integer :: Ntiles 
#endif

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) THEN
        allocate ( Msrc(Ngrids) )
        allocate ( Nsrc(Ngrids) )
      END IF

#ifndef ANA_PSOURCE
!
!  Inquire about the number of point sources.
!
      IF (ng.eq.1) THEN
        DO mg=1,Ngrids
          foundit=.FALSE.
          IF (LuvSrc(mg).or.LwSrc(mg).or.ANY(LtracerSrc(:,mg))) THEN
!!!Modified PWA 20/11/2019 following Hedstrom 14092019 code
            DO i=1,nSSFfiles(ng)
              CALL netcdf_inq_var (ng, iNLM, SSF(i,mg)%name,            &
     &                             MyVarName = TRIM(Vname(1,idRxpo)),   &
     &                             SearchVar = foundit,                 &
     &                             VarID = Vid,                         &
     &                             nVardim = nvdim,                     &
     &                             nVarAtt = nvatt)
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__//", allocate_sources")) RETURN

              IF (foundit) THEN
                Nsrc(mg)=var_Dsize(1)         ! first dimension
                Msrc(mg)=Nsrc(mg)
              END IF
            END DO
!!!End modified PWA 20/11/2019
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
            DO i=1,nISSFfiles(ng)
              CALL netcdf_inq_var (ng, iNLM, ISSF(i,mg)%name,           &
     &                             MyVarName = TRIM(Vname(1,idIxpo)),   &
     &                             SearchVar = foundit,                 &
     &                             VarID = Vid,                         &
     &                             nVardim = nvdim,                     &
     &                             nVarAtt = nvatt)
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__//", allocate_sources")) RETURN

              IF (foundit) THEN
                Nisrc(mg)=var_Dsize(1)         ! first dimension
              END IF
            END DO
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
      END DO
      DO k=1,N(ng)
        DO is=1,Nsrc(ng)
          SOURCES(ng) % Qshape(is,k) = IniVal
          SOURCES(ng) % Qsrc(is,k) = IniVal
        END DO
      END DO
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO is=1,Nsrc(ng)
            SOURCES(ng) % Tsrc(is,k,itrc) = IniVal
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

      RETURN
      END SUBROUTINE allocate_sources
      END MODULE mod_sources
