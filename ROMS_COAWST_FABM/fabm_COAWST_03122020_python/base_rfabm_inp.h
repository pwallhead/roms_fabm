      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: npzd_Franks_inp.h 1001 2020-01-10 22:41:16Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in the tracer diffusivity, nudging, and output   !
!  parameters specified in input script "rfabm.in".                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status

      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc
      logical, dimension(NSAT,Ngrids) :: Lsat
      logical, dimension(NBAT,Ngrids) :: Lbat

      real(r8), dimension(NBT,Ngrids) :: Rbio
      real(r8), dimension(NSAT,Ngrids) :: Rsat
      real(r8), dimension(NBAT,Ngrids) :: Rbat

      real(dp), dimension(800) :: Rval ! Modified PWA 28/07/2019: Increased dimension from nRval (100) to 800 to allow for high complexity

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(800) :: Cval ! Modified PWA 28/07/2019: Increased dimension from nCval (200) to 800 to allow for high complexity
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter


!!!PWA Inserted 08/03/2017
!-----------------------------------------------------------------------
!  Read in tracer diffusivity and nudging parameters.
!-----------------------------------------------------------------------
!
! PWA: When FABM is activated, all BIOLOGICAL parameters have already been
! read into FMODELS(ng)%fmodel(tile) from the .yaml file when roms_fabm_init
! was called. Consequently, only Lbiology, BioIter and the diffusivity coefficients
! need to be read from model_name.in.
! Option 1: Read the biological parameter values again from .yaml file
! Option 2: Extract their values from FMODELS(ng)%fmodel(tile) - SEE BELOW commented out
! Option 3: Do not create biological parameters as internal ROMS variables (CHOSEN 21/07/2015)


!!!!PWA Inserted 21/07/2015
!      property => FMODELS(ng)%fmodel(tile)%root%parameters%first
!      DO WHILE (associated(property))
!        ! Each parameter ("property" object) has the following attributes:
!        ! - name: short name [string]
!        ! - long_name: descriptive name [string]
!        ! - units: units [string]
!        ! - value: parameter value [data type depends on type of "property" object]
!        ! - default: default parameter value provided by the BGC model that registered the parameters [data type depends on type of "property" object]
!        ! - has_default: whether a default value has been provided during registration [logical]
!        SELECT TYPE (property)
!        class is (type_real_property)
!          WRITE (*,*) property%name,property%long_name,property%value,property%units
!        class is (type_integer_property)
!        !  WRITE (*,*) property%name,property%long_name,property%value
!          Npts=load_i(Nval, Rval, Ngrids, BioIter)
!        class is (type_logical_property)
!          WRITE (*,*) property%name,property%long_name,property%value
!        class is (type_string_property)
!          WRITE (*,*) property%name,property%long_name,property%value
!        END SELECT
!        property => property%next
!      END DO
!!!!PWA Inserted 21/07/2015
!!!End insertion PWA 08/03/2017

!!!Deleted PWA 08/03/2017 #ifdef ANA_BIOLOGY ... endif

      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lbiology')
              Npts=load_l(Nval, Cval, Ngrids, Lbiology)
!!!Deleted PWA 08/03/2017 CASE ('BioIter') ... Npts=load_r(Nval, Rval, Ngrids, wDet)
!!!PWA Inserted 13/12/2020
            CASE ('nFABM')
                Npts=load_i(Nval, Rval, Ngrids, nFABM)
!!!End insert PWA 13/12/2020
!!!PWA Inserted 08/03/2017
            CASE ('icheckmax')
              Npts=load_i(Nval, Rval, Ngrids, icheckmax)
!Note: the function load_i expects real-valued input vector Rval (see inp_par.F)
            CASE ('dBdt1max')
              Npts=load_r(Nval, Rval, Ngrids, dBdt1max)
!            write(*,*) "Done input of icheckmax, dBdt1max"
!!!End insertion PWA 08/03/2017
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  nl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  nl_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu2(i,ng)=Rbio(itrc,ng)
                  tl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu4(i,ng)=Rbio(itrc,ng)
                  tl_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('LtracerSponge')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSponge(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Akt_bak(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_AKT_fac')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_Akt_fac(i,ng)=Rbio(itrc,ng)
                  tl_Akt_fac(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Tnudg(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('Hadvection')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idbio(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idbio(iTrcStr), idbio(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Hadvection)
            CASE ('Vadvection')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idbio(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idbio(iTrcStr), idbio(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Vadvection)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_Hadvection')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idbio(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idbio(iTrcStr), idbio(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       ad_Hadvection)
            CASE ('Vadvection')
              IF (itracer.lt.(NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idbio(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idbio(iTrcStr), idbio(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       ad_Vadvection)
#endif
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idbio(iTrcStr), idbio(iTrcEnd),             &
     &                      Vname(1,idTvar(idbio(itracer))), LBC)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idbio(iTrcStr), idbio(iTrcEnd),             &
     &                      Vname(1,idTvar(idbio(itracer))), ad_LBC)
#endif
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
!!!Inserted PWA 07/02/2020
#ifdef ROMS_ISOURCES
            CASE ('LtracerIbio')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  LtracerIbio(itrc,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif  
#ifdef TS_ISOURCE
            CASE ('LtracerISrc')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerISrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
!!!End insert PWA 07/02/2020
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('LnudgeTCLM')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LnudgeTCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTvar)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idTvar(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTsur)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idTsur(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idTvar)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  Qout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idsurT)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idsurT(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idsurT(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Qout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idTsur)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  Qout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
!!!PWA Inserted 25/04/2017
            CASE('Hout(idSAvar)')
              Npts=load_l(Nval, Cval, NSAT, Ngrids, Lsat)
              DO ng=1,Ngrids
                DO i=1,NSAT
                  IF (idSAvar(i).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idSAvar(', i, ')'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(idSAvar(i),ng)=Lsat(i,ng)
                END DO
              END DO
#ifdef DEBUG
              write(*,*) "Done input of surface attached Hout options"
#endif
            CASE('Hout(idBAvar)')
              Npts=load_l(Nval, Cval, NBAT, Ngrids, Lbat)
              DO ng=1,Ngrids
                DO i=1,NBAT
                  IF (idBAvar(i).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idBAvar(', i, ')'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(idBAvar(i),ng)=Lbat(i,ng)
                END DO
              END DO
#ifdef DEBUG
              write(*,*) "Done input of bottom attached Hout options"
#endif
!!!End PWA insert 25/04/2017

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            CASE ('Aout(idTvar)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idTTav)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUTav)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idUTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idVTav)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idVTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHUTav)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=iHUTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHVTav)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=iHVTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
!!!PWA Inserted 25/04/2017
            CASE ('Aout(idSAvar)')
              Npts=load_l(Nval, Cval, NSAT, Ngrids, Lsat)
              DO ng=1,Ngrids
                DO i=1,NSAT
                  Aout(idSAvar(i),ng)=Lsat(i,ng)
                END DO
              END DO
# ifdef DEBUG
              write(*,*) "Done input of surface attached Aout options"
# endif
            CASE ('Aout(idBAvar)')
              Npts=load_l(Nval, Cval, NBAT, Ngrids, Lbat)
              DO ng=1,Ngrids
                DO i=1,NBAT
                  Aout(idBAvar(i),ng)=Lbat(i,ng)
                END DO
              END DO
# ifdef DEBUG
              write(*,*) "Done input of bottom attached Aout options"
# endif
!!!End PWA insert 25/04/2017
#endif
#ifdef DIAGNOSTICS_TS
            CASE ('Dout(iTrate)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTrate),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iThadv)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iThadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTxadv)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTyadv)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTvadv)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltrc(i,ng)
                END DO
              END DO
# if defined TS_DIF2 || defined TS_DIF4
            CASE ('Dout(iThdif)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iThdif),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTxdif)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTydif)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTydif),ng)=Ltrc(i,ng)
                END DO
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            CASE ('Dout(iTsdif)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltrc(i,ng)
                END DO
              END DO
#  endif
# endif
            CASE ('Dout(iTvdif)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltrc(i,ng)
                END DO
              END DO
#endif

!!!FABM-model dependent code to be inserted HERE

          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,50) ng
!!!Deleted PWA 08/03/2017 WRITE (out,60) BioIter ... 'Detrital sinking rate (m/day).'
!!!PWA Inserted 08/03/2017
            WRITE (out,55) nFABM(ng), 'nFABM',                          &
     &            'Number of ROMS time steps (DT) per FABM time step.'
            WRITE (out,60) icheckmax(ng), 'icheckmax',                  &
     &            'Initial iterations for which FABM calls checked.'
            WRITE (out,70) dBdt1max(ng), 'dBdt1max',                    &
     &            'Threshold value of FABM dBio_dt.'
!!!End insertion PWA 08/03/2017
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu2(i,ng), 'nl_tnu2', i,               &
     &              'NLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu2(i,ng), 'ad_tnu2', i,               &
     &              'ADM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu2(i,ng), 'tl_tnu2', i,               &
     &              'TLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu4(i,ng), 'nl_tnu4', i,               &
     &              'NLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu4(i,ng), 'ad_tnu4', i,               &
     &              'ADM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu4(i,ng), 'tl_tnu4', i,               &
     &              'TLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,100) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef FORWARD_MIXING
            DO itrc=1,NBT
              i=idbio(itrc)
# ifdef ADJOINT
              WRITE (out,90) ad_Akt_fac(i,ng), 'ad_Akt_fac', i,         &
     &              'ADM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_Akt_fac(i,ng), 'tl_Akt_fac', i,         &
     &              'TLM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning ON  point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning OFF point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
!!!Inserted PWA 07/02/2020
# ifdef ROMS_ISOURCES
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerIbio(itrc,ng)) THEN
                WRITE (out,100) LtracerIbio(itrc,ng),'LtracerIbio',itrc,&
     &              'Turning ON internal sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerIbio(itrc,ng),'LtracerIbio',itrc,&
     &              'Turning OFF internal sources/Sink on tracer ',i,   &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
# endif
# ifdef TS_ISOURCE
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerISrc(i,ng)) THEN
                WRITE (out,100) LtracerISrc(i,ng), 'LtracerISrc',       &
     &              i, 'Turning ON internal sources/Sink on tracer ', i,&
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerISrc(i,ng), 'LtracerISrc',       &
     &              i, 'Turning OFF internal sources/Sink on tracer ',i,&
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
# endif
!!!End insert PWA 07/02/2020
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,100) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,100) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Hout(idTvar(i),ng)) WRITE (out,110)                 &
     &              Hout(idTvar(i),ng), 'Hout(idTvar)',                 &
     &              'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Hout(idTsur(i),ng)) WRITE (out,110)                 &
     &              Hout(idTsur(i),ng), 'Hout(idTsur)',                 &
     &              'Write out tracer flux ', i,                        &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
            END IF
            IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Qout(idTvar(i),ng)) WRITE (out,110)                 &
     &              Qout(idTvar(i),ng), 'Qout(idTvar)',                 &
     &              'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Qout(idsurT(i),ng)) WRITE (out,110)                 &
     &              Qout(idsurT(i),ng), 'Qout(idsurT)',                 &
     &              'Write out surface tracer ', i,                     &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Qout(idTsur(i),ng)) WRITE (out,110)                 &
     &              Qout(idTsur(i),ng), 'Qout(idTsur)',                 &
     &              'Write out tracer flux ', i,                        &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
            END IF
!!!PWA Inserted 28/04/2017
            DO i=1,NSAT
              IF (Hout(idSAvar(i),ng)) WRITE (out,110)                  &
     &            Hout(idSAvar(i),ng), 'Hout(idSAvar)',                 &
     &            'Write out surface-attached variable ', i,            &
     &            TRIM(Vname(1,idSAvar(i)))
            END DO
            DO i=1,NBAT
              IF (Hout(idBAvar(i),ng)) WRITE (out,110)                  &
     &            Hout(idBAvar(i),ng), 'Hout(idBAvar)',                 &
     &            'Write out bottom-attached variable ', i,             &
     &            TRIM(Vname(1,idBAvar(i)))
            END DO
!!!End insert PWA 28/04/2017
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            IF ((nAVG(ng).gt.0).and.ANY(Aout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Aout(idTvar(i),ng)) WRITE (out,110)                 &
     &              Aout(idTvar(i),ng), 'Aout(idTvar)',                 &
     &              'Write out averaged tracer ', i,                    &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Aout(idTTav(i),ng)) WRITE (out,110)                 &
     &              Aout(idTTav(i),ng), 'Aout(idTTav)',                 &
     &              'Write out averaged <t*t> for tracer ', i,          &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Aout(idUTav(i),ng)) WRITE (out,110)                 &
     &              Aout(idUTav(i),ng), 'Aout(idUTav)',                 &
     &              'Write out averaged <u*t> for tracer ', i,          &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Aout(idVTav(i),ng)) WRITE (out,110)                 &
     &              Aout(idVTav(i),ng), 'Aout(idVTav)',                 &
     &              'Write out averaged <v*t> for tracer ', i,          &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Aout(iHUTav(i),ng)) WRITE (out,110)                 &
     &              Aout(iHUTav(i),ng), 'Aout(iHUTav)',                 &
     &              'Write out averaged <Huon*t> for tracer ', i,       &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Aout(iHVTav(i),ng)) WRITE (out,110)                 &
     &              Aout(iHVTav(i),ng), 'Aout(iHVTav)',                 &
     &              'Write out averaged <Hvom*t> for tracer ', i,       &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
!!!Inserted PWA 28/04/2017
              DO i=1,NSAT
                IF (Aout(idSAvar(i),ng)) WRITE (out,110)                &
     &              Aout(idSAvar(i),ng), 'Aout(idSAvar)',               &
     &              'Write out averaged surface-attached variable ', i, &
     &              TRIM(Vname(1,idSAvar(i)))
              END DO
              DO i=1,NBAT
                IF (Aout(idBAvar(i),ng)) WRITE (out,110)                &
     &              Aout(idBAvar(i),ng), 'Aout(idBAvar)',               &
     &              'Write out averaged bottom-attached variable ', i,  &
     &              TRIM(Vname(1,idBAvar(i)))
              END DO
!!!End PWA insert 28/04/2017
            END IF
#endif
#ifdef DIAGNOSTICS_TS
            IF ((nDIA(ng).gt.0).and.ANY(Dout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTrate),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTrate)',               &
     &              'Write out rate of change of tracer ', itrc,        &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iThadv),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iThadv)',               &
     &              'Write out horizontal advection, tracer ', itrc,    &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTxadv),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTxadv)',               &
     &              'Write out horizontal X-advection, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTyadv),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTyadv)',               &
     &              'Write out horizontal Y-advection, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTvadv),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTvadv)',               &
     &              'Write out vertical advection, tracer ', itrc,      &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
# if defined TS_DIF2 || defined TS_DIF4
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iThdif),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iThdif)',               &
     &              'Write out horizontal diffusion, tracer ', itrc,    &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(i,iTxdif),ng))                          &
     &            WRITE (out,110) .TRUE., 'Dout(iTxdif)',               &
     &              'Write out horizontal X-diffusion, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTydif),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTydif)',               &
     &              'Write out horizontal Y-diffusion, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTsdif),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTsdif)',               &
     &              'Write out horizontal S-diffusion, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
#  endif
# endif
              DO i=1,NBT
                itrc=idbio(i)
                IF (Dout(idDtrc(itrc,iTvdif),ng))                       &
     &            WRITE (out,110) .TRUE., 'Dout(iTvdif)',               &
     &              'Write out vertical diffusion, tracer ', itrc,      &
     &              TRIM(Vname(1,idTvar(itrc)))
              END DO
            END IF
#endif
!!!PWA Inserted 18/04/2017
#ifdef DIAGNOSTICS_BIO
            IF ((nDIA(ng).gt.0).and.ANY(Dout(:,ng))) THEN
              DO i=1,NDbio2d
                IF (Dout(iDbio2(i),ng))                                 &
     &            WRITE (out,120) .TRUE., 'Dout(iDbio2)',               &
     &              'Write out horizontal diagnostic ', i,              &
     &              TRIM(Vname(1,iDbio2(i)))
              END DO
              DO i=1,NDbio3d
                IF (Dout(iDbio3(i),ng))                                 &
     &            WRITE (out,120) .TRUE., 'Dout(iDbio3)',               &
     &              'Write out bulk diagnostic ', i,                    &
     &              TRIM(Vname(1,iDbio3(i)))
              END DO
            END IF
!Note: Vectors iDbio2,3 are already the correct indices for Dout analogous to idTvar (not idbio)
#endif
!!!End PWA insert 18/04/2017
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters.
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          nl_tnu4(i,ng)=SQRT(ABS(nl_tnu4(i,ng)))
#ifdef ADJOINT
          ad_tnu4(i,ng)=SQRT(ABS(ad_tnu4(i,ng)))
#endif
#if defined TANGENT || defined TL_IOMS
          tl_tnu4(i,ng)=SQRT(ABS(tl_tnu4(i,ng)))
#endif
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)
  40  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' RFABM Model Parameters, Grid: ',i2.2,               &
     &        /,  ' ===============================',/)
  55  FORMAT (1x,i10,2x,a,t32,a)
  60  FORMAT (1x,i10,2x,a,t32,a)
  70  FORMAT (1p,e11.4,2x,a,t32,a)
  80  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 100  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 120  FORMAT (10x,l1,2x,a,t32,a,i4.4,':',1x,a)

      RETURN
      END SUBROUTINE read_BioPar
