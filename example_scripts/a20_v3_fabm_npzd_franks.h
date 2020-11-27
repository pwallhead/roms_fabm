/*
** svn $Id: a20.h 172 2015-12-31 01:45:48Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for NS8KM
**
** Application flag:   NS8KM
** Input script:       ocean_ns8km.in
**
** PWA 30/11/2019: Copied from a20_v2e.h
**                 Defined NO_HIS
** PWA 28/10/2020: Name changes for RFABM cpps ROMS_*->FABM_*
**
*/
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */
#define STATIONS

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state tracjectory.
**-----------------------------------------------------------------------------
*/

#if defined NLM_DRIVER
#define NO_HIS
#define UV_ADV
#define TS_C4VADVECTION
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#undef TS_DIF2
#undef TS_MIX_STABILITY
#define MIX_S_UV
#define MIX_S_TS
#define TS_U3HADVECTION
#define SOLVE3D
#undef TCLM_NUDGING
#undef ANA_NUDGCOEF
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define POWER_LAW
#define MASKING
#define AVERAGES
#define SOLAR_SOURCE
#undef SRELAXATION
#define MY25_MIXING
#define DEFLATE
#define CHARNOK /*Charnok Surface Roughness From Wind Stress */

#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# undef K_C4ADVECTION
# define N2S2_HORAVG /*Horizontal Smoothing of Buoyancy/Shea */
#endif

#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif

#ifdef LMD_MIXING
# define N2S2_HORAVG /*Horizontal Smoothing of Buoyancy/Shea */
# define KANTHA_CLAYSON
# define  LMD_RIMIX       /* Add diffusivity due to shear instability */
# define  LMD_CONVEC      /* Add convective mixing due to shear instability */
# define  LMD_DDMIX       /* Add double-diffusive mixing */
#endif

#define ANA_BSFLUX
#define ANA_BTFLUX
#undef ANA_SSFLUX
#undef ANA_STFLUX
#undef ANA_SMFLUX
#undef FORWARD_MIXING
#undef FORWARD_WRITE
#endif


#undef VISC_GRID

/* ICE MODEL */
#define ICE_MODEL         /* Turn on ice model */
#ifdef ICE_MODEL
# define ICE_THERMO
# define ICE_MK
# define CCSM_ICE_SHORTWAVE
# define ICE_I_O
# define ALBEDO_CSIM
# define ICE_SHORTWAVE_R
# undef ICE_ALB_EC92
# define ICE_MOMENTUM
# define ICE_BULK_FLUXES
# undef  ICE_MOM_BULK
# define ICE_EVP
# define ICE_ADVECT
# define ICE_SMOLAR
# define ICE_UPWIND
#endif

/* ATMOSPHERIC FORCING */
#define BULK_FLUXES        /* turn ON or OFF bulk fluxes computation */
#ifdef BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# undef  ANA_RAIN          /* analytical rain fall rate */
# undef  ANA_PAIR          /* analytical surface air pressure */
# undef  ANA_HUMIDITY      /* analytical surface air humidity */
# undef  ANA_CLOUD         /* analytical cloud fraction */
# undef  ANA_TAIR          /* analytical surface air temperature */
# undef  ANA_WINDS         /* analytical surface winds */
# define EMINUSP           /* turn ON internal calculation of E-P */
# define ALBEDO            /* use albedo equation for shortwave radiation */
# define ALBEDO_DIRDIFF    /* use direct vs. diffuse to set shortwave albedo */
# define SHORTWAVE
# define LONGWAVE_OUT
# define DIURNAL_SRFLUX
# define COOL_SKIN         /* turn ON or OFF cool skin correction *//* Ikke def hos Frode*/
#endif

#define ATM_PRESS          /* use to impose atmospheric pressure onto sea surface */
#define SOLAR_SOURCE       /* define solar radiation source term */
#define SPECIFIC_HUMIDITY  /* if input is specific humidity in kg/kg */

/* TIDES */
#define A20TIDES
#if defined A20TIDES
#define SSH_TIDES          /* turn on computation of tidal elevation */
#define UV_TIDES           /* turn on computation of tidal currents */
#define ADD_FSOBC          /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC          /* Add tidal currents  to processed OBC data */
#define RAMP_TIDES         /* Spin up tidal forcing */
#endif
#define WET_DRY

#undef INLINE_2DIO

/* FABM */
#define RFABM
#ifdef RFABM
# define BIOLOGY            /* MUST be defined for use of FABM */
# undef DEBUGFABM           /* use for RFABM-specific debugging output */
# undef DIAGNOSTICS         /* MUST be defined for outputting FABM diagnostics (specified in rfabm.in) */
# undef DIAGNOSTICS_BIO     /* MUST be defined for outputting FABM diagnostics (specified in rfabm.in) */
# define MASKING            /* MUST be defined to provide input to fabm_set_mask */
# define ANA_SPFLUX         /* MUST be defined, or surface bgc fluxes provided in forcing files */
# define ANA_BPFLUX         /* MUST be defined, or bottom bgc fluxes provided in forcing files */
# define SHORTWAVE          /* MAY be required in order to provide light forcing for FABM model */
# undef FABM_PCO2ATM        /* use to provide atmospheric pCO2 forcing via ROMS (need to input xCO2atm) */
# undef FABM_N3ATMDEPO      /* use to provide atmospheric deposition flux of oxidized nitrogen via ROMS (input N3atmd) */
# undef FABM_N4ATMDEPO      /* use to provide atmospheric deposition flux of reduced nitrogen via ROMS (input N4atmd) */
# undef FABM_AICE           /* use to provide fractional ice area from ROMS internal ice model */
# undef FABM_NONNEG_S       /* use to cap salinity input to FABM at zero PSU */
# undef FABM_CHECK_STATE    /* use to cap bgc variable input to FABM */
# define FABM_INITIAL       /* use to set all bgc initial conditions to FABM default values */
# undef FABM_INITIAL_SB     /* use to set initial conditions to FABM defaults only for surface/bottom attached variables */
#endif
/* --------------------------- */
