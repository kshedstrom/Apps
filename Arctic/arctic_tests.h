/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2017 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for ARCTIC simulation
*/

#define ROMS_MODEL
#undef NO_HIS
#define GLOBAL_PERIODIC
#define HDF5
#define DEFLATE
#undef PARALLEL_IN
#undef PARALLEL_OUT
#define PERFECT_RESTART

/* general */

#define MASK_HACK
#define CURVGRID
#define MASKING
#define SOLVE3D
#ifdef SOLVE3D
# define SALINITY
# define NONLIN_EOS
# undef SPLINES_VDIFF
# undef SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#define STATIONS
#define WET_DRY
#define IMPLICIT_NUDGING

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX
# define ANA_PASSIVE
# define TRC_PSOURCE
# define ANA_TRC_PSOURCE
# define AGE_MEAN
#endif

/* ice */

#ifdef SOLVE3D
# undef CICE_MODEL
# ifdef CICE_MODEL
#  define SNOWFALL
#  undef SNOW_FROM_RAIN
#  define INI_GLORYS_ICE
#  undef ICE_LOG_LAYER
# endif

# undef  ICE_MODEL
# ifdef ICE_MODEL
#  define ANA_ICE
#  define INI_GLORYS_ICE
#  define SNOWFALL
#  define OUTFLOW_MASK
#  define ICE_LANDFAST
#  define ICE_SHALLOW_LIMIT
#  define ICE_THERMO
#  define ICE_MK
#  define ICE_MOMENTUM
#  define ICE_MOM_BULK
#  define ICE_EVP
#  define ICE_STRENGTH_QUAD
#  define ICE_ADVECT
#  define ICE_SMOLAR
#  define ICE_UPWIND
#  define ICE_BULK_FLUXES
#  define ICE_CONVSNOW
#  define ICE_I_O
#  define ICE_DIAGS
#  define CCSM_ICE_SHORTWAVE
# endif
#endif

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#ifndef PERFECT_RESTART
# define RST_SINGLE
#endif
#undef AVERAGES
#ifdef SOLVE3D
# undef AVERAGES2
# undef DIAGNOSTICS_TS
#else
# define DIAGNOSTICS_UV
#endif

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

#define UV_ADV
#define UV_COR

#define UV_VIS2
#undef UV_SMAGORINSKY
#undef VISC_3DCOEF
#define MIX_S_UV
#define VISC_GRID

#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif

/* vertical mixing */

#ifdef SOLVE3D
# define WTYPE_GRID

# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  undef LI_FOX_KEMPER
#  undef LMD_BKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  define LMD_DDMIX
#  define LIMIT_VDIFF
# endif

# undef GLS_MIXING
# undef MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
#endif

/* surface forcing */

#ifdef SOLVE3D
# undef CORE_FORCING
# undef BULK_FLUXES
# undef CCSM_FLUXES
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define SOLAR_SOURCE
#  define EMINUSP
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE  /* for water */
#  undef ICE_ALB_EC92   /* for ice */
#  define ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE    /* for both */
#  undef ALBEDO_DIRDIFF
#  define ICE_SHORTWAVE_R
#  undef LONGWAVE
# endif
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# define SCORRECTION
# define NO_SCORRECTION_ICE
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# define ANA_NUDGCOEF
#else
# define ANA_INITIAL
# define ANA_NUDGCOEF
#endif

/* point sources (rivers, line sources) */

#ifdef SOLVE3D
# undef RUNOFF
# define TWO_D_TRACER_SOURCE
#endif

/* tides */

#ifdef SOLVE3D
# undef LTIDES
#endif
#ifdef LTIDES
# if defined AVERAGES && !defined USE_DEBUG
#  define FILTERED
# endif
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES
#endif

#define UV_DRAG_GRID
#define ANA_DRAG
#define LIMIT_BSTRESS
#define UV_QDRAG

/* Boundary conditions...careful with grid orientation */

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
# define ANA_SMFLUX
#else
# undef ANA_SMFLUX
# define BULK_FLUXES2D
# undef CCSM_FLUXES2D
#endif

/*
**  Biological model options.
*/
#ifdef SOLVE3D
# undef BIO_COBALT
#endif

#ifdef BIO_COBALT
# define TS_MPDATA
#elif defined SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
#endif

/* #define DEBUG_COBALT */
/*#define COBALT_CONSERVATION_TEST */
/*#define COBALT_NOSOURCE */
/*#define COBALT_DO_NOTHING  */

#if defined BIO_COBALT
# undef FILTERED
# undef AVERAGES2
# define OPTIC_MANIZZA
# define COBALT_MINERALS
# undef COBALT_PHOSPHORUS
# define COBALT_IRON
# define NO_IRON_COAST
# define COBALT_CARBON
# define DIAGNOSTICS_BIO
/*# define BENTHIC  */
/*# define TIMESERIES */
# undef ANA_ALK_DIC
# undef ANA_BIOLOGY        /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define COASTDIAT
#endif
