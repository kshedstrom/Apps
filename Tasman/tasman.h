/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for Northeast Pacific (NEP5) simulation
*/

#undef NO_HIS
#undef NETCDF4
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

/* general */

#define CURVGRID
#define MASKING
#undef SOLVE3D
#ifdef SOLVE3D
# define NONLIN_EOS
# define SALINITY
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#undef STATIONS
#undef WET_DRY

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#undef RST_SINGLE
#define AVERAGES
#undef AVERAGES2
#ifdef SOLVE3D
# undef AVERAGES_DETIDE
# undef DIAGNOSTICS_TS
#endif
#undef DIAGNOSTICS_UV

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

#define UV_ADV
#define UV_COR
#undef UV_SADVECTION
#define UV_VIS2
#define UV_QDRAG

#ifdef SOLVE3D
# undef UV_SMAGORINSKY
# define VISC_3DCOEF
# define MIX_S_UV
# define VISC_GRID
# define SPONGE
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif


#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif


/* vertical mixing */

#ifdef SOLVE3D
# define SOLAR_SOURCE

# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  undef LMD_BKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  undef LMD_DDMIX
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
# define CORE_FORCING
# define BULK_FLUXES
# define CCSM_FLUXES
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  define DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE
#  undef LONGWAVE
# endif
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# define SRELAXATION
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
#endif

/* point sources (rivers, line sources) */

/* Using Runoff instead now */
#undef ANA_PSOURCE
#ifdef SOLVE3D
# define RUNOFF
#endif

/* tides */

#define LTIDES
#ifdef LTIDES
# undef FILTERED
# define SSH_TIDES
# define UV_TIDES
/* # ifdef SOLVE3D */
#  define ADD_FSOBC
#  define ADD_M2OBC
/* # endif */
# undef RAMP_TIDES
# undef TIDES_ASTRO
# undef POT_TIDES
# undef UV_DRAG_GRID
# undef ANA_DRAG
# define LIMIT_BSTRESS
#endif

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
# undef BULK_FLUXES2D
# undef CCSM_FLUXES2D
# define ANA_INITIAL
# define ANA_FSOBC
# define ANA_M2OBC
#endif
