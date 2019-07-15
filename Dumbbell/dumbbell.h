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
**  Tidal half pipe.
*/

#undef STATIONS
#undef FLOATS
#undef  DIAGNOSTICS_UV

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define UV_LDRAG

#define ANA_GRID
#define ANA_MASK
#define ANA_INITIAL
#define ANA_SMFLUX
#define IMPLICIT_NUDGING

#define SOLVE3D
#define SALINITY
#define PERFECT_RESTART
#define MASKING
#define RI_SPLINES
#define NO_WRITE_GRID
#define AVERAGES
#define UV_LDRAG

#ifdef SOLVE3D
# define DJ_GRADPS
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
#endif

#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif
