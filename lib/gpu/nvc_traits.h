/***************************************************************************
                             nvc_texture_traits.h
                             -------------------
                               W. Michael Brown

  Tricks for templating textures

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Jun 23 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#ifndef NVC_TEXTURE_TRAITS_H
#define NVC_TEXTURE_TRAITS_H

template <class numtyp> class nvc_vec_traits;
template <> class nvc_vec_traits<float> { public: typedef float2 vec2; };
template <> class nvc_vec_traits<double> { public: typedef double2 vec2; };

#endif
