/* DOCUMENT opra_structs
 * OPRA_STRUCTS (OTF-based Phase Retrieval Analysis)
 * Copyright (c) Francois Rigaut 2010.
 * $Id$
 *
 * This file is part of the OPRA package.
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 */


struct opra_data_struct     // struct to describe each input plane
// most member xfered to fixed size arrays for svipc data xfer
// nim, npos, im_dim and otf_dim need to be defined before defining
// this structures (including this file).
{
  pointer  psf;           // pointer to PSF data cube (user data)
  pointer  otf;           // Pointer to OTF of user images
  pointer  psf_model;     // model PSF
  pointer  otf_model;     // 2D model OTF (complex, re & im)
  pointer  position;      // position for all psf in datacube
  pointer  focus;         // defocus coef for all psf in datacube
  pointer  amp;           // amplitude correction for all psf in datacube
  pointer  noise;         // noise in data for all psf in datacube
                          // (just for cosmetic, not used in minim.)
};

struct opra_fit_struct
// fit structure.
// all members booleans
// 0 means don't fit
// 1 means fit
// pointer members point to int vector of 0s and 1s
{
  pointer im;       // [nim elements] Use image in fit?
  pointer modes;    // [nmodes elements] Use mode in fit?
  pointer position; // [nim elements] fit position?
  int amp;          // fit image amplitude?
  int psize;        // fit pixel size?
  int kernel;       // fit gaussian kernel?
  int defocus;      // fit defocus scaling factor? (global scaling not each focus value)
  int cobs;         // fit central obstruction?
}

struct opras
// one structure to hold all important parameters and results.
{
  long    nmodes;            // Number of modes
  string  modes_type;        // type of modes ("kl", "zernikes", "dh" or "yao")
  float   teldiam;           // telescope diameter [m]
  float   cobs;              // telescope central obstruction [fraction of diameter]
  float   lambda;            // image wavelength [m]
  float   psize;             // pixel size [arcsec]
  float   kernd(3);          // kernel parameters (Xfwhm, Yfwhm, Pos angle)
  long    winnum;            // window # (or first one in serie for tomo mode)
  long    dpi                // window dpi for displays
  string  action;            // string for plots, current status/action
  long    ndm;               // number of DM (for yao mode)
  long    ncoef_per_dm(ndm); // vector that contains # of coefs (actuator/modes) per dm
  pointer phase;             // modelled phase map [rd]
  pointer coefs;             // modes coefficients [rd]
  opra_data_struct data;     // data
  opra_fit_struct  fit;      // fit variables (bool)

  // internal parameters
  long    _iter              // iteration from start
  float   _pupd;             // internal, pupil diameter
  long    _nim;              // internal, number of images
  long    _npos;             // number of positions for tomographic mode
  long    _dim;           // internal, images size
  long    _otf_sdim;         // internal, OTF size of returned OTF
  pointer _stfmask;          // internal, source TF mask array
  pointer _kernel ;          // internal, blur kernel array
  pointer _pupi;             // internal, pupil (integer, i.e. 0/1)
  pointer _pupr;             // internal, pupil (real = apodized)
  pointer _modes;            // internal, mode map data cube
};


