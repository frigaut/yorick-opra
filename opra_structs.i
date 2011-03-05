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


struct opra_struct     // struct to describe each input plane
{
  pointer psf_data;    // 2D image (data)
  pointer otf_data;    // 2D OTF (complex, re & im)
  pointer psf;         // model PSF
  pointer otf;         // 2D model OTF (complex, re & im)
  float   delta_tt(2); // tip-tilt coefs for this image
  float   delta_foc;   // defocus coef for this image
  float   amp;         // amplitude correction for this image
  float   noise;       // noise in data (just for cosmetic, not used in minim.)
};

struct oprapar_struct  // one structure to hold all important parameters and results.
{
  long    nim;         // number of images
  pointer coefs;       // modes coefficients
  string  modes_type;  // type of modes ("kl" or "zernikes")
  pointer phase;       // modelled phase map
  long    im_dim;      // images size
  long    otf_dim;     // OTF size for OTF calculations
  long    otf_sdim;    // OTF size of returned OTF
  float   pupd;        // pupil diameter
  pointer pupi;        // pupil (integer, i.e. 0/1)
  pointer pupr;        // pupil (real = apodized)
  float   cobs;        // telescope central obstruction (fraction of diameter)
  pointer modes;       // mode map data cube
  string  action;      // string for plots, current status/action
  pointer iter_v;      // iteration vector
  pointer dist_v;      // distance (data-model rms)
  pointer stfmask;     // source TF mask array
  pointer kernel;      // Blur kernel array
};

struct opra_a_struct   // holds variables passed to opra_foo()
{
  float   pupd; //=-9999
  float   kernd;
  float   stfmaskd;
  float   defoc_scaling;
  float   psize;
  pointer diff_tt;
  pointer amps;
  pointer coefs;
}

