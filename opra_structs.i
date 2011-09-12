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
// most member xfered to fixed size arrays for svipc data xfer
// nim, npos, im_dim and otf_dim need to be defined before defining
// this structures (including this file).
{
  float  psf_data(im_dim,im_dim);   // 2D image (data)
  float  otf_data(otf_sdim,otf_sdim,2); // 2D OTF (complex, re & im)
  float  psf(im_dim,im_dim);        // model PSF
  float  otf(otf_sdim,otf_sdim,2);      // 2D model OTF (complex, re & im)
  float  delta_tt(2);     // tip-tilt coefs for this image
  float  delta_foc;       // defocus coef for this image
  float  amp;             // amplitude correction for this image
  float  noise;           // noise in data (just for cosmetic, not used in minim.)
};

struct oprapar_struct
// one structure to hold all important parameters and results.
// most member xfered to fixed size arrays for svipc data xfer
// ndm, otf_dim and otf_sdim need to be defined before defining
// this structures (including this file).
{
  long    nim;                       // number of images
  long    npos;                      // number of positions
  long    ndm;                       // number of DM (for yao mode)
  long    ncoefs;                    // total number of coefs
  long    im_dim;                    // images size
  long    otf_dim;                   // OTF size for OTF calculations
  long    otf_sdim;                  // OTF size of returned OTF
  float   pupd;                      // pupil diameter
  float   cobs;                      // telescope central obstruction (fraction of diameter)
  float   lambda;                    // image wavelength
  long    ncoef_per_dm(ndm);         // vector that contains # of coefs (actuator/modes) per dm
  float   phase(otf_dim,otf_dim,npos); // modelled phase map
  float   pupi(otf_dim,otf_dim);     // pupil (integer, i.e. 0/1)
  float   pupr(otf_dim,otf_dim);     // pupil (real = apodized)
  float   stfmask(otf_sdim,otf_sdim);// source TF mask array
  float   kernel(otf_sdim,otf_sdim); // Blur kernel array
  pointer coefs;                     // modes coefficients
  pointer modes;                     // mode map data cube
  long    winnum;                    // first window # in serie
  string  modes_type;                // type of modes ("kl", "zernikes", "dh" or "yao")
  string  action;                    // string for plots, current status/action
};

struct opra_a_struct   // holds variables passed to opra_foo()
{
  float   pupd;
  float   kernd(3);
  float   stfmaskd;
  float   defoc_scaling;
  float   psize;
  pointer diff_tt;
  pointer amps;
  pointer coefs;
}

