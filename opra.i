/* OPRA (OTF-based Phase Retrieval Analysis)
 * Copyright (c) Francois Rigaut 2010.
 * $Id: opra.i 9 2010-01-16 19:15:52Z frigaut $
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

//require,"yao.i";          // for zernikes
require,"opra_lmfit.i";        // for ...lmfit, nodep
require,"opra_libkl.i";        // for make_kl, nodep
require,"opra_libdh.i";        // for make_diskharmonic
// require,"yaodh.i";        // for make_diskharmonic
require,"opra_utils.i";   // plots, util functions.
require,"opra_structs.i"; // structures declarations

nmodes_max4printout = 20;


OPRA_VERSION = "1.0";

func opra(images, defocs, lambda, pixsize, teldiam, nmodes=, use_mode=, cobs=,
          noise=, pupd=, otf_dim=, progressive=, niter=, fix_amp=, fix_pix=,
          fix_kern=, fix_defoc=, first_nofit_astig=, winnum=, dpi=, pal=,
          yao_parfile=)
/* DOCUMENT opra(images,defocs,nmodes=,use_mode=)
   images:         Image data cube (data). These must be at least shannon sampled,
                   ideally 2x shannon or more.
   defocs:         Vector of estimate of defoc for each image (in radians)
   lambda:         Image wavelength (meter)
   pixsize:        Image pixel size (arcsec)
   teldiam:        Telescope (optics) diameter (meter)
   nmodes=         Number of modes to include in the phase estimate (default 120)
   use_mode=       Can be set to 'zernike' or 'dh' of 'kl' modes (default DH).
                   This can also be set to 'yao' for use of yao machinery
   cobs=           Optics central obstruction, in unit of optics outer diameter
                   (i.e. 0.1 would mean the central obstruction diameter is
                   10% of the optics/pupil/telescope outer diameter). Default 0.
   noise=          rms of noise in input images (one number, it is assumed the noise
                   is the same in all images). This is just used for display, not
                   used in the phase estimation. Default 0.
   pupd=           Force pupil diameter in pixels
   otf_dim=        ?
   progressive=    Introduce modes slower than regular version (more steps)
   niter=          Set maximum number of iteration
   fix_amp=        1 if image intensity should not be a free parameter
   fix_pix=        1 if pixel size should not be a free parameter
   fix_kern=       1 if gaussian kernel size should not be a free parameter
   fix_defoc=      1 if defoc from image to image should not be a free parameter
                   Note that only the global multiplier to the provided defoc
                   values is ever adjusted.
   first_nofit_astig= 1 if astig should not belong to the initial run
   winnum=         Output graphical window number
   dpi=            Output graphical window dpi
   pal=            Output graphical window palette
   yao_parfile=    yao parameter file if using use_mode='yao'

   Note: lambda, pixsize and teldiam are only used to get a starting value
   for the pupil diameter (in pixels). Lambda is also used to convert modes
   coefficients from radians to nm.

   The input images should be prepared as follow:
     - Background must be subtracted (not super critical but helps)
     - Images must be approximately centered
     - Flux of all images should be approximately the same
     - Dimension of all images *have* to be the same.
     - It is suggested to threshold (clip(im,0.,)).

   SEE ALSO: opra_foo
 */
{
  extern aold, opra_iter, lmfititer;
  extern itvec,distvec;

  //#################################
  // Initialize some general variables
  //#################################
  if (!nmodes) nmodes = 120;
  if (winnum==[]) winnum=0;
  if (!noise) noise=0.;
  itvec = distvec = [];
  opra_iter = lmfititer = 0;
  passn = 1;
  tic;

  //#################################
  // create graphic window
  //#################################
  animate,0;
  if (window_exists(winnum)) winkill,winnum;
  if (!dpi) dpi=130;
  window,winnum,wait=1,style="opra.gs",dpi=dpi,width=0,height=0;
  if (pal) palette,pal;
  animate,0;
  plsys,1; limits,square=1;
  plsys,2; limits,square=1;
  plsys,3; limits,square=1;

  //#################################
  // fill structures and prepare data
  //#################################
  opp = oprapar_struct();

  // number of input images
  opp.nim = dimsof(images)(4);

  // side dimension of input images
  opp.im_dim = dimsof(images)(2);

  // type of modes used to build the  phase
  if (strpart(use_mode,1:2) == "ze") opp.modes_type = "zernike";
  else if (use_mode == "kl")  opp.modes_type = "kl";
  else if (use_mode == "yao") opp.modes_type = "yao";
  else if (use_mode == "dh")  opp.modes_type = "dh";
  else error,"use_mode not defined";

  // fwhm_estimate is used to compute a starting value for pupil
  // diameter (pupd, in pixels)
  fwhm_estimate = lambda/teldiam/4.848e-6/pixsize;
  opp.pupd = (opp.im_dim/2.)*(2./fwhm_estimate);

  // we impose pupd to be integer, and will compensate with the pix size:
  psize_corr = round(opp.pupd)/opp.pupd;
  opp.pupd = psize_corr*opp.pupd;

  if(pupd != []) opp.pupd = pupd;

  opp.otf_sdim = long(2*opp.pupd*1.1)/2*2;
  opp.otf_sdim = long(2*opp.pupd*1.02)/2*2;
  opp.otf_dim  = long(2^ceil(log(opp.otf_sdim)/log(2)));

  if(otf_dim != []) {
    opp.otf_dim = otf_dim;
    opp.otf_sdim = opp.otf_dim;
  }

  if (use_mode=="yao") {
    if (!yao_parfile) error,"yao parfile not defined";
    if (findfiles(yao_parfile)==[]) \
      error,swrite(format="Can't find yao parfile %s\n",yao_parfile);
    // else we are good to go.
    cwd = pwd();
    cd,"/home/frigaut/mcao/myst/yorick";
    require,"yao_mcao.i";
    cd,cwd;
    aoread,yao_parfile;
    sim.pupildiam = opp.pupd;
    aoinit,disp=0,clean=1;
    nmodes = sum(dm._nact);
  } else {
    require,"yao.i";
  }

  opp.coefs = &(array(double,nmodes-1));

  opp.cobs = (cobs?cobs:0.);
  op  = array(opra_struct,opp.nim);

  // normalize
  images = images/sum(images(,,1))(,,-);
  for (i=1;i<=opp.nim;i++) {
    op(i).psf_data  = &(images(,,i));
    op(i).otf_data  = &(get_mtf(images(,,i),opp.otf_sdim));
    op(i).delta_foc = defocs(i);
    op(i).noise     = noise;
  }

  data = (*op(1).otf_data)(,,,-);
  for (i=2;i<=opp.nim;i++) grow,data,(*op(i).otf_data)(,,,-);
  center = opp.otf_sdim/2+1;
  // get rid of (0,0) frequency information
  data(center,center,) = 0;
  // now data has dimension dim * dim * 2 (re,im) * nimages

  // Trick: to avoid having to put in extern, I have to pass
  // both opp and op to opra_foo (via x). I don't want to copy
  // the structure, but pass them as is. I can't not make a vector
  // of structure, so I have to pass the address and do an eq_nocopy
  // on the receiving end (in opra_foo() ).
  x = _(&opp,&op);

  aold = [];

  if (progressive) {
    // set up progression in # modes (geometric, factor 2):
    nmodesv = nm = nmodes;
    while ((nm=nm/2)>15) grow,nmodesv,nm;
    nmodesv = _(6,nmodesv(::-1));
    if (niter) {
      nitv    = array(niter,numberof(nmodesv));// 10iteration max (?)
      nitv(2) = niter*2;
    } else {
      nitv    = array(10,numberof(nmodesv));// 10iteration max (?)
      nitv(2) = 20;//we want more iterations on the low order aberrations
    }
  } else {
    nmodesv = _(6,nmodes);
    nitv    = array(10,numberof(nmodesv));// 10iteration max (?)
    nitv(2) = 30;
  }

  if (use_mode=="yao") {
    nmodesv = [nmodes,nmodes];
    nitv = [10,niter];
  }
  //#################################
  // Let's start !
  //#################################
  // See the parameters descriptions in DOCUMENT section of opra_foo.
  // initialize "a" (the coefficients to find)
  a = fresh_a(opp.nim,nmodesv(1));
  a.kernd = 0.;
  a.stfmaskd = -100.;
  *a.amps = 1.;

  //BN addition:
  //a.psize = pixsize;
  //*a.amps = array(1,dimsof(images)(4));
  //a.defoc_scaling = 1.;


  opp.action = swrite(format="Pass %d: masks + aberrations up to %d",\
                      passn++,nmodesv(1));

  // filter some parameters
  fit = fresh_a(opp.nim,nmodesv(1),val=1);
  fit.pupd = 0;
  fit.kernd = 0;
  fit.stfmaskd = 0;
  (*fit.coefs)(4-1) = 0; // let's not optimize focus at first...
  if (first_nofit_astig) {
    (*fit.coefs)(5-1) = 0; // and not the astigs either...
    (*fit.coefs)(6-1) = 0; //
  }

  if (fix_amp) *fit.amps = 0;
  if (fix_pix) fit.psize = 0;
  if (fix_kern) fit.kernd = 0;
  if (fix_defoc) fit.defoc_scaling = 0;

  // call lmfit
  res = lmfit_wrap(opra_foo,x,a,data,opp.nim,fit=fit,\
                   tol=1e-16,eps=0.01,itmax=nitv(1),aregul=0.);

  // increase nmodes gently
  for (n=2;n<=numberof(nmodesv);n++) {
    aold = sdimold = [];
    if (use_mode!="yao") a.coefs = &( _(*a.coefs,array(0,nmodesv(n)-nmodesv(n-1))) );
    a.psize = clip(a.psize,-10,10);
    //    *a.coefs = clip(*a.coefs,-10,10);
    opp.action = swrite(format="Pass %d: masks + aberrations up to %d",\
                        passn++,nmodesv(n));

    // filter some parameters
    fit = fresh_a(opp.nim,nmodesv(n),val=1);
    fit.pupd = 0;
    fit.stfmaskd = 0;

    if (fix_amp) *fit.amps = 0;
    if (fix_pix) fit.psize = 0;
    if (fix_kern) fit.kernd = 0;
    if (fix_defoc) fit.defoc_scaling = 0;

    //(*fit.coefs)(4-1) = 0; // and let's not optimize focus at all

    // call lmfit
    res = lmfit_wrap(opra_foo,x,a,data,opp.nim,fit=fit,\
                     tol=1e-7*1e-6,eps=0.01,itmax=nitv(n),aregul=0.);
    if (opra_debug)                                                     \
      fits_write,"phase.fits",*opp.phase * (*opp.pupi),overwrite=1;
  }
  /*
  opp.action = swrite(format="Pass %d: Finishing pass (aberrations up to %d)", \
                      passn++,nmodesv(0));
  res = lmfit_wrap(opra_foo,x,a,data,opp.nim,fit=fit,                   \
                   tol=1e-16,eps=0.01,itmax=nitv(0));
  */
  if (opra_debug)                                                       \
    fits_write,"phase.fits",*opp.phase * (*opp.pupi),overwrite=1;

  write,format="Elapsed time = %f sec\n",tac();
  animate,0;
  opp.coefs=&(*a.coefs);
  return opp;
}


func opra_foo(x,b)
/* DOCUMENT
   b(1)        = support diameter ( sup.diam = pupd+atan(b(1))*5 )
   b(2)        = gaussian kernel fwhm ( kernel = abs(atan(b(2))*10) )
   b(3)        = mask diameter in TF plane (due to source)
                 using a makepupil real (~ TF of apodized round source)
                 mask diameter = 2*pupd + atan(b(3))*20;
   b(4)        = foc to defoc defocus (defocus =  atan(b(4))*4 )
   b(5,6,..,n) = foc -> defoc differential TT (2 terms / image additional to
                 base image) i.e. regular 2 image problem -> 2 terms,
                 total of 3 images -> 4 terms, etc...
   b(n+1:n+nim)= sum(data_image)/sum(model_image) for each images
   b(n+nim+1:) = mode coefficients (starting at #2)
   SEE ALSO:
 */
{
  extern aold,sdimold;
  extern lmfititer,newiter,opra_iter;

  opra_iter++;
  if (aold==[]) {
    aold=opra_a_struct();
    aold.pupd = -9999;
    aold.kernd = -9999;
  }
  if (sdimold==[]) sdimold=0;

  // retrieve structures from x. Because it's the same address, modif
  // to opp and op will be passed to main opra() (magic).
  eq_nocopy,opp,*x(1);  // -> opp
  eq_nocopy,op,*x(2);   // -> op

  // base + TT of additional images + amplitudes of all images

  a = a_vec2struct(b,opp.nim);

  //FIXME:
  // limits a's:
  // This is a bit adhoc, although the defaults look fine. FIXME?
  a.pupd          = opp.pupd+atan(a.pupd)*5;
  a.kernd         = atan(a.kernd)*10;
  //if(a.kernd > 1.01) a.kernd = 1.0;
  a.stfmaskd      = 1./(0.05+atan(a.stfmaskd)/pi*0.0999);
  a.defoc_scaling = 1.+atan(a.defoc_scaling)*1.;//was *3 orginally
  //lets try something else:
  //if(a.defoc_scaling > 1.101) a.defoc_scaling = 1.1;
  //if(a.defoc_scaling < 0.899) a.defoc_scaling = 0.9;
  a.psize         = 1.+atan(a.psize)*0.2;
  //if(fix_amp)   *a.amps = 1.;//+atan(*a.amps)*0.6;
  //else  *a.amps = 1.+atan(*a.amps)*0.6;

  *a.amps = 1.+atan(*a.amps)*0.6;
  //if (*a.amps == 0)  *a.amps = 1.0;
  //*a.amps = *a.amps;
  //It could be good to bound properly the parameters.
  //But is that possible with lmfit ?

  nmodes = numberof(*a.coefs)+1;
  az = _(0.,*a.coefs); // az = coef, as noll. az(1) = always piston (kl/zern)
  //FIXME:

  // Damp coefficients
  // not used //for (i=2;i<=nmodes;i++) az(i) = atan(az(i))*2./zernumero(i)(1)^1.5;

  //for (i=2;i<=nmodes;i++) az(i) = atan(az(i))*3./zernumero(i)(1)^1.5;
  //BN commented for test

  //not used //for (i=2;i<=nmodes;i++) az(i) = az(i)/zernumero(i)(1)^2;

  if (aold.pupd!=a.pupd) {
    // pupd has changed, need to recompute pupils and possibly modes.
    // Why possibly? kl are computed on even int dim, so if pupd has not
    // gone above the next even int value, no need to recompute (this is
    // done in opra_gen_modes() ).
    offset   = ((opp.modes_type=="kl")? 0.5: 1);

    //BN 4 tests:
    offset = 0.5; // FR?

    //dh are like zernike, they need offset = 1

    if (use_mode=="yao") {
      write,"Generating yao modes";
      sim.pupildiam = opp.pupd;
      aoinit,disp=0,clean=1;
      opp.modes = &(array(0.0f,[3,opp.otf_dim,opp.otf_dim,nmodes])); // FIXME for multiple DMs
      (*opp.modes)(dm(1)._n1:dm(1)._n2,dm(1)._n1:dm(1)._n2,) = *dm(1)._def; // FIXME for multiple DMs
      opp.pupi = &ipupil;
      opp.pupr = &pupil;
    } else {
      write,"Generating modes";
      opp.modes = &( opra_gen_modes(opp.otf_dim,a.pupd,nmodes,\
                                    mode_type=opp.modes_type));
      opp.pupi = &( make_pupil(opp.otf_dim,a.pupd,                        \
                               xc=opp.otf_dim/2+offset,                   \
                               yc=opp.otf_dim/2+offset,                   \
                               cobs=opp.cobs) );

      opp.pupr = &( make_pupil(opp.otf_dim,a.pupd,real=1,                 \
                               xc=opp.otf_dim/2+offset,                   \
                               yc=opp.otf_dim/2+offset,                   \
                               cobs=opp.cobs) );
      //opp.pupr = opp.pupi;
    }
    prepzernike,opp.otf_dim,a.pupd,opp.otf_dim/2+offset,opp.otf_dim/2+offset;
    // above: needed anyway for added defocus
    // BUG: need to modify xc/yc as using kl???
  }

  // if source TF mask diameter has changed. Recompute:
  if (aold.stfmaskd!=a.stfmaskd) {
    // upgrade with Bessel function (FIXME)
    xy = dist(opp.otf_sdim)/a.stfmaskd;
    opp.stfmask = &( sinc(xy) );
  }

  // blur kernel:
  if (aold.kernd!=a.kernd) \
    //    opp.kernel = &( roll(abs(fft(makegaussian(opp.otf_sdim,a.kernd),1))) );
    opp.kernel = &( makegaussian(opp.otf_sdim,100./(1e-6+a.kernd^2.)) );

  // Build phase:
  if (*opp.phase==[]) opp.phase = &(array(float,[2,opp.otf_dim,opp.otf_dim]));
  else *opp.phase *= 0.;
  if (use_mode=="yao") {
    dm(1)._command = &(float(az));
    extern mircube;
    mircube = array(0.0f,[3,sim._size,sim._size,1]);
    mircube(dm(1)._n1:dm(1)._n2,dm(1)._n1:dm(1)._n2,1) = compDmShape(1,dm(1)._command);
    multWfs,1,disp=0;
    n12 = wfs(1).n12;
    (*opp.phase)(n12(1):n12(2),n12(1):n12(2)) = *wfs(1)._fimage;
  } else {
    for (i=2;i<=nmodes;i++) *opp.phase += az(i)*(*opp.modes)(,,i);
  }

  // Iteration accounting
  if (lmfititer==[]) lmfititer=0;
  if (newiter) lmfititer++;

  // build PSF, OTF
  all_otf = [];
  // w_infocus = where(abs(op.delta_foc)==min(abs(op.delta_foc)))(1);
  // loop on images (defocs)
  for (i=1;i<=opp.nim;i++) {
    // compute defoc to add to phase
    defoc = op(i).delta_foc *  a.defoc_scaling * zernike_ext(4);
    // compute tiptilt to add to phase
    if (i>=2) {
    // if (i!=w_infocus) {
      tt = (*a.diff_tt)(1,i-1) * zernike_ext(2) +       \
           (*a.diff_tt)(2,i-1) * zernike_ext(3);
    } else tt = 0; // no special TT passed for image#1 (TT in "phase")
    // compute complex wavefront
    phi     = array(complex,[2,opp.otf_dim,opp.otf_dim]);
    phi.re  = (*opp.pupr) * cos(*opp.phase + tt + defoc);
    phi.im  = (*opp.pupr) * sin(*opp.phase + tt + defoc);
    // compute PSF
    psf     = roll(abs(fft(phi,1))^2.);
    // normalize
    psf     = psf/sum(psf) * (*a.amps)(i);
    op(i).psf  = &psf;
    // Note that above, PSF does not reflect the true psf, as the source
    // mask effect is not included. But we don't need it at each opra_foo()
    // call, just for display, so I did it in opra_info_and_plots()
    // compute OTF, include attenuation by blur (kernel) and mask (source size)
    op(i).otf  = &( get_mtf(psf,opp.otf_sdim ) * \
                    ((*opp.kernel) * (*opp.stfmask))(,,-) );
    // Pixel size adjustement:
    center = opp.otf_sdim/2+1;
    xy = (indgen(opp.otf_sdim)-center)*a.psize+center;
    (*op(i).otf)(,,1) = bilinear((*op(i).otf)(,,1),xy,xy,grid=1);
    (*op(i).otf)(,,2) = bilinear((*op(i).otf)(,,2),xy,xy,grid=1);
    //    (*op(i).otf)(center,center,) = 0.;
    // append to previous OTF all object:
    grow,all_otf,(*op(i).otf)(,,,-);
  }

  if (newiter) opra_info_and_plots,a,*a.amps,az,nmodes,opp,op;

  aold = a;
  newiter=0;

  // get rid of (0,0) frequency information
  all_otf(center,center,) = 0;

  return all_otf;
}

