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

nmodes_max4printout = 30;

OPRA_VERSION = "1.9.9";

func opra(images, defocs, lambda, pixsize, teldiam,
          nmodes=, use_mode=, cobs=, noise=, pupd=, progressive=, niter=, kernd=,
          fix_amp=, fix_pix=, fix_kern=, fix_defoc=, fix_diff_tt=, fix_cobs=,
          first_nofit_astig=, winnum=, dpi=, pal=, gui=, svipc=, nm=)
/* DOCUMENT opra(images, defocs, lambda, pixsize, teldiam,
                 nmodes=, use_mode=, cobs=, noise=, pupd=, progressive=, niter=,
                 fix_amp=, fix_pix=, fix_kern=, fix_defoc=, fix_diff_tt=, first_nofit_astig=,
                 winnum=, dpi=, pal=, gui=, svipc=, nm=)
   images:         Image data cube (data). These must be at least Shannon (Nyquist)
                   sampled, ideally 2x Shannon or more.
   defocs:         Vector of estimate of defocus for each image (in radians)
   lambda:         Image wavelength (meter)
   pixsize:        Image pixel size (arcsec)
   teldiam:        Telescope (optics) diameter (meter)
   nmodes=         Number of modes to include in the phase estimate (default 120)
   use_mode=       Can be set to 'zernike' or 'dh' of 'kl' modes (default DH).
                   This can also be set to 'yao' for use of yao machinery (kind of
                   advanced mode and poorly documented as of Sep2011)
   cobs=           Optics central obstruction, in unit of optics outer diameter
                   (i.e. 0.1 would mean the central obstruction diameter is
                   10% of the optics/pupil/telescope outer diameter). Default 0.
   noise=          rms of noise in input images (one number, it is assumed the noise
                   is the same in all images). This is just used for display, not
                   used in the phase estimation. Default 0.
   pupd=           Force pupil diameter in pixels
   progressive=    Introduce modes slower than regular version (more steps)
   niter=          Set maximum number of iteration
   fix_amp=        1 if image intensity should not be a free parameter
   fix_pix=        1 if pixel size should not be a free parameter
   fix_kern=       1 if gaussian kernel size should not be a free parameter
   fix_defoc=      1 if defoc from image to image should not be a free parameter
                   Note that only the global multiplier to the provided defoc
                   values is ever adjusted.
   fix_cobs=       1 if cobs should not be a free parameter
   fix_diff_tt=    1 if differential TT between images shoud not be a free parameter
   first_nofit_astig= 1 if astig should not belong to the initial run
   winnum=         Output graphical window number
   dpi=            Output graphical window dpi
   pal=            Output graphical window palette
   svipc=          Force svipc on (1) or off (0)
   nm=             Printout mode coefficients in nm instead of mrd.

   Note: lambda, pixsize and teldiam are only used to get a starting value
   for the pupil diameter (in pixels). Lambda is also used to convert modes
   coefficients from radians to nm if requested (nm=1).

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
  extern opp,op;
  extern has_svipc;
  extern want_nanometer;

  // NO, svipc is too fragile. Let's not impose it by default:
  if (svipc==[]) svipc=has_svipc=0;
  else if (svipc==0) has_svipc=0;
  else {
    if (svipc!=[]) has_svipc=svipc; \
    else if ((find_in_path("svipc.i")!=[])&&(gui)) has_svipc=1;
    if (find_in_path("svipc.i")==[]) has_svipc=0;
  }

  //#################################
  // Initialize some general variables
  //#################################
  if (!nmodes)    nmodes = 120;
  if (winnum==[]) winnum = 0;
  if (!noise)     noise = 0.;
  if (!dpi)       dpi = 130;
  if (nm)         want_nanometer=1;
  if (progressive!=0) progressive=1; // default to progressive
  itvec = distvec = [];
  opra_iter = lmfititer = 0;
  passn = 1;
  tic;

  //#################################
  // fill structures and prepare data
  //#################################
  // we need images to be 4-dimensions:
  // im_dimx x imdimy x number of defocused images x number of positions
  if (dimsof(images)(1)==3) images = images(,,,-);
  if (dimsof(images)(2)!=dimsof(images)(3)) error,"images should be square";

  // Normalize so that regularization and stop criteria are
  // independent of the data:
  images = images/sum(images(*));

  // define parameters for structures
  if (use_mode=="yao") {
    if (sim==[]) error,"aoread() not done?";
    // ndm has been defined by the call to aoread()
    ndm = numberof(dm);
  } else ndm = 1;

  im_dim   = dimsof(images)(2);

  // fwhm_estimate is used to compute a starting value for pupil
  // diameter (pupd, in pixels)
  fwhm_estimate = lambda/teldiam/4.848e-6/pixsize;
  pupdiam = (im_dim/2.)*(2./fwhm_estimate);

  if (!pupd) {
    psize_corr = round(pupdiam)/pupdiam;
    pupd = psize_corr*pupdiam;
  } else psize_corr = 1.;

  otf_sdim = long(2*pupd*1.02)/2*2;
  //  otf_dim  = long(2^ceil(log(otf_sdim)/log(2)));
  otf_dim  = dimsof(images)(2);
  nim      = dimsof(images)(4);  // number of input images
  npos     = dimsof(images)(5);  // number of input positions

  require,"opra_structs.i"; // structures declarations
  opp = oprapar_struct();
  op  = array(opra_struct,[2,nim,npos]);

  opp.pupd     = pupd;
  opp.ndm      = ndm;
  opp.nim      = nim;
  opp.npos     = npos;
  opp.psize    = pixsize;
  opp.otf_dim  = otf_dim;
  opp.otf_sdim = otf_sdim;
  opp.im_dim   = im_dim;
  opp.winnum   = winnum;
  opp.cobs     = (cobs?cobs:0.);
  opp.lambda   = lambda;
  if (kernd!=[]) opp.kernd = kernd;

  // type of modes used to build the  phase
  if (strpart(use_mode,1:2) == "ze") opp.modes_type = "zernike";
  else if (use_mode == "kl")  opp.modes_type = "kl";
  else if (use_mode == "yao") opp.modes_type = "yao";
  else if (use_mode == "dh")  opp.modes_type = "dh";
  else error,"use_mode not defined";


  if (use_mode=="yao") {
    sim.pupildiam = opp.pupd;
    aoinit,disp=0,clean=1;
    opp.ncoefs = sum(dm._nact);
    opp.ncoef_per_dm = dm._nact;
    opp.ndm = ndm;
    dpi2 = long(dpi*ndm/3.);
  } else {
    require,"yao.i";
    opp.ncoefs = nmodes;
    opp.ncoef_per_dm = [opp.ncoefs];
  }

  // svipc and spawn GUI:
  if (has_svipc) {
    require,"opra_svipc.i";
    lmfititer_pass=lmfit_itmax=lmfititer=0;
    status = opra_svipc_init(dpi);
    if (wfs==[]) wfs=dm=[0];
    write_yao_struct_to_shm,wfs,dm;
    write_opra_struct_to_shm,fresh_a(opp),op,opp;
  }

  //##########################################
  // create graphic windows if not in gui mode
  //##########################################
  if (!has_svipc) {  // if has_svipc then display handled by fork()
    write,format="%s\n","Creating graphical windows";
    for (n=1;n<=opp.npos;n++) {
      winn = winnum+n-1;
      if (window_exists(winn)) winkill,winn;
      window,winn,wait=1,style="opra.gs",dpi=dpi,width=0,height=0;
      if (pal) palette,pal;
      for (i=1;i<=3;i++) { plsys,i; limits,square=1; }
    }
    if (use_mode=="yao") {
      if (window_exists(winnum+opp.npos)) winkill,winnum+opp.npos;
      window,winnum+opp.npos,style="nobox.gs",width=long(6.*dpi2/ndm),height=6*dpi2,dpi=dpi2,wait=1;
    }
    window,winnum;
  } else {
    write,format="%s\n","Waiting for semaphore 0 from child";
    sem_take,semkey,0; // wait for graphic window to be realized.
  }

  // normalize
  for (i=1;i<=opp.npos;i++) images(,,,i) = images(,,,i)/sum(images(,,1,i));
  center = opp.otf_sdim/2+1;
  data = [];

  for (n=1;n<=opp.npos;n++) {
    for (i=1;i<=opp.nim;i++) {
      op(i,n).psf_data  = images(,,i,n);
      op(i,n).otf_data  = get_mtf(images(,,i,n),opp.otf_sdim);
      op(i,n).otf_data(center,center,)  = 0;
      // get rid of (0,0) frequency information
      // grow,data,(*op(i,n).otf_data)(,,,-);
      op(i,n).delta_foc = defocs(i);
      op(i,n).noise     = noise;
    }
  }
  disk = dist(opp.otf_sdim)<(opp.otf_sdim/2.+1);
  xy = indices(opp.otf_sdim);
  extern wotf;
  wotf = where( (xy(,,1)>=center) & (disk) );
  data = op.otf_data(*,,,)(wotf,,,);
  // data(center,center,,) = 0;

  data = float(data(*)); // float to save RAM
  if (opp.npos>1) grow,data,array(0.0f,opp.ncoefs);

  // Trick: to avoid having to put in extern, I have to pass
  // both opp and op to opra_foo (via x). I don't want to copy
  // the structure, but pass them as is. I can't not make a vector
  // of structure, so I have to pass the address and do an eq_nocopy
  // on the receiving end (in opra_foo() ).
  x = _(&opp,&op);

  aold = [];

  if (progressive) {
    // set up progression in # modes (geometric, factor 2):
    nmodesv = nm = opp.ncoefs;
    if (use_mode=="yao") nmodesv = nm = max(opp.ncoef_per_dm)
    while ((nm=nm/2)>15) grow,nmodesv,nm;
    if (is_scalar(nmodesv)) nmodesv=[nmodesv];
    nmodesv = _(6,nmodesv(::-1));
    if (niter) {
      nitv    = array(niter,numberof(nmodesv));// 10iteration max (?)
      nitv(2) = niter*2;
    } else {
      nitv    = array(10,numberof(nmodesv));// 10iteration max (?)
      nitv(2) = 20;//we want more iterations on the low order aberrations
    }
  } else {
    nmodesv = _(6,opp.ncoefs);
    nitv    = array(10,numberof(nmodesv));// 10iteration max (?)
    if (niter) nitv(2) = niter; else nitv(2) = 30;
  }

  if (use_mode=="yao") { // Can't do progressive with yao (see comment below)
    if (!(allof(dm.type=="dh")||allof(dm.type=="kl")||allof(dm.type=="zernikes"))) {
      nmodesv = [opp.ncoefs,opp.ncoefs];
      // above: has to be, we don't know "modes" are modal (could be zonal)
      nitv = [6,niter];
    }
  }

  //#################################
  // Let's start !
  //#################################
  // See the parameters descriptions in DOCUMENT section of opra_foo.
  // initialize "a" (the coefficients to find)
  a = fresh_a(opp);
  a.kernd    = opp.kernd;
  a.stfmaskd = -100.;
  opp.action = swrite(format="Pass %d: masks + aberrations up to %d",\
                      passn++,nmodesv(1));

  // filter some parameters. Use fit keyword of lmfit.
  fit = fresh_a(opp,val=1);
  fit.pupd     = 0;
  fit.cobs     = 0;
  fit.kernd    = [1,1,1]*0;
  fit.stfmaskd = 0;
  if (fix_amp)   *fit.amps = 0;
  if (fix_pix)   fit.psize = 0;
  if (fix_kern)  fit.kernd = [0,0,0];
  if (fix_defoc) fit.defoc_scaling = 0;
  if (fix_cobs)  fit.cobs = 0;
  // We want to peg first in-focus image:
  (*fit.diff_tt)(,1,1) = 0;
  if (fix_diff_tt) (*fit.diff_tt) = 0;


  // we want to filter focus at first.
  // Note that focus is a different mode depending to modal basis.
  // we want to filter TT on altitude DM when using yao, and piston on all.
  if (use_mode=="yao") {
    for (nm=1;nm<=ndm;nm++){
      (*(*fit.coefs)(nm))(1) = 0; // filter Piston
      // filter focus:
      if (dm(nm).type="dh") (*(*fit.coefs)(nm))(5) = 0;
      else (*(*fit.coefs)(nm))(4) = 0;
      // filter astig on altitude DMs:
      if (nm>1) {
        if (dm(nm).type="dh") (*(*fit.coefs)(nm))([4,6]) = 0;
        else (*(*fit.coefs)(nm))([5,6]) = 0;
      }
      // filter TT on altitude DMs:
      if (nm>1) (*(*fit.coefs)(nm))(1:3) = 0;
      if (nmodesv(1)<opp.ncoefs) (*(*fit.coefs)(nm))(nmodesv(1)+1:) = 0;
    }
  } else {
    (*(*fit.coefs)(1))(1) = 0; // filter piston
    // filter focus:
    if (use_mode=="dh") (*(*fit.coefs)(1))(5) = 0; // for DH
    else (*(*fit.coefs)(1))(4) = 0; // for zernike & KL
    // Filter astigmatism if requested:
    if (first_nofit_astig) {
      if (use_mode=="dh") (*(*fit.coefs)(1))([4,6]) = 0; // for DH
      else (*(*fit.coefs)(1))([5,6]) = 0; // for zernike & KL
    }
    if (nmodesv(1)<opp.ncoefs) (*(*fit.coefs)(1))(nmodesv(1)+1:) = 0;
  }

  // get reference
  // extern all_otf_ref, a_ref;
  // all_otf_ref = opra_foo(x,a);
  // a_ref = a;
  // call lmfit
  res = lmfit_wrap(opra_foo,x,a,data,opp,fit=fit,\
                   tol=1e-6,eps=0.01,itmax=nitv(1),aregul=0.);
  if (stop_all) goto fin;

  // increase nmodes gently
  for (n=2;n<=numberof(nmodesv);n++) {
    aold = sdimold = [];
    if (use_mode!="yao") (*a.coefs)(1) = &( _(*(*a.coefs)(1),array(0,nmodesv(n)-nmodesv(n-1))) );
    opp.action = swrite(format="Pass %d: masks + aberrations up to %d",\
                        passn++,nmodesv(n));

    // filter some parameters
    fit = fresh_a(opp,val=1);
    fit.pupd     = 0;
    fit.stfmaskd = 0;
    if (fix_amp)   *fit.amps = 0;
    if (fix_pix)   fit.psize = 0;
    if (fix_kern)  fit.kernd = [0,0,0];
    if (fix_defoc) fit.defoc_scaling = 0;
    if (fix_cobs)  fit.cobs = 0;
    (*fit.diff_tt)(,1,1) = 0;
    if (fix_diff_tt) (*fit.diff_tt) = 0;

    // we want to filter focus at first.
    // Note that focus is a different mode depending to modal basis.
    // we want to filter TT on altitude DM when using yao, and piston on all.
    for (nm=1;nm<=opp.ndm;nm++){
      (*(*fit.coefs)(nm))(1) = 0; // filter Piston on all.
      // filter TT and quadratic on altitude DMs:
      if (nm>1) (*(*fit.coefs)(nm))(1:6) = 0;
      if (n<numberof(nmodesv)) {
        im = nmodesv(n)+1;
        if (im<opp.ncoef_per_dm(nm)) (*(*fit.coefs)(nm))(im:) = 0;
      }
    }

    // call lmfit
    res = lmfit_wrap(opra_foo,x,a,data,opp,fit=fit,\
                     tol=1e-14,eps=0.01,itmax=nitv(n),aregul=0.);
    if (stop_all) goto fin;
  }

  // FIXME FIXME, just for test
  write,"FIXME FIXME FIXME !!";
  fit.psize=1;
  res = lmfit_wrap(opra_foo,x,a,data,opp,fit=fit,\
                   tol=1e-14,eps=0.01,itmax=5,aregul=0.);

  fin:

  // fill back opp structure:
  opp.coefs=&(*a.coefs);
  opp.kernd = a.kernd; opp.kernd(1:2) *= 2.35; opp.kernd(3) *= 1e3;
  opp.cobs = opp.cobs+cobs_multiplier*abs(a.cobs);
  opp.psize *= 1.+atan(a.psize)*0.2;
  opp.pupd = opp.pupd+atan(a.pupd)*5;

  write,format="Elapsed time = %f sec\n",tac();
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

  write,format="%s",".";
  opra_iter++;
  if (aold==[]) {
    aold=opra_a_struct();
    aold.pupd = -9999;
    aold.cobs = -9999;
    aold.kernd = [-9999,-9999,1];
  }
  if (sdimold==[]) sdimold=0;

  // retrieve structures from x. Because it's the same address, modif
  // to opp and op will be passed to main opra() (magic).
  eq_nocopy,opp,*x(1);  // -> opp
  eq_nocopy,op,*x(2);   // -> op

  // base + TT of additional images + amplitudes of all images

  a = a_vec2struct(b,opp);

  //FIXME:
  // limits a's:
  // This is a bit adhoc, although the defaults look fine. FIXME?
  a.pupd          = opp.pupd+atan(a.pupd)*5;
  // a.kernd         = atan(a.kernd)*10;
  a.stfmaskd      = 1./(0.05+atan(a.stfmaskd)/pi*0.0999);
  a.defoc_scaling = 1.+atan(a.defoc_scaling)*1.;//was *3 orginally
  a.psize         = 1.+atan(a.psize)*0.2;

  a.amps = &(1.+atan(*a.amps)*0.6);
  //if (*a.amps == 0)  *a.amps = 1.0;
  //*a.amps = *a.amps;
  //It could be good to bound properly the parameters.
  //But is that possible with lmfit ?

  if ((aold.pupd!=a.pupd)||(aold.cobs!=a.cobs)) {
    // pupd has changed, need to recompute pupils and possibly modes.
    // Why possibly? kl are computed on even int dim, so if pupd has not
    // gone above the next even int value, no need to recompute (this is
    // done in opra_gen_modes() ).
    offset   = ((opp.modes_type=="kl")? 0.5: 1); // FIXME: not sure it's good.

    //BN 4 tests:
    offset = 0.5; // FR?

    if (opp.modes_type=="yao") {
      write,"Generating yao modes";
      sim.pupildiam = opp.pupd;
      aoinit,disp=0,clean=1;
      opp.pupi = ipupil;
      opp.pupr = pupil;
      prepzernike,sim._size,sim.pupildiam,sim._cent,sim._cent;
    } else {
      if (aold.pupd!=a.pupd) {
        write,"Generating modes";
        opp.modes = &( opra_gen_modes(opp.otf_dim,a.pupd,nmodes,       \
                                      mode_type=opp.modes_type));
      }
      opp.pupi = make_pupil(opp.otf_dim,a.pupd,                        \
                            xc=opp.otf_dim/2+offset,                   \
                            yc=opp.otf_dim/2+offset,                   \
                            cobs=opp.cobs+cobs_multiplier*abs(a.cobs));

      opp.pupr = make_pupil(opp.otf_dim,a.pupd,real=1,                 \
                            xc=opp.otf_dim/2+offset,                   \
                            yc=opp.otf_dim/2+offset,                   \
                            cobs=opp.cobs+cobs_multiplier*abs(a.cobs));
      nn = 4;
      tmp = make_pupil(opp.otf_dim*nn,a.pupd*nn,                  \
                       xc=opp.otf_dim/2*nn+offset,                \
                       yc=opp.otf_dim/2*nn+offset,                \
                       cobs=opp.cobs+cobs_multiplier*a.cobs);
      opp.pupr = bin2d(tmp,nn);
      // opp.pupr = opp.pupi;
      // opp.pupi = opp.pupr;
      prepzernike,opp.otf_dim,a.pupd,opp.otf_dim/2+offset,opp.otf_dim/2+offset;
      // above: needed anyway for added defocus
      // FIXME: need to modify xc/yc as using kl???
    }
    // cache:
    extern z2ext,z3ext,z4ext;
    z2ext = zernike_ext(2);
    z3ext = zernike_ext(3);
    z4ext = zernike_ext(4);

  }

  // if source TF mask diameter has changed. Recompute:
  if (aold.stfmaskd!=a.stfmaskd) {
    // upgrade with Bessel function (FIXME)
    xy = dist(opp.otf_sdim)/a.stfmaskd;
    opp.stfmask = sinc(xy);
  }

  // blur kernel:
  if (anyof(aold.kernd!=a.kernd)) {
    // opp.kernel = makegaussian(opp.otf_sdim,100./(1e-6+a.kernd^2.));
    xy = indices(opp.otf_sdim)-opp.otf_sdim/2.-0.5;
    // xy(,,1) *= 2e-5/(1e-3+a.kernd(1)^2.);
    // xy(,,2) *= 2e-5/(1e-3+a.kernd(2)^2.);
    // xy(,,1) *= 0.3*abs(atan(a.kernd(1)));
    // xy(,,2) *= 0.3*abs(atan(a.kernd(2)));
    xy =  xy(,,+) * mrot(a.kernd(3)*1000.)(+,);
    opp.kernel = exp( -0.5*((xy(,,1)*a.kernd(1)/sqrt(2.))^2 \
                           +(xy(,,2)*a.kernd(2)/sqrt(2.))^2.));
    // kernd 1 and 2 are thus sigmas
  }

  // Iteration accounting
  if (lmfititer==[]) lmfititer=0;
  if (newiter) lmfititer++;

  // build PSF, OTF
  all_otf = [];

  if (opp.modes_type=="yao") {
    // build mircube. we need that only once, irrespective of positions
    extern mircube;
    mircube = array(0.0f,[3,sim._size,sim._size,ndm]);
    for (nm=1;nm<=ndm;nm++) {
      dm(nm)._command = &(float(*(*a.coefs)(nm)));
      mircube(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2,nm) = compDmShape(nm,dm(nm)._command);
    }
    multWfs,1,disp=0;
  }

  // precompute:
  kern_stfmask = opp.kernel*opp.stfmask;

  // loop on position
  for (n=1;n<=opp.npos;n++) {

    opp.phase(,,n) *= 0;

    if (opp.modes_type=="yao") {
      wdim = dimsof(*wfs(n)._fimage)(2);
      n12 = _(opp.otf_dim/2-wdim/2+1,opp.otf_dim/2+wdim/2);
      opp.phase(n12(1):n12(2),n12(1):n12(2),n) = *wfs(n)._fimage;
    } else {
      // FIXME multiple positions? nope.
      az = *(*a.coefs)(1);
      for (i=2;i<=nmodes;i++) opp.phase(,,n) += az(i)*(*opp.modes)(,,i);
    }

    // loop on images (defocs)
    for (i=1;i<=opp.nim;i++) {
      // same = allof(*(*a.coefs)(1)==*(*aold.coefs)(1));
      // for (nm=2;nm<=opp.ndm;nm++) \
        // same = (same && allof(*(*a.coefs)(1)==*(*aold.coefs)(1)));

      // compute defoc to add to phase
      defoc = (op(i,n).delta_foc *  a.defoc_scaling) * z4ext;
      // compute tiptilt to add to phase
      // if (i>=2) {
      if ((n==1)&&(i==1)) tt = 0; // no special TT passed for image#1 (TT in "phase")
      else {
        tt = (*a.diff_tt)(1,i,n) * z2ext +   \
             (*a.diff_tt)(2,i,n) * z3ext;
      }
      // compute phase
      pha     = opp.phase(,,n) + tt + defoc;
      // compute PSF
      psf = calc_psf_fast(opp.pupr,pha);
      // normalize
      psf     = psf/sum(psf) * abs((*a.amps)(i,n));
      op(i,n).psf  = psf;
      // Note that above, PSF does not reflect the true psf, as the source
      // mask effect is not included. But we don't need it at each opra_foo()
      // call, just for display, so I did it in opra_info_and_plots()
      // compute OTF, include attenuation by blur (kernel) and mask (source size)
      op(i,n).otf  = get_mtf(psf,opp.otf_sdim) * kern_stfmask(,,-);
      // Pixel size adjustement:
      center = opp.otf_sdim/2+1;
      xy = (indgen(opp.otf_sdim)-center)*a.psize+center;
      op(i,n).otf(,,1) = bilinear(op(i,n).otf(,,1),xy,xy,grid=1);
      op(i,n).otf(,,2) = bilinear(op(i,n).otf(,,2),xy,xy,grid=1);
      op(i,n).otf(center,center,) = 0.; // done below.
      // append to previous OTF all object:
      grow,all_otf,op(i,n).otf(,,,-);
    }
  }

  if (newiter) {
    if (has_svipc) {
      write_opra_struct_to_shm,a,op,opp,mircube;
      _yorick_gui_proc,"opra_gui_plots\n";
    } else {
      if (opp.modes_type=="yao") {
        window,opp.winnum+opp.npos;
        tv,mircube(,*),square=1;
      }
    }
    opra_info_and_plots,a,opp,op,noplots=has_svipc;
    if (!has_svipc) pause,100; // give time for graphic update
  }


  aold = a;
  newiter=0;

  // get rid of (0,0) frequency information
  // all_otf(center,center,) = 0;
  all_otf = all_otf(*,,)(wotf,,);
  all_otf = float(all_otf(*)); // to save RAM
  rms1=(all_otf)(rms);
  app = [];
  for (nm=1;nm<=opp.ndm;nm++) grow,app,(*(*a.coefs)(nm));
  if (opra_coef_regul==[]) opra_coef_regul=1e-3;
  app = app*rms1*numberof(all_otf)/numberof(app)*opra_coef_regul;
  if (opp.npos>1) grow,all_otf,app;

  return all_otf;
}

if (original_quit==[]) original_quit = quit;

func opra_quit(void)
{
  write,format="%s\n","Asking fork to quit";
  // notify fork to quit:
  _yorick_gui_proc,"fork_quit\n";
  // waiting for acknowledgment:
  tic,7;
  while ((!shm_read(shmkey,"quit!")(1))&&(tac(7)<1.0)) usleep,20;
  write,format="%s\n","Cleaning up svipc";
  shm_cleanup,shmkey;
  sem_cleanup,semkey;
  original_quit;
}

func poll_quit(void)
{
  if (shm_read(shmkey,"quit?")(1)) opra_quit;
  after,0.1,poll_quit;
}
