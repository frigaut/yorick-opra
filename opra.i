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

nmodes_max4printout = 30;

if (find_in_path("svipc.i")!=[]) has_svipc=1;


OPRA_VERSION = "1.8";

func opra(images, defocs, lambda, pixsize, teldiam, nmodes=, use_mode=, cobs=,
          noise=, pupd=, otf_dim=, progressive=, niter=, fix_amp=, fix_pix=,
          fix_kern=, fix_defoc=, fix_diff_tt=, first_nofit_astig=, winnum=,
          dpi=, pal=, yao_parfile=, gui=)
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
   fix_diff_tt=    1 if differential TT between images shoud not be a free parameter
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
  if (!nmodes)    nmodes = 120;
  if (winnum==[]) winnum = 1;
  if (!noise)     noise = 0.;
  if (!dpi)       dpi = 130;
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

  opp = oprapar_struct();

  opp.nim      = dimsof(images)(4);  // number of input images
  opp.npos     = dimsof(images)(5);  // number of input positions
  opp.im_dim   = dimsof(images)(2);  // side dimension of input images
  opp.winnum   = winnum;

  // svipc and GUI:
  if ((gui)&&(has_svipc)) {
    require,"svipc.i";
    has_svipc=1;
    extern shmkey,semkey;
    shmkey = 0x00abacab;
    semkey = 0x0aabacab;
    write,format="shmkey = %#x\n",shmkey;
    shm_init,shmkey,slots=20;
    sem_init,semkey,nums=20;
    shm_write,shmkey,"stop",&([0]);
    shm_write,shmkey,"next_stage",&([0]);
    shm_write,shmkey,"xids",&([0]);
    shm_write,shmkey,"quit?",&([0]);
    shm_write,shmkey,"quit!",&([0]);
    shm_write,shmkey,"redraw",&([0]);
    if (fork()==0) { // I'm the child
      // Start GUI
      extern fdpi; fdpi=dpi;
      images = [];
      require,"opra_gui.i";
      exit;
    } else {
      set_idler,poll_quit;
    }
  }

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
  if (pupd) opp.pupd = pupd;

  opp.otf_sdim = long(2*opp.pupd*1.1)/2*2;
  opp.otf_sdim = long(2*opp.pupd*1.02)/2*2;
  opp.otf_dim  = long(2^ceil(log(opp.otf_sdim)/log(2)));
  if(otf_dim) opp.otf_dim = opp.otf_sdim = otf_dim;

  if (use_mode=="yao") {
    if (!yao_parfile) error,"yao parfile not defined";
    if (findfiles(yao_parfile)==[]) \
      error,swrite(format="Can't find yao parfile %s\n",yao_parfile);
    // else we are good to go.
    // FIXME : use regular yao.
    cwd = pwd();
    cd,"/home/frigaut/mcao/myst/yorick";
    require,"yao_mcao.i";
    cd,cwd;
    write,"yao aoread() and aoinit()";
    aoread,yao_parfile;
    sim.pupildiam = opp.pupd;
    aoinit,disp=0,clean=1;
    // nmodes = sum(dm._nact-1)+1-2*(ndm-1); // also subtract TT in altitude DMs
    opp.ncoefs = sum(dm._nact);
    opp.ncoef_per_dm = &dm._nact;
    opp.ndm = ndm;
    dpi2 = long(dpi*ndm/3.);

  } else {
    require,"yao.i";
    opp.ncoefs = nmodes;
    opp.ncoef_per_dm = &([opp.ncoefs]);
    opp.ndm = 1;
  }

  //#################################
  // create graphic windows
  //#################################
  if ((gui)&&(has_svipc)) {
    sem_take,semkey,0;
    xids = shm_read(shmkey,"xids");
  }
  for (n=1;n<=opp.npos;n++) {
    winn = winnum+n-1;
    if ((gui)&&(has_svipc)) {
      window,winn,wait=0,style="opra.gs",dpi=dpi,width=0,height=0,parent=xids(n);
      pause,50;
    } else {
      if (window_exists(winn)) winkill,winn;
      window,winn,wait=1,style="opra.gs",dpi=dpi,width=0,height=0;
    }
    if (pal) palette,pal;
    for (i=1;i<=3;i++) { plsys,i; limits,square=1; }
  }
  if ((use_mode=="yao")&&(gui)&&(has_svipc)) {
    window,winnum+opp.npos,style="nobox.gs",width=0,height=0,\
      dpi=dpi,wait=0,parent=xids(opp.npos+1);
  } else {
    if (window_exists(winnum+opp.npos)) winkill,winnum+opp.npos;
    window,winnum+opp.npos,style="nobox.gs",width=long(6.*dpi2/ndm),height=6*dpi2,dpi=dpi2,wait=1;
  }

  opp.cobs = (cobs?cobs:0.);
  op  = array(opra_struct,[2,opp.nim,opp.npos]);

  // normalize
  for (i=1;i<=opp.npos;i++) images(,,,i) = images(,,,i)/sum(images(,,1,i));

  center = opp.otf_sdim/2+1;
  data = [];

  for (n=1;n<=opp.npos;n++) {
    for (i=1;i<=opp.nim;i++) {
      op(i,n).psf_data  = &(images(,,i,n));
      op(i,n).otf_data  = &(get_mtf(images(,,i,n),opp.otf_sdim));
      // get rid of (0,0) frequency information
      grow,data,(*op(i,n).otf_data)(,,,-);
      op(i,n).delta_foc = defocs(i);
      op(i,n).noise     = noise;
    }
  }
  data(center,center,,) = 0;

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
    nmodesv = _(6,opp.ncoefs);
    nitv    = array(10,numberof(nmodesv));// 10iteration max (?)
    nitv(2) = 30;
  }

  if (use_mode=="yao") { // Can't do progressive with yao (see comment below)
    nmodesv = [opp.ncoefs,opp.ncoefs];
    // above: has to be, we don't know "modes" are modal (could be zonal)
    nitv = [6,niter];
  }

  //#################################
  // Let's start !
  //#################################
  // See the parameters descriptions in DOCUMENT section of opra_foo.
  // initialize "a" (the coefficients to find)
  a = fresh_a(opp);
  a.kernd    = 0.;
  a.stfmaskd = -100.;
  opp.action = swrite(format="Pass %d: masks + aberrations up to %d",\
                      passn++,nmodesv(1));

  // filter some parameters. Use fit keyword of lmfit.
  fit = fresh_a(opp,val=1);
  fit.pupd     = 0;
  fit.kernd    = 0;
  fit.stfmaskd = 0;
  if (fix_amp)   *fit.amps = 0;
  if (fix_pix)   fit.psize = 0;
  if (fix_kern)  fit.kernd = 0;
  if (fix_defoc) fit.defoc_scaling = 0;
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
  }

  // call lmfit
  res = lmfit_wrap(opra_foo,x,a,data,opp,fit=fit,\
                   tol=1e-6,eps=(use_mode=="yao"?0.001:0.01),itmax=nitv(1),aregul=0.);
  if (stop_all) goto fin;

  // increase nmodes gently
  for (n=2;n<=numberof(nmodesv);n++) {
    aold = sdimold = [];
    if (use_mode!="yao") (*a.coefs)(1) = &( _(*(*a.coefs)(1),array(0,nmodesv(n)-nmodesv(n-1))) );
    // a.psize = clip(a.psize,-10,10);
    opp.action = swrite(format="Pass %d: masks + aberrations up to %d",\
                        passn++,nmodesv(n));

    // filter some parameters
    fit = fresh_a(opp,val=1);
    fit.pupd     = 0;
    fit.stfmaskd = 0;
    if (fix_amp)   *fit.amps = 0;
    if (fix_pix)   fit.psize = 0;
    if (fix_kern)  fit.kernd = 0;
    if (fix_defoc) fit.defoc_scaling = 0;
    (*fit.diff_tt)(,1,1) = 0;
    if (fix_diff_tt) (*fit.diff_tt) = 0;


    // we want to filter focus at first.
    // Note that focus is a different mode depending to modal basis.
    // we want to filter TT on altitude DM when using yao, and piston on all.
    for (nm=1;nm<=opp.ndm;nm++){
      (*(*fit.coefs)(nm))(1) = 0; // filter Piston on all.
      // filter TT and quadratic on altitude DMs:
      if (nm>1) (*(*fit.coefs)(nm))(1:6) = 0;
    }

    // call lmfit
    res = lmfit_wrap(opra_foo,x,a,data,opp,fit=fit,\
                     tol=1e-12,eps=(use_mode=="yao"?0.001:0.01),itmax=nitv(n),aregul=0.);
    if (stop_all) goto fin;
    if (opra_debug)                                                     \
      fits_write,"phase.fits",*opp.phase * (*opp.pupi),overwrite=1;
  }

  if (opra_debug)                                                       \
    fits_write,"phase.fits",*opp.phase * (*opp.pupi),overwrite=1;

  fin:
  write,format="Elapsed time = %f sec\n",tac();

  opp.coefs=&(*a.coefs);

  // if (has_svipc) shm_cleanup,shmkey;

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

  if (shm_read(shmkey,"redraw")(1)) {
	  pause,50;
	  shm_write,shmkey,"redraw",&([0]);
	}

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

  a = a_vec2struct(b,opp);

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

  a.amps = &(1.+atan(*a.amps)*0.6);
  //if (*a.amps == 0)  *a.amps = 1.0;
  //*a.amps = *a.amps;
  //It could be good to bound properly the parameters.
  //But is that possible with lmfit ?

  if (aold.pupd!=a.pupd) {
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
      opp.pupi = &ipupil;
      opp.pupr = &pupil;
      prepzernike,sim._size,sim.pupildiam,sim._cent,sim._cent;
      // prepzernike,opp.otf_dim,a.pupd,opp.otf_dim/2+offset,opp.otf_dim/2+offset;
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
      prepzernike,opp.otf_dim,a.pupd,opp.otf_dim/2+offset,opp.otf_dim/2+offset;
    }
    // above: needed anyway for added defocus
    // FIXME: need to modify xc/yc as using kl???
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

  // Build phase:
  if (*opp.phase==[]) opp.phase = &(array(float,[2,opp.otf_dim,opp.otf_dim]));

  // loop on position
  for (n=1;n<=opp.npos;n++) {

    *opp.phase *= 0.0f;

    if (opp.modes_type=="yao") {
      n12 = wfs(n).n12;
      (*opp.phase)(n12(1):n12(2),n12(1):n12(2)) = *wfs(n)._fimage;
      // tv,*opp.phase; hitReturn;
    } else {
      // FIXME multiple positions? nope.
      az = *(*a.coefs)(1);
      for (i=2;i<=nmodes;i++) (*opp.phase)(,,1) += az(i)*(*opp.modes)(,,i);
    }

    // cache:
    z2ext = zernike_ext(2);
    z3ext = zernike_ext(3);
    z4ext = zernike_ext(4);

    // loop on images (defocs)
    for (i=1;i<=opp.nim;i++) {
      // compute defoc to add to phase
      defoc = op(i,n).delta_foc *  a.defoc_scaling * z4ext;
      // compute tiptilt to add to phase
      // if (i>=2) {
      if ((n==1)&&(i==1)) tt = 0; // no special TT passed for image#1 (TT in "phase")
      else {
        tt = (*a.diff_tt)(1,i,n) * z2ext +       \
             (*a.diff_tt)(2,i,n) * z3ext;
      }
      // compute complex wavefront
      phi     = array(complex,[2,opp.otf_dim,opp.otf_dim]);
      phi.re  = (*opp.pupr) * cos((*opp.phase) + tt + defoc);
      phi.im  = (*opp.pupr) * sin((*opp.phase) + tt + defoc);
      // compute PSF
      psf     = roll(abs(fft(phi,1))^2.);
      // normalize
      psf     = psf/sum(psf) * (*a.amps)(i,n);
      op(i,n).psf  = &psf;
      // Note that above, PSF does not reflect the true psf, as the source
      // mask effect is not included. But we don't need it at each opra_foo()
      // call, just for display, so I did it in opra_info_and_plots()
      // compute OTF, include attenuation by blur (kernel) and mask (source size)
      op(i,n).otf  = &( get_mtf(psf,opp.otf_sdim ) * \
                        ((*opp.kernel) * (*opp.stfmask))(,,-) );
      // Pixel size adjustement:
      center = opp.otf_sdim/2+1;
      xy = (indgen(opp.otf_sdim)-center)*a.psize+center;
      (*op(i,n).otf)(,,1) = bilinear((*op(i,n).otf)(,,1),xy,xy,grid=1);
      (*op(i,n).otf)(,,2) = bilinear((*op(i,n).otf)(,,2),xy,xy,grid=1);
      //    (*op(i).otf)(center,center,) = 0.;
      // append to previous OTF all object:
      grow,all_otf,(*op(i,n).otf)(,,,-);
    }

    if (newiter) {
      if (opp.modes_type!="yao") opra_info_and_plots,a,*a.amps,az,opp,op;
      else {
        if (n==1) {
          window,opp.winnum+opp.npos;
          mircube(,,1) *= ipupil;
          tv,mircube(,*),square=1;
        }
        window,opp.winnum+n-1;
        opra_info_and_plots,a,*a.amps,*a.coefs,opp,op(,n),yao=1,star=n;
        if (n==opp.npos) pause,100; // give time for graphci update
      }
    }
  }

  aold = a;
  newiter=0;

  // get rid of (0,0) frequency information
  all_otf(center,center,) = 0;

  all_otf = float(all_otf(*)); // to save RAM
  rms1=(all_otf)(rms);
  app = [];
  for (nm=1;nm<=opp.ndm;nm++) grow,app,(*(*a.coefs)(nm));
  if (opra_coef_regul==[]) opra_coef_regul=1e-3;
  app = app*rms1*numberof(all_otf)/numberof(app)*opra_coef_regul;
  // app = app*3e-3*rms1;
  // sum(all_otf); app(rms);
  if (opp.npos>1) grow,all_otf,app;

  return all_otf;
}

if (original_quit==[]) original_quit = quit;

func opra_quit(void)
{
  write,format="%s\n","Asking fork to quit";
  // notify fork to quit:
  shm_write,shmkey,"quit?",&([1]);
  // waiting for acknowledgment:
  while (!shm_read(shmkey,"quit!")(1)) usleep,20;
  write,format="%s\n","Cleaning up svipc";
  shm_cleanup,shmkey;
  sem_cleanup,semkey;
  original_quit;
}
if (has_svipc) quit=opra_quit;

func poll_quit(void)
{
  if (shm_read(shmkey,"quit?")(1)) opra_quit;
  after,0.1,poll_quit;
}
