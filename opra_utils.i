local opra_utils;
/* DOCUMENT opra_utils
 * OPRA_UTILS (OTF-based Phase Retrieval Analysis)
 * Copyright (c) Francois Rigaut 2010.
 * $Id: opra_utils.i 9 2010-01-16 19:15:52Z frigaut $
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
 *
 * Utilitary functions for OPRA
 * - plots and printouts
 * - modes generation
 * - test image generation
 * SEE ALSO:
 */



func get_mtf(psf,dim,zero=)
/* DOCUMENT get_mtf(psf,dim)
   return the mtf (more precisely, the TF of the PSF),
   as a 2 plan data cube:
   plan 1 = real part
   plan 2 = imaginary part
   SEE ALSO:
 */
{
  local mtf;
  mtf = fft(roll(psf),-1);
  if (zero) {
    mtf.re(1,1) = 0;
    mtf.im(1,1) = 0;
  }
  mtf = roll(mtf,[dim/2,dim/2])(1:dim,1:dim);
  return [mtf.re,mtf.im];
}


func lmfit_wrap(f,x,&a,y,nim,fit=,tol=,eps=,itmax=,aregul=)
{
  a = a_struct2vec(a);
  if (fit) {
    fit = a_struct2vec(fit);
    w = where(fit);
    if (numberof(w)) fit = w; else fit=[];
  }
  // indices on which to do the additional regularization (a)
  //  aregul_i = indgen(numberof(a)-(npbc+3))+(npbc+3);
  aregul_i = indgen(numberof(a));
  aregul_i(npbc+1:npbc+2) = 0;
  aregul_i(1:npbc+2) = 0;
  aregul_i = aregul_i(where(aregul_i));

  res = opra_lmfit(f,x,a,y,fit=fit,tol=tol,eps=eps,itmax=itmax,\
                   aregul=aregul,aregul_i=aregul_i,neval_max=300);
  a = a_vec2struct(a,nim);
  return res;
}

func a_struct2vec(a)
{
  return _(a.pupd,a.kernd,a.stfmaskd,a.defoc_scaling,a.psize,   \
           (*a.diff_tt)(*),*a.amps,*a.coefs);
}


func a_vec2struct(b,nim)
{
  extern npbc;
  // Number of Parameter Before (mode) Coefs. (hence the name)
  npbc = 5  + (opp.nim-1)*2 + opp.nim;

  a = opra_a_struct();
  a.pupd          = b(1);
  a.kernd         = b(2);
  a.stfmaskd      = b(3);
  a.defoc_scaling = b(4);
  a.psize         = b(5);
  diff_tt         = b(6:6+2*(nim-1)-1);
  a.diff_tt       = &(reform(diff_tt,[2,2,nim-1]));
  a.amps          = &(b(6+2*(nim-1):6+2*(nim-1)+nim-1));
  a.coefs         = &(b(npbc+1:));
  return a;
}

func fresh_a(nim,ncoefs,val=)
{
  if (val==[]) val=0.;
  np = 5+2*(nim-1)+nim+(ncoefs-1);
  a = array(val,np);
  return a_vec2struct(a,nim);
}


func opra_info_and_plots(a,amps,az,nmodes,&opp,&op)
/* DOCUMENT opra_info_and_plots()
   This routine is called at each lmfit() new iteration to
   print out the current state of the optimizer and plot
   OTFs, PSFs, phase, convergence criteria.
   SEE ALSO:
 */
{
  extern itvec,distvec;
  extern lmfititer,newiter,opraiter;

  mdim = opp.otf_sdim;
  all = array(0.,[2,(opp.nim*2)*mdim,3*mdim]);
  k = 0;
  // loop on all images:
  for (i=1;i<=opp.nim;i++) {
    // stuff in real part of OTF (data and model), and difference
    all(k*mdim+1:(k+1)*mdim,) =                                        \
      _( (*op(i).otf_data)(,,1),                                       \
         (*op(i).otf)(,,1),                                            \
         abs((*op(i).otf_data)(,,1)-(*op(i).otf)(,,1)));
    k++;
    // stuff in imag part of OTF (data and model), and difference
    all(k*mdim+1:(k+1)*mdim,) =                                        \
      _( (*op(i).otf_data)(,,2),                                       \
         (*op(i).otf)(,,2),                                            \
         abs((*op(i).otf_data)(,,2)-(*op(i).otf)(,,2)));
    k++;
  }

  // OTF (data, model and difference)
  plsys,0;
  fma;

  pli,(abs(all));
  plt,"OTFs (re,im,re,im...)",0.145,0.635,tosys=0,height=10;
  // plt,"Data",0.130,0.502,orient=1,tosys=0,height=8,justify="CA";
  plt,"Data<- Model ->Diff",0.130,0.535,orient=1,tosys=0,height=8,justify="CA";
  // plt,"Diff",0.130,0.659,orient=1,tosys=0,height=8,justify="CA";

  plt,opp.action,0.145,0.90,tosys=0,height=12;
  txt = swrite(format="Iterations: %d/%d, total %d",lmfititer_pass,     \
               lmfit_itmax,lmfititer);
  plt,txt,0.655,0.90,tosys=0,height=10,justify="RA";

  // Images: data and model
  plsys,2;
  for (i=1;i<=opp.nim;i++) {
    if (opp.nim<=2) {
      pli,*op(i).psf_data,2*(i-1),0,2*(i-1)+1,1;
      plt,"Data",2*(i-1)+0.5,1,tosys=2,height=8,justify="CT",color="white";
    } else {
      pli,*op(i).psf_data,(i-1),0,(i-1)+1,1;
      plt,"Data",(i-1)+0.5,1,tosys=2,height=6,justify="CT",color="white";
    }
    // Now model image. Remember? We have to apply the object
    // mask (see note in opra.i). For better eye comparison with
    // data, I chose to also rebin and put some noise:
    otf2 = array(complex,dimsof((*op(1).otf)(,,1)));
    otf2.re = roll((*op(i).otf)(,,1));
    //    center = opp.otf_sdim/2+1;
    //    otf2.re(center,center) = otf2.re(center+1,center);
    otf2.im = -roll((*op(i).otf)(,,2));
    psf = roll((fft(otf2,-1)).re)/opp.otf_sdim^2.;
    //    psf -= median(psf(1:20,1:20)(*)); // sky subtraction

    /*
      psf = *op(i).psf;
    psf = psf/maxpsf1/amps(i);
    psf = abs(fft(fft(psf,1)*roll(mask),-1))/opp.otf_dim^2.;
    psf = spline2(psf,opp.im_dim,opp.im_dim);
    */
    psf = spline2(psf,opp.im_dim,opp.im_dim);
    psf = psf/sum(psf);
    psf = psf+gaussdev(dimsof(psf))*op(i).noise;
    psf = clip(psf,0.,);
    if (opp.nim<=2) {
      pli,psf,2*(i-1)+1,0,2*(i-1)+2,1;
      plt,"Model",2*(i-1)+1.5,1,tosys=2,height=8,justify="CT",color="white";
    } else {
      pli,psf,(i-1),1,(i-1)+1,2;
      plt,"Model",(i-1)+0.5,2,tosys=2,height=6,justify="CT",color="white";
    }
    //    if (hitReturn()=="s") error;
  }
  plt,"Images",0.145,0.882,tosys=0,height=10;

  // Model phase
  plsys,3;
  ipupr = long(ceil(a.pupd)/2.)+1;
  // compute phase w/o TT:
  pha = *opp.phase * 0.;
  for (i=4;i<=nmodes;i++) pha += az(i)*(*opp.modes)(,,i);
  iminmax = minmax(pha(where(*opp.pupi)));
  phase_rms = pha(where(*opp.pupi))(rms);
  strehl = exp(-phase_rms^2.);
  dim = opp.otf_dim;
  //  tmp = ((*opp.phase-iminmax(1))*(*opp.pupi))               \
  //    (dim/2-ipupr+1:dim/2+ipupr,dim/2-ipupr+1:dim/2+ipupr);
  tmp = ((pha-iminmax(1))*(*opp.pupi))                   \
    (dim/2-ipupr+1:dim/2+ipupr,dim/2-ipupr+1:dim/2+ipupr);
  pli,tmp;
  txt = swrite(format="Phase (%.2f->%.2f rd). Strehl=%.1f%%",iminmax(1),\
               iminmax(2),strehl*100.);
  plt,txt,0.576,0.635,tosys=0,height=8,justify="CA";

  // convergence plot
  grow,itvec,lmfititer;
  otf_diff = [];
  for (i=1;i<=opp.nim;i++) grow,otf_diff,(*op(i).otf_data-*op(i).otf)(*);
  grow,distvec,(otf_diff)(rms);
  if (numberof(itvec)>=2) {
    plsys,4;
    plh,distvec,itvec;
    range,min(distvec)*0.98;
    logxy,0,1;
    txt = swrite(format="Convergence crit. (current=%.2e)",distvec(0));
    plt,txt,0.50,0.474,tosys=0,height=8;
    plt,"Distance to data",0.471,0.425,tosys=0,orient=1,height=8,justify="CA";
    plt,"Iteration",0.58,0.357,tosys=0,height=8,justify="CA";
  }

  system,"clear";
  write,format="%s\n",pass_action;
  write,format="lmfit Iteration#   : %d\n",lmfititer;
  write,format="Pass through foo2  : %d\n",opraiter;
  write,format="Support diameter   : %.2f\n",a.pupd;
  write,format="Kernel [pixels]    : %.2f\n",a.kernd;
  write,format="Mask diameter      : %.2f\n",a.stfmaskd;
  write,format="foc-defoc defocus  : %.2f\n",a.defoc_scaling;
  write,format="Pixel size scaling : %.2f\n",a.psize;
  //  write,format="foc-defoc tiptilt  : %.2f,%.2f\n",a(5),a(6);
  write,format="Intensity ratios   : %.2f",(*a.amps)(1);
  for (i=2;i<=opp.nim;i++) write,format=",%.2f",(*a.amps)(i);
  write,"";
  for (i=2;i<=min(nmodes,nmodes_max4printout);i++) {
    write,format="a(%2d)              : %+7.0f mrd\n",i,1000*az(i);
  }

  ytop = 0.43; ydelta = 0.010; k = 0; xleft=0.125; hf = 7;
  txt = swrite(format="Support diameter [pixels]          : %.2f",a.pupd);
  plt,txt,xleft,ytop-(k++)*ydelta,tosys=0,height=hf,font="courier";
  txt = swrite(format="Gaussian Kernel FWHM [pixels]      : %.2f",abs(a.kernd));
  plt,txt,xleft,ytop-(k++)*ydelta,tosys=0,height=hf,font="courier";
  //  txt = swrite(format="Mask diameter [fixme]              : %.2f",abs(a.stfmaskd));
  //  plt,txt,0.15,ytop-(k++)*ydelta,tosys=0,height=8,font="courier";
  txt = swrite(format="Defocus scaling fact. [no units]   : %.2f",a.defoc_scaling);
  plt,txt,xleft,ytop-(k++)*ydelta,tosys=0,height=hf,font="courier";
  txt = swrite(format="Pixel size scal. fact. [no units]  : %.2f",a.psize);
  plt,txt,xleft,ytop-(k++)*ydelta,tosys=0,height=hf,font="courier";
  txt = swrite(format="Image ampl. scal. fact. [no units] : %.2f",(*a.amps)(1));
  for (i=2;i<=opp.nim;i++) txt+=swrite(format=",%.2f",(*a.amps)(i));
  plt,txt,xleft,ytop-(k++)*ydelta,tosys=0,height=hf,font="courier";
}


func make_object_otf(dim,pupd,exponent,xc=,yc=,cobs=)
{
  return exp(-(dist(dim,xc=xc,yc=yc)/(pupd/2.))^exponent)^0.69314;
}


func opra_gen_test(dim,defocus,amp=,prefix=)
{
  extern opp;
  if (amp==[]) amp=1.;
  if (prefix==[]) prefix="test";

  opp = oprapar_struct();
  opp.nim = 2;
  opp.im_dim = dim;
  opp.modes_type = "kl";
  lambda = 810e-9;
  teldiam = 7.9;
  pixsize = 1.99e-3;
  fwhm_estimate = lambda/teldiam/4.848e-6/pixsize;
  opp.pupd = (opp.im_dim/2.)*(2./fwhm_estimate);
  opp.otf_sdim = long(2*opp.pupd*1.1)/2*2;
  opp.otf_dim  = long(2^ceil(log(opp.otf_sdim)/log(2)));
  opp.cobs = 0.;
  op  = array(opra_struct,opp.nim);
  nmodes = 600;
  sdimold = [];
  opp.modes = &( opra_gen_modes(opp.otf_dim,opp.pupd,nmodes,\
                                  mode_type=opp.modes_type));

  prepzernike,opp.otf_dim,opp.pupd;

  offset   = ((opp.modes_type=="kl")? 0.5: 1);

  opp.pupi = &( make_pupil(opp.otf_dim,opp.pupd,                        \
                           xc=opp.otf_dim/2+offset,                     \
                           yc=opp.otf_dim/2+offset,                     \
                           cobs=opp.cobs) );

  phase = (*opp.pupi)*0.;
  sdim = opp.otf_dim;
  phi = array(complex,[2,dim,dim]);
  phi.re(1:sdim,1:sdim) = (*opp.pupi)*cos(phase);
  phi.im(1:sdim,1:sdim) = (*opp.pupi)*sin(phase);
  airy = roll(abs(fft(phi,1))^2.);

  coefs = amp*(random(nmodes)-0.5);
  coefs(2:4) *= 0.01;
  for (i=4;i<=nmodes;i++) coefs(i) /= zernumero(i)(1)^2.;
  phase = (*opp.modes)(,,+) * coefs(+);
  fits_write,swrite(format="images/phase_%s.fits",prefix),      \
    phase*(*opp.pupi),overwrite=1;

  // generate images:
  for (n=1;n<=numberof(defocus);n++) {
    pha = phase + defocus(n) * zernike_ext(4);
    sdim = opp.otf_dim;
    phi = array(complex,[2,dim,dim]);
    phi.re(1:sdim,1:sdim) = (*opp.pupi)*cos(pha);
    phi.im(1:sdim,1:sdim) = (*opp.pupi)*sin(pha);
    im = roll(abs(fft(phi,1))^2.);
    tv,im; hitReturn;
    if (n==1) write,format="Strehl=%.1f%%\n",max(im)/max(airy)*100.;
    fits_write,swrite(format="images/%s_%d.fits",prefix,n),im,overwrite=1;
  }
  tv,phase*(*opp.pupi); pause,1000;
}

func plot_modes(mode_type)
{
  nm = 100;
  mc=opra_gen_modes(100,90,nm,mode_type=mode_type);
  mc = mc-mc(1,1,)(-,-,);
  pup = (mc(,,rms)!=0);
  mc = mc-mc(*,)(max,)(-,-,);
  mc = mc*pup(,,-);
  mc(,,1) = -pup;
  winkill,3;
  window,3,dpi=170,style="nobox.gs";
  pli,array(1.,[2,100,100]),-0.5,-0.5,10.5,10.5;
  for (i=1;i<=nm;i++) {
    pli,mc(,,i),(i-1)%10,(i-1)/10,(i-1)%10+1,(i-1)/10+1;
    plt,swrite(format="%d",i),(i-1)%10,(i-1)/10,tosys=1,height=6,opaque=0,color="black";
  }
  limits,-0.5,10.5,-0.5,10.5;
}
func opra_gen_modes(dim,pupd,nmodes,mode_type=)
/* DOCUMENT opra_gen_modes(dim,pupd,nmodes,mode_type=)
   Build a set nmodes of KL (mode_type="kl") or Zernike (mode_type="zernike") or disk harmonic (mode_type="dh")
   Returns modes_cube, a cube of modes.
   modes_cube(,,1) = piston (also for kl=1)
   modes_cube(,,2) = tip (also for kl)
   Plane n = Zernike(n) / formally kl(n-1)
   Input parameters:
   dim: size of the array
   pupd: diameter of pupil on which modes are defined.
   nmodes: # of modes
   mode_type: select type of mode
   Note that I use zernike_ext (for zernike) and nopup for KL, so there
   shouldn't be discontinuities at the edge.
   SEE ALSO:
 */
{
  extern sdimold,modes_cube;

  if (mode_type == []) mode_type = "kl";

  if (mode_type=="kl") {
    sdim = long(ceil(pupd/2)*2);
    if (sdim!=sdimold) {
      modes_cube = array(0.0f,[3,dim,dim,nmodes]);
      tmp = make_kl(nmodes-1,sdim,nopup=1,verbose=0);
      modes_cube(dim/2-sdim/2+1:dim/2+sdim/2,dim/2-sdim/2+1:dim/2+sdim/2,2:)=tmp;
      sdimold = sdim;
    }
  }

  if (mode_type=="zernike") {
    modes_cube = array(0.0f,[3,dim,dim,nmodes]);
    if (dim%2 == 0) prepzernike,dim,pupd,dim/2+0.5,dim/2+0.5;
    else prepzernike,dim,pupd;
    for (i=1;i<=nmodes;i++) modes_cube(,,i) = zernike_ext(i);
    //for (i=1;i<=nmodes;i++) modes_cube(,,i) = zernike_ext(i);

  }

  if (mode_type=="dh") {
    modes_cube = array(0.0f,[3,dim,dim,nmodes]);
    modes_cube = make_diskharmonic(dim,pupd,nmodes,xc=dim/2+0.5,yc=dim/2+0.5);
    // special treatment for tip and tilt like modes.
    // we replace them by the zernike tip and tilt:
    if (dim%2 == 0) prepzernike,dim,pupd,dim/2+0.5,dim/2+0.5;
    else prepzernike,dim,pupd;
    modes_cube(,,2) = zernike(3);
    modes_cube(,,3) = zernike(2);
    //error;
    //tmp = make_kl(nmodes-1,sdim,nopup=1,verbose=0);
    //modes_cube(dim/2-sdim/2+1:dim/2+sdim/2,dim/2-sdim/2+1:dim/2+sdim/2,2:)=tmp;
  }

  return modes_cube;
}
