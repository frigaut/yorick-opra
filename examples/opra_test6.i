require,"opra.i";

modes     = "kl";      // 6.43e-3
modes     = "dh";      // 5.52e-3
modes     = "zernike"; // 4.92e-3
lambda    = 2.12e-6;   // [m]
pixsize   = 0.020;     // [arcsec]
allfocs   = span(-3,4,15)*0.45+0.1825;
nmodesmax = 21;

imcube = fits_read("gemini_niri.fits");
subset = [1,3,5,7,9,11,13,15];
allfocs = allfocs(subset);

opp = opra(imcube,allfocs,lambda,pixsize,7.9,nmodes=nmodesmax,use_mode=modes,\
  noise=0.0,cobs=0.15,progressive=0,first_nofit_astig=0,fix_kern=0,fix_pix=0,\
  niter=50,fix_defoc=0,dpi=140,gui=0,nm=1);
