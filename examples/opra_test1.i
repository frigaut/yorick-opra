require,"opra.i";

nmodesmax = 45;      // number of modes
modes     = "dh";
lambda    = 1.68e-6; // [m]
pixsize   = 0.020;   // [arcsec]
focstep   = 0.833; // this is 200nm of z4
allfocs   = [2,0,-2]*focstep;

require,"yao.i";
aoread,"opra_test3.par";

imcube = fits_read("opra_test1_data.fits");

opp = opra(imcube,allfocs,lambda,pixsize,7.9,nmodes=nmodesmax,use_mode=modes,\
  noise=0.00004,cobs=0.,progressive=0,first_nofit_astig=0,fix_kern=0,fix_pix=0,\
  niter=60,fix_defoc=0,dpi=default_dpi,svipc=1,gui=1);
