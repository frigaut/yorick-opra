require,"opra.i";

modes     = "yao";
lambda    = 1.65e-6; // [m]
pixsize   = 0.020;   // [arcsec]
allfocs   = [0.,1.52,-1.52]*1.5;
parfile   = "opra_test4.par";

require,"yao.i";
aoread,parfile;

imcube = [];
grow,imcube,fits_read("opra_test4.foc0.fits")(,,,-);
grow,imcube,fits_read("opra_test4.focp0.4.fits")(,,,-);
grow,imcube,fits_read("opra_test4.focm0.4.fits")(,,,-);

imcube = transpose(imcube,[3,4]);


opp = opra(imcube,allfocs,lambda,pixsize,7.9,nmodes=nmodesmax,use_mode=modes,\
  noise=0.0,cobs=0.,progressive=0,first_nofit_astig=0,fix_kern=0,fix_pix=0,\
  niter=10,fix_defoc=0,dpi=default_dpi);
