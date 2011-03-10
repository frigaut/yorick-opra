require,"opra.i";

modes     = "yao";
lambda    = 1.65e-6; // [m]
pixsize   = 0.020;   // [arcsec]
allfocs   = [0.,1.52,-1.52]*1.5;
parfile   = "opra_test3.par";

imcube = [];
grow,imcube,fits_read("opra_test3.foc0.fits")(,,,-);
grow,imcube,fits_read("opra_test3.focp0.4.fits")(,,,-);
grow,imcube,fits_read("opra_test3.focm0.4.fits")(,,,-);

imcube = transpose(imcube,[3,4]);

// introduce/simulate some real world quirks:
// imcube(,,2,) = roll(imcube(,,2,),[0,1,0]);
// imcube(,,3,) = roll(imcube(,,3,),[0,-1,0]);
// imcube(,,1,2:) = roll(imcube(,,1,2:),[1,-1,0]);
imamp = array(0.0f,[2,3,5]);
imoff = array(0.0f,[3,2,3,5]);
for (i=1;i<=3;i++) {
  for (n=1;n<=5;n++) {
    imcube(,,i,n) *= (imamp(i,n)=(0.75+0.5*random()));
    imcube(,,i,n) = roll(imcube(,,i,n),(imoff(,i,n)=long([-2+4*random(),-2+4*random()])));
  }
}

fits_write,"opra_test3.imamp.fits",imamp,overwrite=1;
fits_write,"opra_test3.imoff.fits",imoff,overwrite=1;

opp = opra(imcube,allfocs,lambda,pixsize,7.9,nmodes=nmodesmax,use_mode=modes,\
  noise=0.0,cobs=0.,progressive=0,first_nofit_astig=0,fix_kern=0,fix_pix=0,\
  fix_diff_tt=1,niter=10,fix_defoc=0,dpi=100,yao_parfile=parfile);
