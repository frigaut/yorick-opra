// cwd = pwd();
// cd,"/home/frigaut/mcao/myst/yorick";
// require,"yao_mcao.i";
// cd,cwd;
if (findfiles("zeroes.fits")==[]) fits_write,"zeroes.fits",array(0,[2,256,256]),overwrite=1;

require,"yao.i";
aoread,"opra_test3.par";

sim.verbose=0;
aoinit,disp=0,clean=1;
loop.gain=0.;
aoloop,disp=0;
go,20,all=1;
// (*dm(3)._command)(26) = 0.2; // units are in microns. so in rd would be 2*pi/target.lambda= this * 3.801
// (*dm(3)._command)(9) = 0.1; // units are in microns. so in rd would be 2*pi/target.lambda= this * 3.801
random_seed,0.424;
(*dm(1)._command) = random_n(dm(1)._nact)*2.*_(0.,noll(dm(1)._nact-1));
(*dm(2)._command) = random_n(dm(2)._nact)*6.*_(0.,noll(dm(2)._nact-1));
// (*dm(3)._command) = random_n(dm(3)._nact)/2./indgen(dm(3)._nact);
// (*dm(1)._command)(1:6) = 0;
(*dm(1)._command)(5) = 0;
(*dm(1)._command)(1:3) = 0;
//(*dm(2)._command)(15) = 0.1;
(*dm(2)._command)(1:3) = 0;
(*dm(2)._command)(5) = 0;
// (*dm(3)._command)(1:3) = 0; (*dm(3)._command)(5) = 0;
go,10,all=1;
// hitReturn;
fits_write,"opra_test3.foc0.fits",im,overwrite=1;
prepzernike,sim._size,sim.pupildiam,sim._cent,sim._cent;
add_dm0_shape = zernike_ext(4)*0.6; // same, that means 0.4*3.801 = 1.52
go,10,all=1;
// hitReturn;
fits_write,"opra_test3.focp0.4.fits",im,overwrite=1;
add_dm0_shape = -zernike_ext(4)*0.6;
go,10,all=1;
fits_write,"opra_test3.focm0.4.fits",im,overwrite=1;
fits_write,"opra_test3.dm1.fits",*dm(1)._command,overwrite=1;
fits_write,"opra_test3.dm2.fits",*dm(2)._command,overwrite=1;
quit;

/*
c1 = fits_read("opra_test3.dm1.fits");
c2 = fits_read("opra_test3.dm2.fits");
imamp = fits_read("opra_test3.imamp.fits");
imoff = fits_read("opra_test3.imoff.fits");
*/
