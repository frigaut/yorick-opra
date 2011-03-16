// cwd = pwd();
// cd,"/home/frigaut/mcao/myst/yorick";
// require,"yao_mcao.i";
// cd,cwd;
require,"yao.i";
aoread,"opra_test5.par";

sim.verbose=0;
aoinit,disp=1,clean=1;
loop.gain=0.;
aoloop,disp=1;
go,20,all=1;
// (*dm(3)._command)(26) = 0.2; // units are in microns. so in rd would be 2*pi/target.lambda= this * 3.801
// (*dm(3)._command)(9) = 0.1; // units are in microns. so in rd would be 2*pi/target.lambda= this * 3.801
(*dm(1)._command) = random_n(dm(1)._nact)*5.*_(0.,noll(dm(1)._nact-1));
(*dm(2)._command) = random_n(dm(2)._nact)*5.*_(0.,noll(dm(2)._nact-1));
(*dm(3)._command) = random_n(dm(3)._nact)*5.*_(0.,noll(dm(3)._nact-1));
// (*dm(3)._command) = random_n(dm(3)._nact)/2./indgen(dm(3)._nact);
// (*dm(1)._command)(1:6) = 0;
(*dm(1)._command)(5) = 0;
(*dm(1)._command)(1:3) = 0;
(*dm(2)._command)(5) = 0;
(*dm(2)._command)(1:3) = 0;
// (*dm(2)._command)(7:11) = 0.1;
//(*dm(2)._command)(15) = 0.1;
(*dm(3)._command)(1:6) = 0;
(*dm(3)._command)(5) = 0;
// (*dm(3)._command)(1:3) = 0; (*dm(3)._command)(5) = 0;
go,10,all=1;
hitReturn;
norm = 6000./max(im);
imn = poisson(im*norm)+gaussdev(dimsof(im));
fits_write,"opra_test5.foc0.fits",imn,overwrite=1;
prepzernike,sim._size,sim.pupildiam,sim._cent,sim._cent;
add_dm0_shape = zernike_ext(4)*0.6; // same, that means 0.4*3.801 = 1.52
go,10,all=1;
hitReturn;
imn = poisson(im*norm)+gaussdev(dimsof(im));
fits_write,"opra_test5.focp0.4.fits",imn,overwrite=1;
add_dm0_shape = -zernike_ext(4)*0.6;
go,10,all=1;
imn = poisson(im*norm)+gaussdev(dimsof(im));
fits_write,"opra_test5.focm0.4.fits",imn,overwrite=1;
quit;
