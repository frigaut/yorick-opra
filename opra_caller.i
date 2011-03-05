require,"opra.i";

if (case==[]) case = 3;
if (nmodesmax==[]) nmodesmax = 100;

if (case==1) {
  im1=fits_read("images/foc_z5100nm.fits")*1e6-0.;
  sumim1 = sum(im1);
  im1 = im1/sumim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc_z5100nm.fits")*1e6-1.60;
  im2 = im2/sumim1*2.2;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==0) {

  if (ns==[]) ns="1";
  im1=fits_read("images/foc_simu"+ns+".fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1*0.72;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc_simu"+ns+".fits")*1e6-1.36;
  im2 = im2/maxim1*0.72;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==2) {

  im1=fits_read("images/foc_z6100nm.fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc_z6100nm.fits")*1e6-1.36;
  im2 = im2/maxim1*2.2;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==3) {

  im1=fits_read("foc.fits");
  im1sky = 20.;
  im1 -= im1sky;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("defoc.fits");
  im2sky = 18.7;
  im2 -= im2sky;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==4) {

  im1=fits_read("images/foc_simu.fits");
  im1sky = 0.;
  im1 -= im1sky;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  // noise, sigma=1.5e-3
  noise = 20./maxim1;
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);
  
  im2=fits_read("images/defoc_simu.fits");
  im2sky = 0.;
  im2 -= im2sky;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==5) {

  im1=fits_read("images/foc_simu_new.fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc_simu_new.fits")*1e6;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);
  
 } else if (case==6) {

  im1=fits_read("images/foc_simu_new2.fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc_simu_new2.fits")*1e6;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  //  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=120,noise=noise);
  opp = opra([im1,im2],[0,2.25],810e-9,2.02e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==8) {
  // defoc defoc:  2.25,1.8,0.9,-0.6,-1.2,-2.0
  im1=fits_read("images/foc_simu_new4.fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc1_simu_new4.fits")*1e6;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  im3=fits_read("images/defoc2_simu_new4.fits")*1e6;
  im3 = im3/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im3 = clip(im3,0.,);

  im4=fits_read("images/defoc3_simu_new4.fits")*1e6;
  im4 = im4/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im4 = clip(im4,0.,);

  im5=fits_read("images/defoc4_simu_new4.fits")*1e6;
  im5 = im5/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im5 = clip(im5,0.,);
  opp = opra([im1,im2,im3,im4,im5],[0,2.25,1.8,0.9,-0.6],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==9) {
  // defoc defoc:  2.25,1.8,0.9,-0.6,-1.2,-2.0
  im1=fits_read("images/foc_simu_new4.fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc1_simu_new4.fits")*1e6;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==10) {
  files = findfiles("images/test_?.fits");
  files = files(sort(files));
  ima=[];
  
  for (i=1;i<=numberof(files);i++) {
    im = fits_read(files(i));
    if (i==1) maxim = sum(im);
    im = clip(im/maxim,0.,);
    if (i==1) noise = im(1:20,1:20)(*)(rms);
    grow,ima,im(,,-);
  }

  opp = opra(ima,[0,1.5,3.],810e-9,1.99e-3,7.9,nmodes=nmodesmax,noise=noise);

 } else if (case==11) {

  im1=fits_read("images/foc_simu_new3.fits")*1e6-0.;
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  noise = im1(1:20,1:20)(*)(rms);
  //  im1 = im1(225-128:225+127,225-128:225+127);
  im1 = clip(im1,0.,);

  im2=fits_read("images/defoc_simu_new3.fits")*1e6;
  im2 = im2/maxim1;
  //  im2 = im2(225-128:225+127,225-128:225+127);
  im2 = clip(im2,0.,);

  //  opp = opra([im1,im2],[0,2.25],810e-9,1.99e-3,7.9,nmodes=120,noise=noise);
  opp = opra([im1,im2],[0,2.25],810e-9,2.02e-3,7.9,nmodes=nmodesmax,noise=noise);

 }

if (case==7) {
  // NACO PD playkit
  im1=fits_read("images/psf_0mm.fits");
  maxim1 = sum(im1);
  im1 = im1/maxim1;
  // noise, sigma=1.5e-3
  noise =  0.;
  im1 = clip(im1,0.,);
  
  im2=fits_read("images/psf_3mm.fits");
  im2 = im2/maxim1;
  im2 = clip(im2,0.,);
  
  im3=fits_read("images/psf_4mm.fits");
  im3 = im3/maxim1;
  im3 = clip(im3,0.,);
  
  ima1 = ima2 = ima3 = array(0.,[2,256,256]);
  ima1(128-64+1:128+64,128-64+1:128+64) = im1;
  ima2(128-64+1:128+64,128-64+1:128+64) = im2;
  ima3(128-64+1:128+64,128-64+1:128+64) = im3;

  opp = opra([ima1,ima2,ima3],[0,3.,4.]*0.5,2.16e-6,13.25e-3,8.0,\
             nmodes=nmodesmax,cobs=0.14,noise=noise);
 }

/* notes:
- advantages w.r.t focal plane method:
  - all or nearly all pixels are used (!= a few speckles in FP)
  - larger dynamical range (e.g. look at TT) 
*/

/*
  TO DO:
  - thiebault's papers: has he done the same.
  - [doesn't work] finish iteration: go on phase pixels
  - [done] modify lmfit to pass structure -> much easier parsing of
    parameters in opra_foo() (no need to split a). Is it possible
    in yorick with structure functions?
  - [partially done] and, of course, find a way to clean up the phase! unwrap?
    add a regularisation term in lmfit? final iter with phase pixels
    + unwrap will solve? 
  - include effect of pixel averaging (noticeable if fwhm close to 2).
  - pixel size: rebin FTO
  - wrapping when phase too large
 */

/* cases:
 */
