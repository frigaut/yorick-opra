/* libkl

   Available functions:
   func dblindgen(n)
   func polar_coord(r,&mask,&rho,&phi,&pts,occ=,xcent=,ycent=,\
   func radii(nr,np,ri)
   func polang(r)
   func setpincs(ax,ay,px,py,ri,&pincx,&pincy,&pincw)
   func pcgeom (nr,np,ncp,ri,ncmar,ap)    
   func set_pctr(bas, ncp =, ncmar=)
   func pol2car(cpgeom,pol,mask=)
   func kolstf(dvec)
   func karmanstf(dvec,outscl=)
   func gkl_radii(nr,ri)
   func gkl_mkker(ri,nr,rad,funct=,outscl=)
   func piston_orth(nr)
   func gkl_fcom(kers,ri,nf,&evals,&nord,&npo,&ord,&rabas)
   func gkl_mkazi(nord, np)
   func gkl_bas(ri=,nr=,np=,nfunc=,verbose=,funct=,outscl=)
   func gkl_sfi(bas, i)
   func make_kl(nmax,dim,&var,&outpolarbase,&pupil,oc=,nr=,nopup=,\
   func kl_basis_in_dm_space_4extrap(nm, n_rm_modes);
   func kl_basis_in_dm_space(nm,n_rm_modes,&eigen_val);

*/




/*
  Collection of routines for building the KL of a Kolmogorov statistics
  Derived from a collection of routines developped both at ESO and ONERA
  The main routine is the last one ...

  for a Kolmogorov statistics :
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64);
  or
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="kolmo");
  

  for a Von Karman statistics :
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman");
  or
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman",outscl=3);
  default is : an outter scale of 3 times the size of the telescope
  
  res is a 128x128x150 array containing the 150 first KL
  of a kolmogorov or Von Karman stat

  D. Gratadour Feb 2006
 */

//require,mcao_i_dir+"lib/libgenutil.i";
require,"digit2.i";

func dblindgen(n)
  /*DOCUMENT res=dblindgen(size)

  D. Gratadour Feb 2006
  
  This routine returns a size x size array containing increasing indices
  from 1 to size x size.

  SEE ALSO : indgen, indices
   */
{
  n=long(n);
  return reform(indgen(n*n),[2,n,n]);
}

func polar_coord(r,&mask,&rho,&phi,&pts,occ=,xcent=,ycent=,\
                 verbose=,leq=,btw4pix=,sizemin=,dbprec=)
  /* DOCUMENT polar_coord,radius,mask,rho,phi;
           or polar_coord,radius,mask,rho,phi,pts,sizemin=1;

  D. Gratadour Jan 2006

  Calculation of polar coordinates rho and phi and an intensity mask
  (the pupil) of a telescope of radius r
  Derived from an IDL routine (polaire2.pro) written by L. Mugnier
  
  INPUTS :
  r = the radius of the mask

  OUTPUT :
  mask = the intensity mask (2d image)
  rho  = the width coordinate (2d image)
  phi  = the angle coordinate (rad) (2d image)
  pts  = (optional) indices of valid (non null) points of the
         pupil (1d vector) if flag sizemin is set 

  OPTIONAL :
  occ     = the occultation level (<=1)
  xcent   = the x position of the center point of the mask 
  xcent   = the y position of the center point of the mask 
  verbose = flag to get info on the process (0/1)
  leq     = flag to set the limits of the mask (<= radius or < radius)
  btw4pix = flag to set the center of the mask on 1 pixel or
            in between 4 pixels (0/1)
  sizemin = flag to set the outputs in a minimum size arrays (null points
            are eliminated) (0/1)
     
  SEE ALSO : ...
  */
{
  if (r==[]) r=128;
  if (!is_set(occ)) occ=0.0;
  if (!btw4pix) diam=long(2*floor(r)+1);
  else diam=long(2*round(r));
  
  if (!is_set(xcent)) xcent=float((diam-1)/2.);  
  if (!is_set(ycent)) ycent=float((diam-1)/2.);

  if (verbose) {
    write,format="Pupil diameter : %d\n",diam;
    if (xycent!=[]) {
      write,format="Pupil X center : %f\n",xycent(1);
      write,format="Pupil Y center : %f\n",xycent(2);
    }
  }

  x=float(((indices(diam))(,,1)-1.0) % diam);
  y=transpose(x);
  x-=xcent;
  y-=ycent;

  if (dbprec) {
    rho=double(sqrt(x^2+y^2))/double(r);
    phi=double(atan(y,x+(rho==0)));
  }
  else {
    rho=float(sqrt(x^2+y^2))/float(r);
    phi=float(atan(y,x+(rho==0)));
  }
  
  if (leq) mask=((rho<=1) & (rho>=occ));
  else {
    if (dbprec) mask=double((rho<1) & (rho>=occ));
    else mask=float((rho<1) & (rho>=occ));
  }
  
  if (sizemin) {
    if (verbose) write,"We will keep only the non-null points";
    pts=where(mask != 0.);
    if (dimsof(pts) != []) {
      sizepts=(dimsof(pts))(2);
      if (verbose) write,"Number of null points : ",sizepts;
      rho=rho(pts);
      phi=phi(pts);
      mask=mask(pts);
    }
    else {
      if (verbose) write,"All points are non-null !";
    }
  }
}


struct gkl_basis_struct
{
  long   nr;            
  long   ni;  
  long   np;       
  long   nfunc;   
  float   ri;  
  pointer   radp;   
  pointer evals;   
  long   nord;   
  pointer   npo;   
  pointer   ord;   
  pointer   rabas;   
  pointer   azbas;   
};

struct geom_struct
{
  pointer   px;            
  pointer   py;  
  pointer   cr;       
  pointer   cp;   
  pointer   pincx;  
  pointer  pincy;   
  pointer pincw;   
  pointer   ap;   
  long   ncp;   
  long   ncmar;   
};


func radii(nr,np,ri)
  /*DOCUMENT res=radii(NumberOfR,NumberOfPhi,Dim)

  D. Gratadour Feb 2006
  
  This routine generates an nr x np array with np copies of the
  radial coordinate array. Radial coordinate span the range from
  r=ri to r=1 with successive annuli having equal areas (ie, the
  area between ri and 1 is divided into nr equal rings, and the
  points are positioned at the half-area mark on each ring). There
  are no points on the border.     

  SEE ALSO : polang
   */
{
  r2 = ri^2 +(float(indgen(nr)-1)+0.)/nr*(1.0 - ri^2);
  rs = sqrt(r2);
  return rs*array(1.,np)(-,);
}

func polang(r)
  /*DOCUMENT res=polang(RadialCoordArray)

  D. Gratadour Feb 2006
  
  This routine generates an array with the same dimensions as r,
  but containing the azimuthal values for a polar coordinate system.     

  SEE ALSO : radii
   */
{
  s =  dimsof(r);
  nr = s(2);
  np = s(3);
  phi1 = float(indgen(np)-1)/float(np)*2.*pi;
  return phi1(-,)*array(1.,nr);
}

func setpincs(ax,ay,px,py,ri,&pincx,&pincy,&pincw)
  /*DOCUMENT res=polang(RadialCoordArray)

  D. Gratadour Feb 2006
  
  This routine determines a set of squares for interpolating
  from cartesian to polar coordinates, using only those points
  with ri < r < 1     

  SEE ALSO : pcgeom
   */
{
  s = dimsof(ax);
  nc = s(2);
  s = dimsof(px);
  nr = s(2);
  np = s(3);
  dcar = (ax(nc) - ax(1)) / (nc-1);
  ofcar = ax(1,1);

  rlx = (px - ofcar)/dcar;
  rly = (py - ofcar)/dcar;
  lx = long(rlx);
  ly = long(rly);
  shx = rlx - lx;
  shy = rly - ly;

  pincx=[lx,lx+1,lx+1,lx]+1;
  pincy=[ly,ly,ly+1,ly+1]+1;
  
  pincw=[(1-shx)*(1-shy),shx*(1-shy),shx*shy,(1-shx)*shy];
  
  axy = ax^2 + ay^2;
  axyinap = clip(axy,ri^2.+1.e-3,0.999);
  sizeaxyinap=(dimsof(axyinap))(2);
  pincw = pincw*axyinap(pincx+(pincy-1)*sizeaxyinap);
  pincw = pincw*(1.0/pincw(,,sum))(,,-);
}

func pcgeom (nr,np,ncp,ri,ncmar,ap)    
  /*DOCUMENT geom=pcgeom(nr, np, ncp, ri, ncmar,ap)

  D. Gratadour Feb 2006
  
  This routine builds a geom_struct. px and py are the x and y
  coordinates of points in the polar arrays.  cr and cp are the
  r and phi coordinates of points in the cartesian grids. ncmar
  allows the possibility that there is a margin of ncmar points
  in the cartesian arays outside the region of interest
    

  SEE ALSO : setpincs, set_pctr
   */
{
  nused = ncp - 2*ncmar;
  ff = 0.5 * nused;
  hw =  float(ncp-1)/2;
    
  r = radii(nr,np,ri); 
  p = polang(r);

  px0 = r * cos(p);
  py0 = r * sin(p);
  px = ff * px0 + hw;
  py = ff * py0 + hw;
  ax = float(dblindgen(ncp)-1) % ncp - 0.5 * (ncp-1);
  ax = ax / (0.5 * nused); 
  ay = transpose(ax);
        
  setpincs, ax, ay, px0, py0, ri,pincx, pincy, pincw;
  dpi = 2 * pi;
  cr2 = (ax^2 + ay^2); 
  ap = clip(cr2,ri^2+1.e-3,0.999);
  // original:
  //cr = (cr2 - ri^2) / (1 - ri^2) * nr - 0.5; 
  // modified?
  cr = (cr2 - ri^2) / (1 - ri^2) * nr; 
  cp = (atan(ay, ax) + dpi) % dpi;
  cp = (np / dpi) * cp;
    
  cr = clip(cr,1.e-3,nr-1.001);
  //fudge -----, but one of the less bad ones
  cp = clip(cp,1.e-3,np -1.001);
  //fudge -----  this is the line which
  //gives that step in the cartesian grid
  //at phi = 0.
  
  geom = geom_struct();
  geom.px=&px;
  geom.py=&py; 
  geom.cr=&cr;
  geom.cp=&cp;
  geom.pincx=&pincx;
  geom.pincy=&pincy;
  geom.pincw=&pincw;
  geom.ap=&ap;
  geom.ncp=ncp;
  geom.ncmar=ncmar; 
         
  return geom;
}

func set_pctr(bas, ncp =, ncmar=)
  /*DOCUMENT geom=set_pctr(bas, ncp =, ncmar=)

  D. Gratadour Feb 2006
  
  This routine calls pcgeom to build a geom_struct with the
  right initializations. bas is a gkl_basis_struct built with
  the gkl_bas routine.
    
  SEE ALSO : pcgeom, setpincs, gkl_bas
   */
{
  if (!is_set(ncmar)) ncmar = 2;
  if (!is_set(ncp)) ncp = 128;
    
  return pcgeom(bas.nr,bas.np,ncp,bas.ri,ncmar,ap);
}

func pol2car(cpgeom,pol,mask=)
  /*DOCUMENT cart=pol2car(cpgeom, pol, mask=)

  D. Gratadour Feb 2006
  
  This routine is used for polar to cartesian conversion.
  pol is built with gkl_bas and cpgeom with pcgeom.
  However, points not in the aperture are actually treated
  as though they were at the first or last radial polar value
  -- a small fudge, but not serious  ?*******
  
  SEE ALSO : pcgeom, gkl_bas
   */
{
  if (sae) error;
  cd = bilinear(pol, *cpgeom.cr+1, *cpgeom.cp+1);
  if (mask!=[]) cd = cd*(*cpgeom.ap);
  return cd;
} 

func kolstf(dvec)
  /*DOCUMENT var=kolstf(dvec)

  D. Gratadour Feb 2006
  
  This routine returns the kolmogorov phase variance at spatial
  dimension (inverse of the spatial frequency) dvec
  
  SEE ALSO : 
   */
{
  return  6.88 * dvec^(5./3.);
}

func karmanstf(dvec,outscl=)
  /*DOCUMENT var=kolstf(dvec)

  D. Gratadour Feb 2006
  
  This routine returns the Von Karman phase variance at spatial
  dimension (inverse of the spatial frequency) dvec. Same as kolstf
  but with a correcting factor to account for the outter scale.
  The latter should be in units of telescope diameter
  
  SEE ALSO : 
   */
{
  if (dimsof(outscl)==[]) outscl = 3.;
  return 6.88 * dvec^(5./3.)*(1-1.485*(dvec/outscl)^(1./3.)+\
                              5.383*(dvec/outscl)^(2)-6.281*\
                              (dvec/outscl)^(7./3.));
}

func gkl_radii(nr,ri)
  /*DOCUMENT rad=gkl_radii(nr,ri)

  D. Gratadour Feb 2006
  
  This routine generates an array of radial polar coordinates along
  which the KL are generated. nr is the number of elements and ri is
  the maximum radius.
  
  SEE ALSO : 
   */
{
  d = (1.-ri*ri)/nr;
  //    rad2 = ri^2 + d/2. + d * float(indgen(nr)-1);
  //    rad2 = ri^2 + d * float(indgen(nr)-1);
  rad2 = ri^2 +d/16.+ d * float(indgen(nr)-1);
  //  rad2 = ri^2 +d/14.+ d * float(indgen(nr)-1);  // nr=64,128
  //  rad2 = ri^2 +d/10.+ d * float(indgen(nr)-1);
  rad = sqrt(rad2);
  
  return rad;
}

func gkl_mkker(ri,nr,rad,funct=,outscl=)
  /*DOCUMENT 

  D. Gratadour Feb 2006
  
  This routine generates the kernel used to find the KL modes.
  The  kernel constructed here should be simply a discretization
  of the continuous kernel. It needs rescaling before it is treated
  as a matrix for finding  the eigen-values. The outter scale
  should be in units of the diameter of the telescope.

  SEE ALSO : 
   */
{
  nth = 5*nr;
  kers  = array(float,[3,nr, nr, nth]);
  cth = cos(float(indgen(nth)-1)*(2.*pi/nth));
  dth = 2.*pi/nth;
  fnorm = -1./(2*pi*(1.-ri^2))*0.5;
  //the 0.5 is to give  the r^2 kernel, not
  //the r kernel
  for (i =1;i<=nr;i++) { 
    for (j=1;j<=i;j++) {
      te = 0.5*sqrt(rad(i)^2+rad(j)^2-(2*rad(i)*rad(j))*cth);
      //te in units of the diameter, not the radius
      if (funct=="kolmo") te = kolstf(te);
      if (funct=="karman") te = karmanstf(te,outscl=outscl);
      if ((funct!="kolmo") & (funct!="karman")) {
        write,"The statistics is not known !";
        error;
      }
      kelt =  fnorm * dth * float(fft(te,-1));
      kers (i, j,) = kelt;
      kers (j, i,) = kelt;
    }
    if (is_set(verbose)) write, i;
  }
  if (is_set (verbose))  write," ";
  
  return kers;
  
}

func piston_orth(nr)
{
  s = array(float,[2,nr,nr]);
  for (j=1;j<=nr-1;j++) {
    rnm = 1./sqrt (float((j)*(j+1)));
    s(1:j,j) = rnm;
    s(j+1,j)= -1*(j)*rnm;
  }
  rnm = 1./sqrt (nr);
  s(,nr) = rnm;
  return s;
}

func gkl_fcom(kers,ri,nf,&evals,&nord,&npo,&ord,&rabas)
  /*DOCUMENT 

  D. Gratadour Feb 2006
  
  This routine does the work : finding the eigenvalues and
  corresponding eigenvectors. Sort them and select the right
  one. It returns the KL modes : in polar coordinates : rabas
  as well as the associated variance : evals. It also returns
  a bunch of indices used to recover the modes in cartesian
  coordinates (nord, npo and ord).

  SEE ALSO : gkl_bas
   */
{
  s = dimsof(kers);
  nr = s(2);
  nt = s(4);
  nxt = 1;
  fktom =  (1.-ri^2)/nr;
  fevtos = sqrt(2*nr);
  evs = array(float,[2,nr,nt]);
  //ff isnt used - the normalisation for
  //the eigenvectors is straightforward:
  //integral of surface^2 divided by area = 1,
  //and the cos^2 term gives a factor
  //half, so multiply zero order by
  //sqrt(n) and the rest by sqrt (2n)

  //zero order is a special case...
  //need to deflate to eliminate infinite eigenvalue - actually want
  //evals/evecs of zom - b where b is big and negative
  zom = kers(,,1);
  s = piston_orth(nr);
  ts =transpose(s);
  b1 = ((ts(,+)*zom(+,))(,+)*s(+,))(1:nr-1, 1:nr-1);
 
  newev = SVdec(fktom*b1,v0,vt);

  v1 = array(float,[2,nr, nr]);
  v1(1:nr-1,1:nr-1) = v0;
  v1(nr,nr) = 1;

  vs = s(,+)*v1(+,);
  grow,newev,0;
  evs(,nxt) = float(newev);
  kers (,, nxt) = sqrt(nr)*vs;
  // the rest are more straightforward
  nxt = 2;
  do {
  newev = SVdec(fktom*kers(,,nxt),vs,vt);
  evs(,nxt) = float(newev);
  kers (,,nxt) = sqrt(2.*nr)*vs;
  mxn = max(float(newev));
  egtmxn = floor(evs(, 1:nxt)>mxn);
  nxt = nxt + 1;
  } while ((2*sum(egtmxn)-sum(egtmxn(,1))) < nf);
  nus = nxt - 1;

  kers = kers (,,1:nus);
  evs = reform (evs (, 1:nus), nr*nus);
  a = (sort(-1.*evs))(1:nf);
  //every eigenvalue occurs twice except
  //those for the zeroth order mode. This
  //could be done without the loops, but
  //it isn't the stricking point anyway...
  no = 1;
  ni = 1;
  oind = array(long,nf+1);
  do {
       if (a(ni) < nr+1) {
         oind(no) = a(ni);
         no = no + 1;
       } else {
         oind(no) = a(ni);
         oind(no+1) = a(ni);
         no = no + 2;
       }
       ni = ni + 1;
  } while (no < (nf+1));
  
  oind = oind (1:nf);
  tord = (oind-1)/nr+1;
  odd = ((long(indgen(nf)-1) % 2) == 1);
  pio = (oind-1) % nr +1;

  evals = evs(oind);
  ord = 2 *(tord-1) - floor(tord>1 & (odd))+1;

  nord = max(ord);
  rabas = array(float,[2,nr, nf]);
  sizenpo=long(max(ord));
  npo = array(long,sizenpo);
  
  for (i=1;i<=nf;i++) {
    npo(long(ord(i))) = npo(long(ord(i))) + 1;
    rabas(, i) = kers (, pio(i), tord(i));
  }
}

func gkl_mkazi(nord, np)
{
  gklazi = array(float,[2,long(1+nord), np]);
  th = float(indgen(np)-1)*(2.*pi/ np);

  gklazi (1,) = 1.0;
  for (i = 2; i<=nord;i+=2)  gklazi (i,) = cos (((i-1)/2+1) * th);
  for (i = 3; i<=nord;i+=2)  gklazi (i,) = sin (((i-1)/2) * th);
  return gklazi;
}

func gkl_bas(ri=,nr=,np=,nfunc=,verbose=,funct=,outscl=)
  /*DOCUMENT 

  D. Gratadour Feb 2006
  
  This routine uses the output of gkl_fcom to fill the gkl_base_struct.

  SEE ALSO : gkl_fcom
   */
{
  if (!is_set(ri)) ri = 0;
  if (!is_set(nr)) nr = 40;
  if (!is_set(np)) np = long(5*nr);
  if (!is_set(nfunc)) nfunc = 500L;

  nr = long(nr);
  np = long(np);

  if ((nr * np)/ nfunc < 8) {
    if (is_set(verbose)) write,"warning: you may need a finer ",\
                           "radial sampling ";
    if (is_set(verbose)) write, "(ie, increased nr) to generate ",\
                           nfunc, "  functions";
  } else if ((nr * np)/ nfunc > 40) {
    if (is_set(verbose)) write,"note, for this size basis ",\
                           "radial discretization on ", nr;
    if (is_set(verbose)) write, "points is finer than necessary",\
                           "-it should work, but you ";
    if (is_set(verbose)) write, "could take a smaller nr without",\
                           "loss of accuracy";
  }


  radp = gkl_radii(nr, ri);

  kers = gkl_mkker(ri, nr, radp,funct=funct,outscl=outscl);
  
  gkl_fcom,kers,ri,nfunc,evals,nord,npo,ord,rabas;

  azbas = gkl_mkazi(nord, np);

  gklbasis = gkl_basis_struct();
  gklbasis.nr=nr;
  gklbasis.np=np;
  gklbasis.nfunc=nfunc; 
  gklbasis.ri=ri;
  gklbasis.radp=&radp;
  gklbasis.evals=&evals;
  gklbasis.nord=nord;
  gklbasis.npo=&npo;
  gklbasis.ord=&ord;
  gklbasis.rabas=&rabas;
  gklbasis.azbas=&azbas;
  
  return gklbasis;
}

func gkl_sfi(bas, i)
  /*DOCUMENT 

  D. Gratadour Feb 2006
  
  This routine returns the i'th function from the generalised KL
  basis bas. bas must be generated first with gkl_bas.

  SEE ALSO : gkl_bas
   */
{    
  if (i>bas.nfunc) { 
    write, "the basis only contains ", nfunc, "functions";
    return 0;
  }
       
  nr = bas.nr;
  np = bas.np;
  ordp = *bas.ord;
  ord=long(ordp(i));

  rabasp=*bas.rabas;
  rabas=rabasp(,i);

  azbasp=*bas.azbas;
  azbas=azbasp(ord, );

  sf1=array(double,[2,nr,np]);
  sf1(,*)=rabas;

  sf2=array(float,[2,np,nr]);
  sf2(,*)=azbas;  

  sf = sf1*transpose(sf2);
  return sf;
}

      
func make_kl(nmax,dim,&var,&outpolarbase,&pupil,oc=,nr=,nopup=,\
             funct=,outscl=,verbose=)
/* DOCUMENT 
  for a Kolmogorov statistics :
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64);
  or
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="kolmo");
  

  for a Von Karman statistics :
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman");
  or
  res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman",outscl=5);
  
  the outter scale is in units of the telescope diameter
  default is : an outter scale of 3 times the size of the telescope
    
  D. Gratadour Feb 2006
  
  This routine is the main program. It returns nmax generalized
  KL in an array dim x dim x nmax. It also returns the associated
  variance as well as the pupil and the polar base used for their
  calculation. Optional keywords includes any occultation, the
  number of samples for the radial coordinate and a flag to avoid
  pupil multiplication.

  SEE ALSO : polar_coord, gkl_bas, set_pctr
*/
{
  if (pupil==[]) polar_coord,dim/2.,pup,rho,phi,occ=oc,btw4pix=1;
  else pup=pupil;
  
  if (!is_set(nr)) nr=64;

  if (dimsof(funct)==[]) {
    if (verbose!=0) write,format="%s\n","make_kl -> Using the Kolmogorov model";
    funct="kolmo";
  }
    
  polarbase = gkl_bas(ri=oc,nr=nr,np=(2*pi*nr),nfunc=nmax,\
                      funct=funct,outscl=outscl,verbose=verbose);
  
  outpolarbase = polarbase;
  
  pc1 = set_pctr(polarbase, ncp= dim);

  kl = array(float,[3,long(dim),long(dim),nmax]);
  
  if (is_set(nopup)) {
    for (i=1;i<=nmax;i++) kl(,,i)=pol2car(pc1, gkl_sfi(polarbase,i));
  } else {
    for (i=1;i<=nmax;i++) {
      sae=0;
      //      if (i==8) sae=1;
      kl(,,i)=pol2car(pc1, gkl_sfi(polarbase,i))*pup;
    }
  }
    
  pupil =  pup; 
  var =  *polarbase.evals;
  
  return kl;
}


func kl_basis_in_dm_space(nm,n_rm_modes,fudge,&eigen_val,extrap=)
/* DOCUMENT
   calcule la base optimale vis a vis des fonctions d'influence
   du mirroir deformable et des conditions de turbulence.
   To be used then for instance to compute MMSE reconstructors.

   kl2dm = array(pointer,3);
   ei_val = array(pointer,3);
   filt_tab = [3,3,3];
   for(nm=1;nm<=3;nm++) {
     b=kl_basis_in_dm_space(nm, filt_tab(nm),eigen,extrap="extrap_kl.mat");
     kl2dm(nm) = &(b);
     ei_val(nm) = &(eigen);}

   mat_kl2dm = array(float, [2,684,684]);
   mat_kl2dm(1:240,1:240) = *kl2dm(1);
   mat_kl2dm(241:564,241:564) = *kl2dm(2);
   mat_kl2dm(565:,565:) = *kl2dm(3);
   fitsWrite, "mat_kl2dm_filt111.fits",mat_kl2dm;
   tab_ei_val = float(684);
   tab_ei_val = _(*ei_val(1),*ei_val(2),*ei_val(3));
   fitsWrite, "tab_ei_val_filt333.fits",tab_ei_val;
 
     
   b = kl_basis_in_dm_space(1, 3,eigen,extrap="extrap_kl.mat");
   will return a 240x239 array
   plot, b(,1); //first mode = Astig (piston/TT have been filtered)
   to see the modes:
   inf_fun = fitsRead("KLDMmodes_DM1_IF_nrmodes3.fits");
   if_nb   = dimsof(inf_fun)(4);//nb of IF
   dim = dimsof(inf_fun)(2);//size of the support
   dm_mode = b(+,)*(inf_fun(*,))(,+);
   dm_mode2 = reform(dm_mode,if_nb-3,dim,dim);
   tv,dm_mode2(1,,);
   nmodes = 40;
   aff_kldm_basis,dm_mode2,nmodes,puptel=[];

   eigen = eigenvalues associated with each mode
   window,1;nm=1;
   xx = indgen(mcaodm(nm).nvalid);
   plot,eigen,xx;
   plg,xx^(-11./6.)*45.,xx,color="red";
   logxy,1,1;
*/
  
{

  coeffi = [1.,1.24,1.47];//factor to take into account D vs h
  L0 = 25.;
  r0 = 0.4;
  k = 2;//oversampling factor for fft of IF
  //2 should be enough, larger than 2 takes a while

  //inf_fun = 3D array of influence functions
  inf_fun = get_mInfluence_withextrap(nm,extrap=extrap);
  if_nb   = dimsof(inf_fun)(4);//nb of IF
  dim = dimsof(inf_fun)(2);//size of the support

  
  //----Step0 : Modes Filtering-----------------
  
  nnn = dimsof(ipupil)(2);
  if(nm == 1) {
    if(0==dim%2) {puptel = ipupil(nnn/2-dim/2+1:nnn/2+dim/2,nnn/2-dim/2+1:nnn/2+dim/2);}
    if(1==dim%2) {puptel = ipupil(nnn/2-dim/2+1:nnn/2+dim/2+1,nnn/2-dim/2+1:nnn/2+dim/2+1);}
    prepzernike,dim,sim.pupildiam,dim/2+0.5,dim/2+0.5;
    pupkl = zernike(1);
  }

  
  if(nm > 1) {
    pup = abs(inf_fun(,,sum));
    wp = where(pup >= max(pup)*0.75);
    pup2 = pup*0.;pup2(wp)=1.;
    dim2=(where(pup2(,dim/2)==1)(0)-where(pup2(,dim/2)==1)(1)+1);
    prepzernike,dim,dim2,dim/2+0.5,dim/2+0.5;
    pupkl = zernike(1);
    puptel = pupkl;
  }
 
  //define IF on the pupil only
  for(cpt=1;cpt<=if_nb;cpt++){ 
    inf_fun(,,cpt) = inf_fun(,,cpt)*puptel;
  }

  //filter modes on the pupil only:
  if (numberof(n_rm_modes) != 0) {
    if(n_rm_modes > 0) {
      for(cpt=1;cpt<=if_nb;cpt++){ 
        cmpnt = 0.0;
          for(kj=1;kj<=n_rm_modes;kj++){
            polz  = zernike(kj)*puptel;
            coef  = sum(inf_fun(,,cpt)*polz)/sum(polz*polz);
            cmpnt = cmpnt + polz*coef;
          }
        inf_fun(,,cpt) = inf_fun(,,cpt) - cmpnt;
      }
    }
  }
  


  //save influence functions for latter:
  if(n_rm_modes == []) n_rm_modes = 0;
  name=swrite(format="KLDMmodes_DM%d_IF_nrmodes%d.fits",nm,n_rm_modes(1));
  fitsWrite,name,inf_fun;
  
  //-------------------------------------------------------
  //Step1 : Computing geometrique covariance matrix ...
  Delta_IF = array(float,[2,if_nb, if_nb]);
  Spup     = numberof(where(puptel));
  nrm      = Spup;//pupil surface in pixel;
  tmp = inf_fun(*,);
  Delta_IF = (tmp(+,)*tmp(+,))/nrm;
  tmp = [];
  write, "-> Géometrique Covariance, DONE !";

  //---------------------------------------------------------
  //Step2 : Calcul des TF des fction d'influence
  phase_variance = (gamma(11./6.)*gamma(5./6.)/(2.*pi^(8./3.)))*
    (24.*gamma(6./5.)/5.)^(5./6.)*
    (L0/r0)^(5./3);

  if(nm==1) coco = 50.;
  if(nm==2) coco = 25.;
  if(nm==2) coco = 50.;

  kdpix = long(dim/2.+50.)*k;//at least 2 times pupil size !
  if(kdpix <= dim) write, "ATTENTION!"

  
  //-------------------------------------
  D = 8.*coeffi(nm);
  //sp_freq        = dist(kdpix)/(k*D*float(dim2)/float(sim.pupildiam));//FIXME
  //This is not good! should be in real pix/m
  //fudge = 5.;
  sp_freq        = dist(kdpix)/(kdpix*D/dim*fudge);
  //------------------------------------

  
  f = sp_freq;
  cst   = (gamma(11./6.)^2/(2.*pi^(11./3.)))*(24.*gamma(6./5.)/5.)^(5./6.) ;
  phase_spectrum = eclat(cst*r0^(-5./3.)*(f^2+(1/L0)^2)^(-11./6.)) ;


  //error;
  
  //print, "Variance theorique [rd^2]: ", (phase_variance);
  //phase_variance_from_spectrum = sum(double(phase_spectrum))/(k*D)/(k*D);
  //print, "Variance from dsp  [rd^2]: ", (phase_variance_from_spectrum);
  //error;
  support = array(float,[2,kdpix, kdpix]);
  FT_inf_fun = array(complex,[3,kdpix, kdpix, if_nb]);
  for(i=1;i<=if_nb;i++){
    support(1:dim,1:dim) = inf_fun(,,i);
    FT_inf_fun(,,i) = fft(gaussm(support,0.5),1)*sqrt(phase_spectrum);
  }
  write, "-> TF des IF, DONE !";
  //error;

  inf_fun = support = [];
  //---------------------------------------------------------
  //Step3 : correlation statistique des fonctions d'influence
  
  H_IF = array(float,[2,if_nb, if_nb]);
  nrm            = (Spup*Spup)*(k*D)*(k*D);

  tmp = FT_inf_fun(*,);
  FT_inf_fun = [];
  H_IF = float(tmp(+,)*conj(tmp)(+,))/nrm;
  write, "-> Correlation Statistiques, DONE !";
  //error;
  //-------------------------------------------------------------
  //Step4 : Double Diagonalisation

  D2 = SVdec(Delta_IF, Mp);//decomposition en vecteur propre/valeur propre
  Mp = Mp(,:-n_rm_modes);//remove as many modes as we have removed before
  D2 = D2(1:-n_rm_modes);
  print,"conditionning:",max(D2)/min(D2);
  M = transpose(Mp/sqrt(D2)(-,));//on pondere les modes propres par leur valp
  Hp = M(,+)*(H_IF(,+)*M(,+))(+,);//on projete H_IF dans l'espace M
  Lp = SVdec(Hp,A);//on le diagonalise
  print,"conditionning:",max(Lp)/min(Lp);
  Bp = M(+,)*A(+,);
  eigen_val = Lp;
  write, "-> Double Diago, DONE !";
  //error;
  return Bp;

  /*
  //-----------------------------------
  mcao_restore_dm,mysttop+"/mystInitFiles/mcao-prod.ybin";

  nnm = 1;
  spup=ipupil(dm(nnm)._n1:dm(nnm)._n2,dm(nnm)._n1:dm(nnm)._n2);
  i=1;
  tv,compDmShape(nnm,&(float((*kl2dm(nnm))(,i))),extrap=0);
   */
  
 }



func kl_basis_in_dm_space_4extrap(nm, n_rm_modes)
/* DOCUMENT
   calcule la base optimale vis a vis des fonctions d'influence
   du mirroir deformable et des conditions de turbulence.
   To be used then for instance with gen_extrap_klopt

   b = kl_basis_in_dm_space_4extrap(1, 1);
   will return a 293x293 array
*/
  
{
  //-----------------------------------
  //Get Influence Functions WITH extrapolated
  //and on a large array
  inf_fun = get_mInfluence_large(nm);//it includes extrap and valids
  if_nb   = dimsof(inf_fun)(4);//nb of IF
  dim = dimsof(inf_fun)(2);//size of the support

  def1 = def2 = def3 = def4 = def5 = [];
  edef1 = edef2 = edef3 = edef4 = edef5 = [];
  dm(1)._def = &def1;dm(2)._def = &def2;dm(3)._def = &def3;
  dm(4)._def = &def4;dm(5)._def = &def5;dm(1)._edef = &edef1;
  dm(2)._edef = &edef2;dm(3)._edef = &edef3;
  
  //--------------------------------------------
  //Step0 : Modes Filtering

  //First, we have to define a pupil...
  pup = inf_fun(,,sum);
  pup = sign(pup(wheremax(abs(pup))(1)))*pup;
  pup = (pup > 0.8*max(pup));
  dim2=where(pup(,dim/2)==1)(0)-where(pup(,dim/2)==1)(1)+1;
  prepzernike,dim,dim2,dim/2+0.5,dim/2+0.5;
  pupkl = zernike(1);

   //define IF on the pupil only
  for(cpt=1;cpt<=if_nb;cpt++){ 
    inf_fun(,,cpt) = inf_fun(,,cpt)*pupkl;
  }

  
  if (numberof(n_rm_modes) != []) {
    if(n_rm_modes > 0) {
      for(cpt=1;cpt<=if_nb;cpt++){ 
        cmpnt = 0.0
          for(kj=1;kj<=n_rm_modes;kj++){
            polz  = zernike(kj)
            coef  = sum(inf_fun(,,cpt)*polz)/sum(polz*polz)
            cmpnt = cmpnt + polz*coef
          }
        inf_fun(,,cpt) = inf_fun(,,cpt) - cmpnt;
      }
    }
  }

  polz = [];
  //-------------------------------------------------------
  //Step1 : Computing geometrique covariance matrix ...
  Delta_IF = array(float,[2,if_nb, if_nb]);
  Spup     = numberof(where(pupkl));
  nrm      = Spup;//pupil surface in pixel;
  tmp = inf_fun(*,);//(w,);
  Delta_IF = (tmp(+,)*tmp(+,))/nrm;
  write, "-> Géometrique Covariance, DONE !";
  tmp = [];

  //---------------------------------------
  //Step1b: Spectre turbulent:
  k = 2;//2 should be enough
  kdpix = long(dim/2.+10.)*k;//at least 2 times pupil size !
  if(kdpix <= dim) write, "ATTENTION!";

  L0 = 25.;
  r0 = 0.4;

  phase_variance = (gamma(11./6.)*gamma(5./6.)/(2.*pi^(8./3.)))*
    (24.*gamma(6./5.)/5.)^(5./6.)*
    (L0/r0)^(5./3);
  coeffi = [1.,1.24,1.47];//factor to take into account D vs h
  D = 8.*coeffi(nm);
  sp_freq        = dist(kdpix)/(k*D);
  f = sp_freq;
  cst   = (gamma(11./6.)^2/(2.*pi^(11./3.)))*(24.*gamma(6./5.)/5.)^(5./6.) ;
  phase_spectrum = eclat(cst*r0^(-5./3.)*(f^2+(1/L0)^2)^(-11./6.)) ;
  print, "Variance theorique [rd^2]: ", (phase_variance);
  phase_variance_from_spectrum = sum(double(phase_spectrum))/(k*D)/(k*D);
  print, "Variance from dsp  [rd^2]: ", (phase_variance_from_spectrum);



  //---------------------------------------------------------
  //Step2 : Calcul des TF des fonction d'influence

  FT_inf_fun = calc_FT_inf_fun(inf_fun,kdpix,if_nb,phase_spectrum);
  
  /*support = array(float,[2,kdpix, kdpix]);
  FT_inf_fun = array(complex,[3,kdpix, kdpix, if_nb]);
  for(i=1;i<=if_nb;i++){
    support(1:dim,1:dim) = inf_fun(,,i);
    FT_inf_fun(,,i) = fft(support,1)*sqrt(phase_spectrum);
  }
  write, "-> Calcul des TFs, DONE !";*/
  inf_fun = support = phase_spectrum = f = sp_freq = [];
  
  //---------------------------------------------------------
  //Step3 : correlation statistique des fonctions d'influence

  H_IF = array(float,[2,if_nb, if_nb]);
  nrm            = (Spup*Spup)*(k*D)*(k*D);
  tmp = FT_inf_fun(*,);
  FT_inf_fun = [];
  H_IF = float(tmp(+,)*conj(tmp)(+,))/nrm;
  write, "-> Correlation Statistiques, DONE !";

  //-------------------------------------------------------------
  //Step4 : Double Diagonalisation
  D2 = SVdec(Delta_IF, Mp);
  print,"conditionning:",max(D2)/min(D2);
  M = transpose(Mp/sqrt(D2)(-,));
  Hp = M(,+)*(H_IF(,+)*M(,+))(+,);
  Lp = SVdec(Hp,A);
  print,"conditionning:",max(Lp)/min(Lp);
  Bp = M(+,)*A(+,);
  eigen_val = Lp;
  write, "-> Double Diago, DONE !";
  
  return Bp;

  /*To check that the extrapolator is really doing what we want:
    mmm = 2;
    
    nb_modes = 150;
    condy = 5000;

    //and compute the extrapolator:
    bval    = b(*mcaodm(nm).valid_ptr,1:nb_modes);
    bext    = b(*mcaodm(nm).extrap_ptr,1:nb_modes);
    bval_i  = (svd_inverse(bval(+,)*bval(+,),condy))(,+) * bval(,+);
    bval_i  = (LUsolve(bval(+,)*bval(+,)))(,+) * bval(,+);
    extrapbn  = bext(,+)*bval_i(+,);

    //Lets try to test it ?
    com1 = bval(,mmm);
    ecom = extrapbn(,+)*com1(+);
    comtot = array(float,mcaodm(nm).nact);
    comtot(*mcaodm(nm).valid_ptr)=com1;
    comtot(*mcaodm(nm).extrap_ptr)=ecom;
    comref = b(,mmm);

    window,0;
    plot, comref;
    plg, comtot,color = "red";
    window,1;
    plot, abs(comref - comtot);
  */
}



func aff_kldm_basis(dm_mode2,nmodes,puptel=)
/* DOCUMENT
   Display zernike, KL or KLDM modes 
*/
  
{

  window, 49;
  winkill,49;
  window,49,dpi=120,style="nobox.gs";

  if(nmodes <= 5){
    n1=long(nmodes);
    n2=1;}

  if(nmodes > 5){
    n1=5;
    n2=nmodes/n1;}

  if(nmodes > 40){
    n1=8;
    n2=nmodes/n1;}

  if(nmodes > 100){
    n1=12;
    n2=nmodes/n1;}

   if(nmodes > 200){
    n1=16;
    n2=nmodes/n1;}

   if(nmodes > 300){
    n1=18;
    n2=nmodes/n1;}
  
   if(puptel == []) puptel = array(float,[2,dimsof(dm_mode2)(3),dimsof(dm_mode2)(3)])+1.;  
  tab = array(float,[2,dimsof(dm_mode2)(3)*n1,dimsof(dm_mode2)(3)*n2]);

  for(i=1;i<=n2;i++){
    for(j=1;j<=n1;j++){
      if((i-1)*n1+j > dimsof(dm_mode2)(2)) {
        tab((j-1)*dimsof(dm_mode2)(3)+1:j*dimsof(dm_mode2)(3),(i-1)*dimsof(dm_mode2)(3)+1:i*dimsof(dm_mode2)(3)) = 0.*puptel;}
      if((i-1)*n1+j <= dimsof(dm_mode2)(2)) {
      tab((j-1)*dimsof(dm_mode2)(3)+1:j*dimsof(dm_mode2)(3),(i-1)*dimsof(dm_mode2)(3)+1:i*dimsof(dm_mode2)(3)) = dm_mode2((i-1)*n1+j,,)*puptel;
      }
    }
  }
  tv,tab;
  limits, square=1;
  palette,"stern.gp";

    }




func klmodes_4extrap(nm, n_rm_modes,nsamp)
/* DOCUMENT
   Compute the kl-dm modes for DM #nm, following J-P. Veran approach.
   We compute covariance matrices from turbulence based on random
   samples. So nsamp should be large to make the covariance matrices
   to converge.
   It returns directly the Extrap Matrix fo DM #nm:
   
   E = klmodes_4extrap(1,0,10000); 
*/
  
{
  inf_fun = get_mInfluence_large(nm);//it includes extrap and valids
  if_nb   = dimsof(inf_fun)(4);//nb of IF
  dim = dimsof(inf_fun)(2);//size of the support

  //-----------------------
  //Step0 : Modes Filtering
  
  //First, we have to define a pupil...
  pup = inf_fun(,,sum);
  pup = sign(pup(wheremax(abs(pup))(1)))*pup;
  pup = (pup > 0.8*max(pup));
  dim2=where(pup(,dim/2)==1)(0)-where(pup(,dim/2)==1)(1)+1;
  prepzernike,dim,dim2,dim/2+0.5,dim/2+0.5;
  pupkl = zernike(1);

  if (numberof(n_rm_modes) != []) {
    if(n_rm_modes > 0) {
      for(cpt=1;cpt<=if_nb;cpt++){ 
        cmpnt = 0.0
          for(kj=1;kj<=n_rm_modes;kj++){
            polz  = zernike(kj)
            coef  = sum(inf_fun(,,cpt)*polz)/sum(polz*polz)
            cmpnt = cmpnt + polz*coef
          }
        inf_fun(,,cpt) = inf_fun(,,cpt) - cmpnt;
      }
    }
  }


  //-------------------------------------------------------
  //Step1 : Computing the projection matrix:
  Delta_IF = array(float,[2,if_nb, if_nb]);
  w        = where(pupkl);
  Spup     = numberof(w);
  nrm      = Spup;//pupil surface in pixel;

  tmp = inf_fun(*,)(w,);
  Delta_IF = (tmp(+,)*tmp(+,))/nrm;
  tmp2 = svd_inverse(Delta_IF,10000,e1);
  MatProj = tmp2(+,)*tmp(,+);
  write, "-> MatProj, DONE !";


  //-------------------------------------------------------
  //Step2 : Projection of turbulent screens on Influence Functions.
  l0 = 25;//en metre
  l0 = dim*l0/8.;//conversion en PIXELS
  Proj=array(float,[2,if_nb,nsamp]);

  if (nsamp == []) nsamp = 1000;

   if(dim >= 256) dimp = 512 ;
   if(dim < 256) dimp = 256 ;
  
  
  for (i=1;i<=nsamp;i++) {
    phase=generatePhaseWithL0(dimp*2,l0,nalias=1)(1:dim,1:dim,1);
    phaseb=phase-avg(phase); // substraction of the piston
    phase=phaseb(*)(w);
    Proj(,i)=MatProj(,+)*phase(+); // projection on the actuators
    write,format="Still %d samples \r",nsamp-i+1;
  }

 VectValid = Proj(*mcaodm(nm).valid_ptr,);
 VectExtrap = Proj(*mcaodm(nm).extrap_ptr,);

 mat1 = VectValid(,+)*VectExtrap(,+)/nsamp;
 mat2 = VectValid(,+)*VectValid(,+)/nsamp;
 mat2_inv = LUsolve(mat2);

 //-------------------------------------------------------
 //Step4 : Compute the Extrapolator

 E = mat1(+,)*mat2_inv(,+);
 return E;
 
}


func calc_FT_inf_fun(inf_fun,kdpix,if_nb,phase_spectrum)
{
  support = array(float,[2,kdpix, kdpix]);
  FT_inf_fun = array(complex,[3,kdpix, kdpix, if_nb]);
  for(i=1;i<=if_nb;i++){
    support(1:dim,1:dim) = inf_fun(,,i);
    FT_inf_fun(,,i) = fft(support,1)*sqrt(phase_spectrum);
  }
  write, "-> Calcul des TFs, DONE !";
  //inf_fun = support = [];
  return FT_inf_fun;
}


/* GENERATION OF ACTUATORS STATISTICAL COVARIANCE MATRIX USING
   PROJECTION ON RANDOMLY GENERATED PHASE SCREENS.
   
  To generate modes from double diagonalization, using the
  brute force approach (computing the statistics from generated
  phase screens)

  1. compute statistical covariance matrix of actuators w.r.t turbulence
  cv=build_statistics_covar_brute_force(1,3,5000, extrap="extrap_kl.mat")
  - > saves in covar_inffun"+rtd_get_timestamp()+".fits"

  2. compute double diagonalization:
  b = build_modes_2diago(1,3,cv,modes,eigen)
  write,eigen(1:8)

  3. Display modes:
  aff_kldm_basis,transpose(modes,[1,2,3]),50;

  4. feed these modes back into statistical covariance calculator
     to check the stat. covar matrix is indeed diagonal:
  cv2 = build_statistics_covar_brute_force(,,1000,modes)
  */

require,"turbulence.i";


func build_statistics_covar_brute_force(ndm, nzfiltered, nit, modes, extrap=,eigen=)
/* DOCUMENT
   ndm: dm#
   nzfiltered: # of filtered zernike (e.g.: 1 => filter piston)
   nit: # of iterations to compute statistics
   modes: modes cube, dimension [dim,dim,# of mode]
   extrap= extrapolator name (e.g. "extrap_kl.mat")
   eigen= just to plot this input vector w.r.t plotted variances.
   SEE ALSO:
   EXAMPLE:
   cv = build_statistics_covariance(1, 3, 10000, extrap="extrap_kl.mat")
   cv = build_statistics_covariance(,, 10000, modes)
   Francois Rigaut, 2009oct14
 */
{
  if (extrap==[]) extrap="extrap_kl.mat";
  if (!nit) nit=5000;
  if (nzfiltered==[]) nzfiltered=1;
  palette,"earth.gp";
  
  if (modes==[]) { // use actuators
    kl_basis_in_dm_space, ndm, nzfiltered, eigen, extrap=extrap;
    inf_fun = fitsRead(swrite(format="KLDMmodes_DM%d_IF_nrmodes%d.fits",\
                  ndm,nzfiltered));
    if_nb   = dimsof(inf_fun)(4);//nb of IF
    dim     = dimsof(inf_fun)(2);//size of the support
    w       = where(inf_fun(,,rms)>0.);
    p       = inf_fun(*,)(w,);
    pup     = inf_fun(,,1)*0.;
    imax    = 50; 
  } else { // compute stat covar matrix on provided modes.
    w       = where(modes(,,rms)>0);
    p       = modes(*,)(w,); 
    pup     = modes(,,1)*0.;
    imax    = 50; 
  }

  nmodes = dimsof(p)(3);

  pup(w) = 1.;
  // compute effective pupil to save cpu time later
  i1 = min(where2(pup)(1,));
  i2 = max(where2(pup)(1,));
  pupdiam = i2-i1+1;
  pupsmall = pup(i1:i2,i1:i2);
  w = where(pupsmall);

  // To filter TT from input phase (not really neccessary if inf_fun
  // have been orthog to it:
  //tmp  = indices(180)-90.5;
  //tip  = tmp(,,1)(*)(w);
  //tilt = tmp(,,2)(*)(w);
  //tip  = tip /sqrt(sum(tip*tip));
  //tilt = tilt/sqrt(sum(tilt*tilt));
  //pup = inf_fun(,,1)*0.;

  coef = array(1.,[2,nmodes,2*nit]);
  step = 50;
  
  for (i=1;i<=nit;i++) {
    write,format="\rIteration %d",i;

    psc = generatePhaseWithL0(512,220.,silent=1)(1:pupdiam,1:pupdiam,);
    pscl = psc(,,1)(*)(w);
    pscl -= avg(pscl);
    //pscl = pscl - sum(pscl*tip)*tip;
    //pscl = pscl - sum(pscl*tilt)*tilt;
    //pup(w) = pscl; tv,pup; pause,100;
    coef(,2*i-1) = p(+,)*pscl(+);
    // if (i==1) coef=coef(,-);
    pscl = psc(,,2)(*)(w);
    pscl -= avg(pscl);
    //pscl = pscl - sum(pscl*tip)*tip;
    //pscl = pscl - sum(pscl*tilt)*tilt;
    coef(,2*i) = p(+,)*pscl(+);
    //    grow,coef,p(,+)*pscl(+);
    if ((i%step)==0) {
      if (i>1000) step=500;
      window,0;
      tv,abs((covar=coef(,:2*i)(,+)*coef(,:2*i)(,+))(1:imax,1:imax));
      window,1;
      dc = diag(covar);
      dc = dc/dc(50);
      plot,dc; 
      if (eigen!=[]) plg,eigen/eigen(50),color="red";
      logxy,0,1;
      if (modes!=[]) logxy,1,1;
      window,0;
    }
  }
  write,"";
  covar=coef(,+)*coef(,+);
  covar = covar/max(covar);
  if (modes==[]) {
    name = swrite(format="covar_inffun_DM%d_nzf%d_%s.fits",ndm,nzfiltered,\
                  rtd_get_timestamp());
    fits_write,name,covar,overwrite=1;
    write,format="HIF saved in %s\n",name;
  }
  
  return covar;
}



func build_modes_2diago(ndm,nzfiltered,hif,&modes,&eigen)
/* DOCUMENT build_modes_2diago(ndm,nzfiltered,hif,&modes,&eigen)
   Build the B matrix (commands to get modes) and can return
   the modes. It is slightly different from Benoit's version
   as it filters out some modes entirely.
   example: 
   hif =  fits_read("covar_inffun.fits");
   b = build_modes_2diago(1,1,hif,modes,eigen)
   SEE ALSO:
   Francois Rigaut 2009oct14
 */
{
  inf_fun = fitsRead(swrite(format="KLDMmodes_DM%d_IF_nrmodes%d.fits",\
             ndm,nzfiltered));
  w = where(inf_fun(,,rms)>0);
  p = inf_fun(*,)(w,);
  p = transpose(p)*1.;
  // p is array(float,240,12824)
  Delta_IF = p(,+)*p(,+);

  hif = symetrize(hif,ndm);

  D2 = SVdec(Delta_IF, Mp);
  Mp = Mp(,:long(-nzfiltered));
  D2 = D2(1:long(-nzfiltered));
  print,"conditionning:",max(D2)/min(D2);
  M = transpose(Mp/sqrt(D2)(-,));
  Hp = M(,+)*(hif(,+)*M(,+))(+,);
  Lp = SVdec(Hp,A);
//  print,"conditionning:",max(Lp)/min(Lp);
  Bp = M(+,)*A(+,);
  eigen = Lp;
  write, "-> Double Diago, DONE !";

  // compute mode surface:
  modes = inf_fun(,,+) * Bp(+,);
  pup = inf_fun(,,1)*0; 
  pup(w) = 1.;
  modes = modes*pup(,,-);

  return Bp;
}

func symetrize(hif,ndm)
/* DOCUMENT symetrize(hif,ndm)
   This function takes a covariance matrix as input.
   It uses the symetry properties of the mcao dm config to
   average elements to improve symetry and SNR
   hif = statistical covariance matrix for a DM
   ndm: DM#
   SEE ALSO:
   Francois Rigaut, 2009oct14
 */

{
  x = *dm(ndm)._x;
  step = x(2)-x(1);
  // deal with dm4.5 flipped X:
  if (step<0) step*=-1;
  x = long((x-min(x))/step+1);
  y = *dm(ndm)._y;
  y = long((y-min(y))/step+1);
  nax = long(max(x));
  nay = long(max(y));
  nact = dm(ndm)._nact;
  locs = array(long,[2,nax,nay]);
  nv = indgen(nact);
  locs(x+(y-1)*nax) = nv;
  locs_orig = locs;
  
  alllocs = locs(,,-);
  // append symetries:
  // flip X
  grow,alllocs,locs(::-1,);
  
  // now rotate 90 degrees:
  locs = rotby90(locs,1)
  grow,alllocs,locs;
  // flip X
  grow,alllocs,locs(::-1,);
  
  // now rotate 90 degrees:
  locs = rotby90(locs,1)
  grow,alllocs,locs;
  // flip X
  grow,alllocs,locs(::-1,);
  
  // now rotate 90 degrees:
  locs = rotby90(locs,1)
  grow,alllocs,locs;
  // flip X
  grow,alllocs,locs(::-1,);
  
  // ok , we have all the symetries (8).
  // now let's use them (average pairs of points)
  
  w = where(locs_orig!=0);
  alllocs = alllocs(*,)(w,);
  hout = hout2 = hif*0.;
  plsys,1;
  animate,1;  
  
  for (i=1;i<=nact;i++) {
    write,format="\r Processing act%d",i;
    act1s = alllocs(i,);
    // actuator indices (orig)
    for (j=i;j<=nact;j++) {
      act2s = alllocs(j,);
      // actuator indices (targets)
      // extract all pairs
      set = array(0.,8);
      for (k=1;k<=8;k++) set(k) = hif(act1s(k),act2s(k));
      // average
      set = avg(set);
      // fill hout:
      for (k=1;k<=8;k++) hout(act1s(k),act2s(k)) = set;
      // check that stuff are not counted twice:
      // ok they are not
      // for (k=1;k<=8;k++) hout2(act1s(k),act2s(k)) = +set;
    }
    tv,hout;
  }
  animate,0;
  // amazingly, it worked from the first time !!! ???
  return hout;
}

func compare_tf_brute(void)
{
  // TF approach:
  b = kl_basis_in_dm_space(1, 1,eigen,extrap="extrap_kl.mat");
  inf_fun = fitsRead("KLDMmodes_DM1_IF_nrmodes1.fits");
  if_nb   = dimsof(inf_fun)(4);//nb of IF
  dim = dimsof(inf_fun)(2);//size of the support
  dm_mode = b(+,)*(inf_fun(*,))(,+);
  dm_mode = reform(dm_mode,if_nb,dim,dim);
  aff_kldm_basis,dm_mode,50;
  
  // stat aproach / brute force
  dm_mode=[];
  b = comp_modes(dm_mode);
  aff_kldm_basis,transpose(dm_mode,[1,2,3]),50; 
}




