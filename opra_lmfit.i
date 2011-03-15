/*
 * lmfit.i --
 *
 *      Non-linear least-squares fit by Levenberg-Marquardt method.
 *
 * Copyright (c) 1997, Eric THIEBAUT (thiebaut@obs.univ-lyon1.fr, Centre de
 * Recherche Astrophysique de Lyon, 9 avenue Charles  Andre,  F-69561 Saint
 * Genis Laval Cedex).
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
 * History:
 *      $Id: lmfit.i 8 2010-01-16 19:14:14Z frigaut $
 *      $Log: lmfit.i,v $
 *      Revision 1.1.1.1  2007/12/11 23:55:12  frigaut
 *      Initial Import - yorick-yutils
 *
 *      Revision 1.4  2003/06/17 12:17:10  eric
 *       - fix doc.
 *
 *      Revision 1.3  1998/09/08 15:31:06  eric
 *       - Make sure that input parameters have floating-point values.
 *
 *      Revision 1.2  1997/07/28 08:26:25  eric
 *       - Fix the doc.
 *
 *      Revision 1.1  1997/04/21 08:34:04  eric
 *      Initial revision
 *-----------------------------------------------------------------------------
 */

require, "random.i";

struct lmfit_result {
/* DOCUMENT lmfit_result -- structure returned by lmfit
 */
    long        neval;
    long        niter;
    long        nfit;
    long        nfree;
    long        monte_carlo;
    double      chi2_first;
    double      chi2_last;
    double      conv;
    double      sigma;
    double      lambda;
    pointer     stdev;
    pointer     stdev_monte_carlo;
    pointer     correl;
}

func opra_lmfit(f, x, &a, y, w, fit=, correl=, stdev=, gain=, tol=, \
                deriv=, itmax=, lambda=, eps=, monte_carlo=, aregul=,   \
                aregul_i=,neval_max=)
/* DOCUMENT lmfit -- Non-linear least-squares fit by Levenberg-Marquardt
                     method.

   DESCRIPTION:
     Implement Levenberg-Marquardt method  to  perform  a  non-linear least
     squares fit to a function of an arbitrary number of  parameters.   The
     function  may  be  any  non-linear  function.   If  available, partial
     derivatives can be calculated by the user function, else  this routine
     will  estimate  partial   derivatives   with   a   forward  difference
     approximation.

   CATEGORY:
     E2 - Curve and Surface Fitting.

   SYNTAX:
     result= lmfit(f, x, a, y, w, ...);

   INPUTS:
     F:  The model function  to  fit.   The  function  must  be  written as
         described under RESTRICTIONS, below.
     X:  Anything useful for the model function, for  instance: independent
         variables, a complex structure of  data  or  even  nothing!.   The
         LMFIT routine does not manipulate or use values  in  X,  it simply
         passes X to the user-written function F.
     A:  A vector that contains the initial estimate for each parameter.
     Y:  Array of dependent variables (i.e., the  data).   Y  can  have any
         geometry, but it must be the same as the result returned by F.
     W:  Optional weight,  must be conformable  with Y and all  values of W
         must be positive  or null (default = 1.0).   Data points with zero
         weight are not fitted.  Here are some examples:
           - For no weighting (lest square fit): W = 1.0
           - For instrumental weighting: W(i) = 1.0/Y(i)
           - Gaussian noise: W(i) = 1.0/Var(Y(i))

   OUTPUTS:
     A:  The vector of fitted parameters.
     Returns a structure lmfit_result with fields:
       NEVAL:  (long) number of model function evaluations.
       NITER:  (long) number of iteration, i.e. successful CHI2 reductions.
       NFIT:   (long) number of fitted parameters.
       NFREE:  (long) number of degrees of freedom  (i.e.,  number of valid
               data points minus number of fitted parameters).
       MONTE_CARLO: (long) number of Monte Carlo simulations.
       CHI2_FIRST: (double) starting error value: CHI2=sum(W*(F(X,A)-Y)^2).
       CHI2_LAST: (double) last best error value: CHI2=sum(W*(F(X,A)-Y)^2).
       CONV:   (double) relative variation of CHI2.
       SIGMA:  (double) estimated uniform standard deviation of data.  If a
               weight is provided, a  value  of  SIGMA  different  from one
               indicates  that,  if  the  model  is  correct,  W  should be
               multiplied    by     1/SIGMA^2.      Computed     so    that
               sum(W*(F(X,A)-Y)^2)/SIGMA^2=NFREE.
       LAMBDA: (double) last value of LAMBDA.
       STDEV:  (pointer) standard deviation vector of the parameters.
       STDEV_MONTE_CARLO: (pointer)  standard   deviation  vector  of  the
               parameters estimated by Monte Carlo simulations.
       CORREL: (pointer) correlation matrice of the parameters.

   KEYWORDS:
     FIT: List  of  indices  of  parameters  to  fit,  the  others  remaing
          constant.  The default is to tune all parameters.
     CORREL: If  set  to a  non  zero and  non-nil  value, the  correlation
          matrice of the parameters is stored into LMFIT result.
     STDEV: If set to a non zero and non-nil value, the standard deviation
          vector of the parameters is stored into LMFIT result.
     DERIV: When set to a non zero and non-nil  value,  indicates  that the
          model function F is able to compute its derivatives  with respect
          to  the  parameters   (see   RESTRICTIONS).    By   default,  the
          derivatives will be estimated by LMFIT using  forward difference.
          If analytical derivatives are  available  they  should  always be
          used.
     EPS: Small positive value  used  to  estimate  derivatives  by forward
          difference.  Must be such that 1.0+EPS  and  1.0  are numerically
          different  and   should   be   about  sqrt(machine_precision)/100
          (default = 1e-6).
     TOL: Stop criteria for the convergence (default =  1e-7).   Should not
          be smaller  than  sqrt(machine_precision).   The  routine returns
          when  the  relative  decrease of  CHI2 is  less  than  TOL  in an
          interation.
     ITMAX: Maximum number of iterations. Default = 100.
     GAIN: Gain factor for tuning LAMBDA (default = 10.0).
     LAMBDA: Starting value for parameter LAMBDA (default = 1.0e-3).
     MONTE_CARLO: Number of Monte Carlo  simulations to perform to estimate
          standard  deviation  of  parameters  (by default  no Monte  Carlo
          simulations are undergone).  May spend a lot of time if you use a
          large number; but should not be too small!
     AREGUL: specific to opra_lmfit. regularization factor for a(rms) (see
          code).

   GLOBAL VARIABLES:
     None.

   SIDE EFFECTS:
     The values of the vector of parameters A are modified.

   PROCEDURE:
     The function to be fitted must be defined as follow:

       func F(x, a) {....}

     and returns a model with same shape as data Y.  If you want to provide
     analytic derivatives, F should be defined as:

       func F(x, a, &grad, deriv=)
       {
           y= ...;
           if (deriv) {
               grad= ...;
           }
           return y;
       }

     Where X are the independent variables (anything the function  needs to
     compute synthetic data except the model parameters), A  are  the model
     parameters, DERIV is  a  flag  set  to  non-nil  and  non-zero  if the
     gradient is needed and the output gradient GRAD  is  a  numberof(Y) by
     numberof(A) array: GRAD(i,j) = derivative of ith data point model with
     respect to jth parameter.

     LMFIT tune  parameters A so as to minimize:  CHI2=sum(W*(F(X,A)-Y)^2).
     The  Levenberg-Marquardt  method   consists  in  varying  between  the
     inverse-Hessian  method  and the  steepest  descent  method  where the
     quadratic  expansion of  CHI2  does  not  yield  a better  model.  The
     initial  guess of  the parameter  values  should  be as  close to  the
     actual values as possible or the solution may not converge or may give
     a wrong answer.

   RESTRICTIONS:
     Beware that  the result  does depend on your  initial guess A.  In the
     case of  numerous  local  minima,  the only  way to  get  the  correct
     solution is to start with A close enough to this solution.

     The estimates of  standard  deviation of the  parameters are  rescaled
     assuming that, for a correct model  and weights, the expected value of
     CHI2 should  be of the  order of NFREE=numberof(Y)-numberof(A)  (LMFIT
     actually compute NFREE from the number of valid data points and number
     of fitted parameters).   If you don't like this you'll have to rescale
     the returned  standard  deviation  to meet  your needs  (all necessary
     information are in the structure returned by LMFIT).

   EXAMPLE:
     This example is from ODRPACK (version 2.01).  The function to fit is
     of the form:
       f(x) = a1+a2*(exp(a3*x)-1.0)^2
     Starting guess:
       a= [1500.0,  -50.0,   -0.1];
     Independent variables:
       x= [   0.0,    0.0,    5.0,    7.0,    7.5,   10.0,
             16.0,   26.0,   30.0,   34.0,   34.5,  100.0];
     Data:
       y= [1265.0, 1263.6, 1258.0, 1254.0, 1253.0, 1249.8,
           1237.0, 1218.0, 1220.6, 1213.8, 1215.5, 1212.0];
     Function definition (without any optimization):
       func foo(x, a, &grad, deriv=)
       {
           if (deriv)
               grad= [array(1.0, dimsof(x)),
                      (exp(a(3)*x)-1.0)^2,
                      2.0*a(2)*x*exp(a(3)*x)*(exp(a(3)*x)-1.0)];
           return a(1)+a(2)*(exp(a(3)*x)-1.0)^2;
       }

     Fitting this model by:
       r= lmfit(foo, x, a, y, 1., deriv=1, stdev=1, monte_carlo=500, correl=1)
     produces typically the following result:
        a                   = [1264.84, -54.9987, -0.0829835]
        r.neval             = 12
        r.niter             = 6
        r.nfit              = 3
        r.nfree             = 9
        r.monte_carlo       = 500
        r.chi2_first        = 40.4383
        r.chi2_last         = 40.4383
        r.conv              = 3.84967e-09
        r.sigma             = 0.471764
        r.lambda            = 1e-09
       *r.stdev             = [1.23727, 1.78309, 0.00575123]
       *r.stdev_monte_carlo = [1.20222, 1.76120, 0.00494790]
       *r.correl            = [[ 1.000, -0.418, -0.574],
                               [-0.418,  1.000, -0.340],
                               [-0.574, -0.340,  1.000]]

   HISTORY:
     - Basic ideas borrowed from "Numerical Recipes in C", CURVEFIT.PRO (an
       IDL version by DMS, RSI, of the routine "CURFIT: least squares fit to
       a non-linear function", Bevington, Data Reduction and Error Analysis
       for the Physical Sciences) and ODRPACK ("Software for Weigthed
       Orthogonal Distance Regression" freely available at: www.netlib.org).
     - Added: fitting of a subset of the parameters, Monte-Carlo
       simulations...
 */
{
  local grad;
  extern newiter, lmfit_itmax, lmfititer_pass;
  extern stop_all;

  /* Maybe subset of parameters to fit. */
  if (structof(a)!=double) {
    a+= 0.0;
    if (structof(a)!=double)
        error, "bad data type for parameters (complex unsupported)";
  }
  na= numberof(a);
  if (is_void(fit))
    fit= indgen(na);
  else if (dimsof(fit)(1) == 0)
    fit= [fit];
  nfit= numberof(fit);
  if (!nfit)
    error, "no parameters to fit.";

  if (itmax) lmfit_itmax=itmax;
  lmfititer_pass = 0;
  if (aregul==[]) aregul=0.;

  /* Check weights. */
  if (is_void(w)) w= 1.0;
  else if (anyof(w < 0.0))
    error, "bad weights.";
  if (numberof(w) != numberof(y))
    w += array(0.0, dimsof(y));
  nfree= sum(w != 0.0) - nfit;        // Degrees of freedom
  if (nfree <= 0)
    error, "not enough data points.";

  /* Other settings. */
  diag= indgen(1:nfit^2:nfit+1);      // Subscripts of diagonal elements
  if (is_void(lambda)) lambda= 1e-3;
  if (is_void(gain)) gain= 10.0;
  if (is_void(itmax)) itmax= 100;
  if (is_void(eps)) eps= 1e-6;        // sqrt(machine_precision)/100
  if (1.0+eps <= 1.0)
      error, "bad value for EPS.";
  if (is_void(tol)) tol= 1e-7;
  monte_carlo= is_void(monte_carlo) ? 0 : long(monte_carlo);
  warn_zero= 0;
  warn= "*** Warning: LMFIT ";
  neval= 0;
  conv= 0.0;
  niter= 0;

  while (1) {
    if ((has_svipc)&&(shm_read(shmkey,"next_stage")(1))) {
      shm_write,shmkey,"next_stage",&([0]);
      goto done;
    }
    if ((has_svipc)&&(shm_read(shmkey,"stop")(1))) { stop_all=1; goto done; }

    newiter=1;
    lmfititer_pass++;
    if (deriv) {
      m= f(x, a, grad, deriv=1);
      neval++;
      grad= nfit == na ? grad(*,) : grad(*,fit);
    } else {
      if (!niter) {
          m= f(x, a);
          neval++;
      }
      inc= eps * abs(a(fit));
      if (numberof((i= where(inc <= 0.0)))) inc(i)= eps;
      grad= array(double, numberof(y), nfit);
      for (i=1; i<=nfit; i++) {
        if ((has_svipc)&&(shm_read(shmkey,"quit?")(1))) { quit; }
        anew= a;        // Copy current parameters
        anew(fit(i)) += inc(i);
        grad(,i)= (f(x,anew)-m)(*)/inc(i);
      }
      neval += nfit;
    }
    beta= w * (chi2= y-m);
    if (niter) chi2= chi2new;
    else {
      chi2= chi2_first= sum(beta * chi2);
      if (aregul) {
        chi2       *= (1.+aregul*anew(aregul_i)(rms));
        chi2_first *= (1.+aregul*anew(aregul_i)(rms));
      }
    }
    beta= grad(+,) * beta(*)(+);
    alpha= ((w(*)(,-) * grad)(+,) * grad(+,));
    gamma= sqrt(alpha(diag));
    if (anyof(gamma <= 0.0)) {
      /* Some derivatives are null (certainly because of rounding
       * errors). */
      if (!warn_zero) {
          write, warn+"founds zero derivatives.";
          warn_zero= 1;
      }
      gamma(where(gamma <= 0.0))= eps * max(gamma);
      /* goto done; */
    }
    gamma= 1.0 / gamma;
    beta *= gamma;
    alpha *= gamma(,-) * gamma(-,);

    neval_intern = 0;
    while (1) {
      if ((has_svipc)&&(shm_read(shmkey,"stop")(1))) { stop_all=1; goto done; }
      alpha(diag)= 1.0 + lambda;
      anew= a;
      anew(fit) += gamma * LUsolve(alpha, beta);
      m= f(x, anew);
      neval++;
      neval_intern++;
      d= y-m;
      chi2new= sum(w*d*d);
      if (aregul) chi2new *= (1.+aregul*anew(aregul_i)(rms));
      if (chi2new < chi2)
        break;
      lambda *= gain;
      if (allof(anew == a) || neval_intern > neval_max) {
        /* No change in parameters. */
        write, warn+"makes no progress.";
        goto done;
      }
    }
    a= anew;
    lambda /= gain;
    niter++;
    conv= 2.0*(chi2-chi2new)/(chi2+chi2new);
    if (conv <= tol)
        break;
    if (niter >= itmax) {
        write, format=warn+"reached maximum number of iterations (%d).\n",
            itmax;
        break;
    }
  }

done:
    sigma= sqrt(nfree/chi2);
    result= lmfit_result(neval=neval, niter=niter, nfree=nfree, nfit=nfit,
            lambda=lambda, chi2_first=chi2_first, chi2_last=chi2, conv=conv,
            sigma=sigma);
    if (correl || stdev) {
        /* Compute correlation matrice and/or standard deviation vector. */
        alpha(diag)= 1.0;
        alpha= LUsolve(alpha);
        if (anyof((tmp1= alpha(diag)) < 0.0))
            write, format=warn+"%s\n", "found negative variance(s)";
        tmp1= sqrt(abs(tmp1));
        if (stdev) {
            /* Standard deviation is rescaled assuming that statistically
             * chi2 = nfree +/- sqrt(2*nfree). */
            (tmp2= array(double,na))(fit)= gamma * tmp1 / sigma;
            result.stdev= &tmp2;
        }
        if (correl) {
            gamma= 1.0 / tmp1;
            alpha *= gamma(-,) * gamma(,-);
            if (nfit == na) {
                result.correl= &alpha;
            } else {
                (tmp2= array(double, na, na))(fit,fit)= alpha;
                result.correl= &tmp2;
            }
        }
    }
    alpha= beta= gamma= [];     // Free some memory.
    if (monte_carlo >= 1) {
        saa= 0.0*a;
        sig= (w > 0.0) /(sqrt(max(nfree/chi2*w, 0.0)) + (w == 0.0));
        for (i=1; i<=monte_carlo; i++) {
            anew= a;
            ynew= y + sig * random_n(dimsof(y));
            lmfit, f, x, anew, ynew, w, fit=fit, gain=gain, tol=tol,
                deriv=deriv, itmax=itmax, lambda=lambda, eps=eps;
            anew -= a;
            saa += anew * anew;
        }
        result.monte_carlo= monte_carlo;
        result.stdev_monte_carlo= &sqrt(saa / monte_carlo);
    }
    return result;
}

