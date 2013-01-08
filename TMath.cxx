#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include "SpecFuncCephes.h"
#include "TMath.h"

double TMath::Gamma(double a, double x)
{
   // Computation of the normalized lower incomplete gamma function P(a,x) as defined in the
   // Handbook of Mathematical Functions by Abramowitz and Stegun, formula 6.5.1 on page 260 .
   // Its normalization is such that TMath::Gamma(a,+infinity) = 1 .
   //
   //  Begin_Latex
   //  P(a, x) = #frac{1}{#Gamma(a) } #int_{0}^{x} t^{a-1} e^{-t} dt
   //   End_Latex
   //
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   //return ::ROOT::Math::inc_gamma(a, x);

/* This was the hardest part of making TRolke and FeldmanCousins stand-alone

TRolke calls TMath::ChisquareQuantile which calls TMath::Gamma,
which is the incomplete gamma function, which is not in the C or C++ std math lib.
ROOT eventually calls an implementation from Cephes in SpecFuncCephes.h/.cxx

If this doesn't work for some reason, here are some references

http://en.wikipedia.org/wiki/C_mathematical_functions
http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html
http://www.gnu.org/software/gsl/manual/html_node/Incomplete-Gamma-Functions.html
*/
   return ROOT::Math::Cephes::igam(a, x);
}

double TMath::ChisquareQuantile(double p, double ndf)
{
   // Evaluate the quantiles of the chi-squared probability distribution function.
   // Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
   // implemented by Anna Kreshuk.
   // Incorporates the suggested changes in AS R85 (vol.40(1), pp.233-5, 1991)
   // Parameters:
   //   p   - the probability value, at which the quantile is computed
   //   ndf - number of degrees of freedom

   double c[]={0, 0.01, 0.222222, 0.32, 0.4, 1.24, 2.2,
                 4.67, 6.66, 6.73, 13.32, 60.0, 70.0,
                 84.0, 105.0, 120.0, 127.0, 140.0, 175.0,
                 210.0, 252.0, 264.0, 294.0, 346.0, 420.0,
                 462.0, 606.0, 672.0, 707.0, 735.0, 889.0,
                 932.0, 966.0, 1141.0, 1182.0, 1278.0, 1740.0,
                 2520.0, 5040.0};
   double e = 5e-7;
   double aa = 0.6931471806;
   int maxit = 20;
   double ch, p1, p2, q, t, a, b, x;
   double s1, s2, s3, s4, s5, s6;

   if (ndf <= 0) return 0;

   double g = lgamma(0.5*ndf);

   double xx = 0.5 * ndf;
   double cp = xx - 1;
   if (ndf >= log(p)*(-c[5])){
   //starting approximation for ndf less than or equal to 0.32
      if (ndf > c[3]) {
         x = TMath::NormQuantile(p);
         //starting approximation using Wilson and Hilferty estimate
         p1=c[2]/ndf;
         ch = ndf*pow((x*sqrt(p1) + 1 - p1), 3);
         if (ch > c[6]*ndf + 6)
            ch = -2 * (log(1-p) - cp * log(0.5 * ch) + g);
      } else {
         ch = c[4];
         a = log(1-p);
         do{
            q = ch;
            p1 = 1 + ch * (c[7]+ch);
            p2 = ch * (c[9] + ch * (c[8] + ch));
            t = -0.5 + (c[7] + 2 * ch) / p1 - (c[9] + ch * (c[10] + 3 * ch)) / p2;
            ch = ch - (1 - exp(a + g + 0.5 * ch + cp * aa) *p2 / p1) / t;
         }while (std::abs(q/ch - 1) > c[1]);
      }
   } else {
      ch = pow((p * xx * exp(g + xx * aa)),(1./xx));
      if (ch < e) return ch;
   }
//call to algorithm AS 239 and calculation of seven term  Taylor series
   for (int i=0; i<maxit; i++){
      q = ch;
      p1 = 0.5 * ch;
      p2 = p - TMath::Gamma(xx, p1);

      t = p2 * exp(xx * aa + g + p1 - cp * log(ch));
      b = t / ch;
      a = 0.5 * t - b * cp;
      s1 = (c[19] + a * (c[17] + a * (c[14] + a * (c[13] + a * (c[12] +c[11] * a))))) / c[24];
      s2 = (c[24] + a * (c[29] + a * (c[32] + a * (c[33] + c[35] * a)))) / c[37];
      s3 = (c[19] + a * (c[25] + a * (c[28] + c[31] * a))) / c[37];
      s4 = (c[20] + a * (c[27] + c[34] * a) + cp * (c[22] + a * (c[30] + c[36] * a))) / c[38];
      s5 = (c[13] + c[21] * a + cp * (c[18] + c[26] * a)) / c[37];
      s6 = (c[15] + cp * (c[23] + c[16] * cp)) / c[38];
      ch = ch + t * (1 + 0.5 * t * s1 - b * cp * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
      if (std::abs(q / ch - 1) > e) break;
   }
   return ch;
}

double TMath::NormQuantile(double p)
{
   // Computes quantiles for standard normal distribution N(0, 1)
   // at probability p
   // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-484.

   if ((p<=0)||(p>=1)) {
      std::cerr << "TMath::NormQuantile: probability outside (0, 1)" << std::endl;
      return 0;
   }

   double  a0 = 3.3871328727963666080e0;
   double  a1 = 1.3314166789178437745e+2;
   double  a2 = 1.9715909503065514427e+3;
   double  a3 = 1.3731693765509461125e+4;
   double  a4 = 4.5921953931549871457e+4;
   double  a5 = 6.7265770927008700853e+4;
   double  a6 = 3.3430575583588128105e+4;
   double  a7 = 2.5090809287301226727e+3;
   double  b1 = 4.2313330701600911252e+1;
   double  b2 = 6.8718700749205790830e+2;
   double  b3 = 5.3941960214247511077e+3;
   double  b4 = 2.1213794301586595867e+4;
   double  b5 = 3.9307895800092710610e+4;
   double  b6 = 2.8729085735721942674e+4;
   double  b7 = 5.2264952788528545610e+3;
   double  c0 = 1.42343711074968357734e0;
   double  c1 = 4.63033784615654529590e0;
   double  c2 = 5.76949722146069140550e0;
   double  c3 = 3.64784832476320460504e0;
   double  c4 = 1.27045825245236838258e0;
   double  c5 = 2.41780725177450611770e-1;
   double  c6 = 2.27238449892691845833e-2;
   double  c7 = 7.74545014278341407640e-4;
   double  d1 = 2.05319162663775882187e0;
   double  d2 = 1.67638483018380384940e0;
   double  d3 = 6.89767334985100004550e-1;
   double  d4 = 1.48103976427480074590e-1;
   double  d5 = 1.51986665636164571966e-2;
   double  d6 = 5.47593808499534494600e-4;
   double  d7 = 1.05075007164441684324e-9;
   double  e0 = 6.65790464350110377720e0;
   double  e1 = 5.46378491116411436990e0;
   double  e2 = 1.78482653991729133580e0;
   double  e3 = 2.96560571828504891230e-1;
   double  e4 = 2.65321895265761230930e-2;
   double  e5 = 1.24266094738807843860e-3;
   double  e6 = 2.71155556874348757815e-5;
   double  e7 = 2.01033439929228813265e-7;
   double  f1 = 5.99832206555887937690e-1;
   double  f2 = 1.36929880922735805310e-1;
   double  f3 = 1.48753612908506148525e-2;
   double  f4 = 7.86869131145613259100e-4;
   double  f5 = 1.84631831751005468180e-5;
   double  f6 = 1.42151175831644588870e-7;
   double  f7 = 2.04426310338993978564e-15;

   double split1 = 0.425;
   double split2=5.;
   double konst1=0.180625;
   double konst2=1.6;

   double q, r, quantile;
   q=p-0.5;
   if (std::abs(q)<split1) {
      r=konst1-q*q;
      quantile = q* (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3)
                 * r + a2) * r + a1) * r + a0) /
                 (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3)
                 * r + b2) * r + b1) * r + 1.);
   } else {
      if(q<0) r=p;
      else    r=1-p;
      //error case
      if (r<=0)
         quantile=0;
      else {
         r=sqrt(-log(r));
         if (r<=split2) {
            r=r-konst2;
            quantile=(((((((c7 * r + c6) * r + c5) * r + c4) * r + c3)
                     * r + c2) * r + c1) * r + c0) /
                     (((((((d7 * r + d6) * r + d5) * r + d4) * r + d3)
                     * r + d2) * r + d1) * r + 1);
         } else{
            r=r-split2;
            quantile=(((((((e7 * r + e6) * r + e5) * r + e4) * r + e3)
                     * r + e2) * r + e1) * r + e0) /
                     (((((((f7 * r + f6) * r + f5) * r + f4) * r + f3)
                     * r + f2) * r + f1) * r + 1);
         }
         if (q<0) quantile=-quantile;
      }
   }
   return quantile;
}

double TMath::Poisson(double x, double par)
{
  // compute the Poisson distribution function for (x,par)
  // The Poisson PDF is implemented by means of Euler's Gamma-function
  // (for the factorial), so for all integer arguments it is correct.
  // BUT for non-integer values it IS NOT equal to the Poisson distribution.
  // see TMath::PoissonI to get a non-smooth function.
  // Note that for large values of par, it is better to call
  //     TMath::Gaus(x,par,sqrt(par),kTRUE)
//Begin_Html
/*
<img src="gif/Poisson.gif">
*/
//End_Html

   if (x<0)
      return 0;
   else if (x == 0.0)
      return 1./exp(par);
   else {
      double lnpoisson = x*log(par)-par-lgamma(x+1.);
      return exp(lnpoisson);
   }
   // An alternative strategy is to transition to a Gaussian approximation for
   // large values of par ...
   //   else {
   //     return Gaus(x,par,Sqrt(par),kTRUE);
   //   }
}

double TMath::PoissonI(double x, double par) {
	int x_int = x;
	return Poisson(x_int, par);
}

// TODO
bool TMath::RootsCubic(const double coef[4], double &a, double &b, double &c)
{
   // Calculates roots of polynomial of 3rd order a*x^3 + b*x^2 + c*x + d, where
   // a == coef[3], b == coef[2], c == coef[1], d == coef[0]
   //coef[3] must be different from 0
   // If the boolean returned by the method is false:
   //    ==> there are 3 real roots a,b,c
   // If the boolean returned by the method is true:
   //    ==> there is one real root a and 2 complex conjugates roots (b+i*c,b-i*c)
   // Author: Francois-Xavier Gentit

   bool complex = false;
   double r,s,t,p,q,d,ps3,ps33,qs2,u,v,tmp,lnu,lnv,su,sv,y1,y2,y3;
   a    = 0;
   b    = 0;
   c    = 0;
   if (coef[3] == 0) return complex;
   r    = coef[2]/coef[3];
   s    = coef[1]/coef[3];
   t    = coef[0]/coef[3];
   p    = s - (r*r)/3;
   ps3  = p/3;
   q    = (2*r*r*r)/27.0 - (r*s)/3 + t;
   qs2  = q/2;
   ps33 = ps3*ps3*ps3;
   d    = ps33 + qs2*qs2;
   if (d>=0) {
      complex = true;
      d   = sqrt(d);
      u   = -qs2 + d;
      v   = -qs2 - d;
      tmp = 1./3.;
      lnu = log(std::abs(u));
      lnv = log(std::abs(v));
      su  = copysign(1.,u);
      sv  = copysign(1.,v);
      u   = su*exp(tmp*lnu);
      v   = sv*exp(tmp*lnv);
      y1  = u + v;
      y2  = -y1/2;
      y3  = ((u-v)*sqrt(3.))/2;
      tmp = r/3;
      a   = y1 - tmp;
      b   = y2 - tmp;
      c   = y3;
   } else {
      double phi,cphi,phis3,c1,c2,c3,pis3;
      ps3   = -ps3;
      ps33  = -ps33;
      cphi  = -qs2/sqrt(ps33);
      phi   = acos(cphi);
      phis3 = phi/3;
      pis3  = M_PI/3;
      c1    = cos(phis3);
      c2    = cos(pis3 + phis3);
      c3    = cos(pis3 - phis3);
      tmp   = sqrt(ps3);
      y1    = 2*tmp*c1;
      y2    = -2*tmp*c2;
      y3    = -2*tmp*c3;
      tmp = r/3;
      a   = y1 - tmp;
      b   = y2 - tmp;
      c   = y3 - tmp;
   }
   return complex;
}

/* Sort index array using corresponding a values

   There's basically two ways to do this with C++ std::sort:
   http://stackoverflow.com/questions/3909272/sorting-two-corresponding-arrays
   Here we implement the make_pair solution
*/
void TMath::Sort(int n, const double* a, int* index)
{
   // make and init the pairs array
   std::vector<std::pair<double, int> > pairs(n);   
   for ( int i = 0; i < n; ++i )
     pairs[ i ] = std::make_pair( a[i], i );

   // sort
   std::sort( pairs.begin(), pairs.end() );

   for ( int i = 0; i < n; ++i ) {
      // set index array (n - i - 1 because we want descending!)
      // std::cout << i << " " << n - i - 1 << std::endl;
   	  index[i] = pairs[ n - i - 1 ].second;
      }
}
