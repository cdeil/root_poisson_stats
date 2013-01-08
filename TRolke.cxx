
//////////////////////////////////////////////////////////////////////////////
//
//  TRolke 
//
//  This class computes confidence intervals for the rate of a Poisson
//  process in the presence of uncertain background and/or efficiency. 
//
//  The treatment and the resulting limits are fully frequentist. The
//  limit calculations make use of the profile likelihood method.
//
//      Author: Jan Conrad (CERN) 2004
//      Updated: Johan Lundberg (CERN) 2009
//
//      Copyright CERN 2004,2009           Jan.Conrad@cern.ch,
//                                     Johan.Lundberg@cern.ch
//
//  For a full list of methods and their syntax, and build instructions,
//  consult the header file TRolke.h.
//  -------------------------------- 
//
//  Examples/tutorials are found in the separate file Rolke.C 
//  ---------------------------------------------------------
//
//////////////////////////////////////////////////////////////////////////////
//
//
//  TRolke implements the following Models
//  =======================================
//
//  The signal is always assumed to be Poisson, with the following
//  combinations of models of background and detection efficiency:
//
//  If unsure, first consider model 3, 4 or 5.
//
//       1: SetPoissonBkgBinomEff(x,y,z,tau,m)
//
//          Background: Poisson
//          Efficiency: Binomial    
//
//          when the background is simultaneously measured 
//          from sidebands (or MC), and
//          the signal efficiency was determined from Monte Carlo 
//
//       2: SetPoissonBkgGaussEff(x,y,em,sde,tau)
//
//          Background: Poisson
//          Efficiency: Gaussian    
//
//          when the background is simultaneously measured 
//          from sidebands (or MC), and
//          the efficiency is modeled as Gaussian
//
//       3: SetGaussBkgGaussEff(x,bm,em,sde,sdb)
//
//          Background: Gaussian
//          Efficiency: Gaussian    
//
//          when background and efficiency can both be
//          modeled as Gaussian.
//
//       4: SetPoissonBkgKnownEff(x,y,tau,e)
//
//          Background: Poisson
//          Efficiency: Known     
//
//          when the background is simultaneously measured 
//          from sidebands (or MC).
//
//       5: SetGaussBkgKnownEff(x,bm,sdb,e)
//
//          Background: Gaussian
//          Efficiency: Known       
//
//          when background is Gaussian
//
//       6: SetKnownBkgBinomEff(x,z,b,m)
//
//          Background: Known    
//          Efficiency: Binomial    
//
//          when signal efficiency was determined from Monte Carlo 
//
//       7: SetKnownBkgGaussEff(x,em,sde,b)
//
//          Background: Known    
//          Efficiency: Gaussian     
//
//          when background is known and efficiency Gaussian
//
//  Parameters and further explanation
//  ==================================
//
//  For all models:
//  ---------------
//    
//    x = number of observed events in the experiment
//    
//    Efficiency (e or em) is the detection probability for signal. 
//    A low efficiency hence generally means weaker limits.
//    If the efficiency of an experiment (with analysis cuts) is 
//    dealt with elsewhere, em or e can be set to one.
//
//  For Poisson background measurements (sideband or MC):
//  -----------------------------------------------------
//
//    y = number of observed events in background region
//    tau = 
//        Either: the ratio between signal and background region
//        in case background is observed. 
//        Or: the ratio between observed and simulated live-time
//        in case background is determined from MC.
//
//  For Gaussian efficiency or background:
//  --------------------------------------
//
//    bm  = estimate of the background
//    sdb = corresponding standard deviation 
//
//    em  = estimate of the efficiency
//    sde = corresponding standard deviation 
//
//        If the efficiency scale of dealt with elsewhere,
//        set em to 1 and sde to the relative uncertainty.
//
//  For Binomial signal efficiency:
//  -------------------------------
//
//     m = number of MC events generated
//     z = number of MC events observed
//
//  For the case of known background expectation or known efficiency:
//  -----------------------------------------------------------------
//
//     e = true efficiency (considered known)
//     b = background expectation value (considered known)
//
//
////////////////////////////////////////////////////////////////////  
//
//
//  The confidence level (CL) is set either at construction
//  time or with either of SetCL or SetCLSigmas
//
//  The TRolke method is very similar to the one used in MINUIT (MINOS).
//
//  Two options are offered to deal with cases where the maximum likelihood
//  estimate (MLE) is not in the physical region. Version "bounded likelihood"
//  is the one used by MINOS if bounds for the physical region are chosen. 
//  Unbounded likelihood (the default) allows the MLE to be in the
//  unphysical region. It has however better coverage.
//  For more details consult the reference (see below).
//
//  For a description of the method and its properties:
//
//  W.Rolke, A. Lopez, J. Conrad and Fred James
//  "Limits and Confidence Intervals in presence of nuisance parameters"
//   http://lanl.arxiv.org/abs/physics/0403059
//   Nucl.Instrum.Meth.A551:493-503,2005
//
//  Should I use TRolke, TFeldmanCousins, TLimit?
//  ----------------------------------------------
//
//  1. Does TRolke make TFeldmanCousins obsolete?
//
//  Certainly not. TFeldmanCousins is the fully frequentist construction and
//  should be used in case of no (or negligible) uncertainties. It is however
//  not capable of treating uncertainties in nuisance parameters. In other
//  words, it does not handle background expectations or signal efficiencies
//  which are known only with some limited accuracy.
//
//  TRolke is designed for this case and it is shown in the reference above
//  that it has good coverage properties for most cases, and can be used
//  where FeldmannCousins can't.
//
//  2. What are the advantages of TRolke over TLimit?
//
//  TRolke is fully frequentist. TLimit treats nuisance parameters Bayesian.
//
//  For a coverage study of a Bayesian method refer to
//  physics/0408039 (Tegenfeldt & J.C). However, this note studies
//  the coverage of Feldman&Cousins with Bayesian treatment of nuisance
//  parameters. To make a long story short: using the Bayesian method you
//  might introduce a small amount of over-coverage (though I haven't shown it
//  for TLimit). On the other hand, coverage of course is a not so interesting
//  when you consider yourself a Bayesian.
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "TMath.h"
#include "TRolke.h"

//__________________________________________________________________________
TRolke::TRolke(double CL) 
:  fCL(CL),
   fUpperLimit(0.0), 
   fLowerLimit(0.0), 
   fBounding(false),  // true gives bounded likelihood   
   fNumWarningsDeprecated1(0),
   fNumWarningsDeprecated2(0)
{
//constructor with optional Confidence Level argument.

   SetModelParameters();
}

//___________________________________________________________________________
TRolke::~TRolke()
{
}

/* Set the Confidence Level in terms of Sigmas. */
void TRolke::SetCLSigmas(double CLsigmas) {
      fCL = erf(CLsigmas / sqrt(2.0)) ;
}


/* ______________________ new interface _____________________ */
void TRolke::SetPoissonBkgBinomEff(int x, int y, int z, double tau, int m)
{
// Model 1: Background - Poisson, Efficiency - Binomial
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    z   : number of MC events observed
//    tau : ratio parameter (read TRolke.cxx for details)
//    m   : number of MC events generated
   SetModelParameters(
	 x  ,       //   int x,
         y  ,       //   int y,
         z  ,       //   int z,
         0  ,       //   double bm,
         0  ,       //   double em,
         0  ,       //   double e,
         1  ,       //   int mid,
         0  ,       //   double sde,
         0  ,       //   double sdb,
         tau,       //   double tau,
         0  ,       //   double b,
         m);        //   int m          
}

void TRolke::SetPoissonBkgGaussEff(int x, int y, double em, double tau, double sde)
{
// Model 2: Background - Poisson, Efficiency - Gaussian
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    em  : estimate of the efficiency
//    tau : ratio parameter (read TRolke.cxx for details)
//    sde : efficiency estimate's standard deviation 
   SetModelParameters(
         x  ,       //   int x,
         y  ,       //   int y,
         0  ,       //   int z,
         0  ,       //   double bm,
         em ,       //   double em,
         0  ,       //   double e,
         2  ,       //   int mid,
         sde,       //   double sde,
         0  ,       //   double sdb,
         tau,       //   double tau,
         0  ,       //   double b,
         0);        //   int m

}

void TRolke::SetGaussBkgGaussEff(int x, double bm, double em, double sde, double sdb)
{
// Model 3: Background - Gaussian, Efficiency - Gaussian (x,bm,em,sde,sdb)
//    x   : number of observed events in the experiment
//    bm  : estimate of the background
//    em  : estimate of the efficiency
//    sde : efficiency estimate's standard deviation 
//    sdb : background estimate's standard deviation 
   SetModelParameters(
         x  ,       //   int x,
         0  ,       //   int y,
         0  ,       //   int z,
         bm ,       //   double bm,
         em ,       //   double em,
         0  ,       //   double e,
         3  ,       //   int mid,
         sde,       //   double sde,
         sdb,       //   double sdb,
         0  ,       //   double tau,
         0  ,       //   double b,
         0);        //   int m

}

void TRolke::SetPoissonBkgKnownEff(int x, int y, double tau, double e)
{
// Model 4: Background - Poisson, Efficiency - known     (x,y,tau,e)
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    tau : ratio parameter (read TRolke.cxx for details)
//    e   : true efficiency (considered known)
   SetModelParameters(
         x  ,       //   int x,
         y  ,       //   int y,
         0  ,       //   int z,
         0  ,       //   double bm,
         0  ,       //   double em,
         e  ,       //   double e,
         4  ,       //   int mid,
         0  ,       //   double sde,
         0  ,       //   double sdb,
         tau,       //   double tau,
         0  ,       //   double b,
         0);        //   int m

}

void TRolke::SetGaussBkgKnownEff(int x, double bm, double sdb, double e)
{
// Model 5: Background - Gaussian, Efficiency - known    (x,bm,sdb,e
//    x   : number of observed events in the experiment
//    bm  : estimate of the background
//    sdb : background estimate's standard deviation 
//    e   : true efficiency (considered known)
   SetModelParameters(
         x  ,       //   int x,
         0  ,       //   int y,
         0  ,       //   int z,
         bm ,       //   double bm,
         0  ,       //   double em,
         e  ,       //   double e,
         5  ,       //   int mid,
         0  ,       //   double sde,
         sdb,       //   double sdb,
         0  ,       //   double tau,
         0  ,       //   double b,
         0);        //   int m

}

void TRolke::SetKnownBkgBinomEff(int x, int z, int m, double b)
{
// Model 6: Background - known, Efficiency - Binomial    (x,z,m,b)
//    x   : number of observed events in the experiment
//    z   : number of MC events observed
//    m   : number of MC events generated
//    b   : background expectation value (considered known)
   SetModelParameters(
         x  ,       //   int x,
         0  ,       //   int y
         z  ,       //   int z,
         0  ,       //   double bm,
         0  ,       //   double em,
         0  ,       //   double e,
         6  ,       //   int mid,
         0  ,       //   double sde,
         0  ,       //   double sdb,
         0  ,       //   double tau,
         b  ,       //   double b,
         m);        //   int m

}

void TRolke::SetKnownBkgGaussEff(int x, double em, double sde, double b)
{
// Model 7: Background - known, Efficiency - Gaussian    (x,em,sde,b)
//    x   : number of observed events in the experiment
//    em  : estimate of the efficiency
//    sde : efficiency estimate's standard deviation 
//    b   : background expectation value (considered known)
   SetModelParameters(
         x  ,       //   int x,
         0  ,       //   int y
         0  ,       //   int z,
         0  ,       //   double bm,
         em ,       //   double em,
         0  ,       //   double e,
         7  ,       //   int mid,
         sde,       //   double sde,
         0  ,       //   double sdb,
         0  ,       //   double tau,
         b  ,       //   double b,
         0);        //   int m

}

bool TRolke::GetLimits(double& low, double& high)
{
/* Calculate and get the upper and lower limits for the pre-specified model */ 
   if ((f_mid<1)||(f_mid>7)) {
      std::cerr << "TRolke - Error: Model id "<< f_mid<<std::endl;
      if (f_mid<1) {
 	 std::cerr << "TRolke - Please specify a model with e.g. 'SetGaussBkgGaussEff' (read the docs in Rolke.cxx )"<<std::endl; 
      }
      return false;
   }

   ComputeInterval(f_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
   low = fLowerLimit;
   high = fUpperLimit;
   if (low < high) {     
      return true;
   }else{ 
      std::cerr << "TRolke - Warning: no limits found" <<std::endl;
      return false;
   }
}

double TRolke::GetUpperLimit()
{
/* Calculate and get upper limit for the pre-specified model */ 
   double low(0), high(0);
   GetLimits(low,high);
   return fUpperLimit;
}

double TRolke::GetLowerLimit()
{
/* Calculate and get lower limit for the pre-specified model */ 
   double low(0), high(0);
   GetLimits(low,high);
   return fLowerLimit;
}

double TRolke::GetBackground()
{
/* Return a simple background value (estimate/truth) given the pre-specified model */ 
   double background = 0;
   switch (f_mid) {
      case 1:
      case 2:
      case 4:
         background = f_y / f_tau;
         break;
      case 3:
      case 5:
         background = f_bm;
         break;
      case 6:
      case 7:
	 background = f_b;
         break;
      default:
	 std::cerr << "TRolke::GetBackground(): Model NR: " <<
	    f_mid << " unknown"<<std::endl;
	 return 0;
   }
   return background;
}

bool TRolke::GetSensitivity(double& low, double& high, double pPrecision)
{
// get the upper and lower average limits based on the specified model.
// No uncertainties are considered for the Poisson weights in the averaging sum
   double background = GetBackground();

   double weight = 0;
   double weightSum = 0;

   int loop_x = 0;

   while (1) {
      ComputeInterval(loop_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
      weight = TMath::PoissonI(loop_x, background);
      low += fLowerLimit * weight;
      high += fUpperLimit * weight;
      weightSum += weight;
      if (loop_x > (background + 1)) { // don't stop too early
         if ((weightSum > (1 - pPrecision)) || (weight < 1e-12)) break;
      }
      loop_x++;
   }
   low /= weightSum;
   high /= weightSum;

   return (low < high); // could also add more detailed test
}


bool TRolke::GetLimitsQuantile(double& low, double& high, int& out_x, double integral)
{
/* get the upper and lower limits for the outcome corresponding to
   a given quantile.
   For integral=0.5 this gives the median limits
   in repeated experiments. The returned out_x is the corresponding
   (e.g. median) value of x.
   No uncertainties are considered for the Poisson weights when calculating 
   the Poisson integral */
   double background = GetBackground();
   double weight = 0;
   double weightSum = 0;
   int loop_x = 0;

   while (1) {
      weight = TMath::PoissonI(loop_x, background);
      weightSum += weight;
      if (weightSum >= integral) {
         break;
      }
      loop_x++;
   }

   out_x = loop_x;

   ComputeInterval(loop_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
   low = fLowerLimit;
   high = fUpperLimit;
   return (low < high); // could also add more detailed test

}

bool TRolke::GetLimitsML(double& low, double& high, int& out_x)
{
/* get the upper and lower limits for the most likely outcome.
   The returned out_x is the corresponding value of x 
   No uncertainties are considered for the Poisson weights when finding ML */
   double background = GetBackground();

   int loop_x = 0; // this can be optimized if needed. 
   int loop_max = 1000 + (int)background; //     --||--

   double max = TMath::PoissonI(loop_x, background);
   while (loop_x <= loop_max) {
      if (TMath::PoissonI(loop_x + 1, background) < max) {
         break;
      }
      loop_x++;
      max = TMath::PoissonI(loop_x, background);
   }
   if (loop_x >= loop_max) {
      std::cout << "internal error finding maximum of distribution" << std::endl;
      return false;
   }

   out_x = loop_x;

   ComputeInterval(loop_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
   low = fLowerLimit;
   high = fUpperLimit;
   return (low < high); // could also add more detailed test

}

bool TRolke::GetCriticalNumber(int& ncrit, int maxtry)
{
// get the value of x corresponding to rejection of the null hypothesis.
// This means a lower limit >0 with the pre-specified Confidence Level. 
// Optionally give maxtry; the maximum value of x to try. Of not, or if
// maxtry<0 an automatic mode is used.
   double background = GetBackground();

   int j = 0;
   int rolke_ncrit = -1;
   int maxj =maxtry ;
   if(maxtry<1){
     maxj = 1000 + (int)background; // max value to try
   }
   for (j = 0;j < maxj;j++) {
      int rolke_x = j;
      ComputeInterval(rolke_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
      double rolke_ll = fLowerLimit;
      if (rolke_ll > 0) {
         rolke_ncrit = j;
         break;
      }
   }

   if (rolke_ncrit == -1) {
     std::cerr << "TRolke GetCriticalNumber : Error: problem finding rolke inverse. Specify a larger maxtry value. maxtry was: " << maxj << ". highest x considered was j "<< j<< std::endl;
      ncrit = -1;
      return false;
   } else {
      ncrit = rolke_ncrit;
      return true;
   }
}

void TRolke::SetSwitch(bool bnd) {
/* Deprecated name for SetBounding. */
   if(fNumWarningsDeprecated1<2){
      std::cerr << "*******************************************" <<std::endl;
      std::cerr << "TRolke - Warning: 'SetSwitch' is depricated and may be removed from future releases:" <<std::endl;
      std::cerr << " - Use 'SetBounding' instead "<<std::endl;
      std::cerr << "*******************************************" <<std::endl;
      fNumWarningsDeprecated1++;
   }
   SetBounding(bnd);
}

void TRolke::Print(void) const {
/* Dump internals. Print members. */
   std::cout << "*******************************************" <<std::endl;
   std::cout << "* TRolke::Print() - dump of internals:                " <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "* model id, mid = "<<f_mid <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "*             x = "<<f_x   <<std::endl;
   std::cout << "*            bm = "<<f_bm  <<std::endl;
   std::cout << "*            em = "<<f_em  <<std::endl;
   std::cout << "*           sde = "<<f_sde <<std::endl;
   std::cout << "*           sdb = "<<f_sdb <<std::endl;
   std::cout << "*             y = "<<f_y   <<std::endl;      
   std::cout << "*           tau = "<<f_tau <<std::endl;
   std::cout << "*             e = "<<f_e   <<std::endl;
   std::cout << "*             b = "<<f_b   <<std::endl;
   std::cout << "*             m = "<<f_m   <<std::endl;
   std::cout << "*             z = "<<f_z   <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "*            CL = "<<fCL <<std::endl;
   std::cout << "*      Bounding = "<<fBounding <<std::endl;      
   std::cout << "*"<<std::endl;
   std::cout << "* calculated on demand only:"<<std::endl;
   std::cout << "*   fUpperLimit = "<<fUpperLimit<<std::endl; 
   std::cout << "*   fLowerLimit = "<<fLowerLimit<<std::endl;      
   std::cout << "*******************************************" <<std::endl;
}


double TRolke::CalculateInterval(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m){
/* Deprecated and error prone model selection interface. 
   It's use is trongly discouraged. 'mid' is the model ID (1 to 7).
   This method is provided for backwards compatibility/developer use only. */
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    z   : number of MC events observed
//    bm  : estimate of the background
//    em  : estimate of the efficiency
//    e   : true efficiency (considered known)
//    mid : internal model id (really, you should not use this method at all)
//    sde : efficiency estimate's standard deviation 
//    sdb : background estimate's standard deviation 
//    tau : ratio parameter (read TRolke.cxx for details)
//    b   : background expectation value (considered known)
//    m   : number of MC events generated
   if (fNumWarningsDeprecated2<2 ) {
      std::cerr << "*******************************************" <<std::endl;
      std::cerr << "TRolke - Warning: 'CalculateInterval' is depricated and may be removed from future releases:" <<std::endl;
      std::cerr << " - Use e.g. 'SetGaussBkgGaussEff' and 'GetLimits' instead (read the docs in Rolke.cxx )"<<std::endl;
      std::cerr << "*******************************************" <<std::endl;
      fNumWarningsDeprecated2++;
   }
   SetModelParameters(
	 x,
         y,
         z,
         bm,
         em,
         e,
         mid,
         sde,
         sdb,
         tau,
         b,
         m);
   return ComputeInterval(f_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
}


// --------------- Private methods ----------------
void TRolke::SetModelParameters(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m)
{
//___________________________________________________________________________
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    z   : number of MC events observed
//    bm  : estimate of the background
//    em  : estimate of the efficiency
//    e   : true efficiency (considered known)
//    mid : internal model id
//    sde : efficiency estimate's standard deviation 
//    sdb : background estimate's standard deviation 
//    tau : ratio parameter (read TRolke.cxx for details)
//    b   : background expectation value (considered known)
//    m   : number of MC events generated
   f_x   = x   ;
   f_y   = y   ;
   f_z   = z   ;
   f_bm  = bm  ;
   f_em  = em  ;
   f_e   = e   ;
   f_mid = mid ;
   f_sde = sde ;
   f_sdb = sdb ;
   f_tau = tau ;
   f_b   = b   ;
   f_m   = m   ;
}

void TRolke::SetModelParameters()
{
/* Clear internal model */
   SetModelParameters(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
   f_mid=0;
}


double TRolke::ComputeInterval(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m)
{
//___________________________________________________________________________
// ComputeInterval, the internals.
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    z   : number of MC events observed
//    bm  : estimate of the background
//    em  : estimate of the efficiency
//    e   : true efficiency (considered known)
//    mid : internal model id (really, you should not use this method at all)
//    sde : efficiency estimate's standard deviation 
//    sdb : background estimate's standard deviation 
//    tau : ratio parameter (read TRolke.cxx for details)
//    b   : background expectation value (considered known)
//    m   : number of MC events generated

   //calculate interval
   int done = 0;
   double limit[2];

   limit[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);

   if (limit[1] > 0) {
      done = 1;
   }

   if (! fBounding) {

      int trial_x = x;

      while (done == 0) {
         trial_x++;
         limit[1] = Interval(trial_x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
         if (limit[1] > 0) done = 1;
      }
   }

   return limit[1];
}

double TRolke::Interval(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m)
{
//_____________________________________________________________________
// Internal helper function 'Interval'
//
//    x   : number of observed events in the experiment
//    y   : number of observed events in background region
//    z   : number of MC events observed
//    bm  : estimate of the background
//    em  : estimate of the efficiency
//    e   : true efficiency (considered known)
//    mid : internal model id (really, you should not use this method at all)
//    sde : efficiency estimate's standard deviation 
//    sdb : background estimate's standard deviation 
//    tau : ratio parameter (read TRolke.cxx for details)
//    b   : background expectation value (considered known)
//    m   : number of MC events generated

   double dchi2 = TMath::ChisquareQuantile(fCL, 1);
   double tempxy[2], limits[2] = {0, 0};
   double slope, fmid, low, flow, high, fhigh, test, ftest, mu0, maximum, target, l, f0;
   double med = 0;
   double maxiter = 1000, acc = 0.00001;
   int i;
   int bp = 0;

   if ((mid != 3) && (mid != 5)) bm = y;
   if ((mid == 3) || (mid == 5)) {
      if (bm == 0) bm = 0.00001;
   }

   if ((mid == 6) || (mid == 7)) {
      if (bm == 0) bm = 0.00001;
   }

   if ((mid <= 2) || (mid == 4)) bp = 1;


   if (bp == 1 && x == 0 && bm > 0) {
      for (i = 0; i < 2; i++) {
         x++;
         tempxy[i] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      }

      slope = tempxy[1] - tempxy[0];
      limits[1] = tempxy[0] - slope;
      limits[0] = 0.0;
      if (limits[1] < 0) limits[1] = 0.0;
      goto done;
   }

   if (bp != 1 && x == 0) {

      for (i = 0; i < 2; i++) {
         x++;
         tempxy[i] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      }
      slope = tempxy[1] - tempxy[0];
      limits[1] = tempxy[0] - slope;
      limits[0] = 0.0;
      if (limits[1] < 0) limits[1] = 0.0;
      goto done;
   }

   if (bp != 1  && bm == 0) {
      for (i = 0; i < 2; i++) {
         bm++;
         limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
         tempxy[i] = limits[1];
      }
      slope = tempxy[1] - tempxy[0];
      limits[1] = tempxy[0] - slope;
      if (limits[1] < 0) limits[1] = 0;
      goto done;
   }

   if (x == 0 && bm == 0) {
      x++;
      bm++;
      limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      tempxy[0] = limits[1];
      x  = 1;
      bm = 2;
      limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      tempxy[1] = limits[1];
      x  = 2;
      bm = 1;
      limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      limits[1] = 3 * tempxy[0] - tempxy[1] - limits[1];
      if (limits[1] < 0) limits[1] = 0;
      goto done;
   }

   mu0 = Likelihood(0, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 1);
   maximum = Likelihood(0, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 2);
   test = 0;
   f0 = Likelihood(test, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
   if (fBounding) {
      if (mu0 < 0) maximum = f0;
   }

   target = maximum - dchi2;
   if (f0 > target) {
      limits[0] = 0;
   } else {
      if (mu0 < 0) {
         limits[0] = 0;
         limits[1] = 0;
      }

      low   = 0;
      flow  = f0;
      high  = mu0;
      fhigh = maximum;
      for (i = 0; i < maxiter; i++) {
         l = (target - fhigh) / (flow - fhigh);
         if (l < 0.2) l = 0.2;
         if (l > 0.8) l = 0.8;
         med = l * low + (1 - l) * high;
         if (med < 0.01) {
            limits[1] = 0.0;
            goto done;
         }
         fmid = Likelihood(med, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
         if (fmid > target) {
            high  = med;
            fhigh = fmid;
         } else {
            low  = med;
            flow = fmid;
         }
         if ((high - low) < acc*high) break;
      }
      limits[0] = med;
   }

   if (mu0 > 0) {
      low  = mu0;
      flow = maximum;
   } else {
      low  = 0;
      flow = f0;
   }

   test = low + 1 ;
   ftest = Likelihood(test, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
   if (ftest < target) {
      high  = test;
      fhigh = ftest;
   } else {
      slope = (ftest - flow) / (test - low);
      high  = test + (target - ftest) / slope;
      fhigh = Likelihood(high, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
   }

   for (i = 0; i < maxiter; i++) {
      l = (target - fhigh) / (flow - fhigh);
      if (l < 0.2) l = 0.2;
      if (l > 0.8) l = 0.8;
      med  = l * low + (1. - l) * high;
      fmid = Likelihood(med, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);

      if (fmid < target) {
         high  = med;
         fhigh = fmid;
      } else {
         low  = med;
         flow = fmid;
      }

      if (high - low < acc*high) break;
   }

   limits[1] = med;

done:

   // apply known efficiency
   if ((mid == 4) || (mid == 5)) {
      limits[0] /= e;
      limits[1] /= e;
   }

   fUpperLimit = limits[1];
   fLowerLimit = fmax(limits[0], 0.0);

   return limits[1];
}


double TRolke::Likelihood(double mu, int x, int y, int z, double bm, double em, int mid, double sde, double sdb, double tau, double b, int m, int what)
{
//___________________________________________________________________________
//Internal helper function
// Chooses between the different profile likelihood functions to use for the
// different models.
// Returns evaluation of the profile likelihood functions.

   switch (mid) {
      case 1:
         return EvalLikeMod1(mu, x, y, z, tau, m, what);
      case 2:
         return EvalLikeMod2(mu, x, y, em, sde, tau, what);
      case 3:
         return EvalLikeMod3(mu, x, bm, em, sde, sdb, what);
      case 4:
         return EvalLikeMod4(mu, x, y, tau, what);
      case 5:
         return EvalLikeMod5(mu, x, bm, sdb, what);
      case 6:
         return EvalLikeMod6(mu, x, z, b, m, what);
      case 7:
         return EvalLikeMod7(mu, x, em, sde, b, what);
      default:
	 std::cerr << "TRolke::Likelihood(...): Model NR: " <<
	    f_mid << " unknown"<<std::endl;
	 return 0;
   }

   return 0;
}

double TRolke::EvalLikeMod1(double mu, int x, int y, int z, double tau, int m, int what)
{
//_________________________________________________________________________
//
// Calculates the Profile Likelihood for MODEL 1:
//  Poisson background/ Binomial Efficiency
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)
   double f  = 0;
   double zm = double(z) / m;

   if (what == 1) {
      f = (x - y / tau) / zm;
   }

   if (what == 2) {
      mu = (x - y / tau) / zm;
      double b  = y / tau;
      double e = zm;
      f = LikeMod1(mu, b, e, x, y, z, tau, m);
   }

   if (what == 3) {
      if (mu == 0) {
         double b = (x + y) / (1.0 + tau);
         double e = zm;
         f = LikeMod1(mu, b, e, x, y, z, tau, m);
      } else {
         double e = 0;
         double b = 0;
         ProfLikeMod1(mu, b, e, x, y, z, tau, m);
         f = LikeMod1(mu, b, e, x, y, z, tau, m);
      }
   }

   return f;
}


double TRolke::LikeMod1(double mu, double b, double e, int x, int y, int z, double tau, int m)
{
//________________________________________________________________________
// Profile Likelihood function for MODEL 1:
// Poisson background/ Binomial Efficiency

   double s = e*mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 
   double bg = tau*b;
   double llb =  -bg;
   if ( y > 0) llb =  y*log( bg) - bg - LogFactorial(y);

   double lle = 0;  // binomial log-like
   if (z == 0)         lle = m * log(1-e); 
   else if ( z == m)   lle = m * log(e); 
   else                lle =   z * log(e) + (m - z)*log(1 - e) + LogFactorial(m) - LogFactorial(m-z) - LogFactorial(z); 

   double f = 2*( lls + llb + lle); 
   return f;
}


// this code is non-sense - // need to solve using Minuit
struct LikeFunction1 { 
};

void TRolke::ProfLikeMod1(double mu, double &b, double &e, int x, int y, int z, double tau, int m)
{
//________________________________________________________________________
// Helper for calculation of estimates of efficiency and background for model 1
//

   double med = 0.0, fmid;
   int maxiter = 1000;
   double acc = 0.00001;
   double emin = ((m + mu * tau) - sqrt((m + mu * tau) * (m + mu * tau) - 4 * mu * tau * z)) / 2 / mu / tau;

   double low  = fmax(1e-10, emin + 1e-10);
   double high = 1 - 1e-10;

   for (int i = 0; i < maxiter; i++) {
      med = (low + high) / 2.;

      fmid = LikeGradMod1(med, mu, x, y, z, tau, m);

      if (high < 0.5) acc = 0.00001 * high;
      else           acc = 0.00001 * (1 - high);

      if ((high - low) < acc*high) break;

      if (fmid > 0) low  = med;
      else         high = med;
   }

   e = med;
   double eta = double(z) / e - double(m - z) / (1 - e);

   b = double(y) / (tau - eta / mu);
}

//___________________________________________________________________________
double TRolke::LikeGradMod1(double e, double mu, int x, int y, int z, double tau, int m)
{
//gradient model likelihood
   double eta, etaprime, bprime, f;
   eta = static_cast<double>(z) / e - static_cast<double>(m - z) / (1.0 - e);
   etaprime = (-1) * (static_cast<double>(m - z) / ((1.0 - e) * (1.0 - e)) + static_cast<double>(z) / (e * e));
   double b = y / (tau - eta / mu);
   bprime = (b * b * etaprime) / mu / y;
   f = (mu + bprime) * (x / (e * mu + b) - 1) + (y / b - tau) * bprime + eta;
   return f;
}

//___________________________________________________________________________
double TRolke::EvalLikeMod2(double mu, int x, int y, double em, double sde, double tau, int what)
{
// Calculates the Profile Likelihood for MODEL 2:
//  Poisson background/ Gauss Efficiency
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)
   double v =  sde * sde;
   double coef[4], roots[3];
   double f = 0;

   if (what == 1) {
      f = (x - y / tau) / em;
   }

   if (what == 2) {
      mu = (x - y / tau) / em;
      double b = y / tau;
      double e = em;
      f = LikeMod2(mu, b, e, x, y, em, tau, v);
   }

   if (what == 3) {
      if (mu == 0) {
         double b = (x + y) / (1 + tau);
         double e = em ; 
         f = LikeMod2(mu, b, e, x, y, em, tau, v);
      } else {
         coef[3] = mu;
         coef[2] = mu * mu * v - 2 * em * mu - mu * mu * v * tau;
         coef[1] = (- x) * mu * v - mu * mu * mu * v * v * tau - mu * mu * v * em + em * mu * mu * v * tau + em * em * mu - y * mu * v;
         coef[0] = x * mu * mu * v * v * tau + x * em * mu * v - y * mu * mu * v * v + y * em * mu * v;

         TMath::RootsCubic(coef, roots[0], roots[1], roots[2]);

         double e = roots[1];
         double b; 
         if ( v > 0) b = y / (tau + (em - e) / mu / v);
         else b = y/tau; 
         f = LikeMod2(mu, b, e, x, y, em, tau, v);
      }
   }

   return f;
}

//_________________________________________________________________________
double TRolke::LikeMod2(double mu, double b, double e, int x, int y, double em, double tau, double v)
{
// Profile Likelihood function for MODEL 2:
// Poisson background/Gauss Efficiency

   double s = e*mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 
   double bg = tau*b;
   double llb =  -bg;
   if ( y > 0) llb =  y*log( bg) - bg - LogFactorial(y);
   double lle = 0; 
   if ( v > 0) lle = - 0.9189385 - log(v) / 2 - (em - e)*(em - e) / v / 2;

   return 2*( lls + llb + lle);
}

//_____________________________________________________________________
double TRolke::EvalLikeMod3(double mu, int x, double bm, double em, double sde, double sdb, int what)
{
// Calculates the Profile Likelihood for MODEL 3:
// Gauss  background/ Gauss Efficiency
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)

   double f = 0.;
   double  v = sde * sde;
   double  u = sdb * sdb;

   if (what == 1) {
      f = (x - bm) / em;
   }


   if (what == 2) {
      mu = (x - bm) / em;
      double b  = bm;
      double e  = em;
      f  = LikeMod3(mu, b, e, x, bm, em, u, v);
   }


   if (what == 3) {
      if (mu == 0.0) {
         double b = ((bm - u) + sqrt((bm - u) * (bm - u) + 4 * x * u)) / 2.;
         double e = em;
         f = LikeMod3(mu, b, e, x, bm, em, u, v);
      } else {
         double e = em; 
         double b = bm; 
         if ( v > 0) { 
            double temp[3];
            temp[0] = mu * mu * v + u;
            temp[1] = mu * mu * mu * v * v + mu * v * u - mu * mu * v * em + mu * v * bm - 2 * u * em;
            temp[2] = mu * mu * v * v * bm - mu * v * u * em - mu * v * bm * em + u * em * em - mu * mu * v * v * x;
            e = (-temp[1] + sqrt(temp[1] * temp[1] - 4 * temp[0] * temp[2])) / 2 / temp[0];
            b = bm - (u * (em - e)) / v / mu;
         }
         f = LikeMod3(mu, b, e, x, bm, em, u, v);
      }
   }

   return f;
}

//____________________________________________________________________
double TRolke::LikeMod3(double mu, double b, double e, int x, double bm, double em, double u, double v)
{
// Profile Likelihood function for MODEL 3:
// Gauss background/Gauss Efficiency

   double s = e*mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 
   double llb =  0;
   if ( u > 0) llb = - 0.9189385 - log(u) / 2 - (bm - b)*(bm - b) / u / 2;
   double lle = 0; 
   if ( v > 0) lle = - 0.9189385 - log(v) / 2 - (em - e)*(em - e) / v / 2;

   return 2*( lls + llb + lle);

}

//____________________________________________________________________
double TRolke::EvalLikeMod4(double mu, int x, int y, double tau, int what)
{
// Calculates the Profile Likelihood for MODEL 4:
// Poiss  background/Efficiency known
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)
   double f = 0.0;

   if (what == 1) f = x - y / tau;
   if (what == 2) {
      mu = x - y / tau;
      double b  = y / tau;
      f  = LikeMod4(mu, b, x, y, tau);
   }
   if (what == 3) {
      if (mu == 0.0) {
         double b = double(x + y) / (1 + tau);
         f = LikeMod4(mu, b, x, y, tau);
      } else {
         double b = (x + y - (1 + tau) * mu + sqrt((x + y - (1 + tau) * mu) * (x + y - (1 + tau) * mu) + 4 * (1 + tau) * y * mu)) / 2 / (1 + tau);
         f = LikeMod4(mu, b, x, y, tau);
      }
   }
   return f;
}

//___________________________________________________________________
double TRolke::LikeMod4(double mu, double b, int x, int y, double tau)
{
// Profile Likelihood function for MODEL 4:
// Poiss background/Efficiency known

   double s = mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 
   double bg = tau*b;
   double llb =  -bg;
   if ( y > 0) llb =  y*log( bg) - bg - LogFactorial(y);

   return 2*( lls + llb);
}

//___________________________________________________________________
double TRolke::EvalLikeMod5(double mu, int x, double bm, double sdb, int what)
{
// Calculates the Profile Likelihood for MODEL 5:
// Gauss  background/Efficiency known
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)

   double u = sdb * sdb;
   double f = 0;

   if (what == 1) {
      f = x - bm;
   }
   if (what == 2) {
      mu = x - bm;
      double b  = bm;
      f  = LikeMod5(mu, b, x, bm, u);
   }

   if (what == 3) {
      double b = ((bm - u - mu) + sqrt((bm - u - mu) * (bm - u - mu) - 4 * (mu * u - mu * bm - u * x))) / 2;
      f = LikeMod5(mu, b, x, bm, u);
   }
   return f;
}

//_______________________________________________________________________
double TRolke::LikeMod5(double mu, double b, int x, double bm, double u)
{
// Profile Likelihood function for MODEL 5:
// Gauss background/Efficiency known

   double s = mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 
   double llb =  0;
   if ( u > 0) llb = - 0.9189385 - log(u) / 2 - (bm - b)*(bm - b) / u / 2;

   return 2*( lls + llb);
}

//_______________________________________________________________________
double TRolke::EvalLikeMod6(double mu, int x, int z, double b, int m, int what)
{
// Calculates the Profile Likelihood for MODEL 6:
// Background known/Efficiency binomial
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)

   double coef[4], roots[3];
   double f = 0.;
   double zm = double(z) / m;

   if (what == 1) {
      f = (x - b) / zm;
   }

   if (what == 2) {
      mu = (x - b) / zm;
      double e  = zm;
      f  = LikeMod6(mu, b, e, x, z, m);
   }
   if (what == 3) {
      double e;
      if (mu == 0) {
         e = zm;
      } else {
         coef[3] = mu * mu;
         coef[2] = mu * b - mu * x - mu * mu - mu * m;
         coef[1] = mu * x - mu * b + mu * z - m * b;
         coef[0] = b * z;
         TMath::RootsCubic(coef, roots[0], roots[1], roots[2]);
         e = roots[1];
      }
      f = LikeMod6(mu, b, e, x, z, m);
   }
   return f;
}

//_______________________________________________________________________
double TRolke::LikeMod6(double mu, double b, double e, int x, int z, int m)
{
// Profile Likelihood function for MODEL 6:
// background known/ Efficiency binomial

   double s = e*mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 

   double lle = 0; 
   if (z == 0)        lle = m * log(1-e); 
   else if ( z == m)  lle = m * log(e); 
   else               lle =   z * log(e) + (m - z)*log(1 - e) + LogFactorial(m) - LogFactorial(m-z) - LogFactorial(z); 

   return 2* (lls + lle); 
}


//___________________________________________________________________________
double TRolke::EvalLikeMod7(double mu, int x, double em, double sde, double b, int what)
{
// Calculates the Profile Likelihood for MODEL 7:
// background known/Efficiency Gauss
// what = 1: Maximum likelihood estimate is returned
// what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
// what = 3: Profile Likelihood of Test hypothesis is returned
// otherwise parameters as described in the beginning of the class)

   double v = sde * sde;
   double f = 0.;

   if (what ==  1) {
      f = (x - b) / em;
   }

   if (what == 2) {
      mu = (x - b) / em;
      double e  = em;
      f  = LikeMod7(mu, b, e, x, em, v);
   }

   if (what == 3) {
      double e;
      if (mu == 0) {
         e = em;
      } else {
         e = (-(mu * em - b - mu * mu * v) - sqrt((mu * em - b - mu * mu * v) * (mu * em - b - mu * mu * v) + 4 * mu * (x * mu * v - mu * b * v + b * em))) / (- mu) / 2;
      }
      f = LikeMod7(mu, b, e, x, em, v);
   }

   return f;
}

//___________________________________________________________________________
double TRolke::LikeMod7(double mu, double b, double e, int x, double em, double v)
{
// Profile Likelihood function for MODEL 6:
// background known/ Efficiency gaussian

   double s = e*mu+b;
   double lls = - s;
   if (x > 0) lls = x*log(s) - s - LogFactorial(x); 

   double lle = 0; 
   if ( v > 0) lle = - 0.9189385 - log(v) / 2 - (em - e)*(em - e) / v / 2;

   return 2*( lls + lle);
}

//______________________________________________________________________
double TRolke::EvalPolynomial(double x, const int  coef[], int N)
{
// evaluate polynomial
   const int   *p;
   p = coef;
   double ans = *p++;
   int i = N;

   do
      ans = ans * x  +  *p++;
   while (--i);

   return ans;
}

//______________________________________________________________________
double TRolke::EvalMonomial(double x, const int coef[], int N)
{   
// evaluate mononomial
   double ans;
   const int   *p;

   p   = coef;
   ans = x + *p++;
   int i = N - 1;

   do
      ans = ans * x  + *p++;
   while (--i);

   return ans;
}

//______________________________________________________________________
double TRolke::LogFactorial(int n)
{
// LogFactorial function (use the logGamma function via the relation Gamma(n+1) = n!

   return lgamma(n+1);
}


