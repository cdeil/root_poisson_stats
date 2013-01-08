// Example of the usage of the TRolke class 
#include <iostream>
#include "TRolke.h"
      
void Rolke()
{
//////////////////////////////////////////////////////////
//
// The TRolke class computes the profile likelihood
// confidence limits for 7 different model assumptions
// on systematic/statistical uncertainties
//
// Author : Jan Conrad (CERN) <jan.conrad@cern.ch> 2004
//          Johan Lundberg (CERN) <johan.lundberg@cern.ch> 2009
//  
// Please read TRolke.cxx and TRolke.h for more docs.
//             ----------     --------
//
//////////////////////////////////////////////////////


   /* variables used throughout the example */
   double bm;
   double tau;
   //int mid;
   int m;
   int z;
   int y;
   int x;
   double e;
   double em;
   double sde;
   double sdb;
   double b;

   double alpha; //Confidence Level

   // make TRolke objects
   TRolke tr;   //

   double ul ; // upper limit 
   double ll ; // lower limit


/////////////////////////////////////////////////////////////
// Model 1 assumes:
//
// Poisson uncertainty in the background estimate
// Binomial uncertainty in the efficiency estimate
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =1;
   x = 5;     // events in the signal region
   y = 10;    // events observed in the background region
   tau = 2.5; // ratio between size of signal/background region
   m = 100;   // MC events have been produced  (signal)
   z = 50;    // MC events have been observed (signal)          

   alpha=0.9; //Confidence Level

   tr.SetCL(alpha);  

   tr.SetPoissonBkgBinomEff(x,y,z,tau,m); 
   tr.GetLimits(ll,ul);
 
   std::cout << "For model 1: Poisson / Binomial" << std::endl; 
   std::cout << "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;

 
/////////////////////////////////////////////////////////////
// Model 2 assumes:
//
// Poisson uncertainty in the background estimate
// Gaussian  uncertainty in the efficiency estimate
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =2;
   y = 3 ;      // events observed in the background region
   x = 10 ;     // events in the signal region
   tau = 2.5;   // ratio between size of signal/background region
   em = 0.9;    // measured efficiency
   sde = 0.05;  // standard deviation of efficiency
   alpha =0.95; // Confidence L evel

   tr.SetCL(alpha);

   tr.SetPoissonBkgGaussEff(x,y,em,tau,sde);
   tr.GetLimits(ll,ul);
 
   std::cout << "For model 2 : Poisson / Gaussian" << std::endl; 
   std::cout << "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;

  

/////////////////////////////////////////////////////////////
// Model 3 assumes:
//
// Gaussian uncertainty in the background estimate
// Gaussian  uncertainty in the efficiency estimate
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =3;
   bm = 5;      // expected background
   x = 10;      // events in the signal region
   sdb = 0.5;   // standard deviation in background estimate
   em = 0.9;    //  measured efficiency
   sde = 0.05;  // standard deviation of efficiency
   alpha =0.99; // Confidence Level

   tr.SetCL(alpha);

   tr.SetGaussBkgGaussEff(x,bm,em,sde,sdb); 
   tr.GetLimits(ll,ul);
   std::cout << "For model 3 : Gaussian / Gaussian" << std::endl; 
   std::cout << "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;


 
   std::cout << "***************************************" << std::endl;
   std::cout << "* some more example's for gauss/gauss *" << std::endl;
   std::cout << "*                                     *" << std::endl;
   double slow,shigh;
   tr.GetSensitivity(slow,shigh);
   std::cout << "sensitivity:" << std::endl;
   std::cout << "[" << slow << "," << shigh << "]" << std::endl; 

   int outx;
   tr.GetLimitsQuantile(slow,shigh,outx,0.5);
   std::cout << "median limit:" << std::endl;
   std::cout << "[" << slow << "," << shigh << "] @ x =" << outx <<std::endl; 

   tr.GetLimitsML(slow,shigh,outx);
   std::cout << "ML limit:" << std::endl;
   std::cout << "[" << slow << "," << shigh << "] @ x =" << outx <<std::endl; 

   tr.GetSensitivity(slow,shigh);
   std::cout << "sensitivity:" << std::endl;
   std::cout << "[" << slow << "," << shigh << "]" << std::endl; 

   tr.GetLimits(ll,ul);
   std::cout << "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;

   int ncrt;

   tr.GetCriticalNumber(ncrt);
   std::cout << "critical number: " << ncrt << std::endl;

   tr.SetCLSigmas(5);
   tr.GetCriticalNumber(ncrt);
   std::cout << "critical number for 5 sigma: " << ncrt << std::endl;

   std::cout << "***************************************" << std::endl;


/////////////////////////////////////////////////////////////
// Model 4 assumes:
//
// Poisson uncertainty in the background estimate
// known efficiency
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =4;
   y = 7;       // events observed in the background region
   x = 1;       // events in the signal region
   tau = 5;     // ratio between size of signal/background region
   e = 0.25;    // efficiency 

   alpha =0.68; // Confidence L evel

   tr.SetCL(alpha);

   tr.SetPoissonBkgKnownEff(x,y,tau,e);
   tr.GetLimits(ll,ul);
 
   std::cout << "For model 4 : Poissonian / Known" << std::endl; 
   std::cout <<  "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;

   
////////////////////////////////////////////////////////
// Model 5 assumes:
//
// Gaussian uncertainty in the background estimate
// Known efficiency
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =5;
   bm = 0;          // measured background expectation
   x = 1 ;          // events in the signal region
   e = 0.65;        // known eff
   sdb = 1.0;       // standard deviation of background estimate
   alpha =0.799999; // Confidence Level

   tr.SetCL(alpha);

   tr.SetGaussBkgKnownEff(x,bm,sdb,e);
   tr.GetLimits(ll,ul);
 
   std::cout << "For model 5 : Gaussian / Known" << std::endl; 
   std::cout <<  "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;

 

////////////////////////////////////////////////////////
// Model 6 assumes:
//
// Known background 
// Binomial uncertainty in the efficiency estimate
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =6;
   b = 10;      // known background
   x = 25;      // events in the signal region
   z = 500;     // Number of observed signal MC events
   m = 750;     // Number of produced MC signal events
   alpha =0.9;  // Confidence L evel

   tr.SetCL(alpha);

   tr.SetKnownBkgBinomEff(x, z,m,b);
   tr.GetLimits(ll,ul);
 
   std::cout << "For model 6 : Known / Binomial" << std::endl; 
   std::cout <<  "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;

  
////////////////////////////////////////////////////////
// Model 7 assumes:
//
// Known Background
// Gaussian  uncertainty in the efficiency estimate
//
   std::cout << std::endl<<" ======================================================== " <<std::endl;
   //mid =7;
   x = 15;      // events in the signal region
   em = 0.77;   // measured efficiency
   sde = 0.15;  // standard deviation of efficiency estimate
   b = 10;      // known background
   alpha =0.95; // Confidence L evel

   y = 1;

   tr.SetCL(alpha);

   tr.SetKnownBkgGaussEff(x,em,sde,b);
   tr.GetLimits(ll,ul);
  
   std::cout << "For model 7 : Known / Gaussian " << std::endl; 
   std::cout <<  "the Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;


////////////////////////////////////////////////////////
// Example of bounded and unbounded likelihood
// Example for Model 1
///////////////////////////////////////////////////////

   bm = 0.0;
   tau = 5;
   //mid = 1;
   m = 100;
   z = 90;
   y = 15;
   x = 0;
   alpha = 0.90;
   
   tr.SetCL(alpha);
   tr.SetPoissonBkgBinomEff(x,y,z,tau,m); 
   tr.SetBounding(true); //bounded
   tr.GetLimits(ll,ul);   
   
   std::cout << "Example of the effect of bounded vs unbounded, For model 1" << std::endl; 
   std::cout <<  "the BOUNDED Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;


   tr.SetBounding(false); //unbounded
   tr.GetLimits(ll,ul);   
   
   std::cout <<  "the UNBOUNDED Profile Likelihood interval is :" << std::endl;
   std::cout << "[" << ll << "," << ul << "]" << std::endl;
  
}

int main(void) {
   Rolke();
   return 0;
}