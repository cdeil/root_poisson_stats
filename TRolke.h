
//////////////////////////////////////////////////////////////////////////////
//
//  TRolke
//
//  This class computes confidence intervals for the rate of a Poisson
//  in the presence of background and efficiency with a fully frequentist
//  treatment of the uncertainties in the efficiency and background estimate
//  using the profile likelihood method.
//
//      Author: Jan Conrad (CERN) 2004
//      Updated: Johan Lundberg (CERN) 2009
//
//      Copyright CERN 2004,2009           Jan.Conrad@cern.ch,
//                                     Johan.Lundberg@cern.ch
//
//  For information about the statistical meaning of the parameters
//  and the syntax, consult TRolke.cxx
//                  ------------------
//
//  Examples are found in the file Rolke.C
//  --------------------------------------
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TRolke
#define ROOT_TRolke

// Class definition. This class is not intended to be used as a base class.
class TRolke
{

private:
   double fCL;         // confidence level as a fraction [0.9 for 90% ]
   double fUpperLimit; // the calculated upper limit
   double fLowerLimit; // the calculated lower limit
   bool fBounding;       // false for unbounded likelihood
                         // true for bounded likelihood   
   int fNumWarningsDeprecated1;
   int fNumWarningsDeprecated2;

   /* ----------------------------------------------------------------- */
   /* These variables are set by the Set methods for the various models */
   int f_x;
   int f_y;
   int f_z;
   double f_bm;
   double f_em;
   double f_e;
   int f_mid;
   double f_sde;
   double f_sdb;
   double f_tau;
   double f_b;
   int f_m;

   /* ----------------------------------------------------------------- */
   /* Internal helper functions and methods */
   // The Calculator
   double Interval(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m);

   // LIKELIHOOD ROUTINE
   double Likelihood(double mu, int x, int y, int z, double bm, double em, int mid, double sde, double sdb, double tau, double b, int m, int what);

   //MODEL 1
   double EvalLikeMod1(double mu, int x, int y, int z, double tau, int m, int what);
   double LikeMod1(double mu, double b, double e, int x, int y, int z, double tau, int m);
   void     ProfLikeMod1(double mu, double &b, double &e, int x, int y, int z, double tau, int m);
   double LikeGradMod1(double e, double mu, int x, int y, int z, double tau, int m);

   //MODEL 2
   double EvalLikeMod2(double mu, int x, int y, double em, double sde, double tau, int what);

   double LikeMod2(double mu, double b, double e, int x, int y, double em, double tau, double v);

   //MODEL 3
   double EvalLikeMod3(double mu, int x, double bm, double em, double sde, double sdb, int what);
   double LikeMod3(double mu, double b, double e, int x, double bm, double em, double u, double v);

   //MODEL 4
   double EvalLikeMod4(double mu, int x, int y, double tau, int what);
   double LikeMod4(double mu, double b, int x, int y, double tau);

   //MODEL 5
   double EvalLikeMod5(double mu, int x, double bm, double sdb, int what);
   double LikeMod5(double mu, double b, int x, double bm, double u);

   //MODEL 6
   double EvalLikeMod6(double mu, int x, int z, double b, int m, int what);
   double LikeMod6(double mu, double b, double e, int x, int z, int m);

   //MODEL 7
   double EvalLikeMod7(double mu, int x, double em, double sde, double b, int what);
   double LikeMod7(double mu, double b, double e, int x, double em, double v);

   //MISC
   static double EvalPolynomial(double x, const int coef[], int N);
   static double EvalMonomial(double x, const int coef[], int N);
   double LogFactorial(int n);

   double ComputeInterval(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m);

   void SetModelParameters(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m);

   void SetModelParameters();

   double GetBackground();

public:

   /* Constructor */
   TRolke(double CL = 0.9);

   /* Destructor */
   virtual ~TRolke();

   /* Get and set the Confidence Level */
   double GetCL() const         {
      return fCL;
   }
   void     SetCL(double CL)  {
      fCL = CL;
   }

   /* Set the Confidence Level in terms of Sigmas. */
   void SetCLSigmas(double CLsigmas);
   
   // The Set methods for the different models are described in Rolke.cxx
   // model 1
   void SetPoissonBkgBinomEff(int x, int y, int z, double tau, int m);

   // model 2
   void SetPoissonBkgGaussEff(int x, int y, double em, double tau, double sde);

   // model 3
   void SetGaussBkgGaussEff(int x, double bm, double em, double sde, double sdb);

   // model 4
   void SetPoissonBkgKnownEff(int x, int y, double tau, double e);

   // model 5
   void SetGaussBkgKnownEff(int x, double bm, double sdb, double e);

   // model 6
   void SetKnownBkgBinomEff(int x, int z, int m, double b);

   // model 7
   void SetKnownBkgGaussEff(int x, double em, double sde, double b);

   /* Deprecated interface method (read Rolke.cxx). May be removed from future releases */
   double CalculateInterval(int x, int y, int z, double bm, double em, double e, int mid, double sde, double sdb, double tau, double b, int m);

   // get the upper and lower limits based on the specified model
   bool GetLimits(double& low, double& high);
   double GetUpperLimit();
   double GetLowerLimit();

   // get the upper and lower average limits
   bool GetSensitivity(double& low, double& high, double pPrecision = 0.00001);

   // get the upper and lower limits for the outcome corresponding to
   // a given quantile.
   bool GetLimitsQuantile(double& low, double& high, int& out_x, double integral = 0.5);

   // get the upper and lower limits for the most likely outcome.
   bool GetLimitsML(double& low, double& high, int& out_x);

   // get the value of x corresponding to rejection of the null hypothesis.
   bool GetCriticalNumber(int& ncrit,int maxtry=-1);

   /* Get the bounding mode flag. True activates bounded mode. Read
      TRolke.cxx and the references therein for details. */
   bool GetBounding() const {
      return fBounding;
   }

   /* Get the bounding mode flag. True activates bounded mode. Read
      TRolke.cxx and the references therein for details. */
   void SetBounding(const bool bnd) {
      fBounding = bnd;
   }

   /* Deprecated name for SetBounding. */
   void SetSwitch(bool bnd) ;

   /* Dump internals. Option is not used */
   void Print() const;

};

//calculate confidence limits using the Rolke method
#endif

