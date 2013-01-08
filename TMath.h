/* ROOT TMath functions needed for TRolke and TFeldmanCousins implementation */

namespace TMath
{
   double Gamma(double a, double x); // only needed in TMath by ChisquareQuantile
   double ChisquareQuantile(double p, double ndf);
   double NormQuantile(double p); // only needed in TMath by ChisquareQuantile
   double Poisson(double x, double par);
   double PoissonI(double x, double par);
   bool RootsCubic(const double coef[4], double &a, double &b, double &c);
   void Sort(int n, const double* a, int* index);
}
