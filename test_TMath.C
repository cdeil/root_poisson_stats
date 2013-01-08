#include <iostream>
#include "TMath.h"

/* Basic tests of the functions in TMath
   Reference values computed at the ROOT prompt.
*/
int main(void) {
   // double Gamma(double a, double x); // only needed in TMath by ChisquareQuantile
	std::cout << TMath::Gamma(42, 43) << " " << 5.81002199739738767e-01 << std::endl;

   // double ChisquareQuantile(double p, double ndf);
	std::cout << TMath::ChisquareQuantile(0.1, 2) << " " << 2.10721031316083163e-01 << std::endl;

   // double NormQuantile(double p); // only needed in TMath by ChisquareQuantile
	std::cout << TMath::NormQuantile(0.1) << " " << -1.28155156554460081e+00 << std::endl;

   // double Poisson(double x, double par);
	std::cout << TMath::Poisson(2.2, 3.3) << " " << 2.10393562195771128e-01 << std::endl;

   // double PoissonI(double x, double par);
	std::cout << TMath::PoissonI(2.2, 3.3) << " " << 2.00828846499751856e-01 << std::endl;

   // bool RootsCubic(const double coef[4], double &a, double &b, double &c);
	double coef[4] = {1, 1, 1, 1};
	double a=0, b=0, c=0;
	TMath::RootsCubic(coef, a, b, c);
	std::cout << a << " " << b << " " << c << " --- ";
	std::cout << "-1 -1.11022e-16 1" << std::endl;

   // void Sort(int n, const double* a, int* index);
	int n = 5;
	double aa[5] = {1, 5, 4, 2, 3};
	int ii[5] = {0, 0, 0, 0, 0};
	TMath::Sort(n, aa, ii);
	std::cout << ii[0] << " " << ii[1] << " " << ii[2] << " " << ii[3] << " " << ii[4] << " --- ";
	std::cout << "1 2 4 3 0" << std::endl;
	std::cout << aa[ii[0]] << " " << aa[ii[1]] << " " << aa[ii[2]] << " " << aa[ii[3]] << " " << aa[ii[4]] << " --- ";
	std::cout << "5 4 3 2 1" << std::endl;
   return 0;
}