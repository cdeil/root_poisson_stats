#include <iostream>
#include "TFeldmanCousins.h"

void FeldmanCousins()
{
 // Example macro of using the TFeldmanCousins class in root.
 //
 // get a FeldmanCousins calculation object with the default limits
 // of calculating a 90% CL with the minimum signal value scanned 
 // = 0.0 and the maximum signal value scanned of 50.0
 //Author : Adrian John Bevan <bevan@SLAC.Stanford.EDU>
  
 TFeldmanCousins f;

  // calculate either the upper or lower limit for 10 observerd
  // events with an estimated background of 3.  The calculation of
  // either upper or lower limit will return that limit and fill
  // data members with both the upper and lower limit for you.
  double Nobserved   = 10.0;
  double Nbackground = 3.0;

  double ul = f.CalculateUpperLimit(Nobserved, Nbackground);
  double ll = f.GetLowerLimit();

  std::cout << "For " <<  Nobserved << " data observed with and estimated background"<<std::endl;
  std::cout << "of " << Nbackground << " candidates, the Feldman-Cousins method of "<<std::endl;
  std::cout << "calculating confidence limits gives:"<<std::endl;
  std::cout << "\tUpper Limit = " <<  ul << std::endl;
  std::cout << "\tLower Limit = " <<  ll << std::endl;
  std::cout << "at the 90% CL"<< std::endl;
}

int main(void) {
   FeldmanCousins();
   return 0;
}