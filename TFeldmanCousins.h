// @(#)root/physics:$Id: TFeldmanCousins.h 20882 2007-11-19 11:31:26Z rdm $
// Author: Adrian Bevan  2001

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * Copyright (C) 2001, Liverpool University.                             *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFeldmanCousins
#define ROOT_TFeldmanCousins

////////////////////////////////////////////////////////////////////////////
// TFeldmanCousins
//
// class to calculate the CL upper limit using
// the Feldman-Cousins method as described in PRD V57 #7, p3873-3889
//
// The default confidence interval calvculated using this method is 90%
// This is set either by having a default the constructor, or using the
// appropriate fraction when instantiating an object of this class (e.g. 0.9)
//
// The simple extension to a gaussian resolution function bounded at zero
// has not been addressed as yet -> `time is of the essence' as they write
// on the wall of the maze in that classic game ...
//
//    VARIABLES THAT CAN BE ALTERED
//    -----------------------------
// => depending on your desired precision: The intial values of fMuMin,
// fMuMax, fMuStep and fNMax are those used in the PRD:
//   fMuMin = 0.0
//   fMuMax = 50.0
//   fMuStep= 0.005
// but there is total flexibility in changing this should you desire.
//
// Author: Adrian Bevan, Liverpool University
//
// Copyright Liverpool University 2001       bevan@slac.stanford.edu
///////////////////////////////////////////////////////////////////////////


class TFeldmanCousins {
protected:
   double fCL;         // confidence level as a fraction [e.g. 90% = 0.9]
   double fUpperLimit; // the calculated upper limit
   double fLowerLimit; // the calculated lower limit
   double fNobserved;  // input number of observed events
   double fNbackground;// input number of background events
   double fMuMin;      // minimum value of signal to use in calculating the tables
   double fMuMax;      // maximum value of signal to use in calculating the tables
   double fMuStep;     // the step in signal to use when generating tables
   int    fNMuStep;    // = (int)(fMuStep)
   int    fNMax;       // = (int)(fMuMax)
   int    fQUICK;      // take a short cut to speed up the process of generating a
                        // lut.  This scans from Nobserved-Nbackground-fMuMin upwards
                        // assuming that UL > Nobserved-Nbackground.

   ////////////////////////////////////////////////
   // calculate the poissonian probability for   //
   // a mean of mu+B events with a variance of N //
   ////////////////////////////////////////////////
   double Prob(int N, double mu, double B);

   ////////////////////////////////////////////////
   // calculate the probability table and see if //
   // fNObserved is in the 100.0 * fCL %         //
   // interval                                   //
   ////////////////////////////////////////////////
   int FindLimitsFromTable(double mu);

public:
   TFeldmanCousins(double newCL=0.9, bool be_quick=false);
   virtual ~TFeldmanCousins();

   ////////////////////////////////////////////////
   // calculate the upper limit given Nobserved  //
   // and Nbackground events                     //
   // the variables fUpperLimit and fLowerLimit  //
   // are set before returning the upper limit   //
   ////////////////////////////////////////////////
   double CalculateUpperLimit(double Nobserved, double Nbackground);
   double CalculateLowerLimit(double Nobserved, double Nbackground);

   inline double GetUpperLimit(void)  const { return fUpperLimit;  }
   inline double GetLowerLimit(void)  const { return fLowerLimit;  }
   inline double GetNobserved(void)   const { return fNobserved;   }
   inline double GetNbackground(void) const { return fNbackground; }
   inline double GetCL(void)          const { return fCL;          }

   inline double GetMuMin(void)       const { return fMuMin;  }
   inline double GetMuMax(void)       const { return fMuMax;  }
   inline double GetMuStep(void)      const { return fMuStep; }
   inline double GetNMax(void)        const { return fNMax;   }

   inline void     SetNobserved(double NObs)         { fNobserved   = NObs;  }
   inline void     SetNbackground(double Nbg)        { fNbackground = Nbg;   }
   inline void     SetCL(double newCL)               { fCL          = newCL; }

   inline void     SetMuMin(double  newMin    = 0.0)   { fMuMin = newMin;  }
   void            SetMuMax(double  newMax    = 50.0);
   void            SetMuStep(double newMuStep = 0.005);

};

#endif






