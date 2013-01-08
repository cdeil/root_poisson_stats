// @(#)root/physics:$Id: TFeldmanCousins.cxx 44507 2012-06-04 12:30:41Z axel $
// Author: Adrian Bevan  2001

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * Copyright (C) 2001, Liverpool University.                             *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

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
//
// see example of use in $ROOTSYS/tutorials/math/FeldmanCousins.C
//
// see note about: "Should I use TRolke, TFeldmanCousins, TLimit?"
//  in the TRolke class description.
//
// Author: Adrian Bevan, Liverpool University
//
// Copyright Liverpool University 2001       bevan@slac.stanford.edu
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMath.h"
#include "TFeldmanCousins.h"

//______________________________________________________________________________
TFeldmanCousins::TFeldmanCousins(double newFC, bool be_quick)
{
   //constructor
   fCL          = newFC;
   fUpperLimit  = 0.0;
   fLowerLimit  = 0.0;
   fNobserved   = 0.0;
   fNbackground = 0.0;
   if (be_quick) fQUICK = 1;
   else          fQUICK = 0;

   fNMax   = 50;
   fMuStep = 0.005;
   SetMuMin();
   SetMuMax();
   SetMuStep();
}


//______________________________________________________________________________
TFeldmanCousins::~TFeldmanCousins()
{
}


//______________________________________________________________________________
double TFeldmanCousins::CalculateLowerLimit(double Nobserved, double Nbackground)
{
////////////////////////////////////////////////////////////////////////////////////////////
// given Nobserved and Nbackground, try different values of mu that give lower limits that//
// are consistent with Nobserved.  The closed interval (plus any stragglers) corresponds  //
// to the F&C interval                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////

   CalculateUpperLimit(Nobserved, Nbackground);
   return fLowerLimit;
}


//______________________________________________________________________________
double TFeldmanCousins::CalculateUpperLimit(double Nobserved, double Nbackground)
{
////////////////////////////////////////////////////////////////////////////////////////////
// given Nobserved and Nbackground, try different values of mu that give upper limits that//
// are consistent with Nobserved.  The closed interval (plus any stragglers) corresponds  //
// to the F&C interval                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////

   fNobserved   = Nobserved;
   fNbackground = Nbackground;

   double mu = 0.0;

   // for each mu construct the ranked table of probabilities and test the
   // observed number of events with the upper limit
   double min = -999.0;
   double max = 0;
   int iLower = 0;

   int i;
   for(i = 0; i <= fNMuStep; i++) {
      mu = fMuMin + (double)i*fMuStep;
      int goodChoice = FindLimitsFromTable( mu );
      if( goodChoice ) {
         min = mu;
         iLower = i;
         break;
      }
   }

   //==================================================================
   // For quicker evaluation, assume that you get the same results when
   // you expect the uppper limit to be > Nobserved-Nbackground.
   // This is certainly true for all of the published tables in the PRD
   // and is a reasonable assumption in any case.
   //==================================================================

   double quickJump = 0.0;
   if (fQUICK)          quickJump = Nobserved-Nbackground-fMuMin;
   if (quickJump < 0.0) quickJump = 0.0;

   for(i = iLower+1; i <= fNMuStep; i++) {
      mu = fMuMin + (double)i*fMuStep + quickJump;
      int goodChoice = FindLimitsFromTable( mu );
      if( !goodChoice ) {
         max = mu;
         break;
      }
   }

   fUpperLimit = max;
   fLowerLimit = min;

   return max;
}


//______________________________________________________________________________
int TFeldmanCousins::FindLimitsFromTable( double mu )
{
///////////////////////////////////////////////////////////////////
// calculate the probability table for a given mu for n = 0, NMAX//
// and return 1 if the number of observed events is consistent   //
// with the CL bad                                               //
///////////////////////////////////////////////////////////////////

   double *p          = new double[fNMax];   //the array of probabilities in the interval MUMIN-MUMAX
   double *r          = new double[fNMax];   //the ratio of likliehoods = P(Mu|Nobserved)/P(MuBest|Nobserved)
   int    *rank       = new int[fNMax];      //the ranked array corresponding to R (largest first)
   double *muBest     = new double[fNMax];
   double *probMuBest = new double[fNMax];

   //calculate P(i | mu) and P(i | mu)/P(i | mubest)
   int i;
   for(i = 0; i < fNMax; i++) {
      muBest[i] = (double)(i - fNbackground);
      if(muBest[i]<0.0) muBest[i] = 0.0;
      probMuBest[i] = Prob(i, muBest[i],  fNbackground);
      p[i]          = Prob(i, mu,  fNbackground);
      if(probMuBest[i] == 0.0) r[i] = 0.0;
      else                     r[i] = p[i]/probMuBest[i];
   }

   //rank the likelihood ratio
   TMath::Sort(fNMax, r, rank);

   //search through the probability table and get the i for the CL
   double sum = 0.0;
   int iMax = rank[0];
   int iMin = rank[0];
   for(i = 0; i < fNMax; i++) {
      sum += p[rank[i]];
      if(iMax < rank[i]) iMax = rank[i];
      if(iMin > rank[i]) iMin = rank[i];
      if(sum >= fCL) break;
   }

   delete [] p;
   delete [] r;
   delete [] rank;
   delete [] muBest;
   delete [] probMuBest;

   if((fNobserved <= iMax) && (fNobserved >= iMin)) return 1;
   else return 0;
}


//______________________________________________________________________________
double TFeldmanCousins::Prob(int N, double mu, double B)
{
////////////////////////////////////////////////
// calculate the poissonian probability for   //
// a mean of mu+B events with a variance of N //
////////////////////////////////////////////////

   return TMath::Poisson( N, mu+B);
}

//______________________________________________________________________________
void TFeldmanCousins::SetMuMax(double newMax)
{
   //set maximum value of signal to use in calculating the tables
   fMuMax   = newMax;
   fNMax    = (int)newMax;
   SetMuStep(fMuStep);
}

//______________________________________________________________________________
void TFeldmanCousins::SetMuStep(double newMuStep)
{
   //set the step in signal to use when generating tables
   if(newMuStep == 0.0) {
      std::cout << "TFeldmanCousins::SetMuStep ERROR New step size is zero - unable to change value"<< std::endl;
      return;
   } else {
      fMuStep = newMuStep;
      fNMuStep = (int)((fMuMax - fMuMin)/fMuStep);
   }
}
