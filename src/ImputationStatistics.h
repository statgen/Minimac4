#ifndef __IMPUTATIONSTATISTICS_H__
#define __IMPUTATIONSTATISTICS_H__

#include "MathVector.h"
#include "HaplotypeSet.h"
#include "IntArray.h"

class ImputationStatistics
   { public:



       int numRefMarkers, numTarMarkers;


      ImputationStatistics()
      {}
      ~ImputationStatistics()
      {

      }

void Initialize(int nRef, int nTar)
{
    numRefMarkers=nRef;
    numTarMarkers=nTar;

    fill(sum.begin(), sum.end(), 0.0);
    fill(sumSq.begin(), sumSq.end(), 0.0);
    fill(sumCall.begin(), sumCall.end(), 0.0);
    fill(looSum.begin(), looSum.end(), 0.0);
    fill(looSumSq.begin(), looSumSq.end(), 0.0);
    fill(looProduct.begin(), looProduct.end(), 0.0);
    fill(looObserved.begin(), looObserved.end(), 0.0);
    fill(count.begin(), count.end(), 0);
    fill(looCount.begin(), looCount.end(), 0);

}


void PreInitialize(int markers)
   {
       sum.resize(markers, 0.0);
       sumSq.resize(markers, 0.0);
       sumCall.resize(markers, 0.0);
       looSum.resize(markers, 0.0);
       looSumSq.resize(markers, 0.0);
       looProduct.resize(markers, 0.0);
       looObserved.resize(markers, 0.0);
       count.resize(markers,0);
       looCount.resize(markers,0);
   }


      void Update(vector<float> &doses, vector<float> &loo,vector<bool> &observed,vector<bool> &Miss, vector<bool> &major);
     void NewUpdate(HaplotypeSet &rHap, HaplotypeSet &tHap, int hapId,
                                     vector<float> *doses, vector<float> *loo);
      double Rsq(int marker);
      double AlleleFrequency(int marker);
      double AverageCallScore(int marker);
      double LooRsq(int marker);
      double LooMajorDose(int marker);
      double LooMinorDose(int marker);
      double EmpiricalR(int marker);
      double EmpiricalRsq(int marker);
      double GoldenR(int marker);
      double GoldenRsq(int marker);
      double GoldenRgeno(int marker);
      double GoldenRsqgeno(int marker);

  // private:
      vector<double>   sum, sumSq, sumCall, looSum, looSumSq, looProduct, looObserved, looObservedSq;
      vector<int> count, looCount;
   };

#endif
