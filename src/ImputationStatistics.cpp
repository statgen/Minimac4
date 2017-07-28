#include "ImputationStatistics.h"

#include <math.h>
#include <iostream>

using namespace std;



void ImputationStatistics::NewUpdate(HaplotypeSet &rHap,HaplotypeSet &tHap, int hapId, vector<float> *doses, vector<float> *loo)
{
    for (int i = 0; i < numRefMarkers; i++)
    {
        float &DoseVal=(*doses)[i];
        sum[i] += DoseVal;
        sumSq[i] += DoseVal * DoseVal;
        sumCall[i] += DoseVal > 0.5 ? DoseVal : 1.0 - DoseVal;
        count[i] ++;
    }

    vector<AlleleType> &UnScaff = tHap.haplotypesUnscaffolded[hapId];
    vector<AlleleType> &MissUnScaff = tHap.MissingSampleUnscaffolded[hapId];

    int index=0;
    for (int i = 0; i < numRefMarkers; i++)
    {
        if(!rHap.Targetmissing[i])
        {
//            assert(index<tHap.numMarkers);
            if(MissUnScaff[index]=='0')
            {
                AlleleType observed=UnScaff[index];
                float &LooVal = (*loo)[index];
                looSum[i] += LooVal;
                looSumSq[i] += LooVal * LooVal;
                looProduct[i] += (observed == '1') ? LooVal : 0.0;
                looObserved[i] += (observed == '1') ? 1.0 : 0.0;
                looCount[i]++;
            }


//            assert(rHap.MapRefToTar[i]==index);
            index++;
        }
    }
//    assert(index==numTarMarkers);

}

double ImputationStatistics::Rsq(int marker)
   {
   if (count[marker] < 2)
      return 0.0;

   double f = sum[marker] / (count[marker] + 1e-30);
   double evar = f * (1.0 - f);
   double ovar = 0.0;

   if((sumSq[marker] - sum[marker] * sum[marker] / (count[marker] + 1e-30))>0)
        ovar=(sumSq[marker] - sum[marker] * sum[marker] / (count[marker] + 1e-30)) / (count[marker] + 1e-30);

   return ovar / (evar + 1e-30);
   }

double ImputationStatistics::LooRsq(int marker)
   {
    if (looCount[marker] < 2)
        return 0.0;

    double f = looSum[marker] / (looCount[marker] + 1e-30);
    double evar = f * (1.0 - f);
    double ovar = 0.0;

    if(((looSumSq[marker] - looSum[marker] * looSum[marker] / (looCount[marker] + 1e-30)))>0)
        ovar = (looSumSq[marker] - looSum[marker] * looSum[marker] / (looCount[marker] + 1e-30)) / (looCount[marker] + 1e-30);

    return ovar / (evar + 1e-30);
   }

double ImputationStatistics::AlleleFrequency(int marker)
   {
   if (count[marker] < 2)
      return 0.0;

   return sum[marker] / (count[marker] + 1e-30);
   }

double ImputationStatistics::EmpiricalR(int marker)
   {
   if (looCount[marker] < 2)
      return 0.0;

   // n * Sum xy - Sum x * Sum y
    double p = looCount[marker] * looProduct[marker] - looSum[marker] * looObserved[marker];
    double qx=0.0,qy=0.0;
   // sqrt(n*Sum xx - Sum x * Sum x)


    if((looCount[marker] * looSumSq[marker] - looSum[marker] * looSum[marker])>0)
        qx = sqrt((looCount[marker] * looSumSq[marker] - looSum[marker] * looSum[marker]));

    if((looCount[marker] * looObserved[marker] - looObserved[marker] * looObserved[marker])>0)
        qy = sqrt((looCount[marker] * looObserved[marker] - looObserved[marker] * looObserved[marker]));

   if (qx / (qy + 1e-30) < 1e-3)
      return 0.0;

   if (qy / (qx + 1e-30) < 1e-3)
      return 0.0;

   double r = p / (qx * qy + 1e-30);

   return r;
   }

double ImputationStatistics::EmpiricalRsq(int marker)
   {
   double r = EmpiricalR(marker);

   return r * r;
   }

double ImputationStatistics::LooMajorDose(int marker)
   {
   return looProduct[marker] / (looObserved[marker] + 1e-30);
   }

double ImputationStatistics::LooMinorDose(int marker)
   {
   return (looSum[marker] - looProduct[marker]) / (looCount[marker] - looObserved[marker] + 1e-30);
   }

double ImputationStatistics::AverageCallScore(int marker)
   {
   return sumCall[marker] / (count[marker] + 1e-30);
   }

