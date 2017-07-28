#ifndef IMPUTATION_H_INCLUDED
#define IMPUTATION_H_INCLUDED

#include "MyVariables.h"

#include "MarkovParameters.h"

#include "HaplotypeSet.h"
#include "DosageData.h"
#include "MemoryInfo.h"
#include "MemoryAllocators.h"
#include "MarkovModel.h"
#include "MarkovParameters.h"
#include "MyVariables.h"
#include <ctime>
#include "Unique.h"
#include <cstdio>
#include <cmath>


#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;




class Imputation
    {
        public:

            // Haplotype Variables and Other Basic Variables
            HaplotypeSet *THapUnchunked, *rHapChunked;
            AllVariable *MyAllVariables;
            int ChunkNo, TotalNovcfParts;
            vector<MarkovModel> MainMarkovModel;


            // File Output Stream handling
            IFILE dosages, hapdose, haps, vcfdosepartial, info;
            ImputationStatistics *stats;

            // Dosage Date for Output
            DosageData SinglePartialDosageData;
            DosageData* CurrentPartialDosageData;

            // Time Variables
            int TimeToCompress, TimeToImpute, TimeToWrite;


            double Dosagesize()
            {
                double S=0;
                S+=SinglePartialDosageData.size();
                return (S);
            };

            double Probsize()
            {
                double S=0;
                for(int i=0;i<MyAllVariables->myModelVariables.cpus;i++)
                {
                    S+=MainMarkovModel[i].size();
                }
                return (S);
            };


void Minimac3ImputeThisChunk(int ChunkId, HaplotypeSet &FullrHap, HaplotypeSet &tgwasHap, HaplotypeSet &rgwasHap);


        void                            performImputationBasedonTypedSites(HaplotypeSet &tHap,
                                                                    HaplotypeSet &rHap,
                                                                    HaplotypeSet &tHapOrig,
                                                                    HaplotypeSet &FullrHap);
        MarkovParameters*               createEstimates             (HaplotypeSet &rHap,HaplotypeSet &tHap,vector<int> &optStructure,bool NoTargetEstimation);
        void                            splitFoldedProb             (vector<float> &SplitProb,vector<float> &totalProb, vector<float> &noRecomProb);
        void                            performImputation           (HaplotypeSet &tHap,HaplotypeSet &rHap, String Golden);
        void                            ImputeThisChunk             (int ChunkId, HaplotypeSet &rHapOriginal
                                                                    ,HaplotypeSet &tgwasHap, HaplotypeSet &rgwasHap);
        void                            performImputationNew        (HaplotypeSet &tHap,HaplotypeSet &rHap);
        void                            LooOptimalStructure         (int loo,HaplotypeSet &rHap, HaplotypeSet &HapLoo);
        double                          CalculateLikelihood         (HaplotypeSet &rHap,MarkovModel &MM);
        void                            ImputeTraverse              (HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                                                                    MarkovModel &MM,int group, vector<float> &recomProb,
                                                                    vector<float> &PrevRightFoldedProb,
                                                                    vector<float> &CurrentRightProb,
                                                                    vector<float> &CurrentNoRecoRightProb,
                                                                    HaplotypeSet &rHapMain);
        void                            ImputeTraverse              (HaplotypeSet &rHap,
                                                                    int hapID,
                                                                    MarkovModel &MM,int group, vector<float> &recomProb,
                                                                    vector<float> &PrevRightFoldedProb,
                                                                    vector<float> &CurrentRightProb,
                                                                    vector<float> &CurrentNoRecoRightProb);
        void                            LeftTraverse                (HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                                                                    MarkovModel &MM,int group, vector<float> &recomProb);
        void                            EMTraverse                  (HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                                                                    MarkovModel &MM,int group, vector<float> &recomProb,
                                                                    vector<float> &PrevRightFoldedProb,
                                                                    vector<float> &CurrentRightProb,
                                                                    vector<float> &CurrentNoRecoRightProb,
                                                                    HaplotypeSet &rHapMain);
        void                            ConditionJunctionProb       (HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                                                                    double e, double freq, AlleleType observed, double backgroundError,
                                                                    ReducedHaplotypeInfo &Info);
        void                            FlushPartialVcf             (HaplotypeSet &rHap,HaplotypeSet &tHap,HaplotypeSet &PartialDosage, string &filename,int &Index);
        void                            MergeFinalVcf               (HaplotypeSet &rHap,HaplotypeSet &tHap,ImputationStatistics &stats,int MaxIndex);
        void                            MergeFinalVcfAllVariants    (HaplotypeSet &rHap,HaplotypeSet &tHap,ImputationStatistics &stats,int MaxIndex);
        void                            PrintDosageData             (int ThisSampleId, vector<float> &ThisDosage1,vector<float> &ThisDosage2);
        void                            PrintHaplotypeData          (int ThisHapId, int ThisSampleId, vector<float> *ThisimputedHap);
        void                            PrintInfoFile               (HaplotypeSet &rHap,HaplotypeSet &tHap, ImputationStatistics &stats);
        void                            OutputFilesInitialize       (HaplotypeSet &rHap,HaplotypeSet &tHapOrig);
        void                            InitializeOutputFiles       (HaplotypeSet &tarInitializer, int maxSample,
                                                                     int maxRefVar, int maxTarVar);
        void                            FreeMemory                  ();


                                        Imputation                  (AllVariable *MyAllVariable, IFILE Dosages, IFILE Hapdose,IFILE Haps,IFILE Vcfdosepartial,
                                                                    IFILE Info, ImputationStatistics &Stats)
                                                                    {
                                                                        MyAllVariables=MyAllVariable;
                                                                        dosages=Dosages;
                                                                        hapdose=Hapdose;
                                                                        haps=Haps;
                                                                        vcfdosepartial=Vcfdosepartial;
                                                                        info=Info;
                                                                        stats=&Stats;
                                                                        TimeToWrite=0;
                                                                        TimeToImpute=0;
                                                                        TimeToCompress=0;
                                                                    };


};




#endif // IMPUTATION_H_INCLUDED
