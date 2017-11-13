#ifndef __MARKOVMODEL_H__
#define __MARKOVMODEL_H__


#include "MarkovParameters.h"
#include "HaplotypeSet.h"
#include "StringBasics.h"
#include "MathVector.h"
#include "Random.h"
#include "Unique.h"
#include "ImputationStatistics.h"
#include "assert.h"


class MarkovModel : public MarkovParameters
{
    public:

        // Basic Variables
        HaplotypeSet *rHap, *tHap, *rHapFull, *tHapFull;
        int NoRefMarkers,NoChipSites;
        int refCount,tarCount,noReducedStatesCurrent;
        int ThisHapId;
        AllVariable *MyAllVariables;


        // Left Probabilities
        vector<vector<vector<float> > > leftProb;
        vector<vector<float> > ThisBlockLeftNoRecoProb ;
        vector<float> CurrentLeftNoRecoProb;


        // Right Probabilities
        vector<vector<float> > ThisBlockRightProb;
        vector<vector<float> > ThisBlockRightNoRecoProb ;
        vector<float> PrevRightFoldedProb;


        // Final Probability Matrix
        vector<vector<double> > probHapMatrix;
        vector<double> SumOfProb;
        vector<double> probHapMinimac3;


         // Junction Probabilities
        vector<vector<float> > junctionLeftProb;
        vector<float> PrevjunctionRightProb;


        // UnfoldingWithThreshold Variables
        vector<double> LeftAdj_Rec;
        vector<double> LeftAdj_NoRrec;
        vector<double> RightAdj_Rec;
        vector<double> RightAdj_NoRec;


        // Variables for Successive Differences
        vector<double> FirstFoldedValue;
        double FirstDiffValue;
        int KeepMovingLeft;
        vector<int> UnfoldTheseSites;
        int NoSitesToUnfold;


        // Best Match Haplotypes
        vector<int> BestMatchHaps, BestMatchFullRefHaps, FinalBestMatchfHaps;
        int NoBestMatchHaps, NoBestMatchFullRefHaps, noNewReducedStates;
        vector<int> FinalBestMatchfHapsIndicator;

        // Final Best Batch Folded Probabilities for Dosage
        vector<double> FoldedProbValue;
        vector<double> probHapFullAverage;
        vector<float>    *DosageHap, *LooDosageHap;


        // Precision Jump Variables
        vector<bool> LeftPrecisionJump, RightPrecisionJump;
        int NoPrecisionJumps;
        float JumpThreshHold = 1e-20;
        float JumpFix = 1e10;


        // Current Variables for Support
        AlleleType CurrentObsMissing, CurrentObs;
        int CurrentTypedSite;


        // Sum Variables for Support
        double PrevTotalSum,InvPrevTotalSum;


        // Other Variables for Support
        vector<float> Constants;
        int MidPoint;
        bool LowMemory;
        double backgroundError;
        double tempMaxVal;
        int MostProbableTemplate;
        double MostProbableTemplateVal;


        void CeateProbSum                   (int bridgeIndex, int noReference);
        void RightCondition                 (int markerPos,vector<float> &FromProb,vector<float> &ToProb,
                                            vector<float> &FromNoRecomProb, vector<float> &ToNoRecomProb,
                                            AlleleType observed,double e,double freq,ReducedHaplotypeInfo &Info);
        bool RightTranspose                 (vector<float> &fromTo, vector<float> &noRecomProb,
                                            double reco,vector<int> &uniqueCardinality);
        bool LeftTranspose                  (vector<float> &from,vector<float> &to,  vector<float> &noRecomProb,
                                            double reco,vector<int> &uniqueCardinality);
        void LeftCondition                  (int markerPos,vector<float> &Prob, vector<float> &noRecomProb,
                                            AlleleType observed,double e,double freq,ReducedHaplotypeInfo &Info);
        void ImputeSitesMinimac3            (int hapID,int group);
        void ImputeChunkMinimac3            ( int group, int hapID, int position,
                                            vector<float> &Leftprob,vector<float> &rightProb,
                                            vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                            vector<float> &leftEndProb,vector<float> &rightEndProb);
        void ImputeSites                    (int hapID,int group);
        void WalkLeft                       (int &hapID,  int group);
        void WalkLeftMinimac3               (int &hapID,  int group);
        void ImputeChunk                    (int group, int &hapID, int &position,
                                             vector<float> &Leftprob,vector<float> &rightProb,
                                             vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                             vector<float> &leftEndProb,vector<float> &rightEndProb);
        void FindPosteriorProbWithThreshold (int group,int position,
                                             vector<float> &Leftprob,vector<float> &rightProb,
                                             vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                             vector<float> &leftEndProb,vector<float> &rightEndProb);
        void unfoldProbabilitiesWithThreshold(int bridgeIndex,
                                             vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                             vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                             vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb);
        void unfoldProbabilitiesAllProb     (int bridgeIndex,
                                             vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                             vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                             vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb);
        void FoldBackProbabilitiesWithThreshold(ReducedHaplotypeInfo &Info);
        void FoldBackProbabilitiesAllProb   (ReducedHaplotypeInfo &Info, ReducedHaplotypeInfo &TarInfo);
        void ImputeRemainingSitesbyBlock    (int bridgeIndex, int TypedMarkerId, int StartPos, int EndPos);
        void ImputeIntermediateRegionFull   (ReducedHaplotypeInfo &Info, int TypedMarkerId, int StartPos, int EndPos);
        void CreatePRefPAlt                 (ReducedHaplotypeInfo &Info, int StartPos, int EndPos);
        void CreateDosages                  (int TypedMarkerId, int StartPos, int EndPos);
        void CreateLooDosage                (int &StartPos, int &TypedMarkerId,
                                            bool &observedMiss, bool &observed);
        
        void initializeMatricesMinimac3     ();
        void initializeMatrices             ();
        void ReinitializeMatrices           ();
        void foldProbabilities              (vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,
                                            int direction,int noReference);
        void unfoldProbabilities            (int bridgeIndex,vector<float> &recomProb,
                                            vector<float> &noRecomProb,
                                            vector<float> &PrevFoldedProb,int direction,
                                            vector<ReducedHaplotypeInfo> &StructureInfo,int noReference);
        void ReCreateLeftNoRecoProb         (HaplotypeSet &tHap, int &hapID,
                                             int group, ReducedHaplotypeInfo &Info,
                                             vector<double> &alleleFreq);
        void ReCreateLeftNoRecoProbMinimac3 (HaplotypeSet &tHap, int &hapID,
                                             int group, ReducedHaplotypeInfo &Info,
                                             vector<double> &alleleFreq);
        int Myfind                          (vector<int> &MyVector, int Value);
        MarkovModel                         ()
                                            {
                                            };
        ~MarkovModel                        ()
                                            {
                                            };
        MarkovModel                         (const MarkovModel &L)
                                            {
                                            };
        void AssignPanels                   (HaplotypeSet &refFullHap,HaplotypeSet &refChipHap,
                                             HaplotypeSet &tarFullHap,HaplotypeSet &tarChipHap,
                                             AllVariable *MyAllVariable)
                                            {
                                                MyAllVariables=MyAllVariable;
                                                rHapFull=&refFullHap;
                                                rHap=&refChipHap;
                                                tHapFull=&tarFullHap;
                                                tHap=&tarChipHap;
                                                refCount=rHapFull->numHaplotypes;
                                                tarCount=tHapFull->numHaplotypes;
                                                NoRefMarkers=rHapFull->numMarkers;
                                                NoChipSites=rHap->numMarkers;
                                                backgroundError = 1e-5;
                                            };
        void AssignPanels                   (HaplotypeSet &refFullHap,
                                             HaplotypeSet &tarFullHap,
                                             AllVariable *MyAllVariable)
                                            {
        MyAllVariables=MyAllVariable;

        rHap=&refFullHap;
        tHap=&tarFullHap;

        refCount=rHap->numHaplotypes;
        tarCount=tHap->numHaplotypes;
        NoRefMarkers=rHap->numMarkers;
        NoChipSites=tHap->numMarkers;
        backgroundError = 1e-5;
    };
        void AssignPanels                   (HaplotypeSet &refFullHap,
                                             AllVariable *MyAllVariable)
                                            {
                                                MyAllVariables=MyAllVariable;
                                                rHap=&refFullHap;
                                                refCount=rHap->numHaplotypes;
                                                NoRefMarkers=rHap->numMarkers;
                                                NoChipSites=rHap->numMarkers;
                                                backgroundError = 1e-5;
                                            };
        void CountExpected                  (int hapID,int group);
        void CreatePosteriorProb            (vector<float> &Leftprob,vector<float> &rightProb,
                                            vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                            vector<float> &leftEndProb,vector<float> &rightEndProb,
                                             ReducedHaplotypeInfo &Info);
        double CountRecombinants            (vector<float> &from, vector<float> &to,
                                            double r,bool PrecisionMultiply);
        double CountErrors                  (int markerPos,
                                            AlleleType observed,
                                            double e,double freq,
                                            ReducedHaplotypeInfo &Info);







    double size()
    {
        double S=0;

        S+=(2*LeftPrecisionJump.size() * sizeof(bool));

        S+=BestMatchHaps.size() * sizeof(int);
        S+=BestMatchFullRefHaps.size() * sizeof(int);
        S+=FinalBestMatchfHaps.size() * sizeof(int);
        S+=FinalBestMatchfHapsIndicator.size() * sizeof(int);

        S+=FirstFoldedValue.size() * sizeof(double);
        S+=UnfoldTheseSites.size() * sizeof(int);
        S+=PrevRightFoldedProb.size() * sizeof(float);

        S+=Constants.size() * sizeof(float);
        S+=PrevjunctionRightProb.size() * sizeof(float);
        S+=CurrentLeftNoRecoProb.size() * sizeof(float);

        S+=FoldedProbValue.size() * sizeof(double);
        S+=probHapFullAverage.size() * sizeof(double);
        S+=LeftAdj_Rec.size() * sizeof(double);
        S+=LeftAdj_NoRrec.size() * sizeof(double);
        S+=RightAdj_Rec.size() * sizeof(double);
        S+=RightAdj_NoRec.size() * sizeof(double);
        S+=SumOfProb.size() * sizeof(double);

        for(int i=0;i<(int)leftProb.size();i++)
        {
            for(int j=0;j<(int)leftProb[i].size();j++)
            {
                S+=leftProb[i][j].size() * sizeof(float);
            }
        }
        for(int i=0;i<(int)ThisBlockLeftNoRecoProb.size();i++)
        {
            S+=ThisBlockLeftNoRecoProb[i].size() * sizeof(float);
        }
        for(int i=0;i<(int)junctionLeftProb.size();i++)
        {
            S+=junctionLeftProb[i].size() * sizeof(float);
        }
        for(int i=0;i<(int)ThisBlockRightProb.size();i++)
        {
            S+=ThisBlockRightProb[i].size() * sizeof(float);
        }
         for(int i=0;i<(int)ThisBlockRightNoRecoProb.size();i++)
        {
            S+=ThisBlockRightNoRecoProb[i].size() * sizeof(float);
        }
         for(int i=0;i<(int)probHapMatrix.size();i++)
        {
            S+=probHapMatrix[i].size() * sizeof(double);
        }


        return (S);
    };





    };

#endif
