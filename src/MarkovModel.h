#ifndef __MARKOVMODEL_H__
#define __MARKOVMODEL_H__


#include "MarkovParameters.h"
#include "HaplotypeSet.h"
//#include "DosageData.h"
#include "StringBasics.h"
#include "MathVector.h"
#include "Random.h"

#include "Unique.h"

#include "ImputationStatistics.h"



class MarkovModel : public MarkovParameters
{
    public:

        HaplotypeSet *rHap, *tHap, *rHapFull, *tHapFull;
        int NoRefMarkers,NoChipSites;

        int ThisHapId;
        bool LowMemory;
        double backgroundError;
        int refCount,tarCount,noReducedStatesCurrent;


        int NoPrecisionJumps;
        int CurrentTypedSite;
        bool CurrentObsMissing, CurrentObs;
        int BeforeLastUntypedSite;
        double SummedProb;
        double tempMaxVal;
        int NoBestMatchHaps, NoBestMatchFullRefHaps, noNewReducedStates;



        // Left Probabilities

        vector<vector<vector<float> > > leftProb;
        vector<vector<float> > ThisBlockLeftNoRecoProb ;
        vector<float> CurrentLeftNoRecoProb;


        // Junction Probabilities


        vector<vector<float> > junctionLeftProb;
        vector<float> PrevjunctionRightProb;

        // Best Match Haplotypes

        vector<int> BestMatchHaps, BestMatchFullRefHaps, FinalBestMatchfHaps;

        // Final Best Batch Folded Probabilities for Dosage

        vector<float> pREF,pALT;
        vector<double> FoldedProbValue;
        vector<double> probHapFullAverage;
        vector<float>    *DosageHap, *LooDosageHap;


        // Other Variables for Support

        vector<float> Constants;
        vector<bool> PrecisionJump;
        vector<float> tempRightProb;
        vector<float> probHap;






        double size()
        {
            double S=0;
            S+=PrecisionJump.size() * sizeof(bool);

            S+=BestMatchHaps.size() * sizeof(int);
            S+=BestMatchFullRefHaps.size() * sizeof(int);
            S+=FinalBestMatchfHaps.size() * sizeof(int);

            S+=probHap.size() * sizeof(float);
            S+=Constants.size() * sizeof(float);
            S+=tempRightProb.size() * sizeof(float);
            S+=pREF.size() * sizeof(float);
            S+=pALT.size() * sizeof(float);
            S+=PrevjunctionRightProb.size() * sizeof(float);
            S+=CurrentLeftNoRecoProb.size() * sizeof(float);

            S+=FoldedProbValue.size() * sizeof(double);
            S+=probHapFullAverage.size() * sizeof(double);


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

            return (S);
        };



        void AssignPanels        (HaplotypeSet &refFullHap,HaplotypeSet &refChipHap,
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



AllVariable *MyAllVariables;




void ImputeSites( int hapID,int group,
                                   vector<float> &PrevRightFoldedProb,
                                   vector<float> &CurrentRightProb,
                                   vector<float> &CurrentNoRecoRightProb);
void ImputeSitesAllProb( int hapID,int group,
                                   vector<float> &PrevRightFoldedProb,
                                   vector<float> &CurrentRightProb,
                                   vector<float> &CurrentNoRecoRightProb);


void ImputeChunk(int group, int &hapID, int &position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb);

void FindPosteriorProbWithThreshold( int group,int position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb);


void unfoldProbabilitiesWithThreshold(int bridgeIndex,
                                         vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                         vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                         vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb);

void unfoldProbabilitiesAllProb(int bridgeIndex,
                                         vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                         vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                         vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb);



void FoldBackProbabilitiesWithThreshold(ReducedHaplotypeInfo &Info, ReducedHaplotypeInfo &TarInfo);
void FoldBackProbabilitiesAllProb(ReducedHaplotypeInfo &Info, ReducedHaplotypeInfo &TarInfo);



void ImputeRemainingSitesbyBlock(int bridgeIndex, int TypedMarkerId, int StartPos, int EndPos);

void ImputeIntermediateRegionFull(ReducedHaplotypeInfo &Info, int TypedMarkerId, int StartPos, int EndPos);
void CreatePRefPAlt(ReducedHaplotypeInfo &Info, int StartPos, int EndPos);
void CreateDosages(int TypedMarkerId, int StartPos, int EndPos);
void CreateLooDosage(int &StartPos, int &TypedMarkerId,
                                    bool &observedMiss, bool &observed);





        int             Myfind                          (vector<int> &MyVector, int Value);


void ImputeIntermediateRegion(ReducedHaplotypeInfo &Info, int position, int TempEndPos);

        void            initializeMatricesNew           ();

        bool            Transpose                       (vector<float> &from,vector<float> &to,  vector<float> &noRecomProb, double reco,vector<int> &uniqueCardinality);
        void            Condition                       (int markerPos,vector<float> &Prob, vector<float> &noRecomProb,
                                                        bool observed,double e,double freq,ReducedHaplotypeInfo &Info);
        void            WalkLeft                        (HaplotypeSet &tHap, int &hapID,
                                                        int group, ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);
        void            ReCreateBothLeftProb            (HaplotypeSet &tHap, int &hapID,
                                                        int group, ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);
        void            ReCreateLeftNoRecoProb          (HaplotypeSet &tHap, int &hapID,
                                                        int group, ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);

        void            Impute                          (HaplotypeSet &tHap,int hapID,int group,
                                                        vector<float> &PrevRightFoldedProb,
                                                        vector<float> &CurrentRightProb, vector<float> &CurrentNoRecoRightProb,
                                                        ReducedHaplotypeInfo &Info,vector<double> &alleleFreq);



        void            ImputeAllSites                          (int group, int position, bool observed, bool observedMiss,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb);




        void ImputeRemainingSitesbyBlock( int StartPos, int &EndPos);
        void ImputeRemainingSites( int StartPos, int &EndPos);

        void unfoldProbabilitiesAtThisSite(int bridgeIndex,
                                         vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                            vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                         vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb);




        void ImputeAllSites( int hapID,int group,
                                   vector<float> &PrevRightFoldedProb,
                                   vector<float> &CurrentRightProb,
                                   vector<float> &CurrentNoRecoRightProb);


        void            Impute                          (int position, bool observed, bool observedMiss,
                                                        vector<float> &leftProb,vector<float> &rightProb,
                                                        vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                                        vector<float> &leftEndProb,vector<float> &rightEndProb,
                                                        vector<float> &Constants,ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);
        void            initializeMatrices              (HaplotypeSet &rHap,HaplotypeSet &tHap);
        void            ReinitializeMatrices            ();


        double          CountErrors                     (vector<float> &probHap,
                                                        int position, bool observed, double e,double freq, ReducedHaplotypeInfo &Info);
        double          CountRecombinants               (vector<float> &from, vector<float> &to,
                                                        vector<float> &probHap,
                                                        double r,bool PrecisionMultiply);
        void            foldProbabilities               (vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference);
        void            unfoldProbabilities             (int bridgeIndex,vector<float> &recomProb, vector<float> &noRecomProb,
                                                        vector<float> &PrevFoldedProb,int direction,vector<ReducedHaplotypeInfo> &StructureInfo,int noReference);

        void            CountExpected                   (HaplotypeSet &tHap,int hapID,int group,
                                                        vector<float> &PrevRightFoldedProb,
                                                        vector<float> &CurrentRightProb, vector<float> &CurrentNoRecoRightProb,
                                                        ReducedHaplotypeInfo &Info,vector<double> &alleleFreq);

        void            CreatePosteriorProb             (vector<float> &leftProb,vector<float> &rightProb,
                                                        vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                                        vector<float> &leftEndProb,vector<float> &rightEndProb,
                                                        vector<float> &Constants,vector<float> &probHap,ReducedHaplotypeInfo &Info);

        void            CheckSize                       (HaplotypeSet & rHap,HaplotypeSet &tHap);



    MarkovModel()
    {

    };                         // constructor; initialize the list to
                                       // be empty
    ~MarkovModel()
    {

    };                            // destructor

    MarkovModel(const MarkovModel &L)
    {
    };
//    };             // copy constructor
//    MarkovModel & operator=(const MarkovModel &L); // assignment



    };

#endif
