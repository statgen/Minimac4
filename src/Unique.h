#ifndef UNIQUE_H_INCLUDED
#define UNIQUE_H_INCLUDED
#include<cmath>
#include<fstream>
#include "StringBasics.h"
#include<vector>

#include <unordered_set>
#include "assert.h"
using namespace std;

class variant
{
public:

    string name;
    int bp;
    string chr;
    string rsid;
//    char refAllele,altAllele;
    string refAlleleString,altAlleleString;

    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
    void assignValues(string &id,string &CHR,int BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
     void assignRefAlt(string &refe,string &alt)
    {
        refAlleleString=refe;
        altAlleleString=alt;
    };

    ~variant(){};


    double size()
    {
        double S=0;
        S+=name.size();
        S+=chr.size();
        S+=rsid.size();
        S+=refAlleleString.size();
        S+=altAlleleString.size();
        return (S+2);
    };

    variant(const variant &obj)
    {
        name=obj.name;
        rsid=obj.rsid;
        bp=obj.bp;
        chr=obj.chr;
        refAlleleString=obj.refAlleleString;
        altAlleleString=obj.altAlleleString;
//        refAllele=obj.refAllele;
//        altAllele=obj.altAllele;
    }
};



class ReducedHaplotypeInfo
{
    public:

        // Basic Compulsory Parameters

        int startIndex,endIndex;
        vector<int> uniqueCardinality; // has number of representatives for each unique representative
        vector<float> InvuniqueCardinality; // has number of representatives for each unique representative
        vector<int> uniqueIndexMap; // maps to corresponding item in the uniqueRep... required to map Constants and unfold fold probs
//        vector<vector<bool> > uniqueHaps;
        vector<vector<bool> > TransposedUniqueHaps;

        int BlockSize, RepSize;



//        bool returnHapAtPosition(int i,int position)
//        {
//            //assert((position-startIndex)>=0);
//            return uniqueHaps[i][position-startIndex];
//        }


        // Special members for RefAtGWAS Panel

        vector< vector<int> > uniqueIndexReverseMaps;
        vector< unordered_set<int> > uniqueIndexOriginalReducedMaps;

        int size()
        {
            int S=0;
            S+=uniqueCardinality.size() * sizeof(int);
            S+=InvuniqueCardinality.size() * sizeof(float);
            S+=uniqueIndexMap.size() * sizeof(int);
            for(int i=0;i<(int)TransposedUniqueHaps.size();i++)
            {
                S+=TransposedUniqueHaps[i].size() * sizeof(bool);
            }
            return (S+8);
        };

        ReducedHaplotypeInfo()
        {

            startIndex=0;
            endIndex=0;
            BlockSize=0;
            RepSize=0;
            uniqueCardinality.clear();
            InvuniqueCardinality.clear();
            uniqueIndexMap.clear();
            TransposedUniqueHaps.clear();
//            uniqueHaps.clear();
        }

        ~ReducedHaplotypeInfo()
        {
        }

        ReducedHaplotypeInfo(const ReducedHaplotypeInfo &obj)
        {
            startIndex=obj.startIndex;
            endIndex=obj.endIndex;
            BlockSize=obj.BlockSize;
            RepSize=obj.RepSize;
            uniqueCardinality=obj.uniqueCardinality;
            InvuniqueCardinality=obj.InvuniqueCardinality;
            uniqueIndexMap=obj.uniqueIndexMap;
            TransposedUniqueHaps=obj.TransposedUniqueHaps;
//            uniqueHaps=obj.uniqueHaps;
        }

};

class ReducedHaplotypeInfoSummary
{
    public:

        // Basic Summary Parameters
        int startIndex,endIndex;
        int BlockSize, RepSize;

        ReducedHaplotypeInfoSummary()
        {

            startIndex=0;
            endIndex=0;
            BlockSize=0;
            RepSize=0;
        }
        ~ReducedHaplotypeInfoSummary()
        {
        }

};


class CompressedHaplotype
{
    public:

        // Basic Compulsory Parameters

        int Length;
        vector<int> ReducedInfoMapper;
        vector<int> ReducedInfoVariantMapper;
        vector<ReducedHaplotypeInfo> *RhapInfo;
        ReducedHaplotypeInfo *CurrentVariantInfo;
        vector<bool> *CurrentVariantHaps;

        CompressedHaplotype()
        {
            ReducedInfoMapper.clear();
            ReducedInfoVariantMapper.clear();


            Length=0;
        }
        ~CompressedHaplotype()
        {
            RhapInfo=NULL;
        }


        void Size(vector<ReducedHaplotypeInfo> &rhapInfo, int &bufferSize)
        {
            RhapInfo=&rhapInfo;
            ReducedInfoMapper.resize(bufferSize);
            ReducedInfoVariantMapper.resize(bufferSize);
        }


        void Clear()
        {
            Length=0;
        }

        void RetrieveVariant(int position)
        {
            CurrentVariantInfo = & (*RhapInfo)[ReducedInfoMapper[position]];
            CurrentVariantHaps = &(CurrentVariantInfo->TransposedUniqueHaps[ReducedInfoVariantMapper[position]]);
        }
        char GetVal(int &HapId)
        {
            int &Index = CurrentVariantInfo->uniqueIndexMap[HapId];
            return ( (*CurrentVariantHaps)[Index]?'1':'0');
        }

        char GetVal(int &HapId, int position)
        {
            ReducedHaplotypeInfo &tempInfo=(*RhapInfo)[ReducedInfoMapper[position]];
            int &Index = tempInfo.uniqueIndexMap[HapId];

            return (tempInfo.TransposedUniqueHaps[ReducedInfoVariantMapper[position]][Index]?'1':'0');
        }


        void Push(int &Val1,int &Val2)
        {
            ReducedInfoMapper[Length]=Val1;
            ReducedInfoVariantMapper[Length++]=Val2;
        }

};

class CrossReducedHaplotypeInfo
{
    public:
//
//        // Basic Compulsory Parameters
//
//        int startIndex,endIndex;
//        vector<int> uniqueCardinality; // has number of representatives for each unique representative
//        vector<float> InvuniqueCardinality; // has number of representatives for each unique representative
//        vector<int> uniqueIndexMap; // maps to corresponding item in the uniqueRep... required to map Constants and unfold fold probs
//        vector<vector<bool> > uniqueHaps;
//
//        int BlockSize, RepSize;
//
//
//
//        bool returnHapAtPosition(int i,int position)
//        {
//            //assert((position-startIndex)>=0);
//            return uniqueHaps[i][position-startIndex];
//        }
//
//
//        // Special members for RefAtGWAS Panel

        vector< vector<int> > uniqueIndexMaps;

//        vector< vector<int> > uniqueIndexMaps;

        vector< vector<int> > InvuniqueCardinality;


//
//        vector< unordered_set<int> > uniqueIndexOriginalReducedMaps;
//


        CrossReducedHaplotypeInfo()
        {
            uniqueIndexMaps.clear();
            InvuniqueCardinality.clear();
        }

        CrossReducedHaplotypeInfo(int N)
        {
            uniqueIndexMaps.resize(N);
            InvuniqueCardinality.resize(N);
        }

};

class findUnique
{

    public:

        int N; // No. of Samples
        int max_block;  // Maximum Block Length to consider
        int M; // No. of Markers
        vector<vector<int> > uniqueCount; // matrix to store unique number of haplotypes between positions
        vector<int> minCost;
        vector<int> minAllocation;
        vector<int> optimalAllocation;
        int transFactor, cisFactor;

         findUnique()
         {


         }
         findUnique(vector<vector<char> > &hap)
        {
            N=hap.size();
            M=hap[0].size();
        }

    void UpdateDeltaMatrix(CompressedHaplotype &haplotypes, vector<int> & index,
          vector<int> & firstDifference, int length, int blockSize,
          vector<int> & oldIndex,  vector<int> & previousPredecessor,  vector<int> & previousDifference);

          void UpdateDeltaMatrix(vector<String> & haplotypes, vector<int> & index,
          vector<int> & firstDifference, int length, int blockSize,
          vector<int> & oldIndex,  vector<int> & previousPredecessor,  vector<int> & previousDifference);
          void AnalyzeBlocks(vector<int> & index, vector<int> & firstDifference, int length, int blockSize, vector<int> & cost,
         vector<int> & bestSlice, vector<int> & bestComplexity, vector<vector<int> > &bestIndex);

         double FlushBlocks(vector<ReducedHaplotypeInfo> &HapInfo, int LastflushPos, CompressedHaplotype & haplotypes, vector<int> & cost,
                   vector<int> & bestComplexity, vector<int> & bestSlice, vector<vector<int> > &bestIndex);

            double FlushBlocks(vector<ReducedHaplotypeInfo> &HapInfo, int LastflushPos,vector<String> & haplotypes, vector<int> & cost,
                   vector<int> & bestComplexity, vector<int> & bestSlice, vector<vector<int> > &bestIndex);







        void updateCoeffs(int trans,int cis)
        {
            transFactor = trans;
            cisFactor = cis;

        }


};






#endif // UNIQUE_H_INCLUDED
