#ifndef DOSAGEDATA_H_INCLUDED
#define DOSAGEDATA_H_INCLUDED

#include "StringBasics.h"
#include "MyVariables.h"
#include "VcfFileReader.h"
#include "Unique.h"
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include "HaplotypeSet.h"
#include "MyVariables.h"
#include "MarkovModel.h"
using namespace std;



class DosageData
{

	public:

	    // Basic Variables

        int ChunkNo, FirstHapId;
		int         ActualnumHaplotypes, BuffernumHaplotypes;
		int         ActualnumSamples, BuffernumSamples;
		int         numMarkers, noGWASSites;
        HaplotypeSet *tHapFull;
        HaplotypeSet *rHapFull;
		vector<string> individualName;


        // Indexing Samples Variable
        int NoSamplesIndexed;
        vector<int> SampleIndex;
        vector<int> InvertSampleIndex;
        vector<int> BufferSampleNoHaplotypes;


		// Main Data Variable
		vector<vector<float> > hapDosage;
		vector<vector<float> > LoohapDosage;


		// Formatting Variables
        AllVariable *MyAllVariables;

        // Time Variables
        int TimeToWrite;

        // Printing Text Variables
        char *PrintStringPointer;
        char *PrintEmpStringPointer;
        int PrintStringLength;
        int PrintEmpStringLength;

        std::vector<float> savDoseBuf;


        double size()
        {
            double S=0;
            S+=BufferSampleNoHaplotypes.size() * sizeof(int);
            S+=InvertSampleIndex.size() * sizeof(int);
            S+=SampleIndex.size() * sizeof(int);

            for(int i=0;i<(int)hapDosage.size();i++)
            {
                S+=hapDosage[i].size() * sizeof(float);
                S+=LoohapDosage[i].size() * sizeof(float);
            }
            for(int i=0;i<(int)individualName.size();i++)
            {
                S+=individualName[i].size();
            }
            return (S);
        };



void PrintDiploidLooDosage(float &x, float &y, AlleleType a, AlleleType b);
void PrintHaploidLooDosage(float &x, AlleleType a);

        void InitializePartialDosageData(HaplotypeSet &tarInitializer, int MaxNoSamples,
                                             int MaxNoRefVariants, int MaxNoTarVariants,
                                               AllVariable *MyAllVariable);
        void ReParameterizePartialDosageData(int No, HaplotypeSet &refFullHap, HaplotypeSet &tarFullHap);
        void UpdatePartialDosageData(int NewMaxVal, int NewFirstHapId);

        pair <int, int> IndexSample(int SampleId);
        void BindSampleMModel(MarkovModel &MM, int Index, int HapNo);

        void FlushPartialVcf(int NovcfParts);
        void PrintDiploidDosage(float &x, float &y);
        void PrintHaploidDosage(float &x);
        void PrintDosageForVcfOutputForID(int MarkerIndex);
        void PrintDosageForVcfOutputForIDFast(int MarkerIndex);
        void PrintGWASOnlyForVcfOutputForID(int MarkerIndex);

        DosageData () {};



};


#endif // DOSAGEDATA_H_INCLUDED
