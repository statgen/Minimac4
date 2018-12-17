#ifndef ANALYSIS_H_INCLUDED
#define ANALYSIS_H_INCLUDED


#if defined(_OPENMP)
#include <omp.h>
#endif
#include "MyVariables.h"
#include "Unique.h"
#include "Imputation.h"
#include "HaplotypeSet.h"
//#include "StringBasics.h"
//#include "VcfFileReader.h"

#include <savvy/reader.hpp>

//#include <fstream>
//#include <string>
//#include <sstream>
//#include <algorithm>
#include <memory>
#include <unordered_map>

using namespace std;


class Analysis
{


    public:

	    // Basic Variables

        HaplotypeSet targetPanel,referencePanel;
        IFILE RefFileStream;
        VcfFileReader TarFileStream;
        VcfRecord record;
        vector< vector<double> > GeneticMapData;

        // Swapping Variables for actual chunk IMPUTATION
        HaplotypeSet CurrentTarPanel,CurrentRefPanel;
        HaplotypeSet CurrentTarPanelChipOnly,CurrentRefPanelChipOnly;

        int PrevChunkFilledTillRef, PrevChunkStartFromRef;
        int PrevChunkFilledTillTar, PrevChunkStartFromTar;
        int PrevChunkFilledTillTarOnly, PrevChunkStartFromTarOnly;


    void InitializeTargetChipOnlyChunkData(HaplotypeSet &ThisTarPanel);

    void InitializeRefChipOnlyChunkData(HaplotypeSet &ThisRefPanel);

    void readm3vcfFileChunk(int ChunkNo, HaplotypeSet &ThisRefPanel);


    void readVcfFileChunk(int ChunkNo, HaplotypeSet &ThisTargetPanel);



        // Target PANEL Indexing Variables
        int overlapImportIndex, importReadIndex, FileReadIndex;
        int GWASOnlySkipIndexpoint , MainTypedOnlyImportIndex;


        // Output File Handle Streams
        IFILE dosages, hapdose, haps, vcfdosepartial, vcfLoodosepartial, info;

        std::unique_ptr<savvy::sav::writer> savOut;

        ImputationStatistics stats;


        // Chunk Information

        int noChunks;
	    vector< vector<int> > MyChunks,MyChunksInfoNumber;
	    vector< vector<int> > MyRefVariantNumber, MyTargetVariantNumber, MyTypdedOnlyVariantNumber;
	    vector< float > MyRatio;

        int MaxInfoVectorSize, MaxGwasMarkerSize, MaxRefMarkerSize, MaxTypedOnlyMarkerSize;



        // Temporary Variables to be deleted later for checking purposes.

        int OverCount, TypOnlyCount, RefCOUNT;


        // Time Variables

        int TimeToCompress, TimeToRead, TimeToImpute, TimeToWrite;

        // Memory Variables
        double RefMem, TarMem, ComRefMem, DosageMem, ProbMem;


        // Printing Text Variables
        char *VcfPrintStringPointer, *LooVcfPrintStringPointer, *InfoPrintStringPointer, *DosePrintStringPointer;
        int VcfPrintStringLength, InfoPrintStringLength, DosePrintStringLength;




        vector<ReducedHaplotypeInfo> InfoforChunk;


		vector< vector<int> > BPListMapper;

	    HaplotypeSet *CurrentChunk, *NextChunk;




        HaplotypeSet MainReferenceFrame;

        OutputFormatVariable *MyOutFormat;
        ModelVariable *MyModelVariables;
        HaplotypeDataVariables *MyHapDataVariables;
        AllVariable *MyAllVariables;



void PrintInfoFile(int ChunkNo);

Analysis()
{
    TimeToCompress=0;
    TimeToRead=0;
    TimeToImpute=0;
    TimeToWrite=0;
    RefMem =0.0;
    TarMem =0.0;
     ComRefMem =0.0;
     DosageMem =0.0;
    ProbMem = 0.0;


}

void AppendtoMainVcf(int ChunkNo, int MaxIndex);

void AppendtoMainLooVcfFaster(int ChunkNo, int MaxIndex);

void AppendtoMainVcfFaster(int ChunkNo, int MaxIndex);

void AppendtoMainSAVFaster(int ChunkNo, int MaxIndex);


void MemDisplay();

 void GetNumChunks();

                String                  RunEstimation(String &Reffilename, String &Recomfilename, String &Errorfilename);
                String                  AnalyzeExperiment(String &Reffilename, String &Tarfilename, String &Recomfilename,String &Errorfilename, AllVariable& MyAllVariable);
                String                  CheckValidity(String &Reffilename, String &Tarfilename, String &Recomfilename, String &Errorfilename);
                String                  RunAnalysis(String &Reffilename, String &Tarfilename, String &Recomfilename, String &Errorfilename);
                void                    readm3vcfFileChunk(HaplotypeSet &ThisRefPanel, HaplotypeSet &NextRefPanel,
                                        int StartBlock, int EndBlock, int StartNextBlock);
                void                    readVcfFileChunk(HaplotypeSet &ThisTargetPanel, HaplotypeSet &NextTargetPanel,
                                        int StartPos, int EndPos, int StartNextPos, int ChunkNo);
                bool                    OpenStreamOutputFiles();
                void                    CloseStreamOutputFiles();
                void                    InitializeRefFileStream(String &Reffilename);
                void                    InitializeTargetFileStream(String &Tarfilename);
                void                    InitializeRefChunkData(HaplotypeSet &ThisRefPanel);
                void                    InitializeTargetChunkData(HaplotypeSet &ThisTarPanel);

                bool                  CreateRecombinationMap();
                bool                    CreateChunks();
                bool                    CreateChunksForParamEstimation();
                void                    InitializeChunkVariables();
                void                    CreateChunksFromReference();
                void                    CreateChunksFromVCFReference();
                void                    ImportChunksToTarget();
                bool                    CheckChunkValidityandPrintChunk();
                bool                    CheckChunkValidityForParamEstimationandPrintChunk();
                bool                    CheckGeneticMapFile();
                void                    CreatePrintIndices(int ChunkNo, HaplotypeSet &ThisRefPanel, HaplotypeSet &ThisTarPanel);
                void                    GetCurrentPanelReady(int ChunkNo, HaplotypeSet &ThisRefPanel, HaplotypeSet &ThisRefChipOnlyPanel, HaplotypeSet &ThisTarPanel,Imputation &thisDataFast);



};


#endif // ANALYSIS_H_INCLUDED
