#ifndef MINIMAC4_ESTIMATION_H
#define MINIMAC4_ESTIMATION_H


#if defined(_OPENMP)
#include <omp.h>
#endif
#include "MyVariables.h"
#include "Unique.h"
#include "Imputation.h"
#include "HaplotypeSet.h"
#include <unordered_map>

using namespace std;


class Estimation
{


public:

    // Basic Variables

    HaplotypeSet targetPanel,referencePanel;
    IFILE RefFileStream;
    VcfFileReader TarFileStream;
    VcfRecord record;


    // Swapping Variables for actual chunk IMPUTATION
    HaplotypeSet CurrentTarPanel, CurrentRefPanel, CurrentRefPanelLoo, CurrentTarPanelLoo;
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
    IFILE m3vcfpartial, recfile, errfile ;
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
    char *m3vcfPrintStringPointer, *recPrintStringPointer, *errPrintStringPointer;
    int VcfPrintStringLength, InfoPrintStringLength, DosePrintStringLength;
    int m3vcfPrintStringLength, recPrintStringLength, errPrintStringLength;




    vector<ReducedHaplotypeInfo> InfoforChunk;


    vector< vector<int> > BPListMapper;

    HaplotypeSet *CurrentChunk, *NextChunk;




    HaplotypeSet MainReferenceFrame;

    OutputFormatVariable *MyOutFormat;
    ModelVariable *MyModelVariables;
    HaplotypeDataVariables *MyHapDataVariables;
    AllVariable *MyAllVariables;



    void PrintInfoFile(int ChunkNo);

    Estimation()
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


    void MemDisplay();

    void GetNumChunks();
    String                  RunEstimation(String &Reffilename, String &Recomfilename, String &Errorfilename,  AllVariable& MyAllVariable);
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
    bool                    CheckChunkValidityandPrintChunkForEstimation();
    bool                    CreateChunks();
    bool                    CreateChunksForParamEstimation();
    void                    InitializeChunkVariables();
    void                    CreateChunksFromReference();
    void                    CreateChunksFromVCFReference();
    void                    ImportChunksToTarget();
    bool                    CheckChunkValidityandPrintChunk();
    bool                    CheckChunkValidityForParamEstimationandPrintChunk();
    void                    InitializeTargetLooChunkData(HaplotypeSet &ThisTarPanel);
    void                    CreatePrintIndices(int ChunkNo, HaplotypeSet &ThisRefPanel);
    void                    GetCurrentPanelReady(int ChunkNo, HaplotypeSet &ThisRefPanel, Imputation &thisDataFast);



};



#endif //MINIMAC4_ESTIMATION_H
