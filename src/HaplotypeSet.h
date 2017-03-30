#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "Unique.h"
#include "MyVariables.h"
#include "StringBasics.h"
#include "VcfFileReader.h"
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>

using namespace std;



class HaplotypeSet
{

	public:

	    // Basic Variables
        String       inFileName;
		int         numHaplotypes,numSamples;
		int         numMarkers,NoBlocks;
		int         numOverlapMarkers, numTypedOnlyMarkers;
		int         maxBlockSize, maxRepSize;
        int RefTypedTotalCount;
        string finChromosome;
        bool PseudoAutosomal;
        bool  vcfType,m3vcfxType;
        vector<int>        optEndPoints;



		// Reduced Haplotype Information

		vector<ReducedHaplotypeInfo> ReducedStructureInfo;
		vector<ReducedHaplotypeInfoSummary> ReducedStructureInfoSummary;
		vector<vector<ReducedHaplotypeInfo> > ReducedStructureInfoBuffer;
		vector<int> MarkerToReducedInfoMapper;





        // Special Haplotype Variables for GWAS Panel

        vector<vector<bool> >     haplotypesUnscaffolded;
		vector<vector<bool> >     MissingSampleUnscaffolded;
        vector<vector<bool> >     GWASOnlyhaplotypesUnscaffolded;
		vector<vector<bool> >     GWASOnlyMissingSampleUnscaffolded;


        // Variant and Sample Information and Allele Freq

        vector<bool> Targetmissing;
        vector<string> individualName;
		vector<int> SampleNoHaplotypes;
		vector<int> *SampleNoHaplotypesPointer;
		vector<variant> VariantList;
		vector<variant> TypedOnlyVariantList;
		vector<variant> OverlapOnlyVariantList;
		vector<double> AlleleFreq;
		vector<double> GWASOnlyAlleleFreq;
        vector<double>      Recom,Error;



        // Variables for Scaffolding and Import Indexing

        vector<int> knownPosition;
        vector<int> importIndexList;
        vector<bool> RefAlleleSwap;
		vector<int>        MapTarToRef;
		vector<int>        MapRefToTar;
        vector<int> TargetMissingTypedOnly;
        vector<int> RefTypedIndex;





        // Printing Indexing Variables

		int PrintStartIndex,PrintEndIndex;
		int PrintTypedOnlyStartIndex,PrintTypedOnlyEndIndex;
        vector<int> FlankRegionStart,FlankRegionEnd;




        double size()
        {
            double S=0;

            S+=FlankRegionStart.size() * sizeof(int);
            S+=FlankRegionEnd.size() * sizeof(int);
            S+=Targetmissing.size() * sizeof(int);
            S+=knownPosition.size() * sizeof(int);
            S+=importIndexList.size() * sizeof(int);
            S+=MapTarToRef.size() * sizeof(int);
            S+=MapRefToTar.size() * sizeof(int);
            S+=TargetMissingTypedOnly.size() * sizeof(int);
            S+=RefAlleleSwap.size() * sizeof(bool);


            S+=optEndPoints.size() * sizeof(int);
            S+=SampleNoHaplotypes.size() * sizeof(bool);
            S+=Targetmissing.size() * sizeof(int);
            S+=MarkerToReducedInfoMapper.size() * sizeof(int);
            S+=AlleleFreq.size() * sizeof(double);
            S+=GWASOnlyAlleleFreq.size() * sizeof(double);
            S+=Recom.size() * sizeof(double);
            S+=Error.size() * sizeof(double);

            for(int i=0;i<(int)individualName.size();i++)
            {
                S+=individualName[i].size();
            }

            for(int i=0;i<(int)VariantList.size();i++)
            {
                S+=VariantList[i].size();
            }
            for(int i=0;i<(int)VariantList.size();i++)
            {
                S+=VariantList[i].size();
            }
            for(int i=0;i<(int)TypedOnlyVariantList.size();i++)
            {
                S+=TypedOnlyVariantList[i].size();
            }
            for(int i=0;i<(int)OverlapOnlyVariantList.size();i++)
            {
                S+=OverlapOnlyVariantList[i].size();
            }

            for(int i=0;i<(int)ReducedStructureInfo.size();i++)
            {
                S+=ReducedStructureInfo[i].size();
            }

            for(int i=0;i<(int)haplotypesUnscaffolded.size();i++)
            {
                S+=haplotypesUnscaffolded[i].size() * sizeof(bool);
                S+=MissingSampleUnscaffolded[i].size() * sizeof(bool);
            }

            for(int i=0;i<(int)GWASOnlyhaplotypesUnscaffolded.size();i++)
            {
                S+=GWASOnlyhaplotypesUnscaffolded[i].size() * sizeof(bool);
                S+=GWASOnlyMissingSampleUnscaffolded[i].size() * sizeof(bool);
            }
            return (S);
        };





        // Formatting Variables

        OutputFormatVariable *MyOutFormat;
        ModelVariable *MyModelVariables;
        HaplotypeDataVariables *MyHapDataVariables;
        AllVariable *MyAllVariables;



		// temp variables for checking, deleted later

		int importIndexListSize;


        // Minor variables for formatting stuff

        bool StartedThisPanel;
        int NoLinesToDiscardatBeginning;
        bool AlreadyReadMiddle;
        vector<string> BlockPiecesforVarInfo;




void UncompressTypedSitesNew(HaplotypeSet &rHap,HaplotypeSet &tHap);
void CreateScaffoldedParameters(HaplotypeSet &rHap);
void CreateAfterUncompressSummary();
void InvertUniqueIndexMap();
void CreateSiteSummary();
void getm3VCFSampleNames(string line);
void UpdateParameterList();
bool ReadBlockHeaderSummary(string &line, ReducedHaplotypeInfoSummary &tempBlocktoCheck);
void GetVariantInfoFromBlock(IFILE m3vcfxStream, ReducedHaplotypeInfoSummary &tempBlock,int &NoMarkersImported);
bool ReadBlockHeader(string &line, ReducedHaplotypeInfo &tempBlocktoCheck);
void ReadThisBlock(IFILE m3vcfxStream, int blockIndex, ReducedHaplotypeInfo &tempBlock);
bool BasicCheckForReferenceHaplotypes(String &Reffilename,String &Recomfilename, String &Errorfilename, AllVariable& MyAllVariable);
bool BasicCheckForTargetHaplotypes(String &Tarfilename, AllVariable& MyAllVariable);
bool ScaffoldGWAStoReference(HaplotypeSet &rHap, AllVariable& MyAllVariable);
void GetSummary(IFILE m3vcfxStream);
bool ReadM3VCFChunkingInformation(String &Reffilename,string checkChr);
bool CheckValidChrom                             (string chr);
void    writem3vcfFile                              (String filename,bool &gzip);
string  DetectFileType                        (String filename);
void CalculateGWASOnlyAlleleFreq();
void CalculateAlleleFreq();
bool RetrieveMissingScaffoldedHaplotype(int sample,int marker);
bool RetrieveScaffoldedHaplotype(int sample,int marker);
void MyTokenize(vector<string> &result,const char *input,const char *delimiter, int Number);
string FindTokenWithPrefix(const char *input,const char *delimiter,string CheckPrefix);
int CheckBlockPosFlag(string &input,String &CHR,int &START,int &END);
int GetNumVariants(string &input);
int GetNumReps(string &input);
double GetRecom(string &input);
double GetError(string &input);


//
//
//
//
//bool ReadVCFChunkingInformation(String &Tarfilename, HaplotypeSet &rHapforChunk);
//
//
//
//
//
//
//bool readm3vcfFileNew(String m3vcfFile,String CHR,int START,int END,int WINDOW);
//
//
//
//
//
//
//
//
//
//vector<vector<CrossReducedHaplotypeInfo> > CrossPanelReducedHapInfo;
//
//bool FindBasicSummaryM3VCF(String &Reffilename);
//
//
//
//bool ReadChunkingInformation(String &Reffilename);
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//         void UncompressTypedSites(HaplotypeSet &rHap,HaplotypeSet &tHap);
//
//
//
//
//
//        bool CopyVcfTargetHaplotypes(HaplotypeSet &rHap,HaplotypeSet &tHap);
//
//        bool    getScaffoldedHaplotype                      (int sample,int marker);
//        bool    getMissingScaffoldedHaplotype               (int sample,int marker);
//        void    Create                                      (int index, HaplotypeSet &rHap);
//        void    calculateFreq                               ();
//		void    CalculateFreq                               ();
//		void    CalculateGWASOnlyFreq                       ();
//        bool    FasterLoadHaplotypes                        (String filename, int maxIndiv, int maxMarker,String CNO,
//                                                            int START,int END,int WINDOW,bool rsid,bool compressOnly,
//                                                            bool filter,  String &outfile, bool &gz);
//        bool    LoadTargetHaplotypes                        (String filename, String targetSnpfile, vector<string> &refSnpList,HaplotypeSet &rHap,bool typedOnly,bool passOnly);
//		bool    LoadVcfTargetHaplotypes                     (String filename, String snpNames, vector<string> &refSnpList,HaplotypeSet &rHap);
//        void    PrintDosageGWASOnlyForVcfOutputForID        (HaplotypeSet &tHap,
//                                                            IFILE vcfdose,
//                                                            int MarkerIndex);
//        void    PrintDosageGWASOnlyForVcfOutputForIDMaleSamples        (HaplotypeSet &tHap,IFILE vcfdose,
//                                                            int MarkerIndex);
//        void    SaveIndexForGWASOnlyForVcfOutput            (int SamID,int HapId);
//        bool    readm3vcfFile                               (String m3vcfFile,String CHR,int START,int END,int WINDOW);
//        void    reconstructHaplotype                        (vector<bool> &reHaplotypes,int &index);
//        void    SaveDosageForVcfOutput                      (int hapID,vector<float> dose,vector<bool> impAlleles);
//        void    SaveDosageForVcfOutputSampleWise            (int SamID,string &SampleName, vector<float> &dose1,vector<float> &dose2,vector<bool> &impAlleles1,vector<bool> &impAlleles2);
//        void    InitializeDosageForVcfOutput                (int NHaps,int NMarkers);
//        void    InitializePartialDosageForVcfOutput         (int NHaps,int NMarkers, vector<bool> &Format);
//        void    InitializePartialDosageForVcfOutputMaleSamples    (int NHaps,int NMarkers, vector<bool> &Format);
//        void    PrintDosageForVcfOutputForID                (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
//        void    PrintPartialDosageForVcfOutputForID         (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
//        bool    BasicCheckForTargetHaplotypes               (String filename);
//
//
//        string  DetectReferenceFileType                     (String filename);
//        void    SaveDosageForVcfOutputSampleWiseChrX        (int SamID,string &SampleName, vector<float> &dose1,
//                                                            vector<bool> &impAlleles1);
//
//
//        void    CreateSummary                               ();
//        void    PrintDosageForVcfOutputForIDMaleSamples     (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);



};



#endif // HAPLOTYPESET_H_INCLUDED
