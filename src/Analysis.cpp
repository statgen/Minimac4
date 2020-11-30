#include "Analysis.h"
#include "Estimation.h"
#include <iomanip>
#include "assert.h"
# define RECOM_MIN 1e-05


void MyTokenize(vector<string> &result, const char *input, const char *delimiter, int Number)
{

    size_t wordCount = 1;
    result[0].clear();
    std::string *word = &result[0];


    while (*input)
    {
        if (*input==*delimiter)
        {
            // we got a delimeter, and since an empty word following
            // a delimeter still counts as a word, we allocate it here
            wordCount++;

            if((int)wordCount>Number)
                return;

            result[wordCount-1].clear();
            word = &result[wordCount-1];
        }
        else
        {
            word->push_back(*input);
        }
        input++;
    }

}



bool Analysis::CreateRecombinationMap()
{
    int FirstIndex=0;
    while(referencePanel.VariantList[FirstIndex].bp<GeneticMapData[0][0])
    {
        referencePanel.Recom[FirstIndex]=RECOM_MIN;
        FirstIndex++;

    }
    double denom =0.0, num=0.0, Val=0.0;
    int SecondIndex=FirstIndex+1;
    int i=0;
    for(i=1; i<(int)GeneticMapData.size() && SecondIndex<referencePanel.numMarkers; i++)
    {
        while(SecondIndex<referencePanel.numMarkers && referencePanel.VariantList[SecondIndex].bp<GeneticMapData[i][0])
        {
            denom += GeneticMapData[i][0] - GeneticMapData[i-1][0];
            num = referencePanel.VariantList[SecondIndex].bp - referencePanel.VariantList[SecondIndex-1].bp;
            Val += GeneticMapData[i][1];
//            referencePanel.Recom[SecondIndex-1]=(1-exp(-Val*num/denom/50))/2;
            if(num==0  || Val==0)
                referencePanel.Recom[SecondIndex-1]=RECOM_MIN;
            else
                referencePanel.Recom[SecondIndex-1]=max(RECOM_MIN, 0.01*Val*num/denom);

//            cout<<" RECOM = "<<SecondIndex-1<<"\t"<<referencePanel.Recom[SecondIndex-1]<<"\t"<<num/denom<<"\t"<<Val*num/denom<<endl;
            SecondIndex++;
            denom = 0.0;
            Val = 0.0;
        }
        denom+=GeneticMapData[i][0] - GeneticMapData[i-1][0];
        Val += GeneticMapData[i][1];
    }
    while(SecondIndex<referencePanel.numMarkers)
    {
        referencePanel.Recom[SecondIndex-1]=RECOM_MIN;
        SecondIndex++;
    }
    int h=0;

    return true;

}


String Analysis::RunAnalysis(String &Reffilename, String &Tarfilename, String &Recomfilename, String &Errorfilename)
{
    int time_prev=time(0);


    if (!referencePanel.BasicCheckForM3VCFReferenceHaplotypes(Reffilename,*MyAllVariables))
    {
        return "Reference.Panel.Load.Error";
    }
    if (!referencePanel.ReadM3VCFChunkingInformation(Reffilename,targetPanel.finChromosome))
    {
            cout << "\n Program Exiting ... \n\n";
            return "Reference.Panel.Load.Error";
    }
    if (!MyAllVariables->myModelVariables.referenceEstimates && !CreateRecombinationMap())
    {
            cout << "\n Program Exiting ... \n\n";
            return "Genetic.Map.Load.Error";
    }

    if (!targetPanel.ScaffoldGWAStoReference(referencePanel,*MyAllVariables))
    {
            cout << "\n Program Exiting ... \n\n";
            return "Reference.Panel.Load.Error";
    }


    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                           CHUNKING INFORMATION                           "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    if (!CreateChunks())
	{
        cout << "\n Program Exiting ... \n\n";
		return "Chunk.Create.Error";
	}

    InitializeRefFileStream(Reffilename);
    InitializeTargetFileStream(Tarfilename);
    TimeToRead+=( time(0) - time_prev);


    if(!MyAllVariables->myOutFormat.memUsage)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                           MAIN IMPUTATION ANALYSIS                            "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;

        std::cout << "\n Starting imputation analysis of "<<noChunks <<" chunk(s) ... "  << endl;
    }
    else
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                             MEMORY USAGE ANALYSIS                             "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;
    }

    if (!OpenStreamOutputFiles())
	{
		cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
		cout << "\n Program Exiting ... \n\n";
		return "File.Write.Error";
	}

    Imputation thisDataFast(MyAllVariables,dosages, hapdose, haps,vcfdosepartial,info,stats);

    thisDataFast.MainMarkovModel.resize(MyAllVariables->myModelVariables.cpus);

    if(MyAllVariables->myOutFormat.vcfBuffer >= targetPanel.numSamples)
        MyAllVariables->myOutFormat.vcfBuffer=targetPanel.numSamples;

    thisDataFast.InitializeOutputFiles(targetPanel, MyAllVariables->myOutFormat.vcfBuffer, MaxRefMarkerSize, MaxGwasMarkerSize);
    if(MyAllVariables->myOutFormat.memUsage)
    {

        cout<<" Estimating Memory based on a single chunk ..."<< endl;

        readm3vcfFileChunk(0, CurrentRefPanel);
        readVcfFileChunk(0, CurrentTarPanelChipOnly);
        GetCurrentPanelReady(0, CurrentRefPanel, CurrentRefPanelChipOnly, CurrentTarPanelChipOnly, thisDataFast);
        cout<<endl;

        TarMem = CurrentTarPanelChipOnly.size();
        RefMem = CurrentRefPanel.size();
        ComRefMem = CurrentRefPanelChipOnly.size();
        DosageMem = thisDataFast.Dosagesize();
        ProbMem = thisDataFast.Probsize();
        MemDisplay();
        return "Success";
    }

    for(int i=0;i<noChunks;i++)
    {
         int time_prev = time(0), time_load;

        cout<<"\n -------------------------------------------"<<endl;
        cout<<" Analyzing Chunk "<<i+1<<"/"<<noChunks<<" ["<<referencePanel.finChromosome<<":"<< MyChunks[i][0]<<"-"<<MyChunks[i][3]<<"]"<< endl;
        cout<<" -------------------------------------------"<<endl;

        readm3vcfFileChunk(i, CurrentRefPanel);
        readVcfFileChunk(i, CurrentTarPanelChipOnly);
        GetCurrentPanelReady(i, CurrentRefPanel, CurrentRefPanelChipOnly, CurrentTarPanelChipOnly, thisDataFast);




        if(MyAllVariables->myModelVariables.minimac3)
            thisDataFast.Minimac3ImputeThisChunk(i, CurrentRefPanel, CurrentTarPanelChipOnly, CurrentRefPanelChipOnly);
        else
            thisDataFast.ImputeThisChunk(i, CurrentRefPanel, CurrentTarPanelChipOnly, CurrentRefPanelChipOnly);

//abort();
        TimeToImpute+=(thisDataFast.TimeToImpute);
        TimeToWrite+=(thisDataFast.TimeToWrite);
        TimeToCompress+=(thisDataFast.TimeToCompress);


        AppendtoMainVcfFaster(i,thisDataFast.TotalNovcfParts);
        if(MyAllVariables->myOutFormat.meta)
            AppendtoMainLooVcfFaster(i,thisDataFast.TotalNovcfParts);
        PrintInfoFile(i);

        time_load = time(0) - time_prev;
        cout << "\n Time taken for this chunk = " << time_load << " seconds."<<endl;



    }

    thisDataFast.FreeMemory();

//    if(MyAllVariables->myHapDataVariables.end==0 && MyAllVariables->myOutFormat.TypedOnly)
//        assert(FileReadIndex==targetPanel.importIndexListSize);
//    asserassertt(OverCount==targetPanel.numOverlapMarkers);
//    assert(TypOnlyCount==targetPanel.numTypedOnlyMarkers);
//    assert(RefCOUNT==referencePanel.NoBlocks);

    CloseStreamOutputFiles();

    return "Success";


}



void Analysis::PrintInfoFile(int ChunkNo)

{
    InfoPrintStringLength=0;
    int RefStartPos =  MyRefVariantNumber[ChunkNo][0];

    info = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".info", "a");
    int i=0;
    for (int index = 0; index < CurrentRefPanel.RefTypedTotalCount; index++)
    {

        if(CurrentRefPanel.RefTypedIndex[index]==-1)
        {
            variant &thisVariant =  referencePanel.VariantList[i+RefStartPos];

            if(i>=CurrentRefPanel.PrintStartIndex && i <= CurrentRefPanel.PrintEndIndex)
            {

                double TarFreq = stats.AlleleFrequency(i);

                InfoPrintStringLength+=sprintf(InfoPrintStringPointer+InfoPrintStringLength , "%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t",
                MyAllVariables->myOutFormat.RsId ? thisVariant.rsid.c_str(): thisVariant.name.c_str(),
                thisVariant.refAlleleString.c_str(),
                thisVariant.altAlleleString.c_str(),
                TarFreq,
                TarFreq > 0.5 ? 1.0 - TarFreq : TarFreq,
                stats.AverageCallScore(i),
                stats.Rsq(i));

                if (!CurrentRefPanel.Targetmissing[i])
                {
                    InfoPrintStringLength+=sprintf(InfoPrintStringPointer+InfoPrintStringLength , "Genotyped\t%.3f\t%.3f\t%.5f\t%.5f\t%.5f\n",
                      stats.LooRsq(i), stats.EmpiricalR(i), stats.EmpiricalRsq(i),
                      stats.LooMajorDose(i), stats.LooMinorDose(i));
                }
                else
                 InfoPrintStringLength+=sprintf(InfoPrintStringPointer+InfoPrintStringLength , "Imputed\t-\t-\t-\t-\t-\n");
            }
            i++;
        }
        else
        {
//            if(ChunkNo==2 && index==26751)
//    abort();


            int MappingIndex = CurrentRefPanel.RefTypedIndex[index];

            if(MappingIndex>=CurrentTarPanelChipOnly.PrintTypedOnlyStartIndex && MappingIndex<=CurrentTarPanelChipOnly.PrintTypedOnlyEndIndex)
               {

                     variant ThisTypedVariant = CurrentTarPanelChipOnly.TypedOnlyVariantList[MappingIndex];
                    double TarFreq = CurrentTarPanelChipOnly.GWASOnlyAlleleFreq[MappingIndex];
                    InfoPrintStringLength+=sprintf(InfoPrintStringPointer+InfoPrintStringLength , "%s\t%s\t%s\t%.5f\t%.5f\t-\t-\tTyped_Only\t-\t-\t-\t-\t-\n",
                    MyAllVariables->myOutFormat.RsId ? ThisTypedVariant.rsid.c_str(): ThisTypedVariant.name.c_str(),
                    ThisTypedVariant.refAlleleString.c_str(),
                    ThisTypedVariant.altAlleleString.c_str(),
                    TarFreq,
                    TarFreq > 0.5 ?
                                1.0 - TarFreq : TarFreq);

               }


        }

        if(InfoPrintStringLength > 0.9 * (float)(MyAllVariables->myOutFormat.PrintBuffer))
        {
            ifprintf(info,"%s",InfoPrintStringPointer);
            InfoPrintStringLength=0;
        }
    }
    if(InfoPrintStringLength >0)
    {
        ifprintf(info,"%s",InfoPrintStringPointer);
        InfoPrintStringLength=0;
    }

    ifclose(info);

}


void Analysis::AppendtoMainLooVcfFaster(int ChunkNo, int MaxIndex)
{

    VcfPrintStringLength=0;

    int time_prev = time(0);

    int RefStartPos =  MyRefVariantNumber[ChunkNo][0];

    vector<IFILE> vcfLoodosepartialList(MaxIndex);
    vcfLoodosepartial = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".empiricalDose.vcf" + (MyAllVariables->myOutFormat.gzip ? ".gz" : ""), "a", MyAllVariables->myOutFormat.gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);


    for(int i=1;i<=MaxIndex;i++)
    {
        string tempFileIndex(MyAllVariables->myOutFormat.OutPrefix);
        stringstream strs,strs1;
        strs<<(i);
        strs1<<(ChunkNo+1);
        tempFileIndex+=(".chunk."+(string)(strs1.str())+".empiricalDose.part." +
                         (string)(strs.str())+".vcf.gz");
        vcfLoodosepartialList[i-1] = ifopen(tempFileIndex.c_str(), "r");
    }
    string line;
    int i=0;
    for (int index = 0; index < CurrentRefPanel.RefTypedTotalCount; index++)
    {
        if(CurrentRefPanel.RefTypedIndex[index]==-1)
        {

            if(i>=CurrentRefPanel.PrintStartIndex && i <= CurrentRefPanel.PrintEndIndex)
            {

                if(!CurrentRefPanel.Targetmissing[i])
                {
                    variant &tempVariant = referencePanel.VariantList[i+RefStartPos];

                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"%s\t%d\t%s\t%s\t%s\t.\tPASS",
                         tempVariant.chr.c_str(),
                         tempVariant.bp,
                         MyAllVariables->myOutFormat.RsId?tempVariant.rsid.c_str():tempVariant.name.c_str(),
                         tempVariant.refAlleleString.c_str(),
                         tempVariant.altAlleleString.c_str());
                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\tTYPED\tGT:LDS");

                    vector<string> lines(MaxIndex);
                    #pragma omp parallel for
                    for(int j=1;j<=MaxIndex;j++)
                    {
                        lines[j-1].clear();
                        vcfLoodosepartialList[j-1]->readLine(lines[j-1]);
                    }
                    for(int j=1;j<=MaxIndex;j++){
                        VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"%s",lines[j-1].c_str());
                    }

                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\n");

                }

            }
            i++;
        }

        if(VcfPrintStringLength > 0.9 * (float)(MyAllVariables->myOutFormat.PrintBuffer))
        {
            ifprintf(vcfLoodosepartial,"%s",VcfPrintStringPointer);
            VcfPrintStringLength=0;
        }

    }

    if(VcfPrintStringLength > 0)
    {
        ifprintf(vcfLoodosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringLength=0;
    }


    for(int i=1;i<=MaxIndex;i++)
    {
        ifclose(vcfLoodosepartialList[i-1]);
        string tempFileIndex(MyAllVariables->myOutFormat.OutPrefix);
        stringstream strs,strs1;
        strs<<(i);
        strs1<<(ChunkNo+1);
        tempFileIndex+=(".chunk."+(string)(strs1.str())+".empiricalDose.part." +
                         (string)(strs.str())+".vcf.gz");
        remove(tempFileIndex.c_str());
    }

    TimeToWrite+=( time(0) - time_prev);
    ifclose(vcfLoodosepartial);

}



void Analysis::AppendtoMainVcfFaster(int ChunkNo, int MaxIndex)
{
    char bgzf_partial_mode[16];
    char bgzf_final_mode[16];
    // To avoid threadbombing on reads, only use more than one CPU on reads
    // if there is more than one CPU per file to read.
    int readingCpus = (MyAllVariables->myModelVariables.cpus-1)/MaxIndex;
    if(readingCpus > 1){
        snprintf(bgzf_partial_mode, 16-2, "r@%d", readingCpus);
    } else {
        snprintf(bgzf_partial_mode, 16-2, "r@%d", 1);
    }

    // Use all available CPUs for writing
    if(MyAllVariables->myModelVariables.cpus > 1){
        snprintf(bgzf_final_mode, 16-2, "a@%d", MyAllVariables->myModelVariables.cpus);
    } else {
        // Otherwise, do a single threaded read.
        snprintf(bgzf_final_mode, 16-2, "a@%d", 1);
    }

    VcfPrintStringLength=0;

    int time_prev = time(0);

    int RefStartPos =  MyRefVariantNumber[ChunkNo][0];
    printf("\n Appending chunk to final output VCF File :  %s ",(MyAllVariables->myOutFormat.OutPrefix + ".dose.vcf" + (MyAllVariables->myOutFormat.gzip ? ".gz" : "")).c_str() );
    cout<<endl;

    vector<IFILE> vcfdosepartialList(MaxIndex);
    vcfdosepartial = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".dose.vcf" + (MyAllVariables->myOutFormat.gzip ? ".gz" : ""), bgzf_final_mode, MyAllVariables->myOutFormat.gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);


    for(int i=1;i<=MaxIndex;i++)
    {
        string tempFileIndex(MyAllVariables->myOutFormat.OutPrefix);
        stringstream strs,strs1;
        strs<<(i);
        strs1<<(ChunkNo+1);
        tempFileIndex+=(".chunk."+(string)(strs1.str())+".dose.part." +
                         (string)(strs.str())+".vcf.gz");
        vcfdosepartialList[i-1] = ifopen(tempFileIndex.c_str(), bgzf_partial_mode);
        if (!vcfdosepartialList[i-1])
        {
            std::cout << "Error: Could not open temporary file. Ensure that `ulimit -n` is at least greater than " << (MaxIndex + 5) << "." << std::endl;
            exit(-1);
        }
    }
    string line;
    int i=0;
    for (int index = 0; index < CurrentRefPanel.RefTypedTotalCount; index++)
    {
        if(CurrentRefPanel.RefTypedIndex[index]==-1)
        {

            if(i>=CurrentRefPanel.PrintStartIndex && i <= CurrentRefPanel.PrintEndIndex)
            {
                variant &tempVariant = referencePanel.VariantList[i+RefStartPos];

                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"%s\t%d\t%s\t%s\t%s\t.\tPASS",
                         tempVariant.chr.c_str(),
                         tempVariant.bp,
                         MyAllVariables->myOutFormat.RsId?tempVariant.rsid.c_str():tempVariant.name.c_str(),
                         tempVariant.refAlleleString.c_str(),
                         tempVariant.altAlleleString.c_str());



                double freq = stats.AlleleFrequency(i);

                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\tAF=%.5f;MAF=%.5f;R2=%.5f",
                        freq, freq > 0.5 ? 1.0 - freq : freq, stats.Rsq(i));


                if(!CurrentRefPanel.Targetmissing[i])
                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,";ER2=%.5f;TYPED",stats.EmpiricalRsq(i));
                else
                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,";IMPUTED");


                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\t%s",MyAllVariables->myOutFormat.formatStringForVCF.c_str());
                vector<string> lines(MaxIndex);
                #pragma omp parallel for
                for(int j=1;j<=MaxIndex;j++)
                {
                    lines[j-1].clear();
                    vcfdosepartialList[j-1]->readLine(lines[j-1]);
                }
                for(int j=1;j<=MaxIndex;j++){
                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"%s",lines[j-1].c_str());
                }
                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\n");

            }

            i++;
        }
        else
        {

            int MappingIndex = CurrentRefPanel.RefTypedIndex[index];

            if(MappingIndex>=CurrentTarPanelChipOnly.PrintTypedOnlyStartIndex && MappingIndex<=CurrentTarPanelChipOnly.PrintTypedOnlyEndIndex)
            {
                variant &ThisTypedVariant = (CurrentTarPanelChipOnly.TypedOnlyVariantList[CurrentRefPanel.RefTypedIndex[index]]);
                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"%s\t%d\t%s\t%s\t%s\t.\tPASS",
                         ThisTypedVariant.chr.c_str(),
                         ThisTypedVariant.bp,
                         MyAllVariables->myOutFormat.RsId? ThisTypedVariant.rsid.c_str():ThisTypedVariant.name.c_str(),
                         ThisTypedVariant.refAlleleString.c_str(),
                         ThisTypedVariant.altAlleleString.c_str());

                double &freq = (CurrentTarPanelChipOnly.GWASOnlyAlleleFreq[CurrentRefPanel.RefTypedIndex[index]]);

                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\tAF=%.5f;MAF=%.5f;TYPED_ONLY",
                            freq, freq > 0.5 ? 1.0 - freq : freq);

                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\t%s",MyAllVariables->myOutFormat.formatStringForVCF.c_str());
                vector<string> lines(MaxIndex);
                #pragma omp parallel for
                for(int j=1;j<=MaxIndex;j++)
                {
                    lines[j-1].clear();
                    vcfdosepartialList[j-1]->readLine(lines[j-1]);
                }
                for(int j=1;j<=MaxIndex;j++){
                    VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"%s",lines[j-1].c_str());
                }
                VcfPrintStringLength+=sprintf(VcfPrintStringPointer + VcfPrintStringLength,"\n");
            }
        }


      if(VcfPrintStringLength > 0.9 * (float)(MyAllVariables->myOutFormat.PrintBuffer))
      {

            ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
            VcfPrintStringLength=0;
      }



    }

    if(VcfPrintStringLength > 0)
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringLength=0;
    }


    for(int i=1;i<=MaxIndex;i++)
    {
        ifclose(vcfdosepartialList[i-1]);
        string tempFileIndex(MyAllVariables->myOutFormat.OutPrefix);
        stringstream strs,strs1;
        strs<<(i);
        strs1<<(ChunkNo+1);
        tempFileIndex+=(".chunk."+(string)(strs1.str())+".dose.part." +
                         (string)(strs.str())+".vcf.gz");
        remove(tempFileIndex.c_str());
    }

    TimeToWrite+=( time(0) - time_prev);
    cout << " Appending successful (" << time(0) - time_prev << " seconds) !!!"<<endl;

    ifclose(vcfdosepartial);
}




void Analysis::readm3vcfFileChunk(int ChunkNo, HaplotypeSet &ThisRefPanel)
{

    int time_prev = time(0);
    cout << "\n Reading chunk from reference panel ... "<<endl;

    vector<ReducedHaplotypeInfo> &ThisInfoVector = ThisRefPanel.ReducedStructureInfo;

    string line;
    int ThisInfoVectorIndex=0;

    int ThisChunkFilledTillRef = 0, ThisChunkStartFromRef = 0;


    int StartBlock=MyChunksInfoNumber[ChunkNo][0];
    int EndBlock=MyChunksInfoNumber[ChunkNo][1];
    int StartNextBlock=MyChunksInfoNumber[ChunkNo+1][0];


    int blockIndex=StartBlock;
    int FilltheTop=0;


    while(FilltheTop < PrevChunkFilledTillRef)
    {

        ThisInfoVector[blockIndex-StartBlock] = ThisInfoVector[PrevChunkStartFromRef + FilltheTop] ;
        blockIndex++;
        FilltheTop++;
    }
    ThisChunkStartFromRef=FilltheTop;


    for(;blockIndex<=EndBlock;blockIndex++)
    {
        RefCOUNT++;
        line.clear();
        RefFileStream->readLine(line);

//        assert(ThisInfoVectorIndex + PrevChunkFilledTillRef < MaxInfoVectorSize);

        ReducedHaplotypeInfo &tempBlock=ThisInfoVector[ThisInfoVectorIndex + PrevChunkFilledTillRef];
        ThisInfoVectorIndex++;

        referencePanel.ReadBlockHeader(line, tempBlock);
        referencePanel.ReadThisBlock(RefFileStream, blockIndex, tempBlock);


        if(blockIndex>=StartNextBlock)
        {
//            assert(ThisChunkFilledTillRef < MaxInfoVectorSize);
            ThisChunkFilledTillRef++;
        }
        else
        {
            ThisChunkStartFromRef++;
        }

    }

    PrevChunkFilledTillRef=ThisChunkFilledTillRef;
    PrevChunkStartFromRef=ThisChunkStartFromRef;

    TimeToRead+=( time(0) - time_prev);
}


void Analysis::readVcfFileChunk(int ChunkNo, HaplotypeSet &ThisTargetPanel)
{

    int time_prev = time(0);
    cout << " Reading chunk from target/GWAS panel ... "<<endl;

    int ThisnumtoBeWrittenRecords=0;

    int GWASnumtoBeWrittenRecords=0;

    int ThisChunkStartTillTar = 0, ThisChunkStartFromTar = 0;
    int ThisChunkStartTillTarOnly = 0, ThisChunkStartFromTarOnly = 0;



    int StartPos=MyTargetVariantNumber[ChunkNo][0];
    int EndPos=MyTargetVariantNumber[ChunkNo][1];
    int StartNextPos=MyTargetVariantNumber[ChunkNo+1][0];

    int TypedOnlyStartPos=MyTypdedOnlyVariantNumber[ChunkNo][0];
    int TypedOnlyEndPos=MyTypdedOnlyVariantNumber[ChunkNo][1];
    int TypedOnlyNextStartPos=MyTypdedOnlyVariantNumber[ChunkNo+1][0];

    int MainImportIndex=StartPos;
    int MainTypedOnlyImportIndex=TypedOnlyStartPos;

    int FilltheTop=0;
    while(FilltheTop < PrevChunkFilledTillTar)
    {

//        assert( (MainImportIndex-StartPos) == FilltheTop);
        for (int haplotype_index = 0; haplotype_index<(targetPanel.numHaplotypes); haplotype_index++)
        {
            ThisTargetPanel.MissingSampleUnscaffolded[haplotype_index][MainImportIndex-StartPos ] = ThisTargetPanel.MissingSampleUnscaffolded[haplotype_index][PrevChunkStartFromTar + FilltheTop] ;
            ThisTargetPanel.haplotypesUnscaffolded[haplotype_index][MainImportIndex-StartPos ] = ThisTargetPanel.haplotypesUnscaffolded[haplotype_index][PrevChunkStartFromTar + FilltheTop];
        }
        MainImportIndex++;
        FilltheTop++;
    }
    ThisChunkStartFromTar=FilltheTop;



    FilltheTop=0;
    while(FilltheTop < PrevChunkFilledTillTarOnly)
    {


//        assert( (MainTypedOnlyImportIndex-TypedOnlyStartPos) == FilltheTop);
        for (int haplotype_index = 0; haplotype_index<(targetPanel.numHaplotypes); haplotype_index++)
        {
            ThisTargetPanel.GWASOnlyMissingSampleUnscaffolded[haplotype_index][MainTypedOnlyImportIndex-TypedOnlyStartPos ] = ThisTargetPanel.GWASOnlyMissingSampleUnscaffolded[haplotype_index][PrevChunkStartFromTarOnly + FilltheTop] ;
            ThisTargetPanel.GWASOnlyhaplotypesUnscaffolded[haplotype_index][MainTypedOnlyImportIndex-TypedOnlyStartPos ] = ThisTargetPanel.GWASOnlyhaplotypesUnscaffolded[haplotype_index][PrevChunkStartFromTarOnly + FilltheTop];

            ThisTargetPanel.TypedOnlyVariantList[MainTypedOnlyImportIndex-TypedOnlyStartPos ]=ThisTargetPanel.TypedOnlyVariantList[PrevChunkStartFromTarOnly + FilltheTop];


        }
        MainTypedOnlyImportIndex++;
        FilltheTop++;
    }
    ThisChunkStartFromTarOnly=FilltheTop;



//    if(ChunkNo==3)
//        abort();

    while(MainImportIndex<=EndPos || MainTypedOnlyImportIndex <= TypedOnlyEndPos)
    {

        TarFileStream.readRecord(record);

        if (targetPanel.importIndexList[FileReadIndex++] == 0)
        {
//            assert(importReadIndex<(targetPanel.numMarkers));

            if (targetPanel.knownPosition[importReadIndex] != -1)
            {

//                assert( (ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar) < MaxGwasMarkerSize);
                int haplotype_index=0;
                for (int i = 0; i<(targetPanel.numSamples); i++)
                {

                    //haplotype_index=(2*i);

                    if(targetPanel.SampleNoHaplotypes[i]!=record.getNumGTs(i))
                    {
                        cout << "\n ERROR !!! Sample "<< targetPanel.individualName[i]<<" changes its ploidy at variant ";
                        cout << referencePanel.VariantList[targetPanel.knownPosition[importReadIndex]].name << endl;
                        if(targetPanel.finChromosome=="X" || targetPanel.finChromosome=="23" )
                        {
                            cout<<" Please impute PAR and non-PAR regions on chromosome X separately !!! "<<endl<<endl;
                        }
                        else
                            cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        abort();
                    }
                    for (int j = 0; j<targetPanel.SampleNoHaplotypes[i]; j++)
                    {

//                        assert(haplotype_index<targetPanel.numHaplotypes);

                        int alleleIndex = record.getGT(i, j);

//                        bool MissingPointer = ThisTargetPanel.MissingSampleUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar];
//                        bool AllelePointer = ThisTargetPanel.haplotypesUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar];

                        ThisTargetPanel.MissingSampleUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar]='0';
                        ThisTargetPanel.haplotypesUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar]='0';
                        if (alleleIndex<0)
                        {
                             ThisTargetPanel.MissingSampleUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar] = '1';
                        }
                        else
                        {
                             if(!targetPanel.RefAlleleSwap[overlapImportIndex])
                            {
                                if(alleleIndex==1)
                                {
                                    ThisTargetPanel.haplotypesUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar] = '1';
                                }
                            }
                            else
                            {
                                if(alleleIndex==0)
                                {
                                    ThisTargetPanel.haplotypesUnscaffolded[haplotype_index][ThisnumtoBeWrittenRecords + PrevChunkFilledTillTar] = '1';
                                }
                            }
                        }
                        haplotype_index++;
                    }
                }

                if(MainImportIndex>=StartNextPos)
                {
//                    assert( (ThisChunkStartTillTar) < MaxGwasMarkerSize);
                    ThisChunkStartTillTar++;
                }
                else
                {
                    ThisChunkStartFromTar++;
                }

//                assert(overlapImportIndex<targetPanel.numOverlapMarkers);
//                assert(overlapImportIndex==MainImportIndex);

//                assert(OverCount==MainImportIndex);

                OverCount++;

                ++overlapImportIndex;
                ++ThisnumtoBeWrittenRecords;
                ++MainImportIndex;
            }
            else if(MyAllVariables->myOutFormat.TypedOnly)
            {




//                assert( (GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly) < MaxTypedOnlyMarkerSize);

                variant tempVar=targetPanel.TypedOnlyVariantList[MainTypedOnlyImportIndex];
                ThisTargetPanel.TypedOnlyVariantList[GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly]=tempVar;

                int haplotype_index = 0;
                for (int i = 0; i<(targetPanel.numSamples); i++)
                {

                    //haplotype_index=(2*i);

                    if(targetPanel.SampleNoHaplotypes[i]!=record.getNumGTs(i))
                    {
                        cout << "\n ERROR !!! GWAS Sample "<< targetPanel.individualName[i]<<" changes its ploidy at variant ";
                        cout << tempVar.name <<endl;
                        if(targetPanel.finChromosome=="X" || targetPanel.finChromosome=="23")
                        {
                            cout<<" Please impute PAR and non-PAR regions on chromosome X separately !!! "<<endl<<endl;
                        }
                        else
                            cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " <<endl<<endl;
                        abort();
                    }
                    for (int j = 0; j<targetPanel.SampleNoHaplotypes[i]; j++)
                    {

//                        assert(haplotype_index<targetPanel.numHaplotypes);
                        int alleleIndex = record.getGT(i, j);

//                        bool MissingPointer = ThisTargetPanel.GWASOnlyMissingSampleUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly];
//                        bool AllelePointer = ThisTargetPanel.GWASOnlyhaplotypesUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly];

                        ThisTargetPanel.GWASOnlyMissingSampleUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly]='0';
                        ThisTargetPanel.GWASOnlyhaplotypesUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly]='0';
                        if (alleleIndex<0)
                        {
                             ThisTargetPanel.GWASOnlyMissingSampleUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly] = '1';
                        }
                        else
                        {
                            if(alleleIndex==1)
                            {
                                ThisTargetPanel.GWASOnlyhaplotypesUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords + PrevChunkFilledTillTarOnly] = '1';
                            }
                        }
                        haplotype_index++;
                    }
                }
//
//          if(MainTypedOnlyImportIndex==226)
//           {
//
//               cout<<"WELL = "<<record.get1BasedPosition();
//               for(int ll=0;ll<targetPanel.numHaplotypes;ll++)
//               {
//                   cout<<"\t"<<ThisTargetPanel.GWASOnlyhaplotypesUnscaffolded[ll][GWASnumtoBeWrittenRecords+ PrevChunkFilledTillTarOnly];
//               }
//           }


                if(MainTypedOnlyImportIndex>=TypedOnlyNextStartPos)
                {
//                    assert( (ThisChunkStartTillTarOnly) < MaxTypedOnlyMarkerSize);
                    ThisChunkStartTillTarOnly++;
                }
                else
                {
                    ThisChunkStartFromTarOnly++;
                }


//                assert(MainTypedOnlyImportIndex<targetPanel.numTypedOnlyMarkers);
//                assert(MainTypedOnlyImportIndex==TypOnlyCount);


                ++GWASnumtoBeWrittenRecords;
                ++MainTypedOnlyImportIndex;

                TypOnlyCount++;

            }

            importReadIndex++;

        }

    }


    PrevChunkFilledTillTar=ThisChunkStartTillTar;
    PrevChunkStartFromTar=ThisChunkStartFromTar;

    PrevChunkFilledTillTarOnly=ThisChunkStartTillTarOnly;
    PrevChunkStartFromTarOnly=ThisChunkStartFromTarOnly;

    TimeToRead+=( time(0) - time_prev);
}


bool Analysis::OpenStreamOutputFiles()
{


    bool gzip=MyAllVariables->myOutFormat.gzip;


    if(MyAllVariables->myOutFormat.vcfOutput)
    {
        VcfPrintStringPointer = (char*)malloc(sizeof(char) * (MyAllVariables->myOutFormat.PrintBuffer));

        vcfdosepartial = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".dose.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        if(vcfdosepartial==NULL)
        {
            cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< MyAllVariables->myOutFormat.OutPrefix + ".dose.vcf" + (gzip ? ".gz" : "") <<endl;
            return false;
        }

        ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
        time_t t = time(0);
        struct tm * now = localtime( & t );
        ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
        ifprintf(vcfdosepartial,"##source=Minimac4.v%s\n",VERSION);
        ifprintf(vcfdosepartial,"##contig=<ID=%s>\n",referencePanel.finChromosome.c_str());
        ifprintf(vcfdosepartial,"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description=\"Marker was genotyped but NOT imputed\">\n");

        if(MyAllVariables->myOutFormat.GT)
            ifprintf(vcfdosepartial,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        if(MyAllVariables->myOutFormat.DS)
            ifprintf(vcfdosepartial,"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");
        if(MyAllVariables->myOutFormat.HDS)
            ifprintf(vcfdosepartial,"##FORMAT=<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">\n");
        if(MyAllVariables->myOutFormat.GP)
            ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">\n");
        if(MyAllVariables->myOutFormat.SD)
            ifprintf(vcfdosepartial,"##FORMAT=<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">\n");

        ifprintf(vcfdosepartial,"##minimac4_Command=%s\n",MyAllVariables->myOutFormat.CommandLine.c_str());

        ifprintf(vcfdosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

        for(int Id=0;Id<targetPanel.numSamples;Id++)
        {
            ifprintf(vcfdosepartial,"\t%s",targetPanel.individualName[Id].c_str());
        }
        ifprintf(vcfdosepartial,"\n");
        ifclose(vcfdosepartial);

    }

    if(MyAllVariables->myOutFormat.meta)
    {

        vcfLoodosepartial = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".empiricalDose.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

        if(vcfLoodosepartial==NULL)
        {
            cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< MyAllVariables->myOutFormat.OutPrefix + ".empiricalDose.vcf" + (gzip ? ".gz" : "") <<endl;
            return false;
        }

        ifprintf(vcfLoodosepartial,"##fileformat=VCFv4.1\n");
        time_t t = time(0);
        struct tm * now = localtime( & t );
        ifprintf(vcfLoodosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
        ifprintf(vcfLoodosepartial,"##source=Minimac4.v%s\n",VERSION);
        ifprintf(vcfLoodosepartial,"##contig=<ID=%s>\n",referencePanel.finChromosome.c_str());
        ifprintf(vcfLoodosepartial,"##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">\n");
        ifprintf(vcfLoodosepartial,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotyped alleles from Array\">\n");
        ifprintf(vcfLoodosepartial,"##FORMAT=<ID=LDS,Number=1,Type=String,Description=\"Leave-one-out Imputed Dosage : Estimated Haploid Alternate Allele Dosage assuming site was NOT genotyped \">\n");
        ifprintf(vcfLoodosepartial,"##minimac4_Command=%s\n",MyAllVariables->myOutFormat.CommandLine.c_str());
        ifprintf(vcfLoodosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

        for(int Id=0;Id<targetPanel.numSamples;Id++)
        {
            ifprintf(vcfLoodosepartial,"\t%s",targetPanel.individualName[Id].c_str());
        }
        ifprintf(vcfLoodosepartial,"\n");
        ifclose(vcfLoodosepartial);

    }

    if(MyAllVariables->myOutFormat.doseOutput)
    {
        dosages = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".dose" + (gzip ? ".gz" : ""), "wb",(gzip ? InputFile::BGZF:InputFile::UNCOMPRESSED) );
        if(dosages==NULL)
        {
            cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< MyAllVariables->myOutFormat.OutPrefix + ".dose" + (gzip ? ".gz" : "") <<endl;
            return false;
        }
    }

    InfoPrintStringPointer = (char*)malloc(sizeof(char) * (MyAllVariables->myOutFormat.PrintBuffer));
    info = ifopen(MyAllVariables->myOutFormat.OutPrefix + ".info", "wb");
    if(info==NULL)
    {
        cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< MyAllVariables->myOutFormat.OutPrefix + ".info" <<endl;
        return false;
    }
    ifprintf(info, "SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n");
    ifclose(info);

    return true;
}

void Analysis::MemDisplay()
{
    ifclose(info);

    if (MyAllVariables->myOutFormat.hapOutput && !MyAllVariables->myOutFormat.unphasedOutput)
    {
        ifclose(hapdose);
    }

    if(MyAllVariables->myOutFormat.vcfOutput)
    {
        ifclose(vcfdosepartial);
    }

    if(MyAllVariables->myOutFormat.meta)
    {
        ifclose(vcfLoodosepartial);
    }

    if(MyAllVariables->myOutFormat.doseOutput)
    {
        ifclose(dosages);
    }
    ifclose(RefFileStream);
    TarFileStream.close();

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                         PREDICTED MEMORY USAGE SUMMARY                        "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    double Fact;
    string Tag;

    double Total = RefMem+ComRefMem+ TarMem+ ProbMem+ DosageMem;

    if(Total>1073741824)
    {
        Fact=1073741824;
        Tag="Gb";
    }
    else if(Total>1048576)
    {
        Fact=1048576;
        Tag="Mb";
    }
    else if(Total>1024)
    {
        Fact=1024;
        Tag="Kbytes";
    }
    else
    {
        Fact=1;
        Tag="Bytes";
    }

    printf("\n Memory occupied by Full Reference Panel            ~ %.3f %s", (double) RefMem/(Fact),Tag.c_str());
    printf("\n Memory occupied by Compressed Reference Panel      ~ %.3f %s", (double) ComRefMem/(Fact),Tag.c_str());
    printf("\n Memory occupied by Target/GWAS Panel               ~ %.3f %s", (double) TarMem/(Fact),Tag.c_str());
    printf("\n Memory occupied by Likelihood Matrices             ~ %.3f %s", (double) ProbMem/(Fact),Tag.c_str());
    printf("\n Memory occupied by Dosage Data                     ~ %.3f %s", (double) DosageMem/(Fact),Tag.c_str());
    printf("\n TOTAL Memory by the process                        ~ %.3f %s", (double) (Total)/(Fact),Tag.c_str());

    cout<<endl<<endl;


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



}

void Analysis::CloseStreamOutputFiles()
{
    cout<<endl<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              SUMMARY OF ANALYSIS                              "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    bool gzip=MyAllVariables->myOutFormat.gzip;

    cout<<endl<<" Info file written to                    : "<< MyAllVariables->myOutFormat.OutPrefix + ".info" <<endl;
    ifclose(info);
    free(InfoPrintStringPointer);

    if (MyAllVariables->myOutFormat.hapOutput && !MyAllVariables->myOutFormat.unphasedOutput)
    {
        ifclose(hapdose);
        cout<<" Haplotype Dosage information written to : "<<MyAllVariables->myOutFormat.OutPrefix + ".hapDose" + (gzip ? ".gz" : "")<<endl;
    }

    if(MyAllVariables->myOutFormat.vcfOutput)
    {
        ifclose(vcfdosepartial);
        free(VcfPrintStringPointer);
        cout<<" Imputed VCF information written to      : "<< MyAllVariables->myOutFormat.OutPrefix + ".dose.vcf" + (gzip ? ".gz" : "")<<endl;
    }

    if(MyAllVariables->myOutFormat.meta)
    {
        ifclose(vcfLoodosepartial);
        cout<<" Empirical Dosage VCF written to         : "<< MyAllVariables->myOutFormat.OutPrefix + ".empiricalDose.vcf" + (gzip ? ".gz" : "")<<endl;
    }

    if(MyAllVariables->myOutFormat.doseOutput)
    {
        ifclose(dosages);
        cout<<" Dosage information written to           : "<< MyAllVariables->myOutFormat.OutPrefix + ".dose" + (gzip ? ".gz" : "")<<endl;
    }


    cout<<endl<<" Time Taken for Reading File             = "<<TimeToRead<<" seconds";
    cout<<endl<<" Time Taken for Re-compression           = "<<TimeToCompress<<" seconds";
    cout<<endl<<" Time Taken for Imputation               = "<<TimeToImpute<<" seconds";
    cout<<endl<<" Time Taken for Writing File             = "<<TimeToWrite<<" seconds"<<endl<<endl;

    ifclose(RefFileStream);
    TarFileStream.close();

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



}




void Analysis::CreatePrintIndices(int ChunkNo, HaplotypeSet &ThisRefPanel, HaplotypeSet &ThisTarPanel)
{
    if(MyTypdedOnlyVariantNumber[ChunkNo][2]>0)
    {
        int TOnlyStartPos =  MyTypdedOnlyVariantNumber[ChunkNo][0];
        int TOnlyEndPos =  MyTypdedOnlyVariantNumber[ChunkNo][1];
        int i=TOnlyStartPos;

        if(ChunkNo==0)
        {
            if(MyAllVariables->myHapDataVariables.CHR!="")
            {
                while(targetPanel.TypedOnlyVariantList[i].bp <  MyAllVariables->myHapDataVariables.start)
                {
                    i++;
                }
            }
        }
        else
        {
           while(targetPanel.TypedOnlyVariantList[i].bp < MyChunks[ChunkNo][1])
            {
                i++;
            }
        }
        ThisTarPanel.PrintTypedOnlyStartIndex=i-TOnlyStartPos;


        if(ChunkNo==(noChunks-1))
        {
            if(MyAllVariables->myHapDataVariables.CHR!="")
            {
                for(; i <= TOnlyEndPos; i++)
                    if(targetPanel.TypedOnlyVariantList[i].bp > MyAllVariables->myHapDataVariables.end)
                        break;
            }
            else
                i=TOnlyEndPos+1;
        }
        else
        {
            for(; i <= TOnlyEndPos; i++)
                if(targetPanel.TypedOnlyVariantList[i].bp > MyChunks[ChunkNo][2])
                    break;
        }
        ThisTarPanel.PrintTypedOnlyEndIndex=i-TOnlyStartPos-1;

//        assert(ThisTarPanel.PrintTypedOnlyStartIndex <= ThisTarPanel.numTypedOnlyMarkers);
//        assert(ThisTarPanel.PrintTypedOnlyEndIndex <= ThisTarPanel.numTypedOnlyMarkers);



    }

    int RefStartPos =  MyRefVariantNumber[ChunkNo][0];
    int RefEndPos =  MyRefVariantNumber[ChunkNo][1];

    int i=RefStartPos;
    if(ChunkNo==0)
    {
        if(MyAllVariables->myHapDataVariables.CHR!="")
        {
            while(referencePanel.VariantList[i].bp <  MyAllVariables->myHapDataVariables.start)
            {
                i++;
            }
        }
    }
    else
    {
       while(referencePanel.VariantList[i].bp < MyChunks[ChunkNo][1])
        {
            i++;
        }
    }
    ThisRefPanel.PrintStartIndex=i-RefStartPos;


    if(ChunkNo==(noChunks-1))
    {
        if(MyAllVariables->myHapDataVariables.CHR!="")
        {
            for(; i <= RefEndPos; i++)
                if(referencePanel.VariantList[i].bp > MyAllVariables->myHapDataVariables.end)
                    break;
        }
        else
            i=RefEndPos+1;
    }
    else
    {
        for(; i <= RefEndPos; i++)
            if(referencePanel.VariantList[i].bp > MyChunks[ChunkNo][2])
                break;
    }
    ThisRefPanel.PrintEndIndex=i-RefStartPos-1;


//    assert(ThisRefPanel.PrintStartIndex < ThisRefPanel.numMarkers);
//    assert(ThisRefPanel.PrintEndIndex < ThisRefPanel.numMarkers);

}


void Analysis::GetCurrentPanelReady(int ChunkNo, HaplotypeSet &ThisRefPanel,
                                    HaplotypeSet &ThisRefChipOnlyPanel, HaplotypeSet &ThisTarPanel,
                                    Imputation &thisDataFast)
{

    int RefStartPos =  MyRefVariantNumber[ChunkNo][0];
    int RefEndPos =  MyRefVariantNumber[ChunkNo][1];
    int TarStartPos =  MyTargetVariantNumber[ChunkNo][0];
//    int TarEndPos =  MyTargetVariantNumber[ChunkNo][1];
    int RefStartInfo =  MyChunksInfoNumber[ChunkNo][0];
    int RefEndInfo =  MyChunksInfoNumber[ChunkNo][1];


    int NoRefMarkers = RefEndPos - RefStartPos + 1;
    int NoGwasMarkers =MyTargetVariantNumber[ChunkNo][2];
    int NoTOnlyMarkers =MyTypdedOnlyVariantNumber[ChunkNo][2];

    // initialize typical variants

    ThisRefPanel.numMarkers=NoRefMarkers;
    ThisRefChipOnlyPanel.numMarkers=NoGwasMarkers;
    ThisRefPanel.NoBlocks = MyChunksInfoNumber[ChunkNo][2];
    ThisTarPanel.numMarkers=NoGwasMarkers;
    ThisTarPanel.numOverlapMarkers=NoGwasMarkers;
    ThisTarPanel.numTypedOnlyMarkers=NoTOnlyMarkers;
    int i;

     // create new start and end indices for RefInfo

    {
        ThisRefPanel.ReducedStructureInfo[0].startIndex = 0;
        for(i=RefStartInfo; i<RefEndInfo; i++)
        {
            int temp1 = ThisRefPanel.ReducedStructureInfo[i - RefStartInfo].startIndex + referencePanel.ReducedStructureInfoSummary[i].BlockSize - 1;

            ThisRefPanel.ReducedStructureInfo[i - RefStartInfo].endIndex = temp1;
            ThisRefPanel.ReducedStructureInfo[i - RefStartInfo + 1 ].startIndex =  temp1;
        }
        ThisRefPanel.ReducedStructureInfo[RefEndInfo - RefStartInfo].endIndex = (ThisRefPanel.ReducedStructureInfo[RefEndInfo - RefStartInfo].startIndex )
                                                            + (referencePanel.ReducedStructureInfoSummary[RefEndInfo].BlockSize) - 1;

//        assert(ThisRefPanel.ReducedStructureInfo[RefEndInfo - RefStartInfo].endIndex == (NoRefMarkers-1));
//        assert((RefEndInfo - RefStartInfo + 1)==ThisRefPanel.NoBlocks);
    }

    // Create Marker to ReducedInfoMapper for Main Ref Info (NOT REFCHIP ONLY)

    {
        vector<int> &Mapper = ThisRefPanel.MarkerToReducedInfoMapper;
        int j;
        for(i=0;i<ThisRefPanel.NoBlocks;i++)
        {
            ReducedHaplotypeInfo &TempBlock=ThisRefPanel.ReducedStructureInfo[i];

            for(j=TempBlock.startIndex;j<TempBlock.endIndex;j++)
            {
                Mapper[j]=i;
            }

            if(i==(ThisRefPanel.NoBlocks-1))
                 Mapper[j]=i;
        }

    }


    // create the combined vector RefTypedIndex

    {
        ThisRefPanel.RefTypedIndex.clear();
        int ThisIndex=0,ThisMappedIndex=0;
        while(  ThisIndex<(int) targetPanel.TargetMissingTypedOnly.size()
                && RefStartPos>0
                //  && targetPanel.TargetMissingTypedOnly[ThisIndex]!=-1
                && targetPanel.TargetMissingTypedOnly[ThisIndex] < (RefStartPos))
        {
            ThisIndex++;
        }
        int counter=RefStartPos;
        while(counter<=RefEndPos && ThisIndex<(int) targetPanel.TargetMissingTypedOnly.size())
        {
            if(counter<=targetPanel.TargetMissingTypedOnly[ThisIndex])
            {
                ThisRefPanel.RefTypedIndex.push_back(-1);
                counter++;
            }
            else
            {
    //            assert(ThisMappedIndex<NoTOnlyMarkers);
                ThisRefPanel.RefTypedIndex.push_back(ThisMappedIndex++);
                ThisIndex++;
            }
        }
        while(counter<=RefEndPos)
        {
            ThisRefPanel.RefTypedIndex.push_back(-1);
            counter++;
        }
        if(ChunkNo==(noChunks-1))
        {
            while(ThisIndex<(int)targetPanel.TargetMissingTypedOnly.size())
            {
                if(counter>targetPanel.TargetMissingTypedOnly[ThisIndex])
                {
    //                assert(ThisMappedIndex<NoTOnlyMarkers);
                    ThisRefPanel.RefTypedIndex.push_back(ThisMappedIndex++);
                    ThisIndex++;
                }
                else
                    break;
            }
    }
    ThisRefPanel.RefTypedTotalCount=(int)ThisRefPanel.RefTypedIndex.size();
//    assert(ThisRefPanel.RefTypedTotalCount == (NoRefMarkers+NoTOnlyMarkers));

    }


    // Create Print Indices
    CreatePrintIndices(ChunkNo,ThisRefPanel,ThisTarPanel);

    // Create unscaffold And Flank Start and End Region for target and Targetmissing for ref and TARTOREF for REF

    {
        fill(ThisRefPanel.Targetmissing.begin(), ThisRefPanel.Targetmissing.end(), true);
        ThisTarPanel.FlankRegionStart[0]=0;

        for(i=0;i<NoGwasMarkers;i++)
        {

    //        assert((i+TarStartPos)<targetPanel.numOverlapMarkers);
    //        assert( (targetPanel.MapTarToRef[i+TarStartPos] - RefStartPos ) >= 0);
    //        assert( (targetPanel.MapTarToRef[i+TarStartPos] - RefStartPos ) < NoRefMarkers);
            ThisRefPanel.MapTarToRef[i] = targetPanel.MapTarToRef[i+TarStartPos] - RefStartPos;
            ThisRefPanel.Targetmissing[ThisRefPanel.MapTarToRef[i]]=false;

    //        ThisTarPanel.VariantList[i]=targetPanel.OverlapOnlyVariantList[i+TarStartPos];

            if(i>0)
            {
                ThisTarPanel.FlankRegionEnd[i-1]= (ThisRefPanel.MapTarToRef[i-1] + ThisRefPanel.MapTarToRef[i])/2;
                ThisTarPanel.FlankRegionStart[i]=ThisTarPanel.FlankRegionEnd[i-1]+1;
            }
        }
        ThisTarPanel.FlankRegionEnd[i-1]=NoRefMarkers-1;
    }


    // Create PrintStart and End Index and scaffold for target and Recom and Error and MAPTOREF for REF

    {

        fill(ThisRefPanel.MapRefToTar.begin(), ThisRefPanel.MapRefToTar.end(), -1);

        for(i=RefStartPos; i<=RefEndPos; i++)
        {
            if(MyAllVariables->myOutFormat.verbose)
                ThisRefPanel.VariantList[i-RefStartPos]=referencePanel.VariantList[i];

    //        ThisTarPanel.missing[i-RefStartPos]=targetPanel.missing[i];
            if(i<RefEndPos)
                ThisRefPanel.Recom[i-RefStartPos]=referencePanel.Recom[i];
            ThisRefPanel.Error[i-RefStartPos]=referencePanel.Error[i];

            if( targetPanel.MapRefToTar[i] != -1 )
            {
    //            assert( (targetPanel.MapRefToTar[i] - TarStartPos) >= 0);
    //            assert( (targetPanel.MapRefToTar[i] - TarStartPos) < NoGwasMarkers);
                ThisRefPanel.MapRefToTar[i-RefStartPos] = targetPanel.MapRefToTar[i] - TarStartPos;
            }

        }
    }

    ThisRefPanel.CalculateAlleleFreq();
    ThisTarPanel.CalculateGWASOnlyAlleleFreq();

    if(!MyAllVariables->myModelVariables.minimac3)
    {
        int time_prev = time(0);

        cout << "\n Compressing reference panel at GWAS sites ... "<<endl;

        if(MyAllVariables->myModelVariables.intermediate=="")
            ThisRefChipOnlyPanel.UncompressTypedSitesNew(ThisRefPanel,ThisTarPanel,ChunkNo);
        else
        {
            stringstream strs;
            strs<<(ChunkNo+1);
            string ss = (string)MyAllVariables->myModelVariables.intermediate.c_str()+".chunk."+(string)strs.str()+".GWAS.m3vcf.gz";
            ThisRefChipOnlyPanel.readm3vcfFile(ss.c_str(), "", 0, 0, 0);
            ThisRefChipOnlyPanel.CreateScaffoldedParameters(ThisRefPanel);
        }


        thisDataFast.TimeToCompress = time(0) - time_prev;
        time_prev=time(0);
        cout << " Re-compression successful (" << thisDataFast.TimeToCompress<< " seconds) !!!"<<endl;

        for(int i=0;i<MyAllVariables->myModelVariables.cpus;i++)
        {
            thisDataFast.MainMarkovModel[i].AssignPanels(ThisRefPanel,ThisRefChipOnlyPanel,ThisTarPanel,ThisTarPanel,MyAllVariables);
            thisDataFast.MainMarkovModel[i].initializeMatrices();

        }
    }
    else
    {
        int time_prev = time(0);



        for(int i=0;i<MyAllVariables->myModelVariables.cpus;i++)
        {
            thisDataFast.MainMarkovModel[i].AssignPanels(ThisRefPanel,ThisTarPanel,MyAllVariables);
            thisDataFast.MainMarkovModel[i].initializeMatricesMinimac3();

        }
    }


}





 void Analysis::InitializeRefFileStream(String &Reffilename)
{

    RefFileStream = ifopen(Reffilename, "r");

    string line;
    RefFileStream->readLine(line);
    bool Header=true;
    while(Header)
    {
        line.clear();
        RefFileStream->readLine(line);
        if(line.substr(0,6).compare("#CHROM")==0)
            break;
    }

    for(int tempIndex=0;tempIndex<referencePanel.NoLinesToDiscardatBeginning;tempIndex++)
        RefFileStream->discardLine();

    PrevChunkFilledTillRef=0;
    PrevChunkStartFromRef=0;
    PrevChunkFilledTillTar=0;
    PrevChunkStartFromTar=0;
    PrevChunkFilledTillTarOnly=0;
    PrevChunkStartFromTarOnly=0;

    InitializeRefChunkData(CurrentRefPanel);
    InitializeRefChipOnlyChunkData(CurrentRefPanelChipOnly);
    stats.PreInitialize(MaxRefMarkerSize);

//    if(MyAllVariables->myOutFormat.memUsage)
//    {
//        RefMem = CurrentRefPanel.size();
//        ComRefMem = CurrentRefPanelChipOnly.size();
//    }
 }


 void Analysis::InitializeTargetFileStream(String &Tarfilename)
 {
    VcfHeader header;

    TarFileStream.open(Tarfilename, header);
    TarFileStream.setSiteOnly(false);

    FileReadIndex=0;
    overlapImportIndex=0;
    importReadIndex=0;
    GWASOnlySkipIndexpoint = 0;
    MainTypedOnlyImportIndex =0;

    InitializeTargetChipOnlyChunkData(CurrentTarPanelChipOnly);

//    if(MyAllVariables->myOutFormat.memUsage)
//    {
//        TarMem = CurrentTarPanelChipOnly.size();
//    }

 }


 String Analysis::RunEstimation(String &Reffilename, String &Recomfilename, String &Errorfilename)
{
//    int time_prev=time(0);


//	if (!targetPanel.ScaffoldGWAStoReference(referencePanel,*MyAllVariables))
//	{
//		cout << "\n Program Exiting ... \n\n";
//		return "Reference.Panel.Load.Error";
//	}
//
//
    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                           CHUNKING INFORMATION                           "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    if (!CreateChunksForParamEstimation())
	{
		cout << "\n Program Exiting ... \n\n";
		return "Chunk.Create.Error";

	}

//
//    InitializeRefFileStream(Reffilename);
//    InitializeTargetFileStream(Tarfilename);
//    TimeToRead+=( time(0) - time_prev);
//
//
//    if(!MyAllVariables->myOutFormat.memUsage)
//    {
//        cout<<" ------------------------------------------------------------------------------"<<endl;
//        cout<<"                           MAIN IMPUTATION ANALYSIS                            "<<endl;
//        cout<<" ------------------------------------------------------------------------------"<<endl;
//
//        std::cout << "\n Starting imputation analysis of "<<noChunks <<" chunk(s) ... "  << endl;
//    }
//    else
//    {
//        cout<<" ------------------------------------------------------------------------------"<<endl;
//        cout<<"                             MEMORY USAGE ANALYSIS                             "<<endl;
//        cout<<" ------------------------------------------------------------------------------"<<endl;
//    }
//
//
//    OpenStreamOutputFiles();
//
//    Imputation thisDataFast(MyAllVariables,dosages, hapdose, haps,vcfdosepartial,info,stats);
//
//    thisDataFast.MainMarkovModel.resize(MyAllVariables->myModelVariables.cpus);
//
//    if(MyAllVariables->myOutFormat.vcfBuffer >= targetPanel.numSamples)
//        MyAllVariables->myOutFormat.vcfBuffer=targetPanel.numSamples;
//
//    thisDataFast.InitializeOutputFiles(targetPanel, MyAllVariables->myOutFormat.vcfBuffer, MaxRefMarkerSize, MaxGwasMarkerSize);
//    if(MyAllVariables->myOutFormat.memUsage)
//    {
//
//        cout<<" Estimating Memory based on a single chunk ..."<< endl;
//
//        readm3vcfFileChunk(0, CurrentRefPanel);
//        readVcfFileChunk(0, CurrentTarPanelChipOnly);
//        GetCurrentPanelReady(0, CurrentRefPanel, CurrentRefPanelChipOnly, CurrentTarPanelChipOnly, thisDataFast);
//        cout<<endl;
//
//        DosageMem = thisDataFast.Dosagesize();
//        ProbMem = thisDataFast.Probsize();
//        MemDisplay();
//        return "Success";
//    }
//
//    for(int i=0;i<noChunks;i++)
//    {
//         int time_prev = time(0), time_load;
//
//        cout<<"\n -------------------------------------------"<<endl;
//        cout<<" Analyzing Chunk "<<i+1<<"/"<<noChunks<<" ["<<referencePanel.finChromosome<<":"<< MyChunks[i][0]<<"-"<<MyChunks[i][3]<<"]"<< endl;
//        cout<<" -------------------------------------------"<<endl;
//
//        readm3vcfFileChunk(i, CurrentRefPanel);
//        readVcfFileChunk(i, CurrentTarPanelChipOnly);
//        GetCurrentPanelReady(i, CurrentRefPanel, CurrentRefPanelChipOnly, CurrentTarPanelChipOnly, thisDataFast);
//
//
//
//
//
//        thisDataFast.ImputeThisChunk(i, CurrentRefPanel, CurrentTarPanelChipOnly, CurrentRefPanelChipOnly);
//
////abort();
//        TimeToImpute+=(thisDataFast.TimeToImpute);
//        TimeToWrite+=(thisDataFast.TimeToWrite);
//        TimeToCompress+=(thisDataFast.TimeToCompress);
//
//
//        AppendtoMainVcfFaster(i,thisDataFast.TotalNovcfParts);
//        if(MyAllVariables->myOutFormat.meta)
//            AppendtoMainLooVcfFaster(i,thisDataFast.TotalNovcfParts);
//        PrintInfoFile(i);
//
//        time_load = time(0) - time_prev;
//        cout << "\n Time taken for this chunk = " << time_load << " seconds."<<endl;
//
//
//
//    }
//
//
////    if(MyAllVariables->myHapDataVariables.end==0 && MyAllVariables->myOutFormat.TypedOnly)
////        assert(FileReadIndex==targetPanel.importIndexListSize);
////    asserassertt(OverCount==targetPanel.numOverlapMarkers);
////    assert(TypOnlyCount==targetPanel.numTypedOnlyMarkers);
////    assert(RefCOUNT==referencePanel.NoBlocks);
//
//    CloseStreamOutputFiles();

    return "Success";


}



 String Analysis::AnalyzeExperiment(String &Reffilename, String &Tarfilename, String &Recomfilename, String &Errorfilename, AllVariable& MyAllVariable)
{
    String mysuccessresult="Error";
    cout << " Starting Main Imputation/Estimation Analysis ... " << endl;
    cout << "\n Performing preliminary check on input parameters... "  ;

    MyOutFormat=&MyAllVariable.myOutFormat;
    MyModelVariables=&MyAllVariable.myModelVariables;
    MyHapDataVariables=&MyAllVariable.myHapDataVariables;
    MyAllVariables=&MyAllVariable;

    int time_prev=time(0);
    mysuccessresult=CheckValidity(Reffilename, Tarfilename, Recomfilename, Errorfilename);
    TimeToRead+=( time(0) - time_prev);

    if(mysuccessresult!="Success")
        return mysuccessresult;

    if(!MyModelVariables->processReference)
    {
        mysuccessresult=RunAnalysis(Reffilename, Tarfilename, Recomfilename, Errorfilename);
    }
    else
    {
        Estimation MyEstimation;
        mysuccessresult = MyEstimation.RunEstimation(Reffilename, Recomfilename, Errorfilename, MyAllVariable);
    }
    if(mysuccessresult!="Success")
        return mysuccessresult;


    return "Success";


}



bool Analysis::CheckGeneticMapFile()
{
    std::cout << "\n Checking genetic map file : "<<MyAllVariables->myHapDataVariables.mapFile << endl<<endl;
    IFILE fileStream = ifopen(MyAllVariables->myHapDataVariables.mapFile, "r");
    string line;
    GeneticMapData.clear();
    const char* tabSep="\t";
    vector<string> Pieces(4);
    vector<double> AppendPiece(2);
    double PrevSumVal=0.0;

    if(fileStream)
    {
        while(fileStream->readLine(line)!=-1)
        {
            MyTokenize(Pieces, line.c_str(), tabSep,4);
            if(Pieces[0]==targetPanel.finChromosome)
            {
                AppendPiece[0]=atof(Pieces[3].c_str());
                AppendPiece[1]=atof(Pieces[2].c_str())-PrevSumVal;
                PrevSumVal=atof(Pieces[2].c_str());
                GeneticMapData.push_back(AppendPiece);
                //cout<<GeneticMapData.back()[0]<<endl;
            }
//            else
//                cout<<Pieces[0]<<endl;

            line.clear();
        }
    }
    else
    {
        cout << "\n ERROR !!! \n Program could NOT open genetic map file : " << MyAllVariables->myHapDataVariables.mapFile << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }


    if(GeneticMapData.size()<2)
    {
        cout << "\n ERROR !!! Chromosome "<<targetPanel.finChromosome <<" not found in genetic map file : " << MyAllVariables->myHapDataVariables.mapFile << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    ifclose(fileStream);

    cout<<" Successful !!! "<<endl;
    return true;
 }

 String Analysis::CheckValidity(String &Reffilename, String &Tarfilename, String &Recomfilename, String &Errorfilename)
{
    cout<<endl<<endl;
    if (Reffilename == "")
    {
        cout<< " Missing \"--refHaps\", a required parameter.\n";
        cout<< " Type \"--help\" for more help.\n\n";
        return "Command.Line.Error";
    }


    if(!MyAllVariables->myModelVariables.CheckValidity())
    {
        return "Command.Line.Error";
    }
    if(!MyAllVariables->myHapDataVariables.CheckValidity())
    {
        return "Command.Line.Error";
    }
    if(!MyAllVariables->myOutFormat.CheckValidity())
    {
        return "Command.Line.Error";
    }

//    if(MyAllVariables->myOutFormat.verbose)
//    {
//        MyAllVariables->myHapDataVariables.ChunkLengthMb=299;
//    }


    if(!MyAllVariables->myModelVariables.processReference)
    {

        if (Tarfilename == "")
        {
            cout <<" ERROR !!! \n Missing \"--haps\", a required parameter (for imputation).\n";
            cout <<" OR use \"--processReference\" to just process the reference panel.\n";
            cout << "\n Program Exiting ... \n\n";
            return "Command.Line.Error";
        }
    }
    else if (Tarfilename != "")
    {
        cout <<" ERROR !!! \n \"--haps\" is NOT valid when \"--estimate\" is ON.\n";
        cout <<" Parameter estimation only requires a reference panel. \n";
        cout <<" No GWAS panel is required.\n";
        cout << "\n Program Exiting ... \n\n";
        return "Command.Line.Error";
    }


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             PRELIMINARY FILE CHECK                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



    if(!MyAllVariables->myModelVariables.processReference)
    {
        if (!targetPanel.BasicCheckForTargetHaplotypes(Tarfilename, "GWAS", *MyAllVariables))
        {
            return "Target.Panel.Load.Error";
        }
    }

    if(!MyAllVariables->myModelVariables.referenceEstimates)
    {
        if (!CheckGeneticMapFile())
        {
            return "Genetic.Map.File.Load.Error";
        }
    }

    return "Success";

}
