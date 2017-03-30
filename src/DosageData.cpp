#include "DosageData.h"
#include "HaplotypeSet.h"

#include "assert.h"





void DosageData::FlushPartialVcf(int NovcfParts)
{
    PrintStringLength=0;

    int time_Start = time(0);


    string PartialVcfFileName(MyAllVariables->myOutFormat.OutPrefix);
    stringstream strs,strs1;
    strs<<(NovcfParts);
    strs1<<(ChunkNo+1);
    string tempStringforMeta = PartialVcfFileName+(".chunk."+(string)(strs1.str())+ ".maskedDose.vcf.part."+ (string)(strs.str()));
    PartialVcfFileName+=(".chunk."+(string)(strs1.str())+ ".dose.vcf.part."+ (string)(strs.str()));

    IFILE vcfdosepartial = ifopen(PartialVcfFileName.c_str(), "wb", InputFile::UNCOMPRESSED);
    IFILE vcfLoodosepartial = NULL;
    if(MyAllVariables->myOutFormat.meta)
        vcfLoodosepartial=ifopen(tempStringforMeta.c_str(), "wb", InputFile::UNCOMPRESSED);


    for(int Id=0;Id<BuffernumSamples;Id++)
    {
        int indexID = InvertSampleIndex[Id];
        int FullSamID=SampleIndex[indexID];
        BufferSampleNoHaplotypes[indexID]=tHapFull->SampleNoHaplotypes[FullSamID];
    }

    int i=0;
    for (int index = 0; index < rHapFull->RefTypedTotalCount; index++)
    {

        if(rHapFull->RefTypedIndex[index]==-1)
        {

            if(i>=rHapFull->PrintStartIndex && i <= rHapFull->PrintEndIndex)
            {
                PrintDosageForVcfOutputForID(i);
            }
            i++;
        }
        else
        {
            int MappingIndex = rHapFull->RefTypedIndex[index];

            if(MappingIndex>=tHapFull->PrintTypedOnlyStartIndex && MappingIndex<=tHapFull->PrintTypedOnlyEndIndex)
            {
                PrintGWASOnlyForVcfOutputForID(MappingIndex);
            }
        }

        if(PrintStringLength > 0.9 * (float)(MyAllVariables->myOutFormat.PrintBuffer))
        {
            ifprintf(vcfdosepartial,"%s",PrintStringPointer);
            PrintStringLength=0;
        }

    }
    if(PrintStringLength>0)
    {
        ifprintf(vcfdosepartial,"%s",PrintStringPointer);
        PrintStringLength=0;
    }

    ifclose(vcfdosepartial);

    if(MyAllVariables->myOutFormat.meta)
        ifclose(vcfLoodosepartial);

    TimeToWrite=time(0) - time_Start;

}


void DosageData::PrintDiploidLooDosage(bool a, bool b, float &x, float &y)
{
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"\t");
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%d|%d",(a),(b));
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f|%.3f",x , y);
}


void DosageData::PrintDiploidDosage(float &x, float &y)
{

    bool colonIndex=false;
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"\t");

    if(MyAllVariables->myOutFormat.GT)
    {

        if(!MyAllVariables->myOutFormat.unphasedOutput)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%d|%d",(x>0.5),(y>0.5));
        else
        {
            bool a1=(x>0.5);
            bool a2=(y>0.5);

            if((a1^a2)==1)
                PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"0/1");
            else if(a1 && a2)
                PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"1/1");
            else
                PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"0/0");
        }
        colonIndex=true;
    }
    if(MyAllVariables->myOutFormat.DS)
    {

        if(colonIndex)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f",x+ y);
        colonIndex=true;
    }
    if(MyAllVariables->myOutFormat.HDS)
    {

        if(colonIndex)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f|%.3f",x , y);
        colonIndex=true;
    }
    if(MyAllVariables->myOutFormat.GP)
    {

        if(colonIndex)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
        colonIndex=false;
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f,%.3f,%.3f",(1-x)*(1-y),x*(1-y)+y*(1-x),x*y);
    }

}

void DosageData::PrintHaploidDosage(float &x)
{
    bool colonIndex=false;
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"\t");

    if(MyAllVariables->myOutFormat.GT)
    {
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%d",(x>0.5));
        colonIndex=true;
    }
    if(MyAllVariables->myOutFormat.DS)
    {

        if(colonIndex)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f",x);
        colonIndex=true;
    }
    if(MyAllVariables->myOutFormat.HDS)
    {

        if(colonIndex)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f",x );
        colonIndex=true;
    }
    if(MyAllVariables->myOutFormat.GP)
    {

        if(colonIndex)
            PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
        colonIndex=false;
        PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f,%.3f",1-x,x);
    }

}

void DosageData::PrintHaploidLooDosage(bool a, float &x)
{
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"\t");
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%d",(a));
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,":");
    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"%.3f",x);
}



void DosageData::PrintDosageForVcfOutputForID(int MarkerIndex)
{
    for(int Id=0;Id<BuffernumSamples;Id++)
    {
        int IndexId = InvertSampleIndex[Id];
        int NoHaps = BufferSampleNoHaplotypes[IndexId];
//        assert(SampleIndex[IndexId]==(Id+FirstHapId));

        if(NoHaps==2)
            PrintDiploidDosage((hapDosage[2*IndexId])[MarkerIndex] , (hapDosage[2*IndexId+1])[MarkerIndex] );
        else if(NoHaps==1)
            PrintHaploidDosage((hapDosage[2*IndexId])[MarkerIndex] );
        else
            abort();
    }

  PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"\n");

}

void DosageData::PrintGWASOnlyForVcfOutputForID(int MarkerIndex)
{
    float freq=(float)tHapFull->GWASOnlyAlleleFreq[MarkerIndex];
    float x,y;

    for(int Id=0;Id<BuffernumSamples;Id++)
    {
        int IndexId = InvertSampleIndex[Id];
        int FullSamID=SampleIndex[IndexId];
        int NoHaps=BufferSampleNoHaplotypes[IndexId];

//        assert(FullSamID==(Id+FirstHapId));

        if(NoHaps==2)
        {
            bool a1=tHapFull->GWASOnlyMissingSampleUnscaffolded[2*FullSamID][MarkerIndex];
            bool a2=tHapFull->GWASOnlyMissingSampleUnscaffolded[2*FullSamID+1][MarkerIndex];

            if(a1 || a2)
                PrintDiploidDosage(freq,freq);
            else
             {
                 x=(float)tHapFull->GWASOnlyhaplotypesUnscaffolded[2*FullSamID][MarkerIndex];
                 y=(float)tHapFull->GWASOnlyhaplotypesUnscaffolded[2*FullSamID+1][MarkerIndex];
                 PrintDiploidDosage(x, y);
            }
        }

        else if(NoHaps==1)
        {
            bool a1=tHapFull->GWASOnlyMissingSampleUnscaffolded[2*FullSamID][MarkerIndex];
            if(a1)
                PrintHaploidDosage(freq);
            else
             {
                 x=(float)tHapFull->GWASOnlyhaplotypesUnscaffolded[2*FullSamID][MarkerIndex];
                 PrintHaploidDosage(x);
             }

        }
        else
            abort();
    }

    PrintStringLength+=sprintf(PrintStringPointer+PrintStringLength,"\n");
}


void DosageData::BindSampleMModel(MarkovModel &MM, int Index, int HapNo)
{
    MM.DosageHap=&hapDosage[ (2*Index) + HapNo] ;
    MM.LooDosageHap=&LoohapDosage[ (2*Index) + HapNo] ;
}



pair <int, int> DosageData::IndexSample(int SampleId)
{
    pair <int, int> Result;
    InvertSampleIndex[SampleId-FirstHapId]=NoSamplesIndexed;
    SampleIndex[NoSamplesIndexed++]=SampleId;
    Result.second=NoSamplesIndexed-1;
    if(NoSamplesIndexed==BuffernumSamples)
    {
        Result.first=1;
    }

    return Result;
}


void DosageData::InitializePartialDosageData(HaplotypeSet &tarInitializer, int MaxNoSamples,
                                             int MaxNoRefVariants, int MaxNoTarVariants,
                                               AllVariable *MyAllVariable)
{

    MyAllVariables=MyAllVariable;

    ActualnumHaplotypes = tarInitializer.numHaplotypes;
    ActualnumSamples = tarInitializer.numSamples;

    BuffernumSamples = MaxNoSamples;
    BuffernumHaplotypes = 2*MaxNoSamples;

    hapDosage.resize(BuffernumHaplotypes);
    LoohapDosage.resize(BuffernumHaplotypes);
    BufferSampleNoHaplotypes.resize(BuffernumSamples);
    SampleIndex.resize(BuffernumSamples);
    InvertSampleIndex.resize(BuffernumSamples);

    for(int i=0;i<BuffernumHaplotypes;i++)
        {
            hapDosage[i].resize(MaxNoRefVariants,0.0);
            LoohapDosage[i].resize(MaxNoTarVariants,0.0);
        }

    PrintStringPointer = (char*)malloc(sizeof(char) * (MyAllVariables->myOutFormat.PrintBuffer));
    individualName=tarInitializer.individualName;

}



void DosageData::ReParameterizePartialDosageData(int No, HaplotypeSet &refFullHap, HaplotypeSet &tarFullHap)
{
    ChunkNo=No;
    rHapFull=&refFullHap;
    tHapFull=&tarFullHap;
    FirstHapId=0;
    numMarkers = rHapFull->numMarkers;
    noGWASSites= tHapFull->numTypedOnlyMarkers;
    NoSamplesIndexed=0;

}


void DosageData::UpdatePartialDosageData(int NewMaxVal, int NewFirstHapId)
{
    if(NewMaxVal<=0)
        return;

    BuffernumSamples = NewMaxVal;
    BuffernumHaplotypes = 2*NewMaxVal;
    NoSamplesIndexed=0;
    FirstHapId=NewFirstHapId;

}

