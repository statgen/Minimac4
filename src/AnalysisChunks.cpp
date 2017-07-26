#include "Analysis.h"
#include <iomanip>
#include "assert.h"



 void Analysis::InitializeRefChunkData(HaplotypeSet &ThisRefPanel)
 {


    int NoRefHaps=referencePanel.numHaplotypes;

    ThisRefPanel.numHaplotypes=NoRefHaps;
    ThisRefPanel.numSamples=referencePanel.numSamples;
    ThisRefPanel.SampleNoHaplotypes=referencePanel.SampleNoHaplotypes;
    ThisRefPanel.CummulativeSampleNoHaplotypes=referencePanel.CummulativeSampleNoHaplotypes;
    ThisRefPanel.individualName=referencePanel.individualName;
    ThisRefPanel.maxBlockSize=referencePanel.maxBlockSize;
    ThisRefPanel.maxRepSize=referencePanel.maxRepSize;
    ThisRefPanel.ReducedStructureInfo.resize(MaxInfoVectorSize);
    ThisRefPanel.MyAllVariables=MyAllVariables;


	ThisRefPanel.MapRefToTar.resize(MaxRefMarkerSize);
	ThisRefPanel.MapTarToRef.resize(MaxGwasMarkerSize);

	int Mem=0;
    for(int i=0;i<MaxInfoVectorSize;i++)
    {
        ThisRefPanel.ReducedStructureInfo[i].uniqueIndexMap.resize(NoRefHaps);
        ThisRefPanel.ReducedStructureInfo[i].uniqueCardinality.resize(referencePanel.maxRepSize);
        ThisRefPanel.ReducedStructureInfo[i].InvuniqueCardinality.resize(referencePanel.maxRepSize);
        ThisRefPanel.ReducedStructureInfo[i].TransposedUniqueHaps.resize(referencePanel.maxBlockSize);

        for(int j=0;j<referencePanel.maxBlockSize;j++)
            ThisRefPanel.ReducedStructureInfo[i].TransposedUniqueHaps[j].resize(referencePanel.maxRepSize);

        Mem+=ThisRefPanel.ReducedStructureInfo[i].size();
    }

	ThisRefPanel.Targetmissing.resize(MaxRefMarkerSize);

    ThisRefPanel.MarkerToReducedInfoMapper.resize(MaxRefMarkerSize);
    ThisRefPanel.Recom.resize(MaxRefMarkerSize);
    ThisRefPanel.Error.resize(MaxRefMarkerSize);
    ThisRefPanel.AlleleFreq.resize(MaxRefMarkerSize);
    if(MyAllVariables->myOutFormat.verbose)
    {
        ThisRefPanel.VariantList.resize(MaxRefMarkerSize);
    }

 }


 void Analysis::InitializeRefChipOnlyChunkData(HaplotypeSet &ThisRefPanel)
 {

    int NoRefHaps=referencePanel.numHaplotypes;

    ThisRefPanel.numHaplotypes=NoRefHaps;
    ThisRefPanel.numSamples=referencePanel.numSamples;
    ThisRefPanel.MarkerToReducedInfoMapper.resize(MaxGwasMarkerSize);
    ThisRefPanel.Recom.resize(MaxGwasMarkerSize);
    ThisRefPanel.Error.resize(MaxGwasMarkerSize);
    ThisRefPanel.AlleleFreq.resize(MaxGwasMarkerSize);
    ThisRefPanel.MyAllVariables=MyAllVariables;

    if(MyAllVariables->myOutFormat.verbose)
    {
        ThisRefPanel.individualName=referencePanel.individualName;
        ThisRefPanel.SampleNoHaplotypes=referencePanel.SampleNoHaplotypes;
        ThisRefPanel.VariantList.resize(MaxGwasMarkerSize);
    }
 }


 void Analysis::InitializeTargetChipOnlyChunkData(HaplotypeSet &ThisTarPanel)
 {

    int tempnumHaplotypes=targetPanel.numHaplotypes;

    ThisTarPanel.MyAllVariables=MyAllVariables;
    ThisTarPanel.numHaplotypes=targetPanel.numHaplotypes;
    ThisTarPanel.numSamples=targetPanel.numSamples;
    ThisTarPanel.SampleNoHaplotypesPointer=&targetPanel.SampleNoHaplotypes;
    ThisTarPanel.SampleNoHaplotypes=targetPanel.SampleNoHaplotypes;
    ThisTarPanel.CummulativeSampleNoHaplotypes=targetPanel.CummulativeSampleNoHaplotypes;
    ThisTarPanel.individualName=targetPanel.individualName;

    ThisTarPanel.haplotypesUnscaffolded.resize(tempnumHaplotypes);
    ThisTarPanel.MissingSampleUnscaffolded.resize(tempnumHaplotypes);


    if(MyAllVariables->myOutFormat.TypedOnly)
    {
        ThisTarPanel.GWASOnlyhaplotypesUnscaffolded.resize(tempnumHaplotypes);
        ThisTarPanel.GWASOnlyMissingSampleUnscaffolded.resize(tempnumHaplotypes);
    }



	for (int i = 0; i<tempnumHaplotypes; i++)
	{
            ThisTarPanel.haplotypesUnscaffolded[i].resize(MaxGwasMarkerSize, '0');
            ThisTarPanel.MissingSampleUnscaffolded[i].resize(MaxGwasMarkerSize, '0');
            if(MyAllVariables->myOutFormat.TypedOnly)
            {
                ThisTarPanel.TypedOnlyVariantList.resize(MaxTypedOnlyMarkerSize);
                ThisTarPanel.GWASOnlyhaplotypesUnscaffolded[i].resize(MaxTypedOnlyMarkerSize, '0');
                ThisTarPanel.GWASOnlyMissingSampleUnscaffolded[i].resize(MaxTypedOnlyMarkerSize, '0');
            }
	}


	ThisTarPanel.GWASOnlyAlleleFreq.resize(MaxTypedOnlyMarkerSize);
	ThisTarPanel.FlankRegionStart.resize(MaxGwasMarkerSize);
	ThisTarPanel.FlankRegionEnd.resize(MaxGwasMarkerSize);

 }





void Analysis::InitializeChunkVariables()
{
    int i;
    MyChunks.resize(noChunks);
    MyChunksInfoNumber.resize(noChunks+1); // extra element added just to reduce checks later on
    MyRefVariantNumber.resize(noChunks+1); // extra element added just to reduce checks later on
    MyTargetVariantNumber.resize(noChunks+1); // extra element added just to reduce checks later on
    MyTypdedOnlyVariantNumber.resize(noChunks+1); // extra element added just to reduce checks later on
    MyRatio.resize(noChunks+1,-1.0); // extra element added just to reduce checks later on


    for(i=0;i<(noChunks);i++)
    {
        MyChunks[i].resize(4);
        MyChunksInfoNumber[i].resize(3);
        MyRefVariantNumber[i].resize(3);
        MyTargetVariantNumber[i].resize(3,-3);
        MyTypdedOnlyVariantNumber[i].resize(3,-3);
    }

    MyChunksInfoNumber[i].resize(1,99999999); // extra element added just to reduce checks later on
    MyRefVariantNumber[i].resize(1,99999999); // extra element added just to reduce checks later on
    MyTargetVariantNumber[i].resize(1,99999999);  // extra element added just to reduce checks later on
    MyTypdedOnlyVariantNumber[i].resize(1,99999999);  // extra element added just to reduce checks later on




}


void Analysis::CreateChunksFromReference()
{
    vector<ReducedHaplotypeInfoSummary> &RefInfo=referencePanel.ReducedStructureInfoSummary;
    vector<variant> &VarList=referencePanel.VariantList;
    vector<int> &InfoMapper=referencePanel.MarkerToReducedInfoMapper;
    int NoMarkers=referencePanel.numMarkers;

    int i=0,index;
    MyChunksInfoNumber[0][0]=0;
    MyRefVariantNumber[0][0]=0;
    MyChunks[0][0]=VarList[0].bp;
    MyChunks[0][1]=VarList[0].bp;


    if(noChunks>1)
    {

        MyChunksInfoNumber[0][1]=InfoMapper[0 + MyHapDataVariables->ChunkSize + MyHapDataVariables->ChunkWindow];
        MyRefVariantNumber[0][1]= RefInfo[MyChunksInfoNumber[0][1]].endIndex ;

        MyChunks[0][2]=VarList[0 + MyHapDataVariables->ChunkSize].bp;
        MyChunks[0][3]=VarList[MyRefVariantNumber[0][1]].bp;


        for(i=1;i<(noChunks-1);i++)
        {
            index=(i*MyHapDataVariables->ChunkSize);

            MyChunksInfoNumber[i][0]=InfoMapper[index-MyHapDataVariables->ChunkWindow];
            MyChunksInfoNumber[i][1]=InfoMapper[index+MyHapDataVariables->ChunkSize+MyHapDataVariables->ChunkWindow];

            MyRefVariantNumber[i][0]=RefInfo[MyChunksInfoNumber[i][0]].startIndex;
            MyRefVariantNumber[i][1]= RefInfo[MyChunksInfoNumber[i][1]].endIndex;

            MyChunks[i][0]=VarList[MyRefVariantNumber[i][0]].bp;
            MyChunks[i][1]=MyChunks[i-1][2]+1;
            MyChunks[i][2]=VarList[index+MyHapDataVariables->ChunkSize].bp;
            MyChunks[i][3]=VarList[MyRefVariantNumber[i][1]].bp;

        }

    index=(i*MyHapDataVariables->ChunkSize);

    MyChunksInfoNumber[i][0]=InfoMapper[index-MyHapDataVariables->ChunkWindow];
    MyRefVariantNumber[i][0]=RefInfo[MyChunksInfoNumber[i][0]].startIndex;

    MyChunks[i][0]=VarList[MyRefVariantNumber[i][0]].bp;
    MyChunks[i][1]=MyChunks[i-1][2]+1;

    }

    MyChunks[i][2]=VarList[NoMarkers-1].bp;
    MyChunks[i][3]=VarList[NoMarkers-1].bp;
    MyChunksInfoNumber[i][1]=InfoMapper[NoMarkers-1];
    MyRefVariantNumber[i][1]=NoMarkers-1;


    for(int i=0;i<noChunks;i++)
    {
        MyChunksInfoNumber[i][2] = MyChunksInfoNumber[i][1] - MyChunksInfoNumber[i][0] + 1;
        MyRefVariantNumber[i][2] = MyRefVariantNumber[i][1] - MyRefVariantNumber[i][0] + 1;
    }

}

void Analysis::CreateChunksFromVCFReference()
{
    MyChunks.resize(noChunks);
    MyRefVariantNumber.resize(noChunks+1); // extra element added just to reduce checks later on
    int i=0;
    for( i=0;i<(noChunks);i++)
    {
        MyChunks[i].resize(4);
        MyRefVariantNumber[i].resize(3);
    }
    MyRefVariantNumber[i].resize(1,99999999); // extra element added just to reduce checks later on


    vector<variant> &VarList=referencePanel.VariantList;
    int NoMarkers=referencePanel.numMarkers;
    int index;

    MyRefVariantNumber[0][0]=0;
    MyChunks[0][0]=VarList[0].bp;
    MyChunks[0][1]=VarList[0].bp;
    i=0;

    if(noChunks>1)
    {
        MyChunks[0][2]=VarList[0 + MyHapDataVariables->ChunkSize].bp;
        MyChunks[0][3]=VarList[0 + MyHapDataVariables->ChunkSize + MyHapDataVariables->ChunkWindow].bp;
        MyRefVariantNumber[0][1]=MyHapDataVariables->ChunkSize + MyHapDataVariables->ChunkWindow;

        for(i=1;i<(noChunks-1);i++)
        {
            index=(i*MyHapDataVariables->ChunkSize);

            MyRefVariantNumber[i][0]=index-MyHapDataVariables->ChunkWindow;
            MyChunks[i][0]=VarList[index-MyHapDataVariables->ChunkWindow].bp;
            MyChunks[i][1]=MyChunks[i-1][2]+1;
            MyChunks[i][2]=VarList[index+MyHapDataVariables->ChunkSize].bp;
            MyChunks[i][3]=VarList[index+MyHapDataVariables->ChunkSize+MyHapDataVariables->ChunkWindow].bp;
            MyRefVariantNumber[i][1]= index+MyHapDataVariables->ChunkSize+MyHapDataVariables->ChunkWindow;

        }

    index=(i*MyHapDataVariables->ChunkSize);

    MyRefVariantNumber[i][0]=index-MyHapDataVariables->ChunkWindow;
    MyChunks[i][0]=VarList[index-MyHapDataVariables->ChunkWindow].bp;
    MyChunks[i][1]=MyChunks[i-1][2]+1;

    }
    MyChunks[i][2]=VarList[NoMarkers-1].bp;
    MyChunks[i][3]=VarList[NoMarkers-1].bp;
    MyRefVariantNumber[i][1]=NoMarkers-1;

    for(int i=0;i<noChunks;i++)
    {
        MyRefVariantNumber[i][2] = MyRefVariantNumber[i][1] - MyRefVariantNumber[i][0] + 1;
    }

}



void Analysis::ImportChunksToTarget()
{
    OverCount =0;
    TypOnlyCount = 0;
    RefCOUNT = 0;

    int i=0;
    for(int CurrentChunkNo=0; CurrentChunkNo<noChunks; CurrentChunkNo++)
    {

        while(i<targetPanel.numOverlapMarkers && targetPanel.MapTarToRef[i]<MyRefVariantNumber[CurrentChunkNo][0])
        {
            i++;
        }

        MyTargetVariantNumber[CurrentChunkNo][0]=i;

        while(i<targetPanel.numOverlapMarkers && targetPanel.MapTarToRef[i] <= MyRefVariantNumber[CurrentChunkNo][1] )
        {
            i++;
        }
        MyTargetVariantNumber[CurrentChunkNo][1]=i-1;
        i=MyTargetVariantNumber[CurrentChunkNo][0];
    }

//    assert(MyTargetVariantNumber[noChunks-1][1]==(targetPanel.numOverlapMarkers-1));  // Last position of Target Panel indexed should be number of overlapping markers

        i=0;
        for(int CurrentChunkNo=0; CurrentChunkNo<noChunks; CurrentChunkNo++)
        {
            if(CurrentChunkNo>0)
            {
                while(i<targetPanel.numTypedOnlyMarkers && targetPanel.TypedOnlyVariantList[i].bp < MyChunks[CurrentChunkNo][0] )
                {
                    i++;
                }
            }
            MyTypdedOnlyVariantNumber[CurrentChunkNo][0]=i;

            while(i<targetPanel.numTypedOnlyMarkers && targetPanel.TypedOnlyVariantList[i].bp <= MyChunks[CurrentChunkNo][3] )
            {
                i++;
            }
            MyTypdedOnlyVariantNumber[CurrentChunkNo][1]=i-1;
            i=MyTypdedOnlyVariantNumber[CurrentChunkNo][0];

        }
        MyTypdedOnlyVariantNumber[noChunks-1][1]=targetPanel.numTypedOnlyMarkers-1;


    for(int i=0;i<noChunks;i++)
    {
        MyTargetVariantNumber[i][2]=0;
        MyTypdedOnlyVariantNumber[i][2]=0;

        MyTargetVariantNumber[i][2] = ( MyTargetVariantNumber[i][1] + 1 > MyTargetVariantNumber[i][0]? MyTargetVariantNumber[i][1] + 1 - MyTargetVariantNumber[i][0]: 0 );
        MyTypdedOnlyVariantNumber[i][2] = ( MyTypdedOnlyVariantNumber[i][1] + 1 >MyTypdedOnlyVariantNumber[i][0]? MyTypdedOnlyVariantNumber[i][1] + 1 - MyTypdedOnlyVariantNumber[i][0]: 0 );
    }

}


bool Analysis::CheckChunkValidityandPrintChunk()
{

    int MinRefMarkerSize=MyRefVariantNumber[0][2], InIndex=0;
    int MinGwasMarkerSize=MyTargetVariantNumber[0][2], MaIndex=0;
    float minRatio=100.0,raIndex=0;

    MaxInfoVectorSize=0;
    MaxRefMarkerSize=0;
    MaxGwasMarkerSize=0;
    MaxTypedOnlyMarkerSize=0;

    cout<<" No   LeftBuffer      LeftEnd   RightPoint  RightBuffer       #Sites(GWAS/Ref/%)"<<endl;
    cout<<" -------------------------------------------------------------------------------"<<endl;

    for(int i=0;i<noChunks;i++)
    {
        if( MyChunksInfoNumber[i][2]  > MaxInfoVectorSize )
              MaxInfoVectorSize = MyChunksInfoNumber[i][2];

        if( MyRefVariantNumber[i][2]  > MaxRefMarkerSize )
              MaxRefMarkerSize = MyRefVariantNumber[i][2];
        if( MyRefVariantNumber[i][2]  < MinRefMarkerSize )
        {
          InIndex=i;
          MinRefMarkerSize = MyRefVariantNumber[i][2];
        }


        if( MyTargetVariantNumber[i][2]  > MaxGwasMarkerSize )
              MaxGwasMarkerSize = MyTargetVariantNumber[i][2];
        if( MyTargetVariantNumber[i][2]  < MinGwasMarkerSize )
        {
            MaIndex=i;
            MinGwasMarkerSize = MyTargetVariantNumber[i][2];
        }

        if( MyTypdedOnlyVariantNumber[i][2]  > MaxTypedOnlyMarkerSize )
              MaxTypedOnlyMarkerSize = MyTypdedOnlyVariantNumber[i][2];


        if(MyTargetVariantNumber[i][2] >0 && MyRefVariantNumber[i][2] > 0 )
        {
            MyRatio[i]=100.0*(float)MyTargetVariantNumber[i][2]/(float)MyRefVariantNumber[i][2];
            if(minRatio>MyRatio[i])
            {
                minRatio=MyRatio[i];
                raIndex=i;
            }
        }


        cout<<setw(3)<<i+1<<"  "
        <<setw(11)<<MyChunks[i][0]<<"  "
        <<setw(11)<<MyChunks[i][1]<<"  "
        <<setw(11)<<MyChunks[i][2]<<"  "
        <<setw(11)<<MyChunks[i][3]<<"  "
        <<setw(7)<<MyTargetVariantNumber[i][2]<<"/"
        <<setw(8)<<MyRefVariantNumber[i][2]<<"/";

        if(MyRatio[i]==-1.0)
            cout<<"  nan";
        else
            printf(" %2.2f",MyRatio[i]);
        cout<<"%"<< endl;
    }
       cout<<endl<<endl;


    if(MinRefMarkerSize==0)
    {
        cout<<"\n ERROR !!! ERROR !!! ERROR !!! "<<endl;
        cout<<" Chunk "<< InIndex+1<<" has 0 variants from the reference panel in it ... "<<endl;
        cout<<" Please increase the value of \"--ChunkLengthMb\" to analyze larger chunks ..." <<endl<<endl;
        return false;

    }
    if(MinGwasMarkerSize==0)
    {
        cout<<"\n ERROR !!! ERROR !!! ERROR !!! "<<endl;
        cout<<" Chunk "<< MaIndex+1<<" has 0 variants from the GWAS panel in it ... "<<endl;
        cout<<" Please increase the value of \"--ChunkLengthMb\" to analyze larger chunks ..." <<endl<<endl;
        cout<<" Program Exiting ... "<<endl<<endl;
        return false;

    }
    if(minRatio<(MyAllVariables->myHapDataVariables.minRatio))
    {
        cout<<"\n ERROR !!! ERROR !!! ERROR !!! "<<endl;
        cout<<" Chunk "<< raIndex+1<<" has less than "<<MyAllVariables->myHapDataVariables.minRatio <<"\% of variants from the GWAS panel overlapping with the reference panel ... "<<endl;
        cout<<" Please increase the value of \"--ChunkLengthMb\" to analyze larger chunks "<<endl;
        cout<<" or decrease the value of \"--minRatio\" = " <<MyAllVariables->myHapDataVariables.minRatio <<endl<<endl;
        return false;

    }
    return true;

}


bool Analysis::CheckChunkValidityForParamEstimationandPrintChunk()
{
    int MinRefMarkerSize=MyRefVariantNumber[0][2], InIndex=0;
    MaxRefMarkerSize=0;

    cout<<" No   LeftBuffer      LeftEnd   RightPoint  RightBuffer       #Sites"<<endl;
    cout<<" -------------------------------------------------------------------------------"<<endl;

    for(int i=0;i<noChunks;i++)
    {

        if( MyRefVariantNumber[i][2]  > MaxRefMarkerSize )
              MaxRefMarkerSize = MyRefVariantNumber[i][2];

        if( MyRefVariantNumber[i][2]  < MinRefMarkerSize )
        {
          InIndex=i;
          MinRefMarkerSize = MyRefVariantNumber[i][2];
        }

        cout<<setw(3)<<i+1<<"  "
        <<setw(11)<<MyChunks[i][0]<<"  "
        <<setw(11)<<MyChunks[i][1]<<"  "
        <<setw(11)<<MyChunks[i][2]<<"  "
        <<setw(11)<<MyChunks[i][3]<<"  "
        <<setw(8)<<MyRefVariantNumber[i][2]<<"/";

    }
    cout<<endl<<endl;

    if(MinRefMarkerSize==0)
    {
        cout<<"\n ERROR !!! ERROR !!! ERROR !!! "<<endl;
        cout<<" Chunk "<< InIndex+1<<" has 0 variants from the reference panel in it ... "<<endl;
        cout<<" Please increase the value of \"--ChunkLengthMb\" to analyze larger chunks ..." <<endl<<endl;
        return false;

    }
    return true;

}

 bool Analysis::CreateChunks()
{

    GetNumChunks();

    std::cout << "\n Chunking region into "<<noChunks <<" chunk(s) with atleast "<<
    MyHapDataVariables->ChunkSize <<" variants in each chunk ... "  << endl<<endl;
    std::cout << " Details of chunks is given below ..."  << endl<<endl;

    InitializeChunkVariables();
    CreateChunksFromReference();
    ImportChunksToTarget();

    return CheckChunkValidityandPrintChunk();

}



 bool Analysis::CreateChunksForParamEstimation()
{

    GetNumChunks();

    std::cout << "\n Chunking region into "<<noChunks <<" chunk(s) with atleast "<<
    MyHapDataVariables->ChunkSize <<" variants in each chunk ... "  << endl<<endl;
    std::cout << " Details of chunks is given below ..."  << endl<<endl;

    CreateChunksFromVCFReference();

    return CheckChunkValidityForParamEstimationandPrintChunk();

}



 void Analysis::GetNumChunks()
 {

    int NoMarkers=referencePanel.numMarkers;
    int StartPos = referencePanel.VariantList[0].bp;
    int EndPos = referencePanel.VariantList[NoMarkers-1].bp;

    int noOverLaps = ((int)(EndPos-StartPos)/(int)(1000000 *MyAllVariables->myHapDataVariables.ChunkOverlapMb));
    if(noOverLaps==0)
        noOverLaps=1;

    MyHapDataVariables->ChunkWindow=(NoMarkers/noOverLaps);

    noChunks = ((int)(EndPos-StartPos)/(int)(1000000 *MyAllVariables->myHapDataVariables.ChunkLengthMb));
    if(noChunks==0)
        noChunks=1;

    MyHapDataVariables->ChunkSize=(NoMarkers/noChunks);

    return;
 }

