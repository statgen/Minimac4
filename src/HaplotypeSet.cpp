#include "HaplotypeSet.h"
#include "assert.h"
#define ALT_DELIM ','
#define MONOMORPH_INDICATOR '-'

#include "STLUtilities.h"

void HaplotypeSet::UncompressTypedSitesNew(HaplotypeSet &rHap,HaplotypeSet &tHap,int ChunkNo)
{
//    outFile=MyAllVariables->myOutFormat.OutPrefix;
    vcfType=false;

    int CPU=1;

    int    blockSize = 500;
    int    bufferSize = 5000;

    vector<CompressedHaplotype> CompHaplotypes(CPU);

    int currentPiece=0;
    for(currentPiece=0;currentPiece<CPU;currentPiece++)
    {
        CompHaplotypes[currentPiece].Size(rHap.ReducedStructureInfo, bufferSize);

    }
    currentPiece=0;

    vector<int> BufferPosList(CPU);
    BufferPosList[0]=1;
    for(int BufferNo=1;BufferNo<CPU;BufferNo++)
    {
        BufferPosList[BufferNo]= BufferPosList[BufferNo-1]+(bufferSize-1);
    }
    ReducedStructureInfoBuffer.clear();
    ReducedStructureInfoBuffer.resize(CPU);
    ReducedStructureInfo.clear();


    int NoMarkersWritten=1;
    int RefmarkerIndex=-1;
//    int CompressedIndex=0;
    vector<ReducedHaplotypeInfo> &RhapInfo=rHap.ReducedStructureInfo;


//    assert(numMarkers==tHap.numMarkers);

    for(int ThisPanelIndex=0; ThisPanelIndex<numMarkers; ThisPanelIndex++)
    {

        RefmarkerIndex=rHap.MapTarToRef[ThisPanelIndex];
//        assert(RefmarkerIndex<rHap.numMarkers);

        int j=rHap.MarkerToReducedInfoMapper[RefmarkerIndex];
//        assert(j<rHap.NoBlocks);

        ReducedHaplotypeInfo &ThisRhapInfo=RhapInfo[j];

        int ThisMarkerIndex=RefmarkerIndex-ThisRhapInfo.startIndex;
//        assert(ThisMarkerIndex<ThisRhapInfo.BlockSize);

        AlleleFreq[ThisPanelIndex]=rHap.AlleleFreq[RefmarkerIndex];
        if(MyAllVariables->myOutFormat.verbose)
        {
            VariantList[ThisPanelIndex]=rHap.VariantList[RefmarkerIndex];
        }
        CompHaplotypes[currentPiece].Push(j,ThisMarkerIndex);


            int NewPiece=CPU;
            if (CompHaplotypes[currentPiece].Length == bufferSize)
            {
                NewPiece=(currentPiece+1)%CPU;

                vector<String> tempHaplotypes(numHaplotypes);
                if(NewPiece!=0)
                {

                    int temp=CompHaplotypes[currentPiece].ReducedInfoMapper[bufferSize-1];
                    int temp2=CompHaplotypes[currentPiece].ReducedInfoVariantMapper[bufferSize-1];
                    CompHaplotypes[NewPiece].Clear();
                    CompHaplotypes[NewPiece].Push(temp,temp2);
                    currentPiece=NewPiece;
                }
            }


        if(NewPiece==0 || ThisPanelIndex==(numMarkers-1))
        {

            #pragma omp parallel for
            for(int ThisPiece=0;ThisPiece<=currentPiece;ThisPiece++)
            {

                int LastflushPos=BufferPosList[ThisPiece]-1;
                vector<int> index(numHaplotypes),oldIndex;
                vector<int> previousDifference(numHaplotypes);
                vector<int> previousPredecessor(numHaplotypes);
                vector<int> firstDifference(numHaplotypes-1,0);
                vector<int> cost(bufferSize+1,0);
                vector<int> bestSlice(bufferSize+1,0);
                vector<int> bestComplexity(bufferSize+1,0);
                vector<vector<int> > bestIndex(bufferSize+1);

                ReducedStructureInfoBuffer[ThisPiece].clear();
                findUnique RefUnique;
                RefUnique.updateCoeffs(MyAllVariables->myModelVariables.transFactor,MyAllVariables->myModelVariables.cisFactor);
                double blockedCost = 0.0;

                for(int i=0;i<numHaplotypes;i++)
                    index[i]=i;

                for(int length=1;length<=CompHaplotypes[currentPiece].Length;length++)
                {

                    CompHaplotypes[ThisPiece].RetrieveVariant(length-1);

                    vector<int> offsets(3,0);
                    for (int i = 0; i < numHaplotypes; i++)
                    {
                        offsets[CompHaplotypes[ThisPiece].GetVal(i) - '0' + 1]++;
                    }

                    offsets[2]+=offsets[1];
                    oldIndex = index;
                    for (int i = 0; i < numHaplotypes; i++)
                    {
                        index[offsets[CompHaplotypes[ThisPiece].GetVal(oldIndex[i],length - 1) - '0']++] = oldIndex[i];
                    }

                    RefUnique.UpdateDeltaMatrix(CompHaplotypes[ThisPiece], index, firstDifference, length, blockSize,
                           oldIndex, previousPredecessor, previousDifference);

                    RefUnique.AnalyzeBlocks(index, firstDifference, length, blockSize,
                       cost, bestSlice, bestComplexity, bestIndex);

                }

                if(CompHaplotypes[ThisPiece].Length>1)
                    blockedCost += RefUnique.FlushBlocks(ReducedStructureInfoBuffer[ThisPiece],
                                                        LastflushPos,
                                                        CompHaplotypes[ThisPiece], cost,
                                                        bestComplexity, bestSlice, bestIndex);

            }




            NoMarkersWritten+=(CPU*(bufferSize-1));

            BufferPosList[0]=NoMarkersWritten;
            for(int BufferNo=1;BufferNo<CPU;BufferNo++)
            {
                BufferPosList[BufferNo]= BufferPosList[BufferNo-1]+(bufferSize-1);

            }

            int temp=CompHaplotypes[currentPiece].ReducedInfoMapper[bufferSize-1];
            int temp2=CompHaplotypes[currentPiece].ReducedInfoVariantMapper[bufferSize-1];

            CompHaplotypes[0].Clear();
            CompHaplotypes[0].Push(temp,temp2);


            for(int ThisPiece=0;ThisPiece<CPU;ThisPiece++)
            {
                for(int jj=0;jj<(int)ReducedStructureInfoBuffer[ThisPiece].size();jj++)
                    {
                        ReducedStructureInfo.push_back(ReducedStructureInfoBuffer[ThisPiece][jj]);
                    }
                ReducedStructureInfoBuffer[ThisPiece].clear();
            }
            currentPiece=0;
        }

    }

    CreateAfterUncompressSummary();
    CreateScaffoldedParameters(rHap);
    InvertUniqueIndexMap();

    if(MyAllVariables->myOutFormat.verbose)
    {
        stringstream strs;
        strs<<(ChunkNo+1);
        string ss=(string)MyAllVariables->myOutFormat.OutPrefix+".chunk."+(string)(strs.str())+".GWAS";
        String tempString(ss.c_str());

        writem3vcfFile( tempString, MyAllVariables->myOutFormat.gzip);
    }


}

void HaplotypeSet::CreateScaffoldedParameters(HaplotypeSet &rHap)
{

    for(int i=0;i<numMarkers;i++)
    {
        Error[i]=rHap.Error[rHap.MapTarToRef[i]];
        if(i>0)
        {
            int index=rHap.MapTarToRef[i-1];
            double temp=1.0;
            while(index<rHap.MapTarToRef[i])
            {
                temp*=(1.0-(rHap.Recom[index]));
                index++;
            }
            Recom[i-1]=(1.0-temp);
        }
    }
}


void HaplotypeSet::CreateAfterUncompressSummary()
{
    NoBlocks=(int)ReducedStructureInfo.size();

    int i,j;
    maxBlockSize=0;
    maxRepSize=0;
//    optEndPoints.clear();
    for(i=0;i<NoBlocks;i++)
    {
        ReducedHaplotypeInfo &TempBlock=ReducedStructureInfo[i];
//        optEndPoints.push_back(TempBlock.startIndex);

        if(maxBlockSize<TempBlock.BlockSize)
            maxBlockSize=TempBlock.BlockSize;
        if(maxRepSize<TempBlock.RepSize)
            maxRepSize=TempBlock.RepSize;


        for(j=TempBlock.startIndex;j<TempBlock.endIndex;j++)
        {
            MarkerToReducedInfoMapper[j]=i;
        }

        if(i==(NoBlocks-1))
             MarkerToReducedInfoMapper[j]=i;
    }
//    optEndPoints.push_back(ReducedStructureInfo[i-1].endIndex);

}



void printErr(String filename)
{
    cout<<"\n ERROR !!! \n Error in M3VCF File !!! "<<endl;
    cout<<" Please re-construct the following [.m3vcf] file using Minimac3/4 and try again ..."<<endl;
    cout<< " [ "<< filename<<" ] "<<endl;
    cout<<" Contact author if problem still persists : sayantan@umich.edu "<<endl;
    cout<<" Program Exiting ..."<<endl<<endl;
    abort();
}


void HaplotypeSet::InvertUniqueIndexMap()
{

    for(int i=0;i<NoBlocks;i++)
    {
        ReducedHaplotypeInfo &TempBlock=ReducedStructureInfo[i];
        TempBlock.uniqueIndexReverseMaps.resize(TempBlock.RepSize);
        for (int j = 0; j < numHaplotypes; j++)
        {
            TempBlock.uniqueIndexReverseMaps[TempBlock.uniqueIndexMap[j]].push_back(j);
        }

    }
}




void HaplotypeSet::CreateSiteSummary()
 {
    NoBlocks=(int)ReducedStructureInfoSummary.size();
    maxBlockSize=0;
    maxRepSize=0;


    for(int i=0;i<NoBlocks;i++)
    {
        ReducedHaplotypeInfoSummary &TempBlock=ReducedStructureInfoSummary[i];
        if(maxBlockSize<TempBlock.BlockSize)
            maxBlockSize=TempBlock.BlockSize;
        if(maxRepSize<TempBlock.RepSize)
            maxRepSize=TempBlock.RepSize;

    }

    MarkerToReducedInfoMapper.resize(numMarkers, 0);
    int i,j;
    for(i=0;i<NoBlocks;i++)
    {
        ReducedHaplotypeInfoSummary &TempBlock=ReducedStructureInfoSummary[i];
        for(j=TempBlock.startIndex;j<TempBlock.endIndex;j++)
        {
            MarkerToReducedInfoMapper[j]=i;
        }
        if(i==(NoBlocks-1))
             MarkerToReducedInfoMapper[j]=i;
    }


}




void HaplotypeSet::UpdatePloidySummary(string line)
{
    if (m3vcfVERSION == 1)
        return;

    CummulativeSampleNoHaplotypes.resize(numSamples, 0);
    SampleNoHaplotypes.resize(numSamples, 0);

    char *pch;
    char *end_str1;
    string tempString2, tempString, tempString3;
    int colCount = 0, sampleCount = 0;
    pch = strtok_r((char *) line.c_str(), "\t", &end_str1);

    while (pch != NULL)
    {
        colCount++;
        if (colCount > 9)
        {
            tempString = string(pch);
            if(tempString.find('|')==string::npos)
                SampleNoHaplotypes[sampleCount]=1;
            else
                SampleNoHaplotypes[sampleCount]=2;
            sampleCount++;
        }

        pch = strtok_r (NULL,"\t", &end_str1);
    }

    for (int i = 1; i < numSamples; i++)
        CummulativeSampleNoHaplotypes[i] += CummulativeSampleNoHaplotypes[i - 1] + SampleNoHaplotypes[i - 1];

    numHaplotypes = CummulativeSampleNoHaplotypes.back()+SampleNoHaplotypes.back();

}



void HaplotypeSet::getm3VCFSampleNames(string &line)
{

    individualName.clear();
    SampleNoHaplotypes.clear();
    CummulativeSampleNoHaplotypes.clear();
    numSamples=0;
    numHaplotypes=0;

    char * pch;
    char *end_str1;
    string tempString2,tempString,tempString3;
    int colCount=0;
    pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);

    numHaplotypes=0;

    while(pch!=NULL)
    {
        colCount++;
        if(m3vcfVERSION==1 && colCount>9)
        {
            numHaplotypes++;

            tempString=string(pch);
            tempString3=tempString.substr(tempString.size()-1,tempString.size()-1);
            tempString2=tempString.substr(0,tempString.size()-6);

            if(tempString3=="1")
            {
                numSamples++;
                individualName.push_back(tempString2);
                SampleNoHaplotypes.push_back(1);
            }
            else if(tempString3=="2")
            {
                SampleNoHaplotypes.back()=2;
            }
            else
            {
                cout<<endl<<"\n ERROR !!! \n Inconsistent Sample Name !!! "<<endl<<endl;
                cout<<" Haplotype Number suffix cannot be more than 2 ..."<<endl;
                cout<<" Erroneous Sample Name found : "<<tempString<<endl;
                printErr(inFileName);
            }
        }
        else if(m3vcfVERSION==2 && colCount>9)
        {
            numSamples++;
            tempString=string(pch);
            individualName.push_back(tempString);
        }
        pch = strtok_r (NULL,"\t", &end_str1);
    }

    if(m3vcfVERSION==1)
    {
        CummulativeSampleNoHaplotypes.resize(numSamples, 0);
        for (int i = 1; i < numSamples; i++)
            CummulativeSampleNoHaplotypes[i] += CummulativeSampleNoHaplotypes[i - 1] + SampleNoHaplotypes[i - 1];
    }
}




void HaplotypeSet::UpdateParameterList()
{
    MyOutFormat=&(MyAllVariables->myOutFormat);
    MyModelVariables=&(MyAllVariables->myModelVariables);
    MyHapDataVariables=&(MyAllVariables->myHapDataVariables);

}




bool HaplotypeSet::ReadBlockHeaderSummary(string &line, ReducedHaplotypeInfoSummary &tempBlocktoCheck)
{
    const char* tabSep="\t";
    vector<string> BlockPieces(numHaplotypes+9);

    MyTokenize(BlockPieces, line.c_str(), tabSep, 9);

    tempBlocktoCheck.BlockSize=GetNumVariants(BlockPieces[7]);
    tempBlocktoCheck.RepSize=GetNumReps(BlockPieces[7]);

    if(CheckBlockPosFlag(line, MyHapDataVariables->CHR, MyHapDataVariables->START, MyHapDataVariables->END)==1)
        return true;

    return false;

}


void HaplotypeSet::GetVariantInfoFromBlock(IFILE m3vcfxStream, ReducedHaplotypeInfoSummary &tempBlock, int &NoMarkersImported)
{
    string currID, rsID;
    string line;
    int blockEnterFlag=0;
    const char* tabSep="\t";
    int NewBlockSizeCount=0;

    for(int tempIndex=0;tempIndex<tempBlock.BlockSize;tempIndex++)
    {

        line.clear();
        m3vcfxStream->readLine(line);
        MyTokenize(BlockPiecesforVarInfo, line.c_str(), tabSep,9);


        if(StartedThisPanel==false || tempIndex>0)
        {
            StartedThisPanel=true;

            currID=BlockPiecesforVarInfo[0]+":"+BlockPiecesforVarInfo[1]+":"+BlockPiecesforVarInfo[3]+":"+BlockPiecesforVarInfo[4];
            rsID = BlockPiecesforVarInfo[2];
            if(rsID==".")
                rsID=currID;

            variant tempVariant;
            tempVariant.assignValues(currID,rsID,BlockPiecesforVarInfo[0],atoi(BlockPiecesforVarInfo[1].c_str()));
            tempVariant.assignRefAlt(BlockPiecesforVarInfo[3],BlockPiecesforVarInfo[4]);
            VariantList.push_back(tempVariant);
            
            string Info=BlockPiecesforVarInfo[7];

            double tempRecom=GetRecom(Info);
            double tempError=GetError(Info);

            if(MyAllVariables->myModelVariables.constantParam>0.0)
            {
                Recom.push_back(MyAllVariables->myModelVariables.constantParam);
                Error.push_back(0.00999);
            }
            else
            {
                if (tempRecom == -3.0 && MyAllVariables->myModelVariables.referenceEstimates) {
                    cout << "\n ERROR !!! \n NO parameter estimates already found in M3VCF file !!!" << endl;
                    cout << " Please remove handle \"--referenceEstimates\" to ignore this check " << endl;
                    cout << " Else, use M3VCF file with parameter estimates " << endl;
                    cout << "\n Program Exiting ... \n\n";
                    abort();
                }
                Recom.push_back(tempRecom);
                Error.push_back(0.00999);
            }
            NewBlockSizeCount++;
            NoMarkersImported++;
        }

        if(blockEnterFlag==0)
        {
            tempBlock.startIndex=NoMarkersImported-1;
            blockEnterFlag=1;
        }
        tempBlock.endIndex=NoMarkersImported-1;
    }

}


void HaplotypeSet::reconstructHaplotype(vector<AlleleType> &reHaplotypes,int &index)
{
    int markerIndex=0,k;
    for(int j=0;j<NoBlocks;j++)
    {
        int ThisIndex = ReducedStructureInfo[j].uniqueIndexMap[index];
        for(k=0;k<(ReducedStructureInfo[j].BlockSize-1);k++)
        {
            reHaplotypes[markerIndex++]=ReducedStructureInfo[j].TransposedUniqueHaps[k][ThisIndex];
        }

        if(j==(NoBlocks-1))
        {
            reHaplotypes[markerIndex]=ReducedStructureInfo[j].TransposedUniqueHaps[k][ThisIndex];
        }
    }
}


void HaplotypeSet::Create(int index, HaplotypeSet &rHap)
{
    vector<AlleleType> padded(rHap.numMarkers);
    rHap.reconstructHaplotype(padded,index);
    numMarkers=(int)padded.size();
    haplotypesUnscaffolded[0]= padded;
}



bool HaplotypeSet::ReadBlockHeader(string &line, ReducedHaplotypeInfo &tempBlocktoCheck)
{
    const char* tabSep="\t", *dashSep="|";
    vector<string> BlockPieces(numHaplotypes+9);

    MyTokenize(BlockPieces, line.c_str(), tabSep,numHaplotypes+9);
    tempBlocktoCheck.BlockSize=GetNumVariants(BlockPieces[7]);
    tempBlocktoCheck.RepSize=GetNumReps(BlockPieces[7]);

    if(CheckBlockPosFlag(line,MyHapDataVariables->CHR,MyHapDataVariables->START,MyHapDataVariables->END)==1)
        return true;


    fill(tempBlocktoCheck.uniqueCardinality.begin(), tempBlocktoCheck.uniqueCardinality.end(), 0);

    int index=0;

    if(m3vcfVERSION==1) {
        while (index < numHaplotypes) {
            int tempval = atoi(BlockPieces[index + 9].c_str());
            tempBlocktoCheck.uniqueIndexMap[index] = tempval;
            tempBlocktoCheck.uniqueCardinality[tempval]++;
            index++;
        }
    }
    else if (m3vcfVERSION==2)
    {
        int haploIndex = 0;
        while (index < numSamples)
        {
            if(SampleNoHaplotypes[index]==2)
            {
                vector<string> HaploPieces(2);
                MyTokenize(HaploPieces, BlockPieces[index + 9].c_str(), dashSep, 2);

                int tempval = atoi(HaploPieces[0].c_str());
                tempBlocktoCheck.uniqueIndexMap[haploIndex++] = tempval;
                tempBlocktoCheck.uniqueCardinality[tempval]++;

                tempval = atoi(HaploPieces[1].c_str());
                tempBlocktoCheck.uniqueIndexMap[haploIndex++] = tempval;
                tempBlocktoCheck.uniqueCardinality[tempval]++;
            }
            else
            {
                int tempval = atoi(BlockPieces[index + 9].c_str());
                tempBlocktoCheck.uniqueIndexMap[haploIndex++] = tempval;
                tempBlocktoCheck.uniqueCardinality[tempval]++;
            }

            index++;
        }

    }

    for (int i = 0; i < tempBlocktoCheck.RepSize; i++)
    {
        tempBlocktoCheck.InvuniqueCardinality[i]=1.0/(float)tempBlocktoCheck.uniqueCardinality[i];
    }

    return false;

}



void HaplotypeSet::GetTransUniqueHapsVERSION2(int index, ReducedHaplotypeInfo &tempBlock, string &tempString)
{
    vector<AlleleType> &TempHap = tempBlock.TransposedUniqueHaps[index];
    fill(TempHap.begin(), TempHap.end(), '0');
    vector<int> AlternateAlleles;
    AlternateAlleles.clear();
    int prevVal = 0;

    char *input=&tempString[0];
    string word="";
    while (*input)
    {
        if(*input==MONOMORPH_INDICATOR)
        {
            break;
        }

        if (*input==ALT_DELIM)
        {
            word.push_back('\0');
            AlternateAlleles.push_back(prevVal+atoi(word.c_str()));
            prevVal=AlternateAlleles.back();
            word.clear();
        }
        else
        {
            word.push_back(*input);
        }

        input++;
    }

    if(word!="")
    {
        word.push_back('\0');
        AlternateAlleles.push_back(prevVal+atoi(word.c_str()));
    }

    for(int i=0; i<(int)AlternateAlleles.size(); i++) {
        TempHap[AlternateAlleles[i]]='1';
    }
}


void HaplotypeSet::ReadThisBlock(IFILE m3vcfxStream,
                                 int blockIndex, ReducedHaplotypeInfo &tempBlock)
{
    string line;
    const char* tabSep="\t";
    vector<string> BlockPieces(9);
    for(int tempIndex=0;tempIndex<tempBlock.BlockSize;tempIndex++)
    {
        line.clear();
        m3vcfxStream->readLine(line);
        MyTokenize(BlockPieces, line.c_str(), tabSep,9);

        vector<AlleleType> &TempHap = tempBlock.TransposedUniqueHaps[tempIndex];

        string &tempString=BlockPieces[8];

        if(m3vcfVERSION==1) {
            for (int index = 0; index < tempBlock.RepSize; index++) {
                char t = tempString[index];
                TempHap[index] = (t);
            }
        }
        else if(m3vcfVERSION==2) {
            GetTransUniqueHapsVERSION2(tempIndex, tempBlock, tempString);

        }
    }
}




bool HaplotypeSet::BasicCheckForTargetHaplotypes(String &VCFFileName,
                                                 String TypeofFile,
                                                 AllVariable& MyAllVariable)
{
    MyAllVariables=&MyAllVariable;
    std::cout << "\n Checking "<<TypeofFile<<" haplotype file : "<<VCFFileName << endl;

    string FileType=DetectFileType(VCFFileName);

    if(FileType.compare("NA")==0)
    {
        cout << "\n ERROR !!! \n Program could NOT open file : " << VCFFileName << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }
    else if(FileType.compare("vcf")!=0)
    {
        cout << "\n ERROR !!! \n GWAS File provided by \"--haps\" must be a VCF file !!! \n";
        cout << " Please check the following file : "<<VCFFileName<<endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    return GetVariantInfofromVCFFile(VCFFileName, TypeofFile, MyAllVariable);
}


bool HaplotypeSet::BasicCheckForVCFReferenceHaplotypes(String &VCFFileName,
                                                       String TypeofFile,
                                                       AllVariable& MyAllVariable)
{
    MyAllVariables=&MyAllVariable;
    std::cout << "\n Checking Reference haplotype file : "<<VCFFileName << endl;

    string FileType=DetectFileType(VCFFileName);

     if(FileType.compare("NA")==0)
    {
        cout << "\n ERROR !!! \n Program could NOT open file : " << VCFFileName << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    if(FileType.compare("vcf")!=0)
    {
        cout << "\n ERROR !!! \n If \"--processReference\" is ON,";
        cout << " Reference  File provided by \"--refHaps\" must be a VCF file !!! \n";
        cout << " Please check the following file : "<<VCFFileName<<endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    if(MyAllVariables->myModelVariables.rounds==0)
    {

        cout <<"\n NOTE: User has specified \"--rounds\" = 0 !!!\n";
        cout<<"       No parameter estimation will be done on VCF file.\n";
        cout<<"       Program will use default estimates leading to possibly inaccurate estimates."<<endl;
    }

    return GetVariantInfofromVCFFile(VCFFileName, TypeofFile, MyAllVariable);
}


bool HaplotypeSet::BasicCheckForM3VCFReferenceHaplotypes(String &Reffilename,
                                                    AllVariable& MyAllVariable)

{
    MyAllVariables=&MyAllVariable;
    UpdateParameterList();
    std::cout << "\n Checking Reference haplotype file : "<<Reffilename << endl;
    inFileName=Reffilename;

    string refFileType=DetectFileType(Reffilename);

    if(refFileType.compare("NA")==0)
    {
        cout << "\n ERROR !!! \n Program could NOT open file : " << Reffilename << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    if(refFileType.compare("Invalid")==0)
    {
        cout << "\n ERROR !!! \n Reference File provided by \"--refHaps\" must be a M3VCF file !!! \n";
        cout << " Please check the following file : "<<Reffilename<<endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    if(refFileType.compare("vcf")==0)
    {

        cout << "\n ERROR !!! \n VCF Format detected ...";
        cout << "\n The current version of Minimac4 can ONLY handle M3VCF files for imputation "<<endl;
        cout <<   " Please convert the VCF file to M3VCF using Minimac3 "<<endl;
        cout<<    " We will implement this feature in Minimac4 very soon "<<endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }
	return true;

}



bool HaplotypeSet::GetVariantInfofromVCFFile(String &VCFFileName, String TypeofFile, AllVariable& MyAllVariable)
{
	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    inFileName=VCFFileName;

	inFile.setSiteOnly(true);
    inFile.open(VCFFileName, header);


    numSamples=header.getNumSamples();
    if(numSamples==0)
    {
        cout << "\n ERROR !!! \n No samples found in "<< TypeofFile <<" File : "<<VCFFileName<<endl;
		cout << " Please check the file properly..\n";
		cout << "\n Program Exiting ... \n\n";
        return false;
    }

    for (int i = 0; i < numSamples; i++)
	{
		string tempName(header.getSampleName(i));
		individualName.push_back(tempName);
	}


   std::cout << "\n Gathering variant information ..." << endl <<endl;



    int numReadRecords = 0,numActualRecords=0;
    int bp, failFilter=0,notBiallelic=0,inconsistent=0,duplicates=0, outSideRegion = 0;
    string prevID="",currID, refAllele,altAllele,PrefrefAllele,PrevaltAllele,cno,fixCno,id;

    while (inFile.readRecord(record))
    {



        int flag=0;
        cno=record.getChromStr();
        bp=record.get1BasedPosition();
        id=record.getIDStr();
        refAllele = record.getRefStr();
        altAllele = record.getAltStr();

        // Check Valid Chromosome and single chromosome
        {
            if(numActualRecords==0)
            {
                if(!CheckValidChrom(cno))
                {
                    cout << "\n ERROR !!! \n "<< TypeofFile <<" VCF File contains chromosome : "<<cno<<endl;
                    cout << " VCF File can only contain chromosomes 1-22, X(23), Y(24) !!! "<<endl;
                    cout << "\n Program Exiting ... \n\n";
                    return false;
                }
                fixCno=cno;
                finChromosome=fixCno;
            }
            else if(fixCno!=cno)
            {
                cout << "\n ERROR !!! \n "<< TypeofFile <<" VCF File contains multiple chromosomes : "<<cno<<", "<<fixCno<<", ... "<<endl;
                cout << " Please use VCF file with single chromosome !!! "<<endl;
                cout << "\n Program Exiting ... \n\n";
                return false;
            }

        }

        // Check Window Parameters provided by user
        {
            if(MyAllVariable.myHapDataVariables.CHR!="")
            {
                if(cno.compare(MyAllVariable.myHapDataVariables.CHR.c_str())!=0)
                    flag = 1;
                else
                {
                    if(MyAllVariable.myHapDataVariables.END>0)
                    {
                        if(bp>MyAllVariable.myHapDataVariables.END || bp<MyAllVariable.myHapDataVariables.START)
                            flag = 1;
                    }
                    else
                        if(bp<MyAllVariable.myHapDataVariables.START)
                            flag = 1;
                }
            }
            if(flag==1)
                outSideRegion++;
        }


        // Check bi-allelic and FILTER
        {
             if (record.getNumAlts()>1)
            {
                notBiallelic++;
                flag = 1;
            }
            if (MyAllVariable.myHapDataVariables.passOnly && record.getFilter().getString(0).compare("PASS") != 0)
            {
                failFilter++;
                flag = 1;
            }


        }

		// Check ID for duplicate
        {
            stringstream strs3,strs4;
            strs3<<(cno);
            strs4<<(bp);
            currID=(string)strs3.str()+":"+(string)strs4.str()+":"+refAllele+":"+altAllele;
            if(id==".")
                id=currID;
            if(prevID==currID)
            {
                duplicates++;
                if(MyAllVariable.myOutFormat.verbose){cout << " WARNING !!! Duplicate Variant found chr:"<<cno<<":"<<bp<<" with identical REF = "<<refAllele <<" and ALT = "<<altAllele <<"\n";}
                if(!(MyAllVariable.myHapDataVariables.ignoreDuplicates))
                {
                    cout << "\n ERROR !!! \n Duplicate Variant found chr:"<<cno<<":"<<bp<<" with identical REF = "<<refAllele <<" and ALT = "<<altAllele <<"\n";
                    cout<<"\n Use handle \"--ignoreDuplicates\" to ignore duplicate instances ... "<<endl;
                    cout << "\n Program Exiting ... \n\n";
                    return false;
                }
                flag=1;


            }
            prevID=currID;
            PrefrefAllele=refAllele;
            PrevaltAllele=altAllele;
        }

        // Check length of SNPs REF/ALT allele
        {
            // Removed this to allow '-' and '.' in the GWAS Panel
        }

        if(flag==0)
        {
            variant thisVariant(currID,cno,bp);
            VariantList.push_back(thisVariant);
            VariantList[numReadRecords].refAlleleString=refAllele;
            VariantList[numReadRecords].altAlleleString=altAllele;
            VariantList[numReadRecords].rsid=id;
            ++numReadRecords;
        }
        numActualRecords++;
        importIndexList.push_back(flag);

    }

    importIndexListSize=(int)importIndexList.size();


    numMarkers=numReadRecords;

    if(numActualRecords==0)
    {
        cout << "\n ERROR !!! \n No variants found in "<< TypeofFile <<" File : "<<VCFFileName<<endl;
		cout << " Please check the file properly..\n";
		cout << "\n Program Exiting ... \n\n";
        return false;
    }


    inFile.close();
    SampleNoHaplotypes.resize(numSamples,2);
    CummulativeSampleNoHaplotypes.resize(numSamples,0);

    inFile.open(VCFFileName, header);
    inFile.setSiteOnly(false);
    inFile.readRecord(record);
    int tempHapCount=0;
    for (int i = 0; i<(numSamples); i++)
    {
        if(record.getNumGTs(i)==0)
        {
            std::cout << "\n ERROR !!! \n Empty Value for Individual : " << individualName[i] << " at First Marker  " << endl;
            std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
            cout << "\n Program Exiting ... \n\n";
            return false;
        }
        else
        {
            CummulativeSampleNoHaplotypes[i]=tempHapCount;
            SampleNoHaplotypes[i]=(record.getNumGTs(i));
            tempHapCount+=SampleNoHaplotypes[i];
        }
    }
    inFile.close();
    numHaplotypes=tempHapCount;



//    if(outSideRegion + duplicates + notBiallelic + failFilter + inconsistent>0) {cout<<endl;}
    if(MyAllVariable.myOutFormat.verbose)
    {
        std::cout << " NOTE ! "<< numActualRecords << " variants in file ..." << endl;
        if(outSideRegion>0){std::cout << " NOTE ! "<< outSideRegion<< " variants lie outside region ["
                                      <<MyAllVariable.myHapDataVariables.CHR
                                      <<":"<<MyAllVariable.myHapDataVariables.START
                                      <<"-"<<MyAllVariable.myHapDataVariables.END
                                      <<"] ... "<< endl;}
        if(duplicates>0){std::cout << " NOTE ! "<< duplicates<< " instance(s) of duplicated variants discarded ..." << endl;}
        if(notBiallelic>0){std::cout << " NOTE ! "<< notBiallelic<< " multi-allelic variant(s) discarded ..." << endl;}
        if(failFilter>0){std::cout << " NOTE ! "<< failFilter<< " variant(s) failed FILTER = PASS and were discarded ..." << endl;}
        if(inconsistent>0){std::cout << " NOTE ! "<< inconsistent<< " SNP(s) with inconsistent REF/ALT alleles discarded ..." << endl;}
        if(outSideRegion + duplicates + notBiallelic + failFilter + inconsistent>0) {cout<<endl;}
    }
    if(numReadRecords==0)
    {

        if(MyAllVariable.myHapDataVariables.CHR=="")
        {
            cout << "\n ERROR !!! \n No variants left to imported from "<< TypeofFile <<" File "<<endl;
            cout << " Please check the filtering conditions properly ...\n";
        }
        else
        {

            cout << "\n ERROR !!! \n No variants found in region ["<<MyAllVariable.myHapDataVariables.CHR
                 <<":"<<MyAllVariable.myHapDataVariables.START<<"-"<<MyAllVariable.myHapDataVariables.END
                 <<"] in "
                 << TypeofFile << " File : " << VCFFileName << endl;
            cout << " Please check the filtering conditions/input region properly..\n";
        }

        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    cout<<" Successful !!! "<<endl;
    return true;

}



bool HaplotypeSet::ScaffoldGWAStoReference(HaplotypeSet &rHap, AllVariable& MyAllVariable)
{

    int refMarkerCount=(int)rHap.VariantList.size();
    int counter = 0;
    int flag;
	int GWASOnlycounter = 0;
	int OverlapOnlycounter = 0;
	Targetmissing.resize(refMarkerCount, true);
	MapRefToTar.resize(refMarkerCount, -1);
    MapTarToRef.resize(numMarkers, -1);

    TargetMissingTypedOnly.clear();
	knownPosition.resize(numMarkers);
	OverlapOnlyVariantList.resize(numMarkers);
	TypedOnlyVariantList.resize(numMarkers);
	RefAlleleSwap.resize(numMarkers);

	for (int j = 0; j<(int)VariantList.size(); j++)
	{

		int prevCounter = counter;
		flag=0;
		while(counter<refMarkerCount && flag==0 && rHap.VariantList[counter].bp<=VariantList[j].bp)
        {

            if(rHap.VariantList[counter].chr==VariantList[j].chr
             && rHap.VariantList[counter].bp==VariantList[j].bp)
            {
                prevCounter = counter;

                if(rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].altAlleleString)
                    flag=1;
                else if(rHap.VariantList[counter].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else if (VariantList[j].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else if (VariantList[j].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else
                    counter++;
            }
            else
                counter++;
        }


        if(flag==1)
        {
            knownPosition[j]=(counter);
            OverlapOnlyVariantList[OverlapOnlycounter]=(VariantList[j]);

            if(rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString)
                RefAlleleSwap[OverlapOnlycounter]=(false);
            else
                {
                    RefAlleleSwap[OverlapOnlycounter]=(true);
                    string tempa=VariantList[j].refAlleleString;
                    VariantList[j].refAlleleString=VariantList[j].altAlleleString;
                    VariantList[j].altAlleleString=tempa;
                }

            Targetmissing[counter] = false;
            MapRefToTar[counter]=OverlapOnlycounter;
            MapTarToRef[OverlapOnlycounter]=counter;

            OverlapOnlycounter++;
            counter++;
		}
		else
        {
            if(MyAllVariable.myOutFormat.TypedOnly)
            {
                TypedOnlyVariantList[GWASOnlycounter]=(VariantList[j]);
                TargetMissingTypedOnly.push_back(counter-1);
                GWASOnlycounter++;
            }
            knownPosition[j]=(-1);
            counter = prevCounter;
        }

	}

    counter=0;
    rHap.RefTypedIndex.clear();
    int ThisIndex=0;

	while(counter<refMarkerCount && ThisIndex<(int)TargetMissingTypedOnly.size())
    {
        if(counter<=TargetMissingTypedOnly[ThisIndex])
        {
            rHap.RefTypedIndex.push_back(-1);
            counter++;
        }
        else
        {
            rHap.RefTypedIndex.push_back(ThisIndex);
            ThisIndex++;
        }
    }
    while(counter<refMarkerCount)
    {
        rHap.RefTypedIndex.push_back(-1);
        counter++;
    }
    while(ThisIndex<(int)TargetMissingTypedOnly.size())
    {
        rHap.RefTypedIndex.push_back(ThisIndex);
        ThisIndex++;
    }


    rHap.RefTypedTotalCount=GWASOnlycounter+rHap.numMarkers;
    if(rHap.RefTypedTotalCount!=(int)rHap.RefTypedIndex.size())
    {
        cout<<endl<<endl<<" ERROR in Code Construction [ERROR: 007] !!! "<<endl;
        cout<<" Please Contact author with ERROR number urgently : sayantan@umich.edu "<<endl;
        cout<<" Program Exiting ..."<<endl<<endl;
        abort();
    }


    numOverlapMarkers=OverlapOnlycounter;
    numTypedOnlyMarkers=GWASOnlycounter;
    RefAlleleSwap.resize(numOverlapMarkers);
    MapTarToRef.resize(numOverlapMarkers);
    OverlapOnlyVariantList.resize(numOverlapMarkers);
    TypedOnlyVariantList.resize(numTypedOnlyMarkers);


    cout<<"\n Reference Panel   : Found "<<rHap.numSamples<< " samples ("<< rHap.numHaplotypes  <<" haplotypes) and "<< (int)rHap.VariantList.size()<<" variants ..."<<endl;


    cout<<"\n Target/GWAS Panel : Found "<<numSamples<< " samples ("<< numHaplotypes  <<" haplotypes) and "<< (int)VariantList.size()<<" variants ..."<<endl;
    cout<<"                     "<<numOverlapMarkers<<" variants overlap with Reference panel "<<endl;
    cout<<"                     "<< numTypedOnlyMarkers<<" variants imported that exist only in Target/GWAS panel"<<endl;


//    VariantList.clear();
	if (numOverlapMarkers == 0)
	{

		cout << "\n ERROR !!! \n No overlap between Target and Reference markers !!!\n";
		cout << " Please check for consistent marker identifier in reference and target input files..\n";
		cout << "\n Program Exiting ... \n\n";
        return false;

	}
	return true;

}



void HaplotypeSet::GetSummary(IFILE m3vcfxStream)
{
    vector<string> headerTag(2);
    string line;
    const char *equalSep="=";
    m3vcfVERSION=0;

    bool Header=true;

    while(Header)
    {
        line.clear();
        m3vcfxStream->readLine(line);


        if(line.substr(0,2).compare("##")==0)
            Header=true;
        else
            break;

        MyTokenize(headerTag, line.c_str(), equalSep, 2);

        if(headerTag[0].compare("##n_blocks")==0)
        {
            NoBlocks=atoi(headerTag[1].c_str());
        }

        if(m3vcfVERSION==0 && headerTag[0].compare("##fileformat")==0)
        {
            string tempVer = headerTag[1].substr(5,5).c_str();
            if(tempVer.length()==0)
                m3vcfVERSION=1;
            else if (tempVer=="v2.0")
                m3vcfVERSION=2;
            else
            {
                cout << "\n ERROR !!! \n Invalid M3VCF Version  : "<<tempVer<<endl;
                cout << " Please check the file properly..\n";
                cout << "\n Program Exiting ... \n\n";
                abort();
            }
        }

    }

    getm3VCFSampleNames(line);

}






bool HaplotypeSet::ReadM3VCFChunkingInformation(String &Reffilename,string checkChr)
{
    string line;
    int blockIndex=0, NoMarkersImported=0, ReadHeader=0;
    NoLinesToDiscardatBeginning=0;
    finChromosome="NULL";
    StartedThisPanel=false;
    AlreadyReadMiddle=false;
    BlockPiecesforVarInfo.resize(9);

   std::cout << "\n Gathering variant information ..." << endl;

    if(MyHapDataVariables->CHR!="")
    {
        std::cout << "\n Loading markers in region " << MyHapDataVariables->CHR<<":"<<MyHapDataVariables->START<<"-"<<MyHapDataVariables->END<<" ..."<< endl;
    }

    IFILE m3vcfxStream = ifopen(Reffilename, "r");

    if(m3vcfxStream)
    {

        GetSummary(m3vcfxStream);
        if(numSamples==0)
        {
            cout << "\n ERROR !!! \n No samples found in M3VCF Input File  : "<<Reffilename<<endl;
            cout << " Please check the file properly..\n";
            cout << "\n Program Exiting ... \n\n";
            return false;
        }

        ReducedStructureInfoSummary.clear();
        //for(blockIndex=0;blockIndex<NoBlocks;blockIndex++)


        while(m3vcfxStream->readLine(line)!=-1)
        {
            if(ReadHeader==0)
            {
                string tempLine = line.c_str() ;
                UpdatePloidySummary(tempLine);
                ReadHeader=1;
            }

            ReducedHaplotypeInfoSummary tempBlock;
             if(ReadBlockHeaderSummary(line, tempBlock))
            {
                if(finChromosome!=checkChr)
                {
                    cout << "\n ERROR !!! \n Reference Panel is on chromosome ["<<finChromosome<<"] which is ";
                    cout <<" different from chromosome ["<< checkChr<<"] of the GWAS panel  "<<endl;
                    cout << " Please check the file properly..\n";
                    cout << "\n Program Exiting ... \n\n";
                    return false;

                }
                if(!AlreadyReadMiddle)
                    NoLinesToDiscardatBeginning += (tempBlock.BlockSize+1);
                for(int tempIndex=0;tempIndex<tempBlock.BlockSize;tempIndex++)
                    m3vcfxStream->discardLine();
                line.clear();
                continue;
            }

            if(finChromosome!=checkChr)
            {
                cout << "\n ERROR !!! \n Reference Panel is on chromosome ["<<finChromosome<<"] which is ";
                cout <<" different from chromosome ["<< checkChr<<"] of the GWAS panel  "<<endl;
                cout << " Please check the file properly..\n";
                cout << "\n Program Exiting ... \n\n";
                return false;

            }

            AlreadyReadMiddle=true;
            GetVariantInfoFromBlock(m3vcfxStream, tempBlock, NoMarkersImported);
            ReducedStructureInfoSummary.push_back(tempBlock);
            line.clear();
            blockIndex++;
        }



    }
    numMarkers=VariantList.size();
    NoBlocks=ReducedStructureInfoSummary.size();

    if (numMarkers == 0)
    {
        if(MyAllVariables->myHapDataVariables.CHR=="")
        {
            cout << "\n ERROR !!! \n No variants left to imported from reference haplotype file "<<endl;
            cout << " Please check the filtering conditions OR the file properly ..."<<endl;
        }
        else
        {

            cout << "\n ERROR !!! \n No variants found in region ["
                 <<MyAllVariables->myHapDataVariables.CHR<<":"
                 <<MyAllVariables->myHapDataVariables.START<<"-"<<MyAllVariables->myHapDataVariables.END
                 <<"] in reference haplotype file "<<endl;
            cout << " Please check the filtering conditions OR the file properly ..."<<endl;
        }
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

   
    ifclose(m3vcfxStream);

    CreateSiteSummary();

    cout<<"\n Successful !!! "<<endl;

    return true;
}



bool HaplotypeSet::CheckValidChrom(string chr)
{
    bool result=false;

    if(MyAllVariables->myHapDataVariables.MyChromosome!="" && chr==MyAllVariables->myHapDataVariables.MyChromosome.c_str())
        return true;

    string temp[]={"1","2","3","4","5","6","7","8","9","10","11"
            ,"12","13","14","15","16","17","18","19","20","21","22","23","X","Y"
            ,"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11"
            ,"chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20"
            ,"chr21","chr22","chr23","chrX","chrY"};
    std::vector<string> ValidChromList (temp, temp + sizeof(temp) / sizeof(string) );

    for(int counter=0;counter<(int)ValidChromList.size();counter++)
        if(chr==ValidChromList[counter])
            result=true;

    return result;

}

void HaplotypeSet::writem3vcfFile(String filename,bool &gzip)
{

    IFILE m3vcffile = ifopen(filename + ".m3vcf" + (gzip ? ".gz" : ""), "wb",(gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED));
    ifprintf(m3vcffile, "##fileformat=M3VCF\n");
    ifprintf(m3vcffile, "##version=1.2\n");
    ifprintf(m3vcffile, "##compression=block\n");
    ifprintf(m3vcffile, "##n_blocks=%d\n",NoBlocks);
    ifprintf(m3vcffile, "##n_haps=%d\n",numHaplotypes);
    ifprintf(m3vcffile, "##n_markers=%d\n",numMarkers);
    if(finChromosome=="X" || finChromosome=="23")
        ifprintf(m3vcffile, "##chrxRegion=%s\n",PseudoAutosomal?"PseudoAutosomal":"NonPseudoAutosomal");
    ifprintf(m3vcffile, "##<Note=This is NOT a VCF File and cannot be read by vcftools>\n");
    ifprintf(m3vcffile, "#CHROM\t");
    ifprintf(m3vcffile, "POS\t");
    ifprintf(m3vcffile, "ID\t");
    ifprintf(m3vcffile, "REF\t");
    ifprintf(m3vcffile, "ALT\t");
    ifprintf(m3vcffile, "QUAL\t");
    ifprintf(m3vcffile, "FILTER\t");
    ifprintf(m3vcffile, "INFO\t");
    ifprintf(m3vcffile, "FORMAT");
    int i,j,k;

    for(i=0;i<(int)numSamples;i++)
    {
        ifprintf(m3vcffile, "\t%s_HAP_1",individualName[i].c_str());
        if(SampleNoHaplotypes[i]==2)
            ifprintf(m3vcffile, "\t%s_HAP_2",individualName[i].c_str());
    }
    ifprintf(m3vcffile, "\n");

    int length=NoBlocks;
    string cno;

    for(i=0;i<length;i++)
    {

        ReducedHaplotypeInfo &tempInfo = ReducedStructureInfo[i];

        cno=VariantList[tempInfo.startIndex].chr;
        int nvariants=tempInfo.BlockSize;
        int reps=tempInfo.RepSize;


        ifprintf(m3vcffile, "%s\t",cno.c_str());
        ifprintf(m3vcffile, "%d-%d\t",VariantList[tempInfo.startIndex].bp,VariantList[tempInfo.endIndex].bp);
        ifprintf(m3vcffile, "<BLOCK:%d-%d>\t.\t.\t.\t.\t",tempInfo.startIndex,tempInfo.endIndex);

        ifprintf(m3vcffile, "B%d;VARIANTS=%d;REPS=%d\t.",i+1,nvariants,reps);


        for(j=0;j<numHaplotypes;j++)
            ifprintf(m3vcffile, "\t%d",tempInfo.uniqueIndexMap[j]);

        ifprintf(m3vcffile, "\n");

        for(j=0;j<nvariants;j++)
        {
            ifprintf(m3vcffile, "%s\t",cno.c_str());
            ifprintf(m3vcffile, "%d\t",VariantList[j+tempInfo.startIndex].bp);
            ifprintf(m3vcffile, "%s\t",MyAllVariables->myOutFormat.RsId?VariantList[j+tempInfo.startIndex].rsid.c_str():VariantList[j+tempInfo.startIndex].name.c_str());
            ifprintf(m3vcffile, "%s\t%s\t.\t.\t",VariantList[j+tempInfo.startIndex].refAlleleString.c_str(),VariantList[j+tempInfo.startIndex].altAlleleString.c_str());

            ifprintf(m3vcffile, "B%d.M%d",i+1,j+1);
            if(Error.size()>0)
                ifprintf(m3vcffile, ";Err=%.5g;Recom=%.5g",
                     Error[j+tempInfo.startIndex],(j+tempInfo.startIndex)<(int)Recom.size()?Recom[j+tempInfo.startIndex]:0);
            ifprintf(m3vcffile, "\t");

            vector<AlleleType> &TempHap = tempInfo.TransposedUniqueHaps[j];
            for(k=0;k<reps;k++)
            {
                ifprintf(m3vcffile,"%c",TempHap[k]);
            }
            ifprintf(m3vcffile, "\n");
        }
    }

    std::cout << " Successfully written file ... "<<endl;
    ifclose(m3vcffile);

}

string HaplotypeSet::DetectFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");
    string line;

    if(fileStream)
    {

        fileStream->readLine(line);
        if(line.length()<1)
            {
                ifclose(fileStream);
                return "Invalid";
            }
        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }
        if(((string)temp).compare("##fileformat=m3vc")==0)
        {
            ifclose(fileStream);
            return "m3vcf";
        }
        else if(((string)temp).compare("##fileformat=vcfv")==0)
        {
            ifclose(fileStream);
            return "vcf";
        }
        else
        {
            ifclose(fileStream);
            return "Invalid";
        }

    }
    else
    {
        ifclose(fileStream);
        return "NA";
    }

    ifclose(fileStream);
    return "NA";
}


bool HaplotypeSet::readm3vcfFile(String m3vcfFile,String CHR,int START,int END,int WINDOW)
{
    string line,tempString,tempBlockPos,tempString2,tempString3,tempName,tempRsId,tempChr;
    variant tempVariant;
    variant tempVariant2;
    int InitialNMarkers=0,blockIndex,startIndexFlag=0,readIndex=0,writeBlockIndex=0,NoBlocks=0,tempPos,NoMarkersWritten=0,tempVarCount,tempRepCount;
    int OrigStartPos=START;
    int OrigEndPos=END;

    //cout<<"\n Reading Reference Haplotype information from M3VCF files : "<<m3vcfFile<<endl<<endl;
     if (CHR!="" && WINDOW > 0)
    {
        if (START-WINDOW < 0)
            START = 0;
        else
            START -= WINDOW;

        END += WINDOW;
    }
    if(CHR!="")
    {
        std::cout << "\n Loading markers from chromosome " << CHR;
        if(END>0)
            std::cout << " from base position "<<START<<" to base position "<<END<<"."<< endl;
        else
            std::cout << " from base position "<<START<<" till end of M3VCF file."<< endl;
    }


    PrintStartIndex=99999999;
    PrintEndIndex=-1;


    IFILE m3vcfxStream = ifopen(m3vcfFile, "r");

    if(m3vcfxStream)
    {


        m3vcfxStream->readLine(line);
        if(line.compare("##fileformat=M3VCF")!=0 && line.compare("##fileformat=OPTM")!=0)
        {
            cout<<" Incorrect Header Information : "<<line<<endl;
            cout<<" Header line should be : ##fileformat=M3VCF "<<endl<<endl;
            printErr(m3vcfFile);
        }

        bool Header=true;
//        char * pch_split,* pch_split2,* pch_split3;
        char * pch;
        char *end_str1,*end_str2;

        while(Header)
        {
            line.clear();
            m3vcfxStream->readLine(line);
            if(line.substr(0,2).compare("##")==0)
                Header=true;
            else
                break;
            tempString = (string) strtok_r ((char*)line.c_str(),"=", &end_str1);
            pch = strtok_r (NULL, "=", &end_str1);

            if(tempString.compare("##n_blocks")==0)
            {
                NoBlocks=atoi(pch);
                continue;
            }
            else if(tempString.compare("##n_haps")==0)
            {
                numHaplotypes=atoi(pch);
                continue;
            }
            else if(tempString.compare("##n_markers")==0)
            {
                InitialNMarkers=atoi(pch);
//                abort();
                continue;

            }
            else if(tempString.compare("##chrxRegion")==0)
            {
                tempString = (string) pch;
                if(tempString.compare("NonPseudoAutosomal")==0)
                    PseudoAutosomal=false;
                else if(tempString.compare("PseudoAutosomal")==0)
                    PseudoAutosomal=true;
                else
                {
                    cout << "\n Inconsistent Tag for Chr X. "<<endl;
                    cout << " Please check the file properly..\n";
                    cout << " Program Aborting ... "<<endl;
                    return false;

                }
                continue;

            }

        }

        int colCount=0;
        pch = strtok_r ((char*)line.c_str(),"\t", &end_str2);


        while(pch!=NULL)
        {
            colCount++;
            if(colCount>9)
            {
                tempString2=string(pch);
                ///////////////////////////
//                if(colCount%2==0)
                individualName.push_back(tempString2.substr(0,tempString2.size()-6));
            }
            pch = strtok_r (NULL,"\t", &end_str2);
        }


        cout<<" Reading  "<<numHaplotypes<< " haplotypes from data ..."<<endl<<endl;
        ///////////////////////////

        if((int)individualName.size()!=numHaplotypes)
        {
            cout<<endl<<" Error in Data consistency !!! "<<endl<<endl;
            cout<<" Number of Haplotypes should be : "<< numHaplotypes<<", but "<<individualName.size()<<" Haplotypes found in header row !!!"<< endl<<endl;
            printErr(m3vcfFile);
        }


        if(individualName.size()==0)
        {
            cout << "\n No haplotypes recorded from M3VCF Input File : "<<m3vcfFile<<endl;
            cout << " Please check the file properly..\n";
            cout << " Program Aborting ... "<<endl;
            return false;
        }



        for(blockIndex=0;blockIndex<NoBlocks;blockIndex++)
        {


            if (blockIndex % 1000 == 0)
            {
            //   printf("  Loading Block %d out of %d blocks to be loaded... [%.1f%%] "
            //            , blockIndex + 1, NoBlocks, 100*(double)(blockIndex + 1)/(double)NoBlocks);
            //    cout<<endl;
            }

            int flag=0,blockEnterFlag=0,blocktoSave=0;
            line.clear();
            m3vcfxStream->readLine(line);
            char *end_str_new;

            pch = strtok_r ((char*)line.c_str(),"\t", &end_str_new);
            tempChr=string(pch);

            pch = strtok_r (NULL,"\t", &end_str_new);
            tempBlockPos=string(pch);

            char *end_str_new1;
            char *pch1;

            pch1 = strtok_r ((char*)tempBlockPos.c_str(),"-", &end_str_new1);
            int tempStartBlock=atoi(pch1);

            pch1 = strtok_r (NULL,"-", &end_str_new1);
            int tempEndBlock=atoi(pch1);

            if(CHR!="")
            {
                if(tempChr.compare(CHR.c_str())!=0)
                    flag=1;
                else
                {
                    if(END>0)
                    {
                        if(tempStartBlock>END)
                            flag=1;
                    }
                    if(tempEndBlock<START)
                        flag=1;
                }
            }


            pch = strtok_r (NULL,"\t", &end_str_new);
            string blockName=string(pch);

            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);

            pch = strtok_r (NULL,"\t", &end_str_new);
            tempString2=string(pch);

            char *pch2,*end_str_new2;
            pch2 = strtok_r ((char*)tempString2.c_str(),";", &end_str_new2);

            stringstream strs;
            strs<<(blockIndex+1);

            tempString3="B"+ (string)(strs.str());

            if(tempString3.compare((string)pch2)!=0)
            {
                cout<<endl<<" Error in INFO column (Block Identifier) for block : "<<blockName <<endl;
                cout<<" Block Identifier should be : "<< tempString3<<" but is : "<<pch2<<endl<<endl;
                printErr(m3vcfFile);
            }



            tempString3 = (string)strtok_r (NULL,";", &end_str_new2);
            char *pch3,*end_str_new3;
            pch3 = strtok_r ((char*)tempString3.c_str(),"=", &end_str_new3);
            pch3 = strtok_r (NULL,"=", &end_str_new3);
            tempVarCount=atoi(pch3);

            tempString3 = (string)strtok_r (NULL,";", &end_str_new2);
            end_str_new3=NULL;
            pch3 = strtok_r ((char*)tempString3.c_str(),"=", &end_str_new3);
            pch3 = strtok_r (NULL,"=", &end_str_new3);
            tempRepCount=atoi(pch3);

            ReducedHaplotypeInfo tempBlock;

            if(flag==1)
            {
                for(int tempIndex=0;tempIndex<tempVarCount;tempIndex++)
                    m3vcfxStream->discardLine();
                continue;
            }

            tempBlock.TransposedUniqueHaps.resize(tempVarCount);

            tempBlock.RepSize=tempRepCount;
            tempBlock.BlockSize=tempVarCount;

            tempBlock.uniqueCardinality.resize(tempRepCount,0.0);
            tempBlock.InvuniqueCardinality.resize(tempRepCount,0.0);

//            vector<AlleleType> &TempHap = tempBlock.TransposedUniqueHaps[tempIndex];

    //tempBlock.TransposedUniqueHaps  uniqueHaps.resize(tempRepCount);
            tempBlock.uniqueIndexMap.resize(numHaplotypes);

            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);



            int check=0;
            while(pch!=NULL)
            {
                int tempval=atoi(pch);
                tempBlock.uniqueIndexMap[check]=tempval;
                tempBlock.uniqueCardinality[tempval]++;

                check++;
                pch = strtok_r (NULL,"\t", &end_str_new);
            }

            for (int i = 0; i < tempRepCount; i++)
            {
                tempBlock.InvuniqueCardinality[i]=1.0/(float)tempBlock.uniqueCardinality[i];
            }



            if(check!=numHaplotypes)
            {
                cout<<endl<<" Error in Data consistency !!! "<<endl;
                cout<<" Number of Haplotypes should be : "<< numHaplotypes<<", but "<<check<<" indices found in block"<<blockName << endl<<endl;
                printErr(m3vcfFile);
            }




            for(int tempIndex=0;tempIndex<tempVarCount;tempIndex++)
            {
                flag=0;
                line.clear();
                m3vcfxStream->readLine(line);

                end_str_new3=NULL;


                pch3 = strtok_r ((char*)line.c_str(),"\t", &end_str_new3);
                tempChr=(string)pch3;

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                tempPos=atoi(pch3);

                stringstream strs3;
                strs3<<tempPos;


                tempName=tempChr+":"+(string)strs3.str();

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                tempRsId=(string)(pch3);



                if(CHR!="")
                {
                    if(tempChr.compare(CHR.c_str())!=0)
                        flag=1;
                    else
                    {
                        if(END>0)
                        {
                            if(tempPos>END || tempPos<START)
                                flag=1;
                        }
                        else
                        if(tempPos<START)
                            flag=1;
                        if(tempIndex==(tempVarCount-1) && writeBlockIndex==0)
                            flag=1;
                    }

                }
                readIndex++;

                if(flag==1)
                    continue;
                else
                {
                    if(blockEnterFlag==0)
                        tempBlock.startIndex=writeBlockIndex;
                    else
                        tempBlock.endIndex=writeBlockIndex;
                    if(tempIndex<(tempVarCount-1) && tempIndex>0) // to ensure a single marker from a block is NOT read.
                        blocktoSave=1;
                    blockEnterFlag=1;

                }

                if(tempPos>=OrigStartPos && startIndexFlag==0)
                {
                    PrintStartIndex=writeBlockIndex;
                    startIndexFlag=1;
                }

//                if(CHR=="" || tempPos<=OrigEndPos)
//                    PrintEndIndex=writeBlockIndex;
//


                if(CHR=="")
                {
                    PrintEndIndex=writeBlockIndex;
                }
                else
                {
                    if(OrigEndPos==0 || tempPos<=OrigEndPos)
                        PrintEndIndex=writeBlockIndex;
                }


                variant tempVariant;

                tempVariant.assignValues(tempName,tempName,tempChr,tempPos);
                tempVariant.rsid=tempRsId;

                double tempRecom=-3.0,tempError=0.0;

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                string tempString98=(string)(pch3);

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                string tempString99=(string)(pch3);

                tempVariant.assignRefAlt(tempString98,tempString99);

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                pch3 = strtok_r (NULL,"\t", &end_str_new3);

                tempString=(string)(pch3);

                char *pch4, *end_str_new4;

                pch4=strtok_r ((char*)tempString.c_str(),";", &end_str_new4);
                stringstream strs1,strs2;
                strs1<<(blockIndex+1);
                strs2<<(tempIndex+1);
                tempString3="B"+ (string)(strs1.str())+ ".M"+ (string)(strs2.str());
                if(tempString3.compare((string)pch4)!=0)
                {
                    cout<<endl<<" Error in INFO column (Block Identifier) for variant : "<<tempName <<endl;
                    cout<<" Block Identifier should be : "<< tempString3<<" but is : "<<pch4<<endl<<endl;
                    printErr(m3vcfFile);
                }


                pch4=strtok_r (NULL,";", &end_str_new4);

                while(pch4!=NULL)
                {
                    tempString3=(string)(pch4);
                    char *pch5,*pch6,*end_str_new5;
                    pch5=strtok_r ((char*)tempString3.c_str(),"=", &end_str_new5);
                    pch6=strtok_r (NULL,"=", &end_str_new5);

                    if((string)pch5=="Err")
                    {
                        tempError=atof(pch6);
                    }
                    if((string)pch5=="Recom")
                    {
                        tempRecom=atof(pch6);
                    }
                    pch4=strtok_r (NULL,";", &end_str_new4);
                }

                char Rallele=0,Aallele=0;


                string refAllele=tempVariant.refAlleleString;
                string altAllele=tempVariant.altAlleleString;

                if (strlen(refAllele.c_str()) == 1
                    && strlen(altAllele.c_str()) == 1)
                {
                    switch (refAllele[0])
                    {
                        case 'A': case 'a': Rallele = 1; break;
                        case 'C': case 'c': Rallele = 2; break;
                        case 'G': case 'g': Rallele = 3; break;
                        case 'T': case 't': Rallele = 4; break;
                        case 'D': case 'd': Rallele = 5; break;
                        case 'I': case 'i': Rallele = 6; break;
                        case 'R': case 'r': Rallele = 7; break;
                        default:
                        {
                            cout << "\n\n Data Inconsistency !!! \n";
                            cout << " Error with reference allele for marker : " << tempVariant.rsid<< " in M3VCF File : " << m3vcfFile;
                            cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                            cout << " " << tempVariant.rsid << " has " << refAllele << endl;
                            cout << "\n Program Aborting ... \n\n";
                            abort();
//                                       return false;
                        }
                    }

                    switch (altAllele[0])
                    {
                        case 'A': case 'a': Aallele = 1; break;
                        case 'C': case 'c': Aallele = 2; break;
                        case 'G': case 'g': Aallele = 3; break;
                        case 'T': case 't': Aallele = 4; break;
                        case 'D': case 'd': Aallele = 5; break;
                        case 'I': case 'i': Aallele = 6; break;
                        case 'R': case 'r': Aallele = 7; break;
                        default:
                        {
                            cout << "\n\n Data Inconsistency !!! \n";
                            cout << " Error with alternate allele for marker : " <<tempVariant.rsid<< " in M3VCF File : " << m3vcfFile;
                            cout << "\n VCF alternate alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                            cout << " " << tempVariant.rsid << " has " << altAllele << endl;
                            cout << "\n Program Aborting ... \n\n";
                            abort();
                        }
                    }
                }
                else
                {
                    Rallele = 7;
                    if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
                        Aallele=6;
                    else
                        Aallele=5;
                }
                //tempVariant.refAllele=Rallele;
                //tempVariant.altAllele=Aallele;



                if(tempIndex<(tempVarCount-1) || blockIndex==(NoBlocks-1))
                {

                    VariantList.push_back(tempVariant);
                    //refAlleleList.push_back(tempVariant.refAllele);
                    //markerName.push_back(tempName);
//                    if(tempRecom!=-3.0)
//                    {
//                        Recom.push_back(tempRecom);
//                        Error.push_back(tempError);
//                    }
                    writeBlockIndex++;
                }

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                tempString=(string)(pch3);
                vector<AlleleType> &TempHap = tempBlock.TransposedUniqueHaps[tempIndex];
                TempHap.resize(tempRepCount);

                for(check=0;check<tempRepCount;check++)
                {

                    char testVal=tempString[check];

                    TempHap[check] = testVal;

                    //tempBlock.uniqueHaps[check].push_back(t=='0'? false:true);

                }

            }

            NoMarkersWritten+=(tempVarCount-1);
            if(blocktoSave==1)
            {
                //optEndPoints.push_back(tempBlock.startIndex);

                ReducedStructureInfo.push_back(tempBlock);
            }



        }
//        if(ReducedStructureInfo.size()>0)
//            optEndPoints.push_back(ReducedStructureInfo[ReducedStructureInfo.size()-1].endIndex);

        if(Recom.size()>0)
            Recom.pop_back();

        finChromosome=tempChr;

    }
    else
    {
        cout<<" Following M3VCF File Not Available : "<<m3vcfFile<<endl;
        cout<<" Program Exiting ... "<<endl<<endl;
        abort();
    }


    numMarkers=writeBlockIndex;

   // cout<<endl<<" Reference Haplotype information succesfully recorded. "<<endl;
    CreateAfterUncompressSummary();
    InvertUniqueIndexMap();
    CalculateAlleleFreq();

//    if(finChromosome=="X")
//    {
//        cout<<"\n Chromosome X Detected !!! \n";
//    }

//    std::cout << "\n Number of Markers in File                           : " << InitialNMarkers << endl;
//    std::cout << "\n Number of Markers Recorded                          : " << numMarkers << endl;
//    std::cout << " Number of Haplotypes Recorded                       : " << numHaplotypes << endl;

    ifclose(m3vcfxStream);

    if(numMarkers<2)
    {
        cout << "\n None/Single marker left after filtering from Input File : "<<m3vcfFile<<endl;
        cout << " Please check the file or the filtering options properly ...\n";
        cout << " Program Aborting ... "<<endl;
        return false;
    }



//    std::cout << "\n Haplotype Set successfully loaded from M3VCF File     : " << m3vcfFile << endl;

    return true;




}



void HaplotypeSet::CalculateAlleleFreq()
{
    fill(AlleleFreq.begin(), AlleleFreq.end(), 0.0);

    int i,j,k;
    for(k=0;k<NoBlocks;k++)
    {

        ReducedHaplotypeInfo &TempBlock=ReducedStructureInfo[k];
        for(j=TempBlock.startIndex;j<TempBlock.endIndex + (k==(NoBlocks-1)? 1:0) ;j++)
        {
            vector<AlleleType> &TempHap = TempBlock.TransposedUniqueHaps[j-TempBlock.startIndex];
            for (i = 0; i<TempBlock.RepSize; i++)
            {
                if(TempHap[i]=='1')
                {
                    AlleleFreq[j]+=TempBlock.uniqueCardinality[i];
                }
            }
        }
    }

	for (int i = 0; i<numMarkers; i++)
	{
		AlleleFreq[i] /= (double)numHaplotypes;
	}

}



void HaplotypeSet::CalculateGWASOnlyAlleleFreq()
{

    fill(GWASOnlyAlleleFreq.begin(), GWASOnlyAlleleFreq.end(), 0.0);
    vector<int> TotalSample(numTypedOnlyMarkers, 0);

    int i,K,j;
    int haplotypeIndex=0;
    for(K=0;K<numSamples;K++)
    {
        for(j=0; j<(*SampleNoHaplotypesPointer)[K];j++)
        {
            for (i = 0; i<numTypedOnlyMarkers; i++)
            {
                if(GWASOnlyMissingSampleUnscaffolded[haplotypeIndex][i]=='0')
                {
                    TotalSample[i]++;
                    if(GWASOnlyhaplotypesUnscaffolded[haplotypeIndex][i]=='1')
                        GWASOnlyAlleleFreq[i]++;
                }
            }
            haplotypeIndex++;
        }
    }

    for (int i = 0; i<numTypedOnlyMarkers; i++)
	{
		GWASOnlyAlleleFreq[i]/=(double)TotalSample[i];
	}
}


AlleleType HaplotypeSet::RetrieveMissingScaffoldedHaplotype(int sample,int marker)
{
    return MissingSampleUnscaffolded[sample][marker];
}

AlleleType HaplotypeSet::RetrieveScaffoldedHaplotype(int sample,int marker)
{
    return haplotypesUnscaffolded[sample][marker];
}



void HaplotypeSet::MyTokenize(vector<string> &result, const char *input, const char *delimiter, int Number)
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

string HaplotypeSet::FindTokenWithPrefix(const char *input,const char *delimiter, string CheckPrefix)
{

    std::string word = "";
    int Size = (int)CheckPrefix.size();
    int FirstChar = 0;
    while (*input)
    {
        if ( FirstChar==0 || *input==*delimiter)
        {
            if(FirstChar==0)
                FirstChar=1;
            else
                ++input;

            int Index=0;
            while(*input)
            {
                word=word + (*input);
                if(Index<Size && *input!=CheckPrefix[Index++])
                    break;

                ++input;
                if(*input==*delimiter || *input=='\0')
                    return word;
            }
        }
        input++;
        word="";
    }
    return word;

}

int HaplotypeSet::CheckBlockPosFlag(string &input, String &CHR, int &START, int &END)
{

    const char *dashSep="-", *tabSep="\t";
    vector<string> tokens(2);
    vector<string> PosTokens(2);

    MyTokenize(tokens, input.c_str(), tabSep,2);
    string tempChr=tokens[0];

    MyTokenize(PosTokens, tokens[1].c_str(), dashSep, 2);
    int tempStartBlock=atoi(PosTokens[0].c_str());
    int tempEndBlock=atoi(PosTokens[1].c_str());
    if(finChromosome=="NULL")
        finChromosome=tempChr;

    if(CHR!="")
    {
        if(tempChr.compare(CHR.c_str())!=0)
            return 1;
        else
        {
            if(tempStartBlock>END)
                    return 1;
            if(tempEndBlock<START)
                return 1;
        }
    }


    return 0;

}


int HaplotypeSet::GetNumVariants(string &input)
{
    const char *equalSep="=",*semicolSep=";";


    string PrefixString = FindTokenWithPrefix(input.c_str(),semicolSep, "VARIANTS=");
    if(PrefixString!="")
    {
        vector<string> NoVariantsToken(2);
        MyTokenize(NoVariantsToken, PrefixString.c_str(), equalSep, 2);
        return atof(NoVariantsToken[1].c_str());
    }
    else
    {
        abort();
        cout<<endl<<endl<<" ERROR !!! \n ERROR in M3VCF File [ERROR Code : 3278] !!! "<<endl;
        cout<<" Please Contact author with ERROR number : sayantan@umich.edu "<<endl;
        cout<<" Program Exiting ..."<<endl<<endl;
        abort();
    }
}

int HaplotypeSet::GetNumReps(string &input)
{
    const char *equalSep="=",*semicolSep=";";

    string PrefixString = FindTokenWithPrefix(input.c_str(),semicolSep, "REPS=");
    if(PrefixString!="")
    {
        vector<string> NoVariantsToken(2);
        MyTokenize(NoVariantsToken, PrefixString.c_str(), equalSep, 2);
        return atof(NoVariantsToken[1].c_str());
    }
    else
    {
        cout<<endl<<endl<<" ERROR !!! \n ERROR in M3VCF File [ERROR Code : 1472] !!! "<<endl;
        cout<<" Please Contact author with ERROR number : sayantan@umich.edu "<<endl;
        cout<<" Program Exiting ..."<<endl<<endl;
        abort();
    }
}

double HaplotypeSet::GetRecom(string &input)
{

    const char *equalSep="=",*semicolSep=";";


    string PrefixString = FindTokenWithPrefix(input.c_str(),semicolSep, "Recom=");
    if(PrefixString!="")
    {
        vector<string> NoVariantsToken(2);
        MyTokenize(NoVariantsToken, PrefixString.c_str(), equalSep, 2);
        return atof(NoVariantsToken[1].c_str());
    }
    else
        return -3.0;

}


double HaplotypeSet::GetError(string &input)
{
    const char *equalSep="=",*semicolSep=";";


    string PrefixString = FindTokenWithPrefix(input.c_str(),semicolSep, "Err=");
    if(PrefixString!="")
    {
        vector<string> NoVariantsToken(2);
        MyTokenize(NoVariantsToken, PrefixString.c_str(), equalSep, 2);
        return atof(NoVariantsToken[1].c_str());
    }
    else
        return -3.0;

}




