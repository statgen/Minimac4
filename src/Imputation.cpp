#include "Imputation.h"



void Imputation::ImputeThisChunk(int ChunkId, HaplotypeSet &FullrHap, HaplotypeSet &tgwasHap, HaplotypeSet &rgwasHap)
{
    ChunkNo=ChunkId;
    THapUnchunked=&tgwasHap;
    rHapChunked=&FullrHap;

    int time_prev = time(0);

    MarkovParameters *MP=new MarkovParameters(rgwasHap.numMarkers);
    MP->Recom=rgwasHap.Recom;
    MP->Error=rgwasHap.Error;
    stats->Initialize(FullrHap.numMarkers, tgwasHap.numMarkers);

    cout << "\n Starting Imputation ..."<<endl<<endl;

    TotalNovcfParts=0;

    int maxVcfSample= MyAllVariables->myOutFormat.vcfBuffer,NumVcfWritten=0,NumVcfToBeWritten=0;
    if((maxVcfSample)>=tgwasHap.numSamples)
        maxVcfSample=tgwasHap.numSamples;
    int TotalNumSamples=tgwasHap.numSamples;

    SinglePartialDosageData.ReParameterizePartialDosageData(ChunkNo, FullrHap, tgwasHap);
    CurrentPartialDosageData=&SinglePartialDosageData;

    CurrentPartialDosageData->UpdatePartialDosageData(maxVcfSample, NumVcfToBeWritten);

    int StartSamId=0, EndSamId;
    TimeToWrite=0;



    for(int batchNo=0; ; batchNo++)
    {
        EndSamId = StartSamId + (maxVcfSample) <  TotalNumSamples ? StartSamId + (maxVcfSample) : TotalNumSamples;

        printf("  Imputing Samples %d-%d [%0.0f%%] out of %d samples ...", StartSamId + 1, EndSamId, 100*(float)EndSamId/TotalNumSamples, TotalNumSamples);
        cout<<endl;


        #pragma omp parallel for
        for(int SampleId=StartSamId;SampleId<EndSamId;SampleId ++)
        {

            int MMIndexId=0;

            pair <int, int> SwapDosageData;
            int hapId= tgwasHap.CummulativeSampleNoHaplotypes[SampleId];
            //(2*SampleId);

            vector<float> PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,recomProb;
            DosageData *ThisSamplePartialDosageData;

            #pragma omp critical (BindSample)
            {
                 #ifdef _OPENMP
                    MMIndexId = omp_get_thread_num();
                #endif
                ThisSamplePartialDosageData=CurrentPartialDosageData;
                SwapDosageData=CurrentPartialDosageData->IndexSample(SampleId);
            }


            MarkovModel &MM=MainMarkovModel[MMIndexId];
            MM.CopyParametersNew(MP);

            int hapIdIndiv=hapId;
            do{

                MM.ThisHapId=hapIdIndiv;
                MM.ReinitializeMatrices();
                ThisSamplePartialDosageData->BindSampleMModel(MM,SwapDosageData.second,hapIdIndiv%2);

                if(!tgwasHap.RetrieveMissingScaffoldedHaplotype(hapIdIndiv,0))
                {
                    ConditionJunctionProb(rgwasHap,0,MM.junctionLeftProb[0],MM.Error[0],tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,0)? rgwasHap.AlleleFreq[0] : 1-rgwasHap.AlleleFreq[0],
                                          tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,0),MM.backgroundError,rgwasHap.ReducedStructureInfo[0]);
                }

                for(int group=1;group<=rgwasHap.NoBlocks;group++)
                {
                    LeftTraverse(rgwasHap,tgwasHap,hapIdIndiv,MM,group,recomProb,rgwasHap);
                }

                for(int group=rgwasHap.NoBlocks;group>0;group--)
                {
                    ImputeTraverse(rgwasHap, hapIdIndiv,MM,group, recomProb,PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb);
                }

                #pragma omp critical (StatUpdate)
                {
                    stats->NewUpdate(FullrHap, tgwasHap, hapIdIndiv, MM.DosageHap, MM.LooDosageHap);
                }


                if(tgwasHap.SampleNoHaplotypes[SampleId]==1)
                    break;
                hapIdIndiv++;
            }while(hapIdIndiv%2==1);

            if (MyAllVariables->myOutFormat.hapOutput && !MyAllVariables->myOutFormat.unphasedOutput)
            {
                #pragma omp critical (PrintHapData)
                PrintHaplotypeData( hapIdIndiv, SampleId, MM.DosageHap);
            }

            if(MyAllVariables->myOutFormat.doseOutput)
            {
                #pragma omp critical  (PrintDoseData)
                PrintDosageData(SampleId, ThisSamplePartialDosageData->hapDosage[(2*SwapDosageData.second)], ThisSamplePartialDosageData->hapDosage[(2*SwapDosageData.second)+1] );
            }

            if(MyAllVariables->myOutFormat.vcfOutput)
            {
                if(SwapDosageData.first==1)
                {
                    #pragma omp critical (printVCF)
                    {
                        TotalNovcfParts++;
                        printf("       Saving Samples in temporary VCF file ... ");
                        cout<<endl;
                        ThisSamplePartialDosageData->FlushPartialVcf(TotalNovcfParts);

                        TimeToWrite+=ThisSamplePartialDosageData->TimeToWrite;
                        NumVcfToBeWritten = NumVcfWritten + maxVcfSample;
                        ThisSamplePartialDosageData->UpdatePartialDosageData(maxVcfSample < TotalNumSamples-NumVcfToBeWritten ? maxVcfSample : TotalNumSamples-NumVcfToBeWritten, NumVcfToBeWritten);
                        NumVcfWritten+=maxVcfSample;
                    }
                }
            }



        }

        StartSamId=EndSamId ;

        if(StartSamId>=TotalNumSamples)
            break;
    }

    delete MP;



    TimeToImpute = time(0) - time_prev;
    cout <<endl<< " Imputation successful (" << TimeToImpute << " seconds) !!!"<<endl;
    TimeToImpute-=(TimeToWrite);

}


void Imputation::ImputeTraverse(HaplotypeSet &rHap,
                                int hapID,
                              MarkovModel &MM,int group, vector<float> &recomProb,
                              vector<float> &PrevRightFoldedProb,
                              vector<float> &CurrentRightProb,
                              vector<float> &CurrentNoRecoRightProb)
{

    PrevRightFoldedProb.resize(rHap.ReducedStructureInfo[group-1].RepSize);

    MM.foldProbabilities(PrevRightFoldedProb,group-1,
                            rHap.ReducedStructureInfo[group-1],1,rHap.numHaplotypes);


    MM.ImputeSites(hapID,group-1,
                      PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb);

    splitFoldedProb(recomProb,CurrentRightProb,CurrentNoRecoRightProb);
    MM.unfoldProbabilities(group-1,recomProb,CurrentNoRecoRightProb,
                           PrevRightFoldedProb,1,
                           rHap.ReducedStructureInfo,rHap.numHaplotypes);

}


void Imputation::LeftTraverse(HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                              MarkovModel &MM,int group, vector<float> &recomProb,
                              HaplotypeSet &rHapMain)
{
    vector<int> &optStructure=rHapMain.optEndPoints;

    int Start=optStructure[group-1];
    int End=optStructure[group];
    vector<vector<float> > &ThisBlockLeftProb= MM.leftProb[group-1];

    MM.foldProbabilities(ThisBlockLeftProb[0],group-1,rHap.ReducedStructureInfo[group-1],
                            0,rHap.numHaplotypes);
    MM.CurrentLeftNoRecoProb=ThisBlockLeftProb[0];

    MM.WalkLeft(tHap,hapID,group-1,
                   rHap.ReducedStructureInfo[group-1],rHapMain.AlleleFreq);

    splitFoldedProb(recomProb,ThisBlockLeftProb[End-Start],MM.CurrentLeftNoRecoProb);
    MM.unfoldProbabilities(group-1,recomProb,MM.CurrentLeftNoRecoProb,
                                       ThisBlockLeftProb[0],0,
                                       rHap.ReducedStructureInfo,rHap.numHaplotypes);
}





void Imputation::InitializeOutputFiles(HaplotypeSet &tarInitializer, int maxSample, int maxRefVar, int maxTarVar)
{
    SinglePartialDosageData.InitializePartialDosageData(tarInitializer, maxSample, maxRefVar, maxTarVar,  MyAllVariables);
//    TwoPartialDosageData.resize(2);
//    TwoPartialDosageData[0].InitializePartialDosageData(tarInitializer, maxSample, maxVar, MyAllVariables);
//    TwoPartialDosageData[1].InitializePartialDosageData(tarInitializer, maxSample, maxVar, MyAllVariables);
}


void Imputation::PrintDosageData(int ThisSampleId, vector<float> &ThisDosage1,vector<float> &ThisDosage2)
{

//    printf("    Outputting Individual %s for Dosage file...",  THapUnchunked->individualName[ThisSampleId].c_str());
//    cout<<endl;
//    ifprintf(dosages, "%s\tDOSE",THapUnchunked->individualName[ThisSampleId].c_str());
//    int i=0;
//
//    vector<bool> *GWASDosage1, *GWASDosage2;
//
//    if(MyAllVariables.myOutFormat.TypedOnly)
//    {
//        GWASDosage1 = &(THapUnchunked->GWASOnlyMissingSampleUnscaffolded[2*ThisSampleId]);
//        GWASDosage2 = THapUnchunked->GWASOnlyMissingSampleUnscaffolded[2*ThisSampleId + 1];
//    }
//
//    int NoHaps=THapUnchunked->SampleNoHaplotypes[ThisSampleId];
//
//    for (int index =0; index < rHapChunked->RefTypedTotalCount; index++)
//    {
//
//        if(rHapChunked->RefTypedIndex[index]==-1)
//        {
//
//            if(i>=rHapChunked->PrintStartIndex && i <= rHapChunked->PrintEndIndex)
//            {
//                if(NoHaps==1)
//                    ifprintf(dosages, "\t%.3f", ThisDosage1[i]+ThisDosage2[i]);
//                 else
//                    ifprintf(dosages, "\t%.3f", ThisDosage1[i]);
//            }
//            i++;
//
//        }
//        else
//        {
//            int MarkerIndex=rHapChunked->RefTypedIndex[index];
//
//            if(MarkerIndex>=THapUnchunked->PrintTypedOnlyStartIndex && MappingIndex<=THapUnchunked->PrintTypedOnlyEndIndex)
//            {
//                if(NoHaps==1)
//                {
//                    ifprintf(dosages, "\t%.3f", (float)GWASDosage1[MarkerIndex]);
//                }
//                else
//                {
//                     ifprintf(dosages, "\t%.3f", (float)((*GWASDosage1)[MarkerIndex]) + (float)((*GWASDosage2)[MarkerIndex]));
//                }
//            }
//        }
//
//        ifprintf(dosages,"\n");
//    }

}



void Imputation::PrintHaplotypeData(int ThisHapId, int ThisSampleId, vector<float> *ThisimputedHap)
{

    printf("    Outputting HAPLO%d of Individual %s for Haplotype File...",
           ThisHapId%2+1 ,THapUnchunked->individualName[ThisSampleId].c_str());
    cout<<endl;

    ifprintf(hapdose, "%s\tHAPLO%d",  THapUnchunked->individualName[ThisSampleId].c_str(), ThisHapId%2+1  );
    int i=0;
    vector<bool> &GWASDosage1 = THapUnchunked->GWASOnlyMissingSampleUnscaffolded[ThisHapId];

    for (int index =0; index < rHapChunked->RefTypedTotalCount; index++)
    {
        if(rHapChunked->RefTypedIndex[index]==-1)
        {
            if(i>=rHapChunked->PrintStartIndex && i <= rHapChunked->PrintEndIndex)
            {
                ifprintf(hapdose, "\t%.5f", (*ThisimputedHap)[i]);
            }
            i++;
        }
        else
        {
            int MarkerIndex=rHapChunked->RefTypedIndex[index];
            ifprintf(hapdose, "\t%.5f", (float)GWASDosage1[MarkerIndex]);
        }
    }

    ifprintf(hapdose, "\n");

}


void Imputation::splitFoldedProb(vector<float> &SplitProb, vector<float> &totalProb, vector<float> &noRecomProb)
{
    SplitProb.resize(totalProb.size());

    for(int i=0;i<(int)totalProb.size();i++)
    {
        SplitProb[i]=totalProb[i]-noRecomProb[i];
    }
}




void Imputation::Condition(HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                            vector<float> &noRecomProb,
                            double e, double freq, bool observed, double backgroundError, int NoRedStates, ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    vector<bool> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];

    for (int i = 0; i<NoRedStates; i++)
    {

        bool allele= TempHap[i];
        if(allele==observed)
        {
            Prob[i]*=pmatch;
            noRecomProb[i]*=pmatch;
        }
        else
        {
            Prob[i]*=prandom;
            noRecomProb[i]*=prandom;
        }
    }
}



void Imputation::ConditionJunctionProb(HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                            double e, double freq, bool observed, double backgroundError,
                            ReducedHaplotypeInfo &Info)
{
    vector<bool> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];


    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;
    int NoStates=rHap.numHaplotypes;

    for (int i = 0; i<NoStates; i++)
    {
        Prob[i]*= TempHap[Info.uniqueIndexMap[i]] ==observed ? pmatch : prandom;
    }
}





