#include "Imputation.h"


void Imputation::ParameterEstimateThisChunk(int ChunkId, HaplotypeSet &FullrHap, HaplotypeSet &rHap_loo)
{
    ChunkNo=ChunkId;
    rHapChunked=&FullrHap;

    int time_prev = time(0);
    int EstimationStates = MyAllVariables->myModelVariables.states;
    int EstimationRounds = MyAllVariables->myModelVariables.rounds;
    MarkovParameters *MP=new MarkovParameters(FullrHap.numMarkers);


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             PARAMETER ESTIMATION                              "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    if(FullrHap.Recom.size()>0)
    {
        printf(" Reading pre-calculated estimates located in M3VCF file ...\n");
        cout<<endl;
        MP->Recom=FullrHap.Recom;
        MP->Error=FullrHap.Error;
    }


    printf(" Initializing Model Parameters (using EM and up to %d haplotypes) ...", EstimationStates);
    cout<<endl;


    double LogLikelihoodValue=0.0;
    for (int round = 0; round < EstimationRounds; round++)
    {
        cout<<"\n Round "<<round+1<<" of Parameter Refinement ..."<<endl;
        int iterations = EstimationStates < FullrHap.numHaplotypes ? EstimationStates :  FullrHap.numHaplotypes;
        LogLikelihoodValue=0.0;
        #pragma omp parallel for
        for (int i = 0; i < iterations; i++)
        {
            cout<<i<<endl;
            // Reference leave one out (loo) panel
            LooOptimalStructure(i,FullrHap,rHap_loo);

            // Target with one panel at position i
            HaplotypeSet tHap_loo;
            tHap_loo.Create(i,FullrHap);
//
//                // Create Markov Model
//                vector<float> PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,recomProb;
//                PrevRightFoldedProb.reserve(rHap.maxRepSize);
//                CurrentRightProb.reserve(rHap.maxRepSize);
//                CurrentNoRecoRightProb.reserve(rHap.maxRepSize);
//                recomProb.reserve(rHap.maxRepSize);
//
//
//                vector<bool> missing_loo(markerCount,false);
//                MarkovModel MM(tHap_loo,rHap_loo,missing_loo,rHap.major,LowMemory);
//                MM.CopyParameters(MP);
//                MM.initializeMatrices(rHap_loo,tHap_loo);
//                MM.ReinitializeMatrices();
//
//
//                if(!missing_loo[0])
//                {
//                    if(!tHap_loo.getMissingScaffoldedHaplotype(0,0))
//                    {
//                        ConditionJunctionProb(rHap_loo,0,MM.junctionLeftProb[0],
//                                              MM.Error[0],
//                                              tHap_loo.getScaffoldedHaplotype(0,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
//                                              tHap_loo.getScaffoldedHaplotype(0,0),
//                                              MM.backgroundError,
//                                              rHap_loo.ReducedStructureInfo[0]);
//                    }
//                }
//
//                for(int group=1;group<=rHap.NoBlocks;group++)
//                {
//                    LeftTraverse(rHap_loo,tHap_loo,0,MM,group,recomProb,rHap);
//                }
//
//                double tempLogLikelihoodValue=CalculateLikelihood(rHap_loo,MM);
//
//                if(em)
//                {
//                    for(int group=rHap.NoBlocks;group>0;group--)
//                    {
//                        EMTraverse(rHap_loo,tHap_loo,0,MM,group,
//                                   recomProb,PrevRightFoldedProb,
//                                   CurrentRightProb,CurrentNoRecoRightProb,rHap);
//                    }
//                    MM.empiricalCount++;
//                }
//
//#pragma omp critical
//                {
//                    (*MP)+=MM;
//                    LogLikelihoodValue+=tempLogLikelihoodValue;
//                }

            }

//
//
//
//            MP->UpdateModel();
//
//            double crossovers = 0;
//            for (int i = 0; i < rHap.numMarkers - 1; i++)
//                crossovers += MP->Recom[i];
//
//            double errors = 0;
//            for (int i = 0; i <  rHap.numMarkers ; i++)
//            {
//                double heterozygosity = 1.0 - pow(rHap.AlleleFreq[i],2)
//                                        - pow(1-rHap.AlleleFreq[i],2);
//                errors += MP->Error[i] * heterozygosity;
//            }
//            errors /= (double) rHap.numMarkers  + 1e-30;
//
//            printf("      %.0f mosaic crossovers expected per haplotype\n", crossovers);
//            printf("      %.3g errors in mosaic expected per marker\n", errors);
//            printf("    Log-Likelihood of this Iteration : %.5g", LogLikelihoodValue);
//            cout<<endl;
//
//            rHap.Recom=MP->Recom;
//            rHap.Error=MP->Error;
//
        }
//
//
//
//        printf("\n Saving estimated parameters for future use/reference to ...");
//        cout<<endl;
//        MP->WriteParameters(rHap.markerName, outFile, false);
//        int time_load = time(0) - time_prev;
//
//        cout << "\n Time taken for parameter estimation = " << time_load << " seconds. "<<endl;
//
//
//        if(rHap.vcfType || EstimationRounds>0)
//        {
//            std::cout << "\n Writing final reduced haplotype information to [.m3vcf] file  : " <<outFile+".m3vcf" + (gzip ? ".gz" : "")<< endl<<endl;
//            rHap.writem3vcfFile(outFile,gzip);
//
//            if(rHap.vcfType)
//            {
//                remove(outFile+".draft"+".m3vcf" + (gzip ? ".gz" : ""));
//                std::cout << "\n Temporary Draft [.m3vcf] file "<<outFile+".draft"+".m3vcf" + (gzip ? ".gz" : "")<<" deleted ... "<< endl<<endl;
//            }
//
//        }
//
//
//        return MP;
//



//
//    MarkovParameters *MP=new MarkovParameters(FullrHap.numMarkers);
//    MP->Recom=FullrHap.Recom;
//    MP->Error=FullrHap.Error;
//    stats->Initialize(FullrHap.numMarkers, tgwasHap.numMarkers);
//
//    cout << "\n Starting Imputation ..."<<endl<<endl;
//
//    TotalNovcfParts=0;
//
//    int maxVcfSample= MyAllVariables->myOutFormat.vcfBuffer,NumVcfWritten=0,NumVcfToBeWritten=0;
//    if((maxVcfSample)>=tgwasHap.numSamples)
//        maxVcfSample=tgwasHap.numSamples;
//    int TotalNumSamples=tgwasHap.numSamples;
//
//    SinglePartialDosageData.ReParameterizePartialDosageData(ChunkNo, FullrHap, tgwasHap);
//    CurrentPartialDosageData=&SinglePartialDosageData;
//
//    CurrentPartialDosageData->UpdatePartialDosageData(maxVcfSample, NumVcfToBeWritten);
//
//    int StartSamId=0, EndSamId;
//    TimeToWrite=0;
//
//
//
//    for(int batchNo=0; ; batchNo++)
//    {
//        EndSamId = StartSamId + (maxVcfSample) <  TotalNumSamples ? StartSamId + (maxVcfSample) : TotalNumSamples;
//
//        printf("  Imputing Samples %d-%d [%0.0f%%] out of %d samples ...", StartSamId + 1, EndSamId, 100*(float)EndSamId/TotalNumSamples, TotalNumSamples);
//        cout<<endl;
//
//
//#pragma omp parallel for
//        for(int SampleId=StartSamId;SampleId<EndSamId;SampleId ++)
//        {
//
//            int MMIndexId=0;
//
//            pair <int, int> SwapDosageData;
//            int hapId= tgwasHap.CummulativeSampleNoHaplotypes[SampleId];
//            //(2*SampleId);
//
//            vector<float> PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,recomProb;
//            DosageData *ThisSamplePartialDosageData;
//
//#pragma omp critical (BindSample)
//            {
//#ifdef _OPENMP
//                MMIndexId = omp_get_thread_num();
//#endif
//                ThisSamplePartialDosageData=CurrentPartialDosageData;
//                SwapDosageData=CurrentPartialDosageData->IndexSample(SampleId);
//            }
//
//
//            MarkovModel &MM=MainMarkovModel[MMIndexId];
//            MM.CopyParametersNew(MP);
//
//            int hapIdIndiv=hapId;
//            do{
//
//                MM.ThisHapId=hapIdIndiv;
//                MM.ReinitializeMatrices();
//                ThisSamplePartialDosageData->BindSampleMModel(MM,SwapDosageData.second,hapIdIndiv-hapId);
//
//                if(FullrHap.MapRefToTar[0]!=-1 && tgwasHap.RetrieveMissingScaffoldedHaplotype(hapIdIndiv,0)=='0')
//                {
//
//                    int TargetMarkerPosition = FullrHap.MapRefToTar[0];
//                    ConditionJunctionProb(FullrHap,0,MM.junctionLeftProb[0],MM.Error[0],
//                                          tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,TargetMarkerPosition)=='1'? FullrHap.AlleleFreq[0] : 1-FullrHap.AlleleFreq[0],
//                                          tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,TargetMarkerPosition),MM.backgroundError,FullrHap.ReducedStructureInfo[0]);
//                }
//
//                for(int group=1;group<=FullrHap.NoBlocks;group++)
//                {
//                    LeftTraverse(FullrHap,tgwasHap,hapIdIndiv,MM,group,recomProb);
//                }
//
//                MM.CeateProbSum(FullrHap.NoBlocks,FullrHap.numHaplotypes);
//                for(int group=FullrHap.NoBlocks;group>0;group--)
//                {
//                    ImputeTraverse(FullrHap, hapIdIndiv,MM,group, recomProb,PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb);
//                }
//
//#pragma omp critical (StatUpdate)
//                {
//                    stats->NewUpdate(FullrHap, tgwasHap, hapIdIndiv, MM.DosageHap, MM.LooDosageHap);
//                }
//
//                if(tgwasHap.SampleNoHaplotypes[SampleId]==1)
//                {
//                    break;
//                }
//            }while(hapId == hapIdIndiv++);
//
//
//
//            if (MyAllVariables->myOutFormat.hapOutput && !MyAllVariables->myOutFormat.unphasedOutput)
//            {
//#pragma omp critical (PrintHapData)
//                PrintHaplotypeData( hapIdIndiv, SampleId, MM.DosageHap);
//            }
//
//            if(MyAllVariables->myOutFormat.doseOutput)
//            {
//#pragma omp critical  (PrintDoseData)
//                PrintDosageData(SampleId, ThisSamplePartialDosageData->hapDosage[(2*SwapDosageData.second)], ThisSamplePartialDosageData->hapDosage[(2*SwapDosageData.second)+1] );
//            }
//
//            if(MyAllVariables->myOutFormat.vcfOutput)
//            {
//                if(SwapDosageData.first==1)
//                {
//#pragma omp critical (printVCF)
//                    {
//                        TotalNovcfParts++;
//                        printf("       Saving Samples in temporary VCF file ... ");
//                        cout<<endl;
//                        ThisSamplePartialDosageData->FlushPartialVcf(TotalNovcfParts);
//
//                        TimeToWrite+=ThisSamplePartialDosageData->TimeToWrite;
//                        NumVcfToBeWritten = NumVcfWritten + maxVcfSample;
//                        ThisSamplePartialDosageData->UpdatePartialDosageData(maxVcfSample < TotalNumSamples-NumVcfToBeWritten ? maxVcfSample : TotalNumSamples-NumVcfToBeWritten, NumVcfToBeWritten);
//                        NumVcfWritten+=maxVcfSample;
//                    }
//                }
//            }
//
//
//
//        }
//
//        StartSamId=EndSamId ;
//
//        if(StartSamId>=TotalNumSamples)
//            break;
//    }
//
//    delete MP;
//


    TimeToImpute = time(0) - time_prev;
    cout <<endl<< " Imputation successful (" << TimeToImpute << " seconds) !!!"<<endl;
    TimeToImpute-=(TimeToWrite);

}


void Imputation::Minimac3ImputeThisChunk(int ChunkId, HaplotypeSet &FullrHap, HaplotypeSet &tgwasHap, HaplotypeSet &rgwasHap)
{
    ChunkNo=ChunkId;
    THapUnchunked=&tgwasHap;
    rHapChunked=&FullrHap;

    int time_prev = time(0);

    MarkovParameters *MP=new MarkovParameters(FullrHap.numMarkers);
    MP->Recom=FullrHap.Recom;
    MP->Error=FullrHap.Error;
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

        int counter_sample_imputation_finished = 0;

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
                ThisSamplePartialDosageData->BindSampleMModel(MM,SwapDosageData.second,hapIdIndiv-hapId);

                if(FullrHap.MapRefToTar[0]!=-1 && tgwasHap.RetrieveMissingScaffoldedHaplotype(hapIdIndiv,0)=='0')
                {
                    
                    int TargetMarkerPosition = FullrHap.MapRefToTar[0]; 
                    ConditionJunctionProb(FullrHap,0,MM.junctionLeftProb[0],MM.Error[0],
                                            tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,TargetMarkerPosition)=='1'? FullrHap.AlleleFreq[0] : 1-FullrHap.AlleleFreq[0],
                                            tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,TargetMarkerPosition),MM.backgroundError,FullrHap.ReducedStructureInfo[0]);
                }

                for(int group=1;group<=FullrHap.NoBlocks;group++)
                {
                    LeftTraverse(FullrHap,tgwasHap,hapIdIndiv,MM,group,recomProb);
                }

                MM.CeateProbSum(FullrHap.NoBlocks,FullrHap.numHaplotypes);
                for(int group=FullrHap.NoBlocks;group>0;group--)
                {
                    ImputeTraverse(FullrHap, hapIdIndiv,MM,group, recomProb,PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb);
                }

                #pragma omp critical (StatUpdate)
                {
                    stats->NewUpdate(FullrHap, tgwasHap, hapIdIndiv, MM.DosageHap, MM.LooDosageHap);
                }

                if(tgwasHap.SampleNoHaplotypes[SampleId]==1)
                    {
                        break;
                    }
            }while(hapId == hapIdIndiv++);

            #pragma omp critical (addCounter)
            {
                counter_sample_imputation_finished++;
            }


            // if (MyAllVariables->myOutFormat.hapOutput && !MyAllVariables->myOutFormat.unphasedOutput)
            // {
            //     #pragma omp critical (PrintHapData)
            //     PrintHaplotypeData( hapIdIndiv, SampleId, MM.DosageHap);
            // }

            // if(MyAllVariables->myOutFormat.doseOutput)
            // {
            //     #pragma omp critical  (PrintDoseData)
            //     PrintDosageData(SampleId, ThisSamplePartialDosageData->hapDosage[(2*SwapDosageData.second)], ThisSamplePartialDosageData->hapDosage[(2*SwapDosageData.second)+1] );
            // }

            if(MyAllVariables->myOutFormat.vcfOutput)
            {
                if(counter_sample_imputation_finished == EndSamId-StartSamId)
                //if(SwapDosageData.first==1)
                {
                    #pragma omp critical (printVCF)
                    {
                        TotalNovcfParts++;
                        printf("       Saving samples to temporary file ... ");
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


void Imputation::LooOptimalStructure(int loo,HaplotypeSet &rHap, HaplotypeSet &HapLoo)
{

    HapLoo.numHaplotypes=rHap.numHaplotypes-1;
    HapLoo.numSamples=rHap.numSamples;
    HapLoo.numMarkers=rHap.numMarkers;
    HapLoo.maxBlockSize=rHap.maxBlockSize;
    HapLoo.NoBlocks=rHap.NoBlocks;
    HapLoo.maxRepSize=rHap.maxRepSize;
    for(int i=0; i<HapLoo.numMarkers; i++)
    {
        HapLoo.MapTarToRef[i]=i;
        HapLoo.MapRefToTar[i]=i;
    }
    HapLoo.AlleleFreq=rHap.AlleleFreq;


    for(int i=0;i<rHap.NoBlocks;i++)
    {

        ReducedHaplotypeInfo &ToBlock=HapLoo.ReducedStructureInfo[i];
        ReducedHaplotypeInfo &FromBlock=rHap.ReducedStructureInfo[i];


        int looIndex=FromBlock.uniqueIndexMap[loo];

        ToBlock.uniqueCardinality=FromBlock.uniqueCardinality;
        ToBlock.startIndex=FromBlock.startIndex;
        ToBlock.endIndex=FromBlock.endIndex;
        ToBlock.BlockSize=FromBlock.BlockSize;
        ToBlock.RepSize=FromBlock.RepSize;

        ToBlock.uniqueCardinality[looIndex]--;
        if(ToBlock.uniqueCardinality[looIndex]==0)
        {
            ToBlock.RepSize--;
            for(int marker=0; marker<(FromBlock.BlockSize); marker++)
            {
                for (int in = 0, out = 0; in < FromBlock.RepSize; in++)
                {
                    if (in != looIndex)
                    {
                        ToBlock.TransposedUniqueHaps[marker][out] = FromBlock.TransposedUniqueHaps[marker][in];
                        out++;
                    }
                }

            }
            for (int in = 0, out = 0; in < FromBlock.RepSize; in++)
            {
                if (in != looIndex)
                {
                    ToBlock.uniqueCardinality[out++] = FromBlock.uniqueCardinality[in];
                }
            }


            for (int in = 0, out = 0; in <(rHap.numHaplotypes); in++)
            {

                if(FromBlock.uniqueIndexMap[in]<looIndex)
                    ToBlock.uniqueIndexMap[out++] = FromBlock.uniqueIndexMap[in];
                else if(FromBlock.uniqueIndexMap[in]>looIndex)
                    ToBlock.uniqueIndexMap[out++] = FromBlock.uniqueIndexMap[in]-1;
            }
        }
        else
        {
            ToBlock.TransposedUniqueHaps=FromBlock.TransposedUniqueHaps;
            for (int in = 0, out = 0; in <(rHap.numHaplotypes); in++)
            {
                if(in!=loo)
                    ToBlock.uniqueIndexMap[out++] = FromBlock.uniqueIndexMap[in];

            }

        }

        for (int in = 0; in < ToBlock.RepSize; in++)
        {
            ToBlock.InvuniqueCardinality[in]=1.0/(float)ToBlock.uniqueCardinality[in];
        }

    }

}


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

        int counter_sample_imputation_finished = 0;

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
                ThisSamplePartialDosageData->BindSampleMModel(MM,SwapDosageData.second,hapIdIndiv-hapId);

                if(tgwasHap.RetrieveMissingScaffoldedHaplotype(hapIdIndiv,0)=='0')
                {
                    ConditionJunctionProb(rgwasHap,0,MM.junctionLeftProb[0],MM.Error[0],tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,0)=='1'? rgwasHap.AlleleFreq[0] : 1-rgwasHap.AlleleFreq[0],
                                          tgwasHap.RetrieveScaffoldedHaplotype(hapIdIndiv,0),MM.backgroundError,rgwasHap.ReducedStructureInfo[0]);
                }

                for(int group=1;group<=rgwasHap.NoBlocks;group++)
                {
                    LeftTraverse(rgwasHap,tgwasHap,hapIdIndiv,MM,group,recomProb);
                }

                MM.CeateProbSum(rgwasHap.NoBlocks,rgwasHap.numHaplotypes);
                for(int group=rgwasHap.NoBlocks;group>0;group--)
                {
                    ImputeTraverse(rgwasHap, hapIdIndiv,MM,group, recomProb,PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb);
                }

                #pragma omp critical (StatUpdate)
                {
                    stats->NewUpdate(FullrHap, tgwasHap, hapIdIndiv, MM.DosageHap, MM.LooDosageHap);
                }

                if(tgwasHap.SampleNoHaplotypes[SampleId]==1)
                    {
                        break;
                    }
            }while(hapId == hapIdIndiv++);



        }

        if(MyAllVariables->myOutFormat.vcfOutput)
        {
            TotalNovcfParts++;
            printf("       Saving samples to temporary file ... ");
            cout<<endl;
            CurrentPartialDosageData->FlushPartialVcf(TotalNovcfParts);

            TimeToWrite+=CurrentPartialDosageData->TimeToWrite;
            NumVcfToBeWritten = NumVcfWritten + maxVcfSample;
            CurrentPartialDosageData->UpdatePartialDosageData(maxVcfSample < TotalNumSamples-NumVcfToBeWritten ? maxVcfSample : TotalNumSamples-NumVcfToBeWritten, NumVcfToBeWritten);
            NumVcfWritten+=maxVcfSample;
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
//    vector<int> &optStructure=rHap.optEndPoints;
//
//    int Start=optStructure[group-1];
//    int End=optStructure[group];
    
    int Start=rHap.ReducedStructureInfo[group-1].startIndex;
    int End=rHap.ReducedStructureInfo[group-1].endIndex;
  
    MM.foldProbabilities(MM.ThisBlockRightProb[End-Start],group-1,
                            rHap.ReducedStructureInfo[group-1],1,rHap.numHaplotypes);

   if(MyAllVariables->myModelVariables.minimac3)
        MM.ImputeSitesMinimac3(hapID, group-1);
    else
        MM.ImputeSites(hapID, group-1);


    splitFoldedProb(recomProb,MM.ThisBlockRightProb[0],MM.ThisBlockRightNoRecoProb[0]);


    MM.unfoldProbabilities(group-1,recomProb,MM.ThisBlockRightNoRecoProb[0],
                           MM.ThisBlockRightProb[End-Start],1,
                           rHap.ReducedStructureInfo,rHap.numHaplotypes);



}


void Imputation::LeftTraverse(HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                              MarkovModel &MM,int group, vector<float> &recomProb)
{
//    vector<int> &optStructure=rHap.optEndPoints;

    int Start=rHap.ReducedStructureInfo[group-1].startIndex;
    int End=rHap.ReducedStructureInfo[group-1].endIndex;
    vector<vector<float> > &ThisBlockLeftProb= MM.leftProb[group-1];

    MM.foldProbabilities(ThisBlockLeftProb[0],group-1,rHap.ReducedStructureInfo[group-1],
                            0,rHap.numHaplotypes);
    MM.CurrentLeftNoRecoProb=ThisBlockLeftProb[0];
    
    if(MyAllVariables->myModelVariables.minimac3)
        MM.WalkLeftMinimac3(hapID, group-1);
    else
        MM.WalkLeft(hapID, group-1);

    splitFoldedProb(recomProb,ThisBlockLeftProb[End-Start],MM.CurrentLeftNoRecoProb);
    MM.unfoldProbabilities(group-1,recomProb,MM.CurrentLeftNoRecoProb,
                                       ThisBlockLeftProb[0],0,
                                       rHap.ReducedStructureInfo,rHap.numHaplotypes);
}





void Imputation::FreeMemory()
{
    free(SinglePartialDosageData.PrintStringPointer);
    if(MyAllVariables->myOutFormat.meta)
        free(SinglePartialDosageData.PrintEmpStringPointer);
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

//    printf("    Outputting HAPLO%d of Individual %s for Haplotype File...",
//           ThisHapId%2+1 ,THapUnchunked->individualName[ThisSampleId].c_str());
//    cout<<endl;
//
//    ifprintf(hapdose, "%s\tHAPLO%d",  THapUnchunked->individualName[ThisSampleId].c_str(), ThisHapId%2+1  );
//    int i=0;
//    vector<bool> &GWASDosage1 = THapUnchunked->GWASOnlyMissingSampleUnscaffolded[ThisHapId];
//
//    for (int index =0; index < rHapChunked->RefTypedTotalCount; index++)
//    {
//        if(rHapChunked->RefTypedIndex[index]==-1)
//        {
//            if(i>=rHapChunked->PrintStartIndex && i <= rHapChunked->PrintEndIndex)
//            {
//                ifprintf(hapdose, "\t%.5f", (*ThisimputedHap)[i]);
//            }
//            i++;
//        }
//        else
//        {
//            int MarkerIndex=rHapChunked->RefTypedIndex[index];
//            ifprintf(hapdose, "\t%.5f", (float)GWASDosage1[MarkerIndex]);
//        }
//    }
//
//    ifprintf(hapdose, "\n");

}


void Imputation::splitFoldedProb(vector<float> &SplitProb, vector<float> &totalProb, vector<float> &noRecomProb)
{
    SplitProb.resize(totalProb.size());


    if(MyAllVariables->myOutFormat.verbose)
    {
        for(int i=0;i<(int)totalProb.size();i++)
        {
            if(noRecomProb[i]>totalProb[i])
                cout<<" DIFFERENCE = "<<(noRecomProb[i]-totalProb[i])<<"\t"<<totalProb[i]<<"\t"<<(noRecomProb[i]-totalProb[i])/totalProb[i] <<endl;
            SplitProb[i]=max(totalProb[i]-noRecomProb[i],0.0f);
        }
    }
    else {
        for (int i = 0; i < (int) totalProb.size(); i++) {

            SplitProb[i] = max(totalProb[i] - noRecomProb[i], 0.0f);
            assert(SplitProb[i] >= 0.0);
        }
    }
}




void Imputation::ConditionJunctionProb(HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                            double e, double freq, AlleleType observed, double backgroundError,
                            ReducedHaplotypeInfo &Info)
{
    vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];


    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;
    int NoStates=rHap.numHaplotypes;

    for (int i = 0; i<NoStates; i++)
    {
        Prob[i]*= TempHap[Info.uniqueIndexMap[i]] ==observed ? pmatch : prandom;
    }

}





