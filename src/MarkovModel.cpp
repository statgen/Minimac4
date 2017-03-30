
#include "MarkovModel.h"




void MarkovModel::WalkLeft(HaplotypeSet &ThistHap, int &hapID,
                           int group, ReducedHaplotypeInfo &Info,
                           vector<double> &alleleFreq)
{
    vector<vector<float> > &Leftprob =leftProb[group];

    int &Start=Info.startIndex;
    int &End=Info.endIndex;

    noReducedStatesCurrent=Info.RepSize;

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {
        PrecisionJump[markerPos]=Transpose(Leftprob[markerPos-Start-1],
                  Leftprob[markerPos-Start],CurrentLeftNoRecoProb,
                  Recom[markerPos-1],Info.uniqueCardinality);

        if (!ThistHap.RetrieveMissingScaffoldedHaplotype(hapID,markerPos))
        {
                Condition(markerPos,Leftprob[markerPos-Start],
                     CurrentLeftNoRecoProb,
                     ThistHap.RetrieveScaffoldedHaplotype(hapID,markerPos),
                     Error[markerPos],
                     ThistHap.RetrieveScaffoldedHaplotype(hapID,markerPos)?
                          alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);
        }
    }
}





void MarkovModel::ImputeSites(int hapID,int group,
                                   vector<float> &PrevRightFoldedProb,
                                   vector<float> &CurrentRightProb,
                                   vector<float> &CurrentNoRecoRightProb)
{

    vector<float> &juncLeftprob = junctionLeftProb[group];
    vector<float> &juncRightProb = PrevjunctionRightProb;
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double> &alleleFreq=rHap->AlleleFreq;


    vector<double> value(0);


    ReCreateLeftNoRecoProb(*tHap,hapID,group,Info,alleleFreq);


    vector<vector<float> > &Leftprob = leftProb[group];
    vector<vector<float> > &leftNoRecomProb= ThisBlockLeftNoRecoProb;
    int startIndex=Info.startIndex;
    int endIndex=Info.endIndex;
    CurrentRightProb=PrevRightFoldedProb;
    CurrentNoRecoRightProb=PrevRightFoldedProb;
    noReducedStatesCurrent=Info.RepSize;


    fill(Constants.begin(), Constants.end(), 0.0);
    for(int i=0;i<refCount;i++)
            Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);


    ImputeChunk(group, hapID, endIndex,
              Leftprob[endIndex-startIndex], CurrentRightProb,
              leftNoRecomProb[endIndex-startIndex],CurrentNoRecoRightProb,
              Leftprob[0],PrevRightFoldedProb);

    for (int markerPos=endIndex-1; markerPos>startIndex; markerPos--)
    {

        if (!tHap->RetrieveMissingScaffoldedHaplotype(hapID,markerPos+1))
        {
              Condition(markerPos+1,CurrentRightProb,CurrentNoRecoRightProb,
                      tHap->RetrieveScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],
                       tHap->RetrieveScaffoldedHaplotype(hapID,markerPos+1)?
                          alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }


        tempRightProb=CurrentRightProb;
        Transpose(tempRightProb,CurrentRightProb,CurrentNoRecoRightProb,Recom[markerPos],
                            Info.uniqueCardinality);

        ImputeChunk(group, hapID, markerPos,
              Leftprob[markerPos-startIndex], CurrentRightProb,
              leftNoRecomProb[markerPos-startIndex],CurrentNoRecoRightProb,
              Leftprob[0],PrevRightFoldedProb);

    }

    if (!tHap->RetrieveMissingScaffoldedHaplotype(hapID,startIndex+1))
    {

          Condition(startIndex+1,CurrentRightProb,CurrentNoRecoRightProb,
                  tHap->RetrieveScaffoldedHaplotype(hapID,startIndex+1),
                  Error[startIndex+1],
                   tHap->RetrieveScaffoldedHaplotype(hapID,startIndex+1)?
                      alleleFreq[startIndex+1] : 1-alleleFreq[startIndex+1],Info);
    }

    tempRightProb=CurrentRightProb;
    Transpose(tempRightProb,CurrentRightProb,CurrentNoRecoRightProb,Recom[startIndex],Info.uniqueCardinality);


    if(startIndex==0)
    {
        ImputeChunk(group, hapID, startIndex,
              Leftprob[startIndex-startIndex], CurrentRightProb,
              leftNoRecomProb[startIndex-startIndex],CurrentNoRecoRightProb,
              Leftprob[0],PrevRightFoldedProb);
    }
}


void MarkovModel::ImputeChunk(int group, int &hapID, int &position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb)
{

    CurrentTypedSite=rHapFull->MapTarToRef[position];
    CurrentObsMissing=tHap->RetrieveMissingScaffoldedHaplotype(hapID,position);
    CurrentObs=tHap->RetrieveScaffoldedHaplotype(hapID,position);

    FindPosteriorProbWithThreshold( group, position,
                    Leftprob, rightProb,
                    leftNoRecoProb,rightNoRecoProb,
                    leftEndProb,rightEndProb);

    if(MyAllVariables->myModelVariables.probThreshold>0.0)
    {
        unfoldProbabilitiesWithThreshold(group, leftNoRecoProb, Leftprob,
                                  rightNoRecoProb, rightProb, leftEndProb, rightEndProb);
    }

    else
    {
          unfoldProbabilitiesAllProb(group, leftNoRecoProb, Leftprob,
                                  rightNoRecoProb, rightProb, leftEndProb, rightEndProb);
    }

    ImputeRemainingSitesbyBlock( group, position   ,
                                            tHapFull->FlankRegionStart[position],
                                            CurrentTypedSite-1);
    ImputeRemainingSitesbyBlock( group, position   ,
                                            CurrentTypedSite+1,
                                            tHapFull->FlankRegionEnd[position]);
}





void MarkovModel::FindPosteriorProbWithThreshold( int group, int position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb)
{


    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double> &alleleFreq=rHap->AlleleFreq;
    float Pref=0.0,Palt=0.0,ptotal;
    double sum=0.0;
    SummedProb=0.0;

    float *value = (float *)alloca(noReducedStatesCurrent*sizeof(float));
    for(int i=0; i<noReducedStatesCurrent; i++)
    {
        // careful: order of operations is important to avoid overflows
        value[i] = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);

        sum+=value[i];

    }

    vector<bool> &TempHap = Info.TransposedUniqueHaps[position-Info.startIndex];

    for (int i=0; i<noReducedStatesCurrent;)
    {
        bool hp = TempHap[i];
        float pp=0.0;
        pp= value[i] + (hp? Palt:Pref);

        i++;
        while ((i < noReducedStatesCurrent) && (hp == TempHap[i]))
        {
            pp += value[i];
            i++;
        }
        if(hp)
            Palt = pp ;
        else
            Pref = pp;
    }


    ptotal=Pref+Palt;
    (*DosageHap)[CurrentTypedSite] =  (Palt / ptotal);


    if(!CurrentObsMissing)
    {
        double Err=Error[position];
        double freq=CurrentObs? alleleFreq[position] : 1-alleleFreq[position];
        double fmatch = 1.0 / ( (1.0 - Err) + Err*freq + backgroundError);
        double fmismatch = 1.0 / (Err * freq + backgroundError);


        if(CurrentObs)
        {
            Palt *= fmatch;
            Pref *= fmismatch;
        }
        else
        {
            Pref *= fmatch;
            Palt *= fmismatch;
        }

        ptotal =Pref+Palt;
        (*LooDosageHap)[position] =  (Palt / ptotal);

    }


    if(MyAllVariables->myModelVariables.probThreshold>0.0)
    {
        sum=1.0/sum;
        int Index=0;
        NoBestMatchHaps=0;
        for (int i=0; i<noReducedStatesCurrent;i++)
        {
            double tempVal=value[i]*sum;
            if(tempVal >= MyAllVariables->myModelVariables.probThreshold)
            {
                BestMatchHaps[Index++]=i;
                NoBestMatchHaps++;
            }
        }
    }

}





void MarkovModel::unfoldProbabilitiesWithThreshold(int bridgeIndex,
                                         vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                         vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                         vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb)
{
    ReducedHaplotypeInfo &thisInfo = rHap->ReducedStructureInfo[bridgeIndex];
    int N = NoBestMatchHaps;

    float *Leftadj_rec = (float *)alloca(N*sizeof(float));
    float *Leftadj_norec = (float *)alloca(N*sizeof(float));
    float *Rightadj_rec = (float *)alloca(N*sizeof(float));
    float *Rightadj_norec = (float *)alloca(N*sizeof(float));

    NoBestMatchFullRefHaps=0;

    for (int index=0; index<NoBestMatchHaps; index++)
    {

        int i=BestMatchHaps[index];
        double temp=LeftNoRecomProb[i];
        double tempInvCardinality=thisInfo.InvuniqueCardinality[i];

        Leftadj_rec[index] = (LeftTotalProb[i] - temp ) * tempInvCardinality;
        Leftadj_norec[index] = temp / PrevLeftFoldedProb[i];

        temp=RightNoRecomProb[i];
        Rightadj_rec[index] = (RightTotalProb[i] - temp )  * tempInvCardinality;
        Rightadj_norec[index] = RightNoRecomProb[i] / PrevRightFoldedProb[i];
    }


    vector<float> &LeftPrev=junctionLeftProb[bridgeIndex];
    vector<float> &RightPrev=PrevjunctionRightProb;
    double ThisLeftProb,ThisRightProb;

    for (int index=0; index<NoBestMatchHaps; index++)
    {


        int i=BestMatchHaps[index];

        for(int K=0;K< thisInfo.uniqueCardinality[i];K++)
        {

            int UnMappedIndex=thisInfo.uniqueIndexReverseMaps[i][K];
            ThisLeftProb = Leftadj_rec[index] + Leftadj_norec[index]*LeftPrev[UnMappedIndex];
            ThisRightProb = Rightadj_rec[index] + Rightadj_norec[index]*RightPrev[UnMappedIndex];
            probHapFullAverage[NoBestMatchFullRefHaps]=ThisLeftProb*ThisRightProb;

            BestMatchFullRefHaps[NoBestMatchFullRefHaps]=UnMappedIndex;
            NoBestMatchFullRefHaps++;
        }
    }

}
void MarkovModel::unfoldProbabilitiesAllProb(int bridgeIndex,
                                         vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                         vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                         vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb)
{
    ReducedHaplotypeInfo &thisInfo = rHap->ReducedStructureInfo[bridgeIndex];
    int N = thisInfo.RepSize;;
    int noReference=rHap->numHaplotypes;

    float *Leftadj_rec = (float *)alloca(N*sizeof(float));
    float *Leftadj_norec = (float *)alloca(N*sizeof(float));
    float *Rightadj_rec = (float *)alloca(N*sizeof(float));
    float *Rightadj_norec = (float *)alloca(N*sizeof(float));

    NoBestMatchFullRefHaps=0;

    for (int index=0; index<N; index++)
    {

        int i=index;
        double temp=LeftNoRecomProb[i];
        double tempInvCardinality=thisInfo.InvuniqueCardinality[i];

        Leftadj_rec[index] = (LeftTotalProb[i] - temp ) * tempInvCardinality;
        Leftadj_norec[index] = temp / PrevLeftFoldedProb[i];

        temp=RightNoRecomProb[i];
        Rightadj_rec[index] = (RightTotalProb[i] - temp )  * tempInvCardinality;
        Rightadj_norec[index] = RightNoRecomProb[i] / PrevRightFoldedProb[i];
    }


    vector<float> &LeftPrev=junctionLeftProb[bridgeIndex];
    vector<float> &RightPrev=PrevjunctionRightProb;
    double ThisLeftProb,ThisRightProb;


    for (int index=0; index<noReference; index++)
    {
        int UnMappedIndex=thisInfo.uniqueIndexMap[index];
        ThisLeftProb = Leftadj_rec[UnMappedIndex] + Leftadj_norec[UnMappedIndex]*LeftPrev[index];
        ThisRightProb = Rightadj_rec[UnMappedIndex] + Rightadj_norec[UnMappedIndex]*RightPrev[index];
        probHapFullAverage[index]=ThisLeftProb*ThisRightProb;
    }

//cout<<" POSI\t"<<N<<"\t"<<noReference<<"\t"<<thisInfo.RepSize<<"\t"<<refCount<<endl;

// cout<<" JAAA\t"<<sum2 <<"\t"<<tempMaxVal<<"\t"<< MaxVal <<"\t"<<NoBestMatchHaps<<"\t"<<COUNT<<"\t"<<thisInfo.RepSize<<endl;

}






void MarkovModel::FoldBackProbabilitiesAllProb(ReducedHaplotypeInfo &Info, ReducedHaplotypeInfo &TarInfo)
{

    fill(FoldedProbValue.begin(), FoldedProbValue.end(), 0.0);
    fill(pREF.begin(), pREF.end(), 0.0 );
    fill(pALT.begin(), pALT.end(), 0.0 );
    noNewReducedStates=Info.RepSize;


//    int FullRefIndex=0;

    for (int index=0; index<refCount; index++)
    {
        int MappedIndex=Info.uniqueIndexMap[index];
        FoldedProbValue[MappedIndex]+=probHapFullAverage[index];
    }
//
//    for(int i=0;i<noNewReducedStates ;i++)
//    {
//        FinalBestMatchfHaps[i]=i;
//    }

}
void MarkovModel::FoldBackProbabilitiesWithThreshold(ReducedHaplotypeInfo &Info, ReducedHaplotypeInfo &TarInfo)
{

    fill(FoldedProbValue.begin(), FoldedProbValue.end(), 0.0);
    fill(pREF.begin(), pREF.end(), 0.0 );
    fill(pALT.begin(), pALT.end(), 0.0 );
    noNewReducedStates = 0;

    int FullRefIndex=0;

    for (int index=0; index<NoBestMatchFullRefHaps; index++)
    {

        int MappedIndex=Info.uniqueIndexMap[BestMatchFullRefHaps[index]];
        int tempIndex= Myfind (FinalBestMatchfHaps, MappedIndex);

        if (tempIndex<0)
        {
            FinalBestMatchfHaps[noNewReducedStates]=MappedIndex;
            tempIndex=noNewReducedStates++;
        }

        FoldedProbValue[tempIndex]+=probHapFullAverage[FullRefIndex];
        FullRefIndex++;
    }


}





void MarkovModel::ImputeRemainingSitesbyBlock(int bridgeIndex, int TypedMarkerId,
                                              int StartPos, int EndPos)
{



    for(int position=StartPos;position<=EndPos;)
    {

        int NewBlockIndex=rHapFull->MarkerToReducedInfoMapper[position];
        ReducedHaplotypeInfo &Info=rHapFull->ReducedStructureInfo[NewBlockIndex];
        ReducedHaplotypeInfo &TarInfo=rHap->ReducedStructureInfo[bridgeIndex];
        int TempEndPos= EndPos<Info.endIndex ? EndPos:Info.endIndex;


        if(MyAllVariables->myModelVariables.probThreshold>0.0)
            FoldBackProbabilitiesWithThreshold (Info,TarInfo);
        else
            FoldBackProbabilitiesAllProb (Info,TarInfo);


        ImputeIntermediateRegionFull(Info, TypedMarkerId, position, TempEndPos);


//        if(CurrentTypedSite<=TempEndPos && CurrentTypedSite>=position)
//            CreateLooDosage(position,TypedMarkerId, CurrentObsMissing, CurrentObs );


        position=TempEndPos+1;
    }

}







void MarkovModel::ImputeIntermediateRegionFull(ReducedHaplotypeInfo &Info, int TypedMarkerId, int StartPos, int EndPos)
{

    CreatePRefPAlt(Info, StartPos,EndPos);
    CreateDosages(TypedMarkerId, StartPos,EndPos);

}
void MarkovModel::CreatePRefPAlt(ReducedHaplotypeInfo &Info, int StartPos, int EndPos)
{
    int tempPosition;
    int ThisStartIndex=Info.startIndex;

    for(tempPosition = StartPos; tempPosition<=EndPos; tempPosition++)
    {
        vector<bool> &ThisMapped=Info.TransposedUniqueHaps[tempPosition-ThisStartIndex];

        float &Palt = pALT[tempPosition-StartPos];
        float &Pref = pREF[tempPosition-StartPos];

        for (int i=0; i< noNewReducedStates;i++)
        {
            int MappedIndex=FinalBestMatchfHaps[i];

            if(ThisMapped[MappedIndex])
                Palt+=FoldedProbValue[i];
            else
                Pref+=FoldedProbValue[i];
        }
    }

//
//    for (int i=0; i< noNewReducedStates;i++)
//    {
//        tempPosition=StartPos;
//        int MappedIndex=FinalBestMatchfHaps[i];
//
//        vector<bool> &ThisMapped=Info.uniqueHaps[MappedIndex];
//        int ThisStartIndex=Info.startIndex;
//
//        while(tempPosition <= EndPos)
//        {
////            assert((tempPosition-StartPos)<rHapFull->maxBlockSize);
//            if(ThisMapped[tempPosition-ThisStartIndex])
//                pALT[tempPosition-StartPos]+=FoldedProbValue[i];
//            else
//                pREF[tempPosition-StartPos]+=FoldedProbValue[i];
//
//            tempPosition++;
//        }
//    }

}
 void MarkovModel::CreateDosages(int TypedMarkerId, int StartPos, int EndPos)
{
    float Pref,Palt,ptotal;
//    bool mle;


//    if(noNewReducedStates==0)
//    {
//
//
//        for(int tempPosition=StartPos;tempPosition<=EndPos;tempPosition++)
//        {
//
//            Palt=rHapFull->AlleleFreq[tempPosition];
//            mle=false;
//            if(Palt>0.5)
//            {
//                mle=true;
//            }
//
//
//            imputedDose[tempPosition] += imputedHap[tempPosition]= (Palt);
//            imputedAlleleNumber[tempPosition] = mle;
//
//        }
//        return;
//
//
//    }
//


    for(int tempPosition=StartPos;tempPosition<=EndPos;tempPosition++)
    {
        Pref=pREF[tempPosition-StartPos];
        Palt=pALT[tempPosition-StartPos];
        ptotal=Pref+Palt;
//
//        mle = false;
//        if(Pref<Palt)
//        {
//            mle=true;
//        }


//        imputedDose[tempPosition] += imputedHap[tempPosition]= (Palt / ptotal);
        (*DosageHap)[tempPosition] =  (Palt / ptotal);

//        if(tempPosition==0)
//        cout<<" NOW = "<<(*DosageHap)[0]<<endl;
//        imputedAlleleNumber[tempPosition] = mle;




    }
}
int MarkovModel::Myfind(vector<int> &MyVector, int Value)
{
    for(int i=0;i<noNewReducedStates;i++)
    {
        if(Value==MyVector[i])
            return i;
    }
    return -1;
}





void MarkovModel::initializeMatricesNew()
{

    // Left Probabilities Initialize (NO NEED TO RE-INITIALIZE)
    leftProb.resize(rHap->NoBlocks);
    for(int i=0;i<rHap->NoBlocks;i++)
    {
        vector<vector<float> > &TempLeft=leftProb[i];
        TempLeft.resize(rHap->maxBlockSize);
        for(int j=0;j<rHap->maxBlockSize;j++)
        {
            TempLeft[j].resize(rHap->maxRepSize);
        }
    }


    // Left No Recom Probabilities Initialize (NO NEED TO RE-INITIALIZE)
    CurrentLeftNoRecoProb.resize(rHap->maxRepSize);
    ThisBlockLeftNoRecoProb.resize(rHap->maxBlockSize);
    for(int i=0;i<rHap->maxBlockSize;i++)
    {
        ThisBlockLeftNoRecoProb[i].resize(rHap->maxRepSize);
    }


    // Junction Probabilities Initialize (ONLY junctionLeftProb[0] and
    // PrevjunctionRightProb nneds to be re-initialized. DONE

    junctionLeftProb.resize(rHap->NoBlocks+1);
    PrevjunctionRightProb.resize(rHap->numHaplotypes);

    for(int i=0;i<=rHap->NoBlocks;i++)
    {
        junctionLeftProb[i].resize(rHap->numHaplotypes);
    }


    // Other Variables Initialize (NO NEED TO RE-INITIALIZE
    Constants.resize(rHap->maxRepSize);
    tempRightProb.reserve(rHap->maxRepSize);

    PrecisionJump.resize(rHap->numMarkers); // RE-INTIALIZE DONE

    // Full Dimension Variable Initialize (NO NEED TO REINITLAIZE)

    FoldedProbValue.resize(rHapFull->maxRepSize);
    pREF.resize(rHapFull->maxBlockSize);
    pALT.resize(rHapFull->maxBlockSize);
    probHapFullAverage.resize(rHapFull->numHaplotypes);



    FinalBestMatchfHaps.resize(rHapFull->maxRepSize);
    if(MyAllVariables->myModelVariables.probThreshold>0.0)
    {
        BestMatchHaps.resize(rHap->maxRepSize);
        BestMatchFullRefHaps.resize(refCount);
    }
    else
    {
        for(int i=0;i<rHapFull->maxRepSize ;i++)
        {
            FinalBestMatchfHaps[i]=i;
        }
    }




}


void MarkovModel::ReinitializeMatrices()
{
    BeforeLastUntypedSite=rHapFull->numMarkers-1;
    NoPrecisionJumps=0;
    fill(PrecisionJump.begin(), PrecisionJump.end(), false);
    fill(junctionLeftProb[0].begin(), junctionLeftProb[0].end(), 1.0);
    fill(PrevjunctionRightProb.begin(), PrevjunctionRightProb.end(), 1.0);

}







void MarkovModel::ReCreateLeftNoRecoProb(HaplotypeSet &tHap, int &hapID,
                           int group, ReducedHaplotypeInfo &Info,
                           vector<double> &alleleFreq)
{

    int &Start=Info.startIndex;
    int &End=Info.endIndex;
    noReducedStatesCurrent=Info.RepSize;
    ThisBlockLeftNoRecoProb[0]=leftProb[group][0];

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {
        vector<float> &NextnoRecomProb = ThisBlockLeftNoRecoProb[markerPos-Start];
        double complement = 1. - Recom[markerPos-1];
        NextnoRecomProb=ThisBlockLeftNoRecoProb[markerPos-Start-1];
        double freq=tHap.RetrieveScaffoldedHaplotype(hapID,markerPos)? alleleFreq[markerPos] : 1-alleleFreq[markerPos];
        double e=Error[markerPos];
        bool observed=tHap.RetrieveScaffoldedHaplotype(hapID,markerPos);
        double prandom = e*freq+backgroundError;
        double pmatch = (1.0 - e)+e*freq+backgroundError;

        for (int i = 0; i <noReducedStatesCurrent; i++)
        {
            NextnoRecomProb[i]*=(complement);
        }

        if (PrecisionJump[markerPos])
        {
            for (int i = 0; i <noReducedStatesCurrent; i++)
            {
                NextnoRecomProb[i]*=(1e15);
            }
        }

        if (!tHap.RetrieveMissingScaffoldedHaplotype(hapID,markerPos))
        {
            vector<bool> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];

            for (int i = 0; i<noReducedStatesCurrent; i++)
            {
                NextnoRecomProb[i]*=(TempHap[i]==observed)?pmatch:prandom;
            }
        }
    }
}








void MarkovModel::foldProbabilities(vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference) //0 - left; 1 - right
{
    vector<int> *TempuniqueIndexMap=&Info.uniqueIndexMap;
    fill(foldProb.begin(), foldProb.end(), 0.0);
    if(direction==0)
    {
        vector<float> *PrevjunctionLeftProb=&junctionLeftProb[bridgeIndex];
        for(int i=0;i<noReference;i++)
        {
            foldProb[(*TempuniqueIndexMap)[i]]+=(*PrevjunctionLeftProb)[i];
        }
    }
    else if(direction==1)
    {

        for(int i=0;i<noReference;i++)
        {
            foldProb[(*TempuniqueIndexMap)[i]]+=PrevjunctionRightProb[i];
        }
    }
}


void MarkovModel::unfoldProbabilities(int bridgeIndex,vector<float> &recomProb,
                                       vector<float> &noRecomProb,vector<float> &PrevFoldedProb,
                                     int direction,vector<ReducedHaplotypeInfo> &StructureInfo,
                                     int noReference)
{
    ReducedHaplotypeInfo &thisInfo = StructureInfo[bridgeIndex];
    int N = thisInfo.RepSize;

    float *adj_rec = (float *)alloca(N*sizeof(float));
    float *adj_norec = (float *)alloca(N*sizeof(float));

    for (int i=0; i<N; i++)
    {
        adj_rec[i] = recomProb[i] * thisInfo.InvuniqueCardinality[i];
    }
    for (int i=0; i<N; i++)
    {
        adj_norec[i] = noRecomProb[i] / PrevFoldedProb[i];
    }
    vector<float> &prev = direction ? PrevjunctionRightProb : junctionLeftProb[bridgeIndex];
    vector<float> &next = direction ? PrevjunctionRightProb : junctionLeftProb[bridgeIndex+1];

    if(direction)
    {
        for (int i=0; i<noReference; i++)
        {
            int m = thisInfo.uniqueIndexMap[i];
            prev[i]*=adj_norec[m];
            prev[i]+=adj_rec[m];
        }
    }
    else
    {
         for (int i=0; i<noReference; i++)
        {
            int m = thisInfo.uniqueIndexMap[i];
            next[i] = adj_rec[m] + adj_norec[m]*prev[i];
        }
    }

}




void MarkovModel::Condition(int markerPos,vector<float> &Prob,
                            vector<float> &noRecomProb, bool observed,double e,double freq,ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    vector<bool> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];


    for (int i = 0; i<noReducedStatesCurrent; i++)
    {

        bool allele=TempHap[i];
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


bool MarkovModel::Transpose(vector<float> &from,
                            vector<float> &to, vector<float> &noRecomProb,
                            double reco,vector<int> &uniqueCardinality)
{
    bool tempPrecisionJumpFlag=false;
    if (reco == 0)
    {
        to=from;
        return false;
    }

    double sum = 0.0;
    for (int i = 0; i <noReducedStatesCurrent; i++)
    {
        sum += from[i];
        noRecomProb[i]*=(1.-reco);
    }

    sum*=(reco/(double)refCount);
    double complement = 1. - reco;

    // avoid underflows
    if (sum < 1e-10)
    {
        tempPrecisionJumpFlag=true;
        sum*= 1e15;
        complement *= 1e15;
        for(int i=0;i<noReducedStatesCurrent;i++)
            noRecomProb[i]*=1e15;
        NoPrecisionJumps++;
    }

    for (int i = 0; i <noReducedStatesCurrent; i++)
    {
        to[i]=from[i]*complement+(uniqueCardinality[i]*sum);
    }

    return tempPrecisionJumpFlag;


 }









