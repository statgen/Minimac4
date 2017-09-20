#include "MarkovModel.h"



void MarkovModel::CeateProbSum(int bridgeIndex, int noReference)
{
    PrevTotalSum = 0.0;
    for (int i=0; i<noReference; i++)
    {
        PrevTotalSum+=junctionLeftProb[bridgeIndex][i];
    }
    InvPrevTotalSum=1.0/PrevTotalSum;
}

void MarkovModel::WalkLeft(int &hapID, int group)
{
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
     vector<double> &alleleFreq=rHap->AlleleFreq;
    vector<vector<float> > &Leftprob =leftProb[group];

    int &Start=Info.startIndex;
    int &End=Info.endIndex;

    noReducedStatesCurrent=Info.RepSize;

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {
        LeftPrecisionJump[markerPos]=LeftTranspose(Leftprob[markerPos-Start-1],
                  Leftprob[markerPos-Start],CurrentLeftNoRecoProb,
                  Recom[markerPos-1],Info.uniqueCardinality);

        if (tHap->RetrieveMissingScaffoldedHaplotype(hapID,markerPos)=='0')
        {
                LeftCondition(markerPos,Leftprob[markerPos-Start],
                     CurrentLeftNoRecoProb,
                     tHap->RetrieveScaffoldedHaplotype(hapID,markerPos),
                     Error[markerPos],
                     tHap->RetrieveScaffoldedHaplotype(hapID,markerPos)=='1'?
                          alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);
        }
    }
}

void MarkovModel::WalkLeftMinimac3(int &hapID, int group)
{
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
     vector<double> &alleleFreq=rHap->AlleleFreq;
    vector<vector<float> > &Leftprob =leftProb[group];

    int &Start=Info.startIndex;
    int &End=Info.endIndex;

    noReducedStatesCurrent=Info.RepSize;

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {
        LeftPrecisionJump[markerPos]=LeftTranspose(Leftprob[markerPos-Start-1],
                  Leftprob[markerPos-Start],CurrentLeftNoRecoProb,
                  Recom[markerPos-1],Info.uniqueCardinality);
        
        int TargetMarkerPosition = rHap->MapRefToTar[markerPos];
        if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID,TargetMarkerPosition)=='0')
        {
          
            LeftCondition(markerPos,Leftprob[markerPos-Start],
                     CurrentLeftNoRecoProb,
                     tHap->RetrieveScaffoldedHaplotype(hapID,TargetMarkerPosition),
                     Error[markerPos],
                     tHap->RetrieveScaffoldedHaplotype(hapID,TargetMarkerPosition)=='1'?
                          alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);
        
        }
    }
}





void MarkovModel::ImputeSitesMinimac3(int hapID,int group)
{

    vector<float> &juncLeftprob = junctionLeftProb[group];
    vector<float> &juncRightProb = PrevjunctionRightProb;
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double> &alleleFreq=rHap->AlleleFreq;


    ReCreateLeftNoRecoProbMinimac3(*tHap,hapID,group,Info,alleleFreq);


    vector<vector<float> > &Leftprob = leftProb[group];
    vector<vector<float> > &leftNoRecomProb= ThisBlockLeftNoRecoProb;
    int startIndex=Info.startIndex;
    int endIndex=Info.endIndex;
    PrevRightFoldedProb = ThisBlockRightProb[endIndex-startIndex];
    ThisBlockRightNoRecoProb[endIndex-startIndex] = ThisBlockRightProb[endIndex-startIndex];
    noReducedStatesCurrent=Info.RepSize;

    fill(Constants.begin(), Constants.end(), 0.0);
    for(int i=0;i<refCount;i++)
            Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);


    ImputeChunkMinimac3(group, hapID, endIndex,
              Leftprob[endIndex-startIndex], ThisBlockRightProb[endIndex-startIndex],
              leftNoRecomProb[endIndex-startIndex],ThisBlockRightProb[endIndex-startIndex],
              Leftprob[0],PrevRightFoldedProb);

    for (int markerPos=endIndex-1; markerPos>startIndex; markerPos--)
    {
        int TargetMarkerPosition = rHap->MapRefToTar[markerPos+1];
        if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID, TargetMarkerPosition)=='0')
        {
            RightCondition(markerPos+1,ThisBlockRightProb[markerPos-startIndex+1],
                           ThisBlockRightProb[markerPos-startIndex],
                           ThisBlockRightNoRecoProb[markerPos-startIndex+1],
                           ThisBlockRightNoRecoProb[markerPos-startIndex],
                           tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                           Error[markerPos+1],
                           tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                           alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }
        else
        {
             ThisBlockRightProb[markerPos-startIndex]= ThisBlockRightProb[markerPos-startIndex+1];
             ThisBlockRightNoRecoProb[markerPos-startIndex]= ThisBlockRightNoRecoProb[markerPos-startIndex+1];
            
        }


        RightPrecisionJump[markerPos] = RightTranspose(ThisBlockRightProb[markerPos-startIndex],
                                                  ThisBlockRightNoRecoProb[markerPos-startIndex],
                                                  Recom[markerPos],
                                                  Info.uniqueCardinality);

        ImputeChunkMinimac3(group, hapID, markerPos,
              Leftprob[markerPos-startIndex], ThisBlockRightProb[markerPos-startIndex],
              leftNoRecomProb[markerPos-startIndex], ThisBlockRightNoRecoProb[markerPos-startIndex],
              Leftprob[0],PrevRightFoldedProb);

    }
    
    int TargetMarkerPosition = rHap->MapRefToTar[startIndex+1];
    if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID,TargetMarkerPosition)=='0')
    {

          RightCondition(startIndex+1, ThisBlockRightProb[1],
                         ThisBlockRightProb[0],
                         ThisBlockRightNoRecoProb[1],
                         ThisBlockRightNoRecoProb[0],
                         tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                         Error[startIndex+1],
                         tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                         alleleFreq[startIndex+1] : 1-alleleFreq[startIndex+1],Info);

    }
    else
    {
         ThisBlockRightProb[0]= ThisBlockRightProb[1];
         ThisBlockRightNoRecoProb[0]= ThisBlockRightNoRecoProb[1];

    }

    RightPrecisionJump[startIndex] = RightTranspose(ThisBlockRightProb[0],
                                                ThisBlockRightNoRecoProb[0],Recom[startIndex],
                                               Info.uniqueCardinality);

    if(startIndex==0)
    {
        ImputeChunkMinimac3(group, hapID, startIndex,
              Leftprob[startIndex-startIndex], ThisBlockRightProb[startIndex-startIndex],
              leftNoRecomProb[startIndex-startIndex],ThisBlockRightNoRecoProb[0],
              Leftprob[0],PrevRightFoldedProb);
    }
}



void MarkovModel::ImputeChunkMinimac3( int group, int hapID, int position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb)
{
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double> &alleleFreq=rHap->AlleleFreq;
    float Pref=0.0,Palt=0.0,ptotal;
    double tempVal;
    vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[position-Info.startIndex];
    CurrentObsMissing=tHap->RetrieveMissingScaffoldedHaplotype(hapID,rHap->MapRefToTar[position]);
    CurrentObs=tHap->RetrieveScaffoldedHaplotype(hapID,rHap->MapRefToTar[position]);

    for(int i=0; i<noReducedStatesCurrent; i++)
    {
        // careful: order of operations is important to avoid overflows
        probHapMinimac3[i] = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);


    }


    for (int i=0; i<noReducedStatesCurrent;)
    {
        AlleleType hp = TempHap[i];
        float pp=0.0;
        pp = probHapMinimac3[i] + (hp == '1' ? Palt:Pref);

        i++;
        while ((i < noReducedStatesCurrent) && (hp == TempHap[i]))
        {
            pp += probHapMinimac3[i];
            i++;
        }
        if(hp == '1')
            Palt = pp ;
        else
            Pref = pp;
    }

    ptotal=Pref+Palt;
    (*DosageHap)[position] =  min(1.0,(double)((Palt / ptotal)));
 

    
    if(CurrentObsMissing=='0')
    {
        double Err=Error[position];
        double freq=CurrentObs == '1' ? alleleFreq[position] : 1-alleleFreq[position];
        double fmatch = 1.0 / ( (1.0 - Err) + Err*freq + backgroundError);
        double fmismatch = 1.0 / (Err * freq + backgroundError);

        if(CurrentObs == '1')
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
        (*LooDosageHap)[rHap->MapRefToTar[position]] =  (Palt / ptotal);

    }
}




void MarkovModel::ImputeSites(int hapID,int group)
{

    vector<float> &juncLeftprob = junctionLeftProb[group];
    vector<float> &juncRightProb = PrevjunctionRightProb;
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double> &alleleFreq=rHap->AlleleFreq;


    ReCreateLeftNoRecoProb(*tHap,hapID,group,Info,alleleFreq);


    vector<vector<float> > &Leftprob = leftProb[group];
    vector<vector<float> > &leftNoRecomProb= ThisBlockLeftNoRecoProb;
    int startIndex=Info.startIndex;
    int endIndex=Info.endIndex;
    PrevRightFoldedProb = ThisBlockRightProb[endIndex-startIndex];
    ThisBlockRightNoRecoProb[endIndex-startIndex] = ThisBlockRightProb[endIndex-startIndex];
    noReducedStatesCurrent=Info.RepSize;
    NoSitesToUnfold=0;
    MostProbableTemplate=-1;

    fill(Constants.begin(), Constants.end(), 0.0);
    for(int i=0;i<refCount;i++)
            Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);


    ImputeChunk(group, hapID, endIndex,
              Leftprob[endIndex-startIndex], ThisBlockRightProb[endIndex-startIndex],
              leftNoRecomProb[endIndex-startIndex],ThisBlockRightProb[endIndex-startIndex],
              Leftprob[0],PrevRightFoldedProb);

    for (int markerPos=endIndex-1; markerPos>startIndex; markerPos--)
    {

        if (tHap->RetrieveMissingScaffoldedHaplotype(hapID,markerPos+1)=='0')
        {
            RightCondition(markerPos+1,ThisBlockRightProb[markerPos-startIndex+1],
                           ThisBlockRightProb[markerPos-startIndex],
                           ThisBlockRightNoRecoProb[markerPos-startIndex+1],
                           ThisBlockRightNoRecoProb[markerPos-startIndex],
                           tHap->RetrieveScaffoldedHaplotype(hapID,markerPos+1),
                           Error[markerPos+1],
                           tHap->RetrieveScaffoldedHaplotype(hapID,markerPos+1)=='1'?
                           alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }
        else
        {
             ThisBlockRightProb[markerPos-startIndex]= ThisBlockRightProb[markerPos-startIndex+1];
             ThisBlockRightNoRecoProb[markerPos-startIndex]= ThisBlockRightNoRecoProb[markerPos-startIndex+1];
            
        }


        RightPrecisionJump[markerPos] = RightTranspose(ThisBlockRightProb[markerPos-startIndex],
                                                  ThisBlockRightNoRecoProb[markerPos-startIndex],
                                                  Recom[markerPos],
                                                  Info.uniqueCardinality);

        ImputeChunk(group, hapID, markerPos,
              Leftprob[markerPos-startIndex], ThisBlockRightProb[markerPos-startIndex],
              leftNoRecomProb[markerPos-startIndex], ThisBlockRightNoRecoProb[markerPos-startIndex],
              Leftprob[0],PrevRightFoldedProb);

    }

    if (tHap->RetrieveMissingScaffoldedHaplotype(hapID,startIndex+1)=='0')
    {

          RightCondition(startIndex+1, ThisBlockRightProb[1],
                         ThisBlockRightProb[0],
                         ThisBlockRightNoRecoProb[1],
                         ThisBlockRightNoRecoProb[0],
                         tHap->RetrieveScaffoldedHaplotype(hapID,startIndex+1),
                         Error[startIndex+1],
                         tHap->RetrieveScaffoldedHaplotype(hapID,startIndex+1)=='1'?
                         alleleFreq[startIndex+1] : 1-alleleFreq[startIndex+1],Info);

    }
    else
    {
         ThisBlockRightProb[0]= ThisBlockRightProb[1];
         ThisBlockRightNoRecoProb[0]= ThisBlockRightNoRecoProb[1];

    }

    RightPrecisionJump[startIndex] = RightTranspose(ThisBlockRightProb[0],
                                                ThisBlockRightNoRecoProb[0],Recom[startIndex],
                                               Info.uniqueCardinality);

    if(startIndex==0)
    {
        ImputeChunk(group, hapID, startIndex,
              Leftprob[startIndex-startIndex], ThisBlockRightProb[startIndex-startIndex],
              leftNoRecomProb[startIndex-startIndex],ThisBlockRightNoRecoProb[0],
              Leftprob[0],PrevRightFoldedProb);
    }
}


void MarkovModel::ImputeChunk(int group, int &hapID, int &position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb)
{
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    int startIndex=Info.startIndex;


    CurrentTypedSite=rHapFull->MapTarToRef[position];
    CurrentObsMissing=tHap->RetrieveMissingScaffoldedHaplotype(hapID,position);
    CurrentObs=tHap->RetrieveScaffoldedHaplotype(hapID,position);

    FindPosteriorProbWithThreshold( group, position,
                    Leftprob, rightProb,
                    leftNoRecoProb,rightNoRecoProb,
                    leftEndProb,rightEndProb);

    UnfoldTheseSites[NoSitesToUnfold++]=position;
    vector<vector<float> > &ThisLeftprob = leftProb[group];

    if(KeepMovingLeft==0)
    {

        int StartPoint = 0;
        int EndPoint = NoSitesToUnfold-2;
        MidPoint = (UnfoldTheseSites[StartPoint] + UnfoldTheseSites[EndPoint])/2 - startIndex ;

        if(MyAllVariables->myModelVariables.probThreshold>0.0)
        {
            unfoldProbabilitiesWithThreshold(group, ThisBlockLeftNoRecoProb[MidPoint],
                                             ThisLeftprob[MidPoint],
                                             ThisBlockRightNoRecoProb[MidPoint],
                                             ThisBlockRightProb[MidPoint],
                                             leftEndProb, rightEndProb);
        }


        for(int i=EndPoint;i>=StartPoint;i--)
        {

            int tempPosition = UnfoldTheseSites[i];
            CurrentTypedSite=rHapFull->MapTarToRef[tempPosition];
            ImputeRemainingSitesbyBlock( group, tempPosition,
                                            tHapFull->FlankRegionStart[tempPosition],
                                            CurrentTypedSite-1);
            ImputeRemainingSitesbyBlock( group, tempPosition,
                                            CurrentTypedSite+1,
                                            tHapFull->FlankRegionEnd[tempPosition]);

        }

        UnfoldTheseSites[StartPoint]=UnfoldTheseSites[EndPoint+1];
        NoSitesToUnfold=1;
    }
    else if(MyAllVariables->myOutFormat.verbose)
    {
        cout<<" SKIPPED A SITE "<<endl;
    }
    if(  (position==(rHap->ReducedStructureInfo[group].startIndex+1) && position!=1) || (position == 0)  )
    {
        int StartPoint = 0;
        int EndPoint = NoSitesToUnfold-1;
        MidPoint = (UnfoldTheseSites[StartPoint] + UnfoldTheseSites[EndPoint])/2 - startIndex ;

        if(MyAllVariables->myModelVariables.probThreshold>0.0)
        {
            unfoldProbabilitiesWithThreshold(group, ThisBlockLeftNoRecoProb[MidPoint],
                                             ThisLeftprob[MidPoint],
                                             ThisBlockRightNoRecoProb[MidPoint],
                                             ThisBlockRightProb[MidPoint],
                                             leftEndProb, rightEndProb);
        }


        for(int i=EndPoint;i>=StartPoint;i--)
        {
            int tempPosition = UnfoldTheseSites[i];
            CurrentTypedSite=rHapFull->MapTarToRef[tempPosition];
            ImputeRemainingSitesbyBlock( group, tempPosition,
                                            tHapFull->FlankRegionStart[tempPosition],
                                            CurrentTypedSite-1);
            ImputeRemainingSitesbyBlock( group, tempPosition,
                                            CurrentTypedSite+1,
                                            tHapFull->FlankRegionEnd[tempPosition]);

        }

    }

}





void MarkovModel::FindPosteriorProbWithThreshold( int group, int position,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb)
{
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double>  &probHap = probHapMatrix[position-Info.startIndex];
    vector<double> &alleleFreq=rHap->AlleleFreq;
    float Pref=0.0,Palt=0.0,ptotal;
    double tempVal;
    vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[position-Info.startIndex];

    if(LeftPrecisionJump[position+1])
        InvPrevTotalSum*=1e15;

    if(RightPrecisionJump[position])
        InvPrevTotalSum*=1e-15;
    SumOfProb[position-Info.startIndex]=InvPrevTotalSum;


     if(MostProbableTemplate!=-1 && MostProbableTemplateVal > (1-MyAllVariables->myModelVariables.topThreshold) )
    {
        int i=MostProbableTemplate;
        fill(probHap.begin(), probHap.end(), 0);

        probHap[i] = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);

        tempVal=probHap[i]*InvPrevTotalSum;

        if(tempVal>(1-MyAllVariables->myModelVariables.topThreshold))
        {
            AlleleType hp = TempHap[i];
            if(hp=='1')
            {
                (*DosageHap)[CurrentTypedSite] =  1.0;
                if(CurrentObsMissing=='0')
                    (*LooDosageHap)[position] =  1.0;
            }
            else
            {
                (*DosageHap)[CurrentTypedSite] =  0.0;
                if(CurrentObsMissing=='0')
                    (*LooDosageHap)[position] =  0.0;
            }

            if(position>(Info.startIndex+1))
                KeepMovingLeft=1;
            else
                {
                    KeepMovingLeft=0;
                }
            return ;
        }
    }


    FirstDiffValue=0.0;
    KeepMovingLeft=0;
    for(int i=0; i<noReducedStatesCurrent; i++)
    {
        // careful: order of operations is important to avoid overflows
        probHap[i] = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);

        tempVal = (probHap[i] - FirstFoldedValue[i]) * InvPrevTotalSum;
//        FirstDiffValue+= tempVal * tempVal ;
        FirstDiffValue+= abs(tempVal) ;
    }

    if(position==Info.endIndex)
    {
        FirstDiffValue=0.0;
        FirstFoldedValue = probHap;
    }

    if(FirstDiffValue< MyAllVariables->myModelVariables.diffThreshold )
    {
        if(position>(Info.startIndex+1))
            KeepMovingLeft=1;
        else
        {
            KeepMovingLeft=0;
        }
    }
    if(KeepMovingLeft==0)
    {
        FirstFoldedValue = probHap;
    }


    for (int i=0; i<noReducedStatesCurrent;)
    {
        AlleleType hp = TempHap[i];
        float pp=0.0;
        pp= probHap[i] + (hp == '1' ? Palt:Pref);

        i++;
        while ((i < noReducedStatesCurrent) && (hp == TempHap[i]))
        {
            pp += probHap[i];
            i++;
        }
        if(hp == '1')
            Palt = pp ;
        else
            Pref = pp;
    }

    ptotal=Pref+Palt;
    (*DosageHap)[CurrentTypedSite] =  (Palt / ptotal);

    if(CurrentObsMissing=='0')
    {
        double Err=Error[position];
        double freq=CurrentObs == '1' ? alleleFreq[position] : 1-alleleFreq[position];
        double fmatch = 1.0 / ( (1.0 - Err) + Err*freq + backgroundError);
        double fmismatch = 1.0 / (Err * freq + backgroundError);

        if(CurrentObs == '1')
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

    MostProbableTemplateVal=0.0;
    for (int i=0; i<noReducedStatesCurrent;i++)
    {
        tempVal=probHap[i];
        if(tempVal>MostProbableTemplateVal)
        {
            MostProbableTemplateVal=tempVal;
            MostProbableTemplate=i;
        }
    }
    MostProbableTemplateVal*=InvPrevTotalSum;


//    if(MyAllVariables->myModelVariables.probThreshold>0.0)
//    {
////      sum=1.0/sum;
//        int Index=0;
//        NoBestMatchHaps=0;
//        for (int i=0; i<noReducedStatesCurrent;i++)
//        {
//            double tempVal=probHap[i]*InvPrevTotalSum;
//            if(tempVal >= MyAllVariables->myModelVariables.probThreshold)
//            {
//                BestMatchHaps[Index++]=i;
//                NoBestMatchHaps++;
//            }
//        }
//    }

}





void MarkovModel::unfoldProbabilitiesWithThreshold(int bridgeIndex,
                                         vector<float> &LeftNoRecomProb, vector<float> &LeftTotalProb,
                                         vector<float> &RightNoRecomProb, vector<float> &RightTotalProb,
                                         vector<float> &PrevLeftFoldedProb, vector<float> &PrevRightFoldedProb)
{
    vector<double>  &probHap = probHapMatrix[MidPoint];
    double tempInvSum = SumOfProb[MidPoint];
    double maxVal = 0.0, sum = 0.0;

    if(!MyAllVariables->myOutFormat.verbose && MyAllVariables->myModelVariables.probThreshold>0.0)
    {
        int Index=0;
        NoBestMatchHaps=0;
        for (int i=0; i<noReducedStatesCurrent;i++)
        {
            if(probHap[i]*tempInvSum >= MyAllVariables->myModelVariables.probThreshold)
            {
                BestMatchHaps[Index++]=i;
                NoBestMatchHaps++;
            }
        }
        if(NoBestMatchHaps==0)
        {
            NoBestMatchHaps=noReducedStatesCurrent;
            for (int i=0; i<noReducedStatesCurrent;i++)
            {
                BestMatchHaps[i]=i;
            }
        }
    }

    if(MyAllVariables->myOutFormat.verbose && MyAllVariables->myModelVariables.probThreshold>0.0)
    {
        int Index=0;
        NoBestMatchHaps=0;
        sum = 0.0;

        for (int i=0; i<noReducedStatesCurrent;i++)
        {
            double tempVal=probHap[i]*tempInvSum;
            if(maxVal<tempVal)
                maxVal=tempVal;
            if(tempVal >= MyAllVariables->myModelVariables.probThreshold)
            {
                sum+=tempVal;
                BestMatchHaps[Index++]=i;
                NoBestMatchHaps++;
            }
        }
        cout<<" No_Best_MatchHaps = "<<NoBestMatchHaps<<"\t "<<maxVal<<"\t"<<sum<<endl;
        if(NoBestMatchHaps==0)
        {
            cout<<" VALUE "<<noReducedStatesCurrent<<endl;
        }
    }

    ReducedHaplotypeInfo &thisInfo = rHap->ReducedStructureInfo[bridgeIndex];
    NoBestMatchFullRefHaps=0;
    double temp,tempInvCardinality;
    int i;

    for (int index=0; index<NoBestMatchHaps; index++)
    {

        i=BestMatchHaps[index];
        tempInvCardinality=thisInfo.InvuniqueCardinality[i];

        temp=LeftNoRecomProb[i];
        LeftAdj_Rec[i] = (LeftTotalProb[i] - temp ) * tempInvCardinality;
        LeftAdj_NoRrec[i] = temp / PrevLeftFoldedProb[i];

        temp=RightNoRecomProb[i];
        RightAdj_Rec[i] = (RightTotalProb[i] - temp )  * tempInvCardinality;
        RightAdj_NoRec[i] = RightNoRecomProb[i] / PrevRightFoldedProb[i];
    }


    vector<float> &LeftPrev=junctionLeftProb[bridgeIndex];
    vector<float> &RightPrev=PrevjunctionRightProb;
    int UnMappedIndex;

    for (int index=0; index<NoBestMatchHaps; index++)
    {
        i=BestMatchHaps[index];
        for(int K=0;K< thisInfo.uniqueCardinality[i];K++)
        {
            UnMappedIndex=thisInfo.uniqueIndexReverseMaps[i][K];
            probHapFullAverage[UnMappedIndex]=(LeftAdj_Rec[i] + LeftAdj_NoRrec[i]*LeftPrev[UnMappedIndex])
                                               * (RightAdj_Rec[i] + RightAdj_NoRec[i]*RightPrev[UnMappedIndex]);
            BestMatchFullRefHaps[NoBestMatchFullRefHaps++]=UnMappedIndex;
        }
    }

}



void MarkovModel::FoldBackProbabilitiesWithThreshold(ReducedHaplotypeInfo &Info)
{

    fill(FinalBestMatchfHapsIndicator.begin(), FinalBestMatchfHapsIndicator.end(), 0);
    fill(FoldedProbValue.begin(), FoldedProbValue.end(), 0.0);
    noNewReducedStates = 0;
    int UnMappedIndex, MappedIndex;

    for (int index=0; index<NoBestMatchFullRefHaps; index++)
    {
        UnMappedIndex=BestMatchFullRefHaps[index];
        MappedIndex=Info.uniqueIndexMap[UnMappedIndex];
        FinalBestMatchfHapsIndicator[MappedIndex]=1;
        FoldedProbValue[MappedIndex]+=probHapFullAverage[UnMappedIndex];
    }

    for (int i=0; i<rHapFull->maxRepSize; i++)
    {
        if(FinalBestMatchfHapsIndicator[i]==1)
            FinalBestMatchfHaps[noNewReducedStates++]=i;
    }
}





void MarkovModel::ImputeRemainingSitesbyBlock(int bridgeIndex, int TypedMarkerId,
                                              int StartPos, int EndPos)
{

    for(int position=StartPos;position<=EndPos;)
    {

        int NewBlockIndex=rHapFull->MarkerToReducedInfoMapper[position];
        ReducedHaplotypeInfo &Info=rHapFull->ReducedStructureInfo[NewBlockIndex];
        int TempEndPos= EndPos<Info.endIndex ? EndPos:Info.endIndex;

        FoldBackProbabilitiesWithThreshold (Info);
        CreatePRefPAlt(Info, position, TempEndPos);
        position=TempEndPos+1;
    }

}



void MarkovModel::CreatePRefPAlt(ReducedHaplotypeInfo &Info, int StartPos, int EndPos)
{
    int tempPosition;
    int ThisStartIndex=Info.startIndex;
    double tempInvSum = SumOfProb[MidPoint];
    int MappedIndex;

    for(tempPosition = StartPos; tempPosition<=EndPos; ++tempPosition)
    {
        vector<AlleleType> &ThisMapped=Info.TransposedUniqueHaps[tempPosition-ThisStartIndex];
        double tempPALT = 0.0;

        for (int i=0; i< noNewReducedStates;i++)
        {
            MappedIndex=FinalBestMatchfHaps[i];
            if(ThisMapped[MappedIndex] == '1')
                tempPALT+=FoldedProbValue[MappedIndex];
        }

        (*DosageHap)[tempPosition] =  min(1.0,(tempPALT * tempInvSum));

    }
}



void MarkovModel::initializeMatricesMinimac3()
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
    ThisBlockRightNoRecoProb.resize(rHap->maxBlockSize);
    probHapMatrix.resize(rHap->maxBlockSize);
    ThisBlockRightProb.resize(rHap->maxBlockSize);
    PrevRightFoldedProb.resize(rHap->maxRepSize);
    SumOfProb.resize(rHap->maxBlockSize);

    probHapMinimac3.resize(rHap->maxRepSize);


    for(int i=0;i<rHap->maxBlockSize;i++)
    {
        ThisBlockLeftNoRecoProb[i].resize(rHap->maxRepSize);
        ThisBlockRightNoRecoProb[i].resize(rHap->maxRepSize);
        ThisBlockRightProb[i].resize(rHap->maxRepSize);
        probHapMatrix[i].resize(rHap->maxRepSize);

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

    LeftPrecisionJump.resize(rHap->numMarkers+1); // RE-INTIALIZE DONE
    RightPrecisionJump.resize(rHap->numMarkers); // RE-INTIALIZE DONE

    // Full Dimension Variable Initialize (NO NEED TO REINITLAIZE)







}



void MarkovModel::initializeMatrices()
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
    ThisBlockRightNoRecoProb.resize(rHap->maxBlockSize);
    probHapMatrix.resize(rHap->maxBlockSize);
    ThisBlockRightProb.resize(rHap->maxBlockSize);
    PrevRightFoldedProb.resize(rHap->maxRepSize);
    SumOfProb.resize(rHap->maxBlockSize);


    for(int i=0;i<rHap->maxBlockSize;i++)
    {
        ThisBlockLeftNoRecoProb[i].resize(rHap->maxRepSize);
        ThisBlockRightNoRecoProb[i].resize(rHap->maxRepSize);
        ThisBlockRightProb[i].resize(rHap->maxRepSize);
        probHapMatrix[i].resize(rHap->maxRepSize);

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

    LeftPrecisionJump.resize(rHap->numMarkers+1); // RE-INTIALIZE DONE
    RightPrecisionJump.resize(rHap->numMarkers); // RE-INTIALIZE DONE

    // Full Dimension Variable Initialize (NO NEED TO REINITLAIZE)

    FoldedProbValue.resize(rHapFull->maxRepSize);
    probHapFullAverage.resize(rHapFull->numHaplotypes);
    FirstFoldedValue.resize(rHap->maxRepSize);
    UnfoldTheseSites.resize(rHapFull->maxBlockSize);

    LeftAdj_Rec.resize(rHap->maxRepSize);
    LeftAdj_NoRrec.resize(rHap->maxRepSize);
    RightAdj_Rec.resize(rHap->maxRepSize);
    RightAdj_NoRec.resize(rHap->maxRepSize);



    FinalBestMatchfHaps.resize(rHapFull->maxRepSize);
    FinalBestMatchfHapsIndicator.resize(rHapFull->maxRepSize);
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
    NoPrecisionJumps=0;
    PrevTotalSum=-1.0;
    fill(LeftPrecisionJump.begin(), LeftPrecisionJump.end(), false);
    fill(RightPrecisionJump.begin(), RightPrecisionJump.end(), false);
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

        for (int i = 0; i <noReducedStatesCurrent; i++)
        {
            NextnoRecomProb[i]*=(complement);
        }

        if (LeftPrecisionJump[markerPos])
        {
            for (int i = 0; i <noReducedStatesCurrent; i++)
            {
                NextnoRecomProb[i]*=(1e15);
            }
        }

        if (tHap.RetrieveMissingScaffoldedHaplotype(hapID,markerPos)=='0')
        {
            double freq=tHap.RetrieveScaffoldedHaplotype(hapID,markerPos)=='1'? alleleFreq[markerPos] : 1-alleleFreq[markerPos];
            double e=Error[markerPos];
            AlleleType observed=tHap.RetrieveScaffoldedHaplotype(hapID,markerPos);
            double prandom = e*freq+backgroundError;
            double pmatch = (1.0 - e)+e*freq+backgroundError;
            vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];

            for (int i = 0; i<noReducedStatesCurrent; i++)
            {
                NextnoRecomProb[i]*=(TempHap[i]==observed)?pmatch:prandom;
            }
        }
    }


}



void MarkovModel::ReCreateLeftNoRecoProbMinimac3(HaplotypeSet &tHap, int &hapID, int group,
                                                 ReducedHaplotypeInfo &Info,
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
        int TargetMarkerPosition = rHap->MapRefToTar[markerPos];

        for (int i = 0; i <noReducedStatesCurrent; i++)
        {
            NextnoRecomProb[i]*=(complement);
        }

        if (LeftPrecisionJump[markerPos])
        {
            for (int i = 0; i <noReducedStatesCurrent; i++)
            {
                NextnoRecomProb[i]*=(1e15);
            }
        }

        if (TargetMarkerPosition!=-1 && tHap.RetrieveMissingScaffoldedHaplotype(hapID,TargetMarkerPosition)=='0')
        {
            double freq=tHap.RetrieveScaffoldedHaplotype(hapID,TargetMarkerPosition)=='1'? alleleFreq[markerPos] : 1-alleleFreq[markerPos];
            double e=Error[markerPos];
            AlleleType observed=tHap.RetrieveScaffoldedHaplotype(hapID,TargetMarkerPosition);
            double prandom = e*freq+backgroundError;
            double pmatch = (1.0 - e)+e*freq+backgroundError;

            vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];

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



void MarkovModel::RightCondition(int markerPos,vector<float> &FromProb,vector<float> &ToProb,
                            vector<float> &FromNoRecomProb, vector<float> &ToNoRecomProb,
                            AlleleType observed,double e,double freq,ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];


    for (int i = 0; i<noReducedStatesCurrent; i++)
    {

        AlleleType allele=TempHap[i];
        if(allele==observed)
        {
            ToProb[i] = FromProb[i] * pmatch;
            ToNoRecomProb[i] = FromNoRecomProb[i] * pmatch;
        }
        else
        {
            ToProb[i] = FromProb[i] * prandom;
            ToNoRecomProb[i] = FromNoRecomProb[i] * prandom;
        }
    }


}



void MarkovModel::LeftCondition(int markerPos,vector<float> &Prob,
                            vector<float> &noRecomProb, AlleleType observed,double e,double freq,ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];


    for (int i = 0; i<noReducedStatesCurrent; i++)
    {

        AlleleType allele=TempHap[i];
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


bool MarkovModel::RightTranspose(vector<float> &fromTo, vector<float> &noRecomProb,
                            double reco,vector<int> &uniqueCardinality)
{
    bool tempPrecisionJumpFlag=false;
    if (reco == 0)
    {
        return false;
    }

    double sum = 0.0;
    for (int i = 0; i <noReducedStatesCurrent; i++)
    {
        sum += fromTo[i];
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
        fromTo[i]= (fromTo[i]*complement+(uniqueCardinality[i]*sum));
    }

    return tempPrecisionJumpFlag;


 }



bool MarkovModel::LeftTranspose(vector<float> &from,
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




double MarkovModel::CountErrors(int markerPos,
                                AlleleType observed,
                                double e,double freq,
                                ReducedHaplotypeInfo &Info)
{

    double match = 0;
    double mismatch = 0;
    double background = 0;
    vector<AlleleType> &TempHap = Info.TransposedUniqueHaps[markerPos-Info.startIndex];

    for (int i = 0; i < noReducedStatesCurrent; i++)
    {

        if(TempHap[i]==observed)
            match += probHapMinimac3[i];
        else
            mismatch += probHapMinimac3[i];
    }


    background = (match + mismatch) * backgroundError;
    mismatch = (match + mismatch) * e *freq;
    match *= 1.0 - e;

    return mismatch / (mismatch + match + background);
}


double MarkovModel::CountRecombinants(vector<float> &from, vector<float> &to,
                                      double r,bool PrecisionMultiply)
{
    if (r == 0)
        return 0.0;

    double fromSum = 0.0,toSum=0.0,totalSum=0.0;

    for (int i = 0; i < noReducedStatesCurrent; i++)
    {
        fromSum += from[i];
        toSum += to[i];
        totalSum+=probHapMinimac3[i];
    }

    double rsum = fromSum*r*toSum/(double)refCount;

    if(PrecisionMultiply)
        return (1e15*rsum / totalSum);
    else
        return (rsum / totalSum);
}


void MarkovModel::CreatePosteriorProb( vector<float> &Leftprob,vector<float> &rightProb,
                                       vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                       vector<float> &leftEndProb,vector<float> &rightEndProb,
                                       ReducedHaplotypeInfo &Info)
{

    for(int i=0;i<noReducedStatesCurrent;i++)
    {
        probHapMinimac3[i] = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
                             +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);
    }
}



void MarkovModel::CountExpected(int hapID,int group)
{

    vector<float> &juncLeftprob = junctionLeftProb[group];
    vector<float> &juncRightProb = PrevjunctionRightProb;
    ReducedHaplotypeInfo &Info=rHap->ReducedStructureInfo[group];
    vector<double> &alleleFreq=rHap->AlleleFreq;
    ReCreateLeftNoRecoProbMinimac3(*tHap,hapID,group,Info,alleleFreq);
    vector<vector<float> > &Leftprob = leftProb[group];
    vector<vector<float> > &leftNoRecomProb= ThisBlockLeftNoRecoProb;
    int startIndex=Info.startIndex;
    int endIndex=Info.endIndex;
    PrevRightFoldedProb = ThisBlockRightProb[endIndex-startIndex];
    ThisBlockRightNoRecoProb[endIndex-startIndex] = ThisBlockRightProb[endIndex-startIndex];
    noReducedStatesCurrent=Info.RepSize;
    fill(Constants.begin(), Constants.end(), 0.0);
    for(int i=0;i<refCount;i++)
        Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);


    CreatePosteriorProb( Leftprob[endIndex-startIndex], ThisBlockRightProb[endIndex-startIndex],
                         leftNoRecomProb[endIndex-startIndex],ThisBlockRightProb[endIndex-startIndex],
                         Leftprob[0],PrevRightFoldedProb,Info);

    int TargetMarkerPosition = rHap->MapRefToTar[endIndex];
    if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID, TargetMarkerPosition)=='0')
    {
        empError[endIndex]+=CountErrors(endIndex,
                                        tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                                        Error[endIndex],
                                        tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                                        alleleFreq[endIndex] : 1-alleleFreq[endIndex],
                                        Info);

    }
    else
        empError[endIndex]+=Error[endIndex];


    for (int markerPos=endIndex-1; markerPos>startIndex; markerPos--)
    {


        int TargetMarkerPosition = rHap->MapRefToTar[markerPos+1];
        if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID, TargetMarkerPosition)=='0')
        {
            RightCondition(markerPos+1,ThisBlockRightProb[markerPos-startIndex+1],
                           ThisBlockRightProb[markerPos-startIndex],
                           ThisBlockRightNoRecoProb[markerPos-startIndex+1],
                           ThisBlockRightNoRecoProb[markerPos-startIndex],
                           tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                           Error[markerPos+1],
                           tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                           alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }
        else
        {
            ThisBlockRightProb[markerPos-startIndex]= ThisBlockRightProb[markerPos-startIndex+1];
            ThisBlockRightNoRecoProb[markerPos-startIndex]= ThisBlockRightNoRecoProb[markerPos-startIndex+1];

        }


        empRecom[markerPos]+=CountRecombinants(Leftprob[markerPos-startIndex],
                                               ThisBlockRightProb[markerPos-startIndex],
                                               Recom[markerPos],
                                               LeftPrecisionJump[markerPos+1]);

        RightTranspose(ThisBlockRightProb[markerPos-startIndex],
                       ThisBlockRightNoRecoProb[markerPos-startIndex],
                       Recom[markerPos],
                       Info.uniqueCardinality);

        CreatePosteriorProb(Leftprob[markerPos-startIndex], ThisBlockRightProb[markerPos-startIndex],
                            leftNoRecomProb[markerPos-startIndex], ThisBlockRightNoRecoProb[markerPos-startIndex],
                            Leftprob[0],PrevRightFoldedProb,Info);

        TargetMarkerPosition = rHap->MapRefToTar[markerPos];
        if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID, TargetMarkerPosition)=='0')
        {
            empError[markerPos]+=CountErrors(markerPos,
                                            tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                                            Error[markerPos],
                                            tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                                            alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);

        }
        else
            empError[markerPos]+=Error[markerPos];

    }


     TargetMarkerPosition = rHap->MapRefToTar[startIndex+1];
    if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID,TargetMarkerPosition)=='0')
    {

        RightCondition(startIndex+1, ThisBlockRightProb[1],
                       ThisBlockRightProb[0],
                       ThisBlockRightNoRecoProb[1],
                       ThisBlockRightNoRecoProb[0],
                       tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                       Error[startIndex+1],
                       tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                       alleleFreq[startIndex+1] : 1-alleleFreq[startIndex+1],Info);

    }
    else
    {
        ThisBlockRightProb[0]= ThisBlockRightProb[1];
        ThisBlockRightNoRecoProb[0]= ThisBlockRightNoRecoProb[1];

    }

    empRecom[startIndex]+=CountRecombinants(Leftprob[0],
                                           ThisBlockRightProb[0],
                                           Recom[startIndex],
                                           LeftPrecisionJump[startIndex+1]);

    RightTranspose(ThisBlockRightProb[0],
                   ThisBlockRightNoRecoProb[0],
                   Recom[startIndex],
                   Info.uniqueCardinality);

    if(startIndex==0)
    {
        CreatePosteriorProb(Leftprob[0],
                            ThisBlockRightProb[startIndex-startIndex],
                            leftNoRecomProb[startIndex-startIndex],ThisBlockRightNoRecoProb[0],
                            Leftprob[0],PrevRightFoldedProb,Info);


        TargetMarkerPosition = rHap->MapRefToTar[startIndex];
        if (TargetMarkerPosition!=-1 && tHap->RetrieveMissingScaffoldedHaplotype(hapID, TargetMarkerPosition)=='0')
        {
            empError[startIndex]+=CountErrors(startIndex,
                                             tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition),
                                            Error[startIndex],
                                            tHap->RetrieveScaffoldedHaplotype(hapID, TargetMarkerPosition)=='1'?
                                            alleleFreq[startIndex] : 1-alleleFreq[startIndex],Info);

        }
        else
            empError[startIndex]+=Error[startIndex];
    }

}







//
//
// void MarkovModel::CreateDosages(int TypedMarkerId, int StartPos, int EndPos)
//{
//    //float Pref,Palt,ptotal;
//    double tempInvSum = SumOfProb[MidPoint];
//
//    for(int tempPosition=StartPos;tempPosition<=EndPos;++tempPosition)
//    {
//
//        //Palt=pALT[tempPosition-StartPos];
//        //Pref=pREF[tempPosition-StartPos];
//        //ptotal=Pref+Palt;
//        //(*DosageHap)[tempPosition] =  (Palt / ptotal);
//
//        (*DosageHap)[tempPosition] =  (pALT[tempPosition-StartPos] * tempInvSum);
//
//    }
//}
//

//void MarkovModel::ImputeIntermediateRegionFull(ReducedHaplotypeInfo &Info, int TypedMarkerId, int StartPos, int EndPos)
//{
//
//    CreatePRefPAlt(Info, StartPos,EndPos);
//    //CreateDosages(TypedMarkerId, StartPos,EndPos);
//
//}
