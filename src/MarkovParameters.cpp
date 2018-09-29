#include "MarkovParameters.h"

using namespace std;

//
//void MarkovParameters::ScaffoldParameters(HaplotypeSet &rHap, HaplotypeSet &tHap,MarkovParameters FromParameters)
//{
//
//
//
//
//    for(int i=0;i<(int)tHap.numMarkers;i++)
//    {
//
//        FromParameters.Error[i]=Error[tHap.UnScaffoldIndex[i]];
//
//
//
//        if(i>0)
//        {
//            int index=tHap.UnScaffoldIndex[i-1];
//
//            double temp=1.0;
//            while(index<tHap.UnScaffoldIndex[i])
//            {
//                temp*=(1.0-(Recom[index]));
//                index++;
//
//            }
//            FromParameters.Recom[i-1]=(1.0-temp);
//        }
//
//
//    }
//
//    rHap.Recom=FromParameters.Recom;
//    rHap.Recom=FromParameters.Recom;
//
//}
//


void MarkovParameters::ReadErrorRates(String &filename)
{
    cout<<"\n Reading Genotyping Error values from Input File "<<filename<<endl;
    Error.clear();
    ifstream ifs(filename);
    string line;
    int count=-1;
    while(getline(ifs,line))
    {

        count++;
        if(count>0)
        {
            char *end_str1;
            char * pch=strtok_r ((char*)line.c_str(),"\t", &end_str1);;
            pch = strtok_r (NULL, "\t", &end_str1);

            double val=atof(pch);

            if(val<0.0 || val > 1.0)
            {
                cout<<" Genotpying Error is outside boundary values [0.0 - 1.0] at position "<<count<<" in Input File "<<filename<<endl;
                cout<<" Try --support for more information.\n\n";
                cout<<" Program Exiting ... \n\n";
                abort();
            }

            Error.push_back(val);
        }
    }

    if((int)Error.size()!=noMarker)
    {
        cout<<" Input File for Genotyping Error ["<<filename<<"] should have "<<noMarker<<" entries but has "<<Error.size()<<" entries."<<endl;
        cout<<" Try --support for more information.\n\n";
        cout<<" Program Exiting ... \n\n";
        abort();
    }

    ifs.close();

}

void MarkovParameters::ReadCrossoverRates(String &filename)
{
    cout<<"\n Reading Recombination Fraction values from Input File "<<filename<<endl<<endl;
    Recom.clear();
    ifstream ifs(filename);
    string line;
    int count=-1;
    while(getline(ifs,line))
    {

        count++;
        if(count>0)
        {
            char *end_str1;
            char * pch=strtok_r ((char*)line.c_str(),"\t", &end_str1);;
            pch = strtok_r (NULL, "\t", &end_str1);

            double val=atof(pch);
            if(val<0.0 || val > 1)
            {
                cout<<" Recombination Fraction is outside boundary values [0.0 - 1.0] at position "<<count<<" in Input File "<<filename<<endl;
                cout<<" Try --support for more information.\n\n";
                cout<<" Program Exiting ... \n\n";
                abort();
            }

            Recom.push_back(val);
        }
    }

    if((int)Recom.size()!=noMarker-1)
    {
        cout<<" Input File for Recombination Fraction ["<<filename<<"] should have "<<noMarker-1<<" entries but has "<<(int)Recom.size()<<" entries."<<endl;
        cout<<" Try --support for more information.\n\n";
        cout<<" Program Exiting ... \n\n";
        abort();
    }


    ifs.close();
}


void MarkovParameters::CopyParameters(MarkovParameters & rhs)
{

    if(noMarker!=rhs.noMarker)
    {
        cout<<"\n\n Run Time Error !!! Markov Parameters cannot be copied. Inconsistent number of markers ...\n";
        cout<<" Contact Author for more details ...\n";
        cout<<" Program Exiting ...\n\n";
        abort();
    }

    empiricalCount=0;
    noMarker = rhs.noMarker;
    Recom = rhs.Recom;
    Error = rhs.Error;
    empError.resize(noMarker,0.0);
    empRecom.resize(noMarker-1,0.0);



}


void MarkovParameters::CopyParameters(MarkovParameters * rhs)
{

    if(noMarker!=rhs->noMarker)
    {
        cout<<"\n\n Run Time Error !!! Markov Parameters cannot be copied. Inconsistent number of markers ...\n";
        cout<<" Contact Author for more details ...\n";
        cout<<" Program Exiting ...\n\n";
        abort();
    }

    empiricalCount=0;
    noMarker = rhs->noMarker;
    Recom = rhs->Recom;
    Error = rhs->Error;
    empError.resize(noMarker,0.0);
    empRecom.resize(noMarker-1,0.0);
}

void MarkovParameters::CopyParametersNew(MarkovParameters * rhs)
{

//    if(noMarker!=rhs->noMarker)
//    {
//        cout<<"\n\n Run Time Error !!! Markov Parameters cannot be copied. Inconsistent number of markers ...\n";
//        cout<<" Contact Author for more details ...\n";
//        cout<<" Program Exiting ...\n\n";
//        abort();
//    }
//
//    empiricalCount=0;
//    noMarker = rhs->noMarker;
    Recom = rhs->Recom;
    Error = rhs->Error;
//    empError.resize(noMarker,0.0);
//    empRecom.resize(noMarker-1,0.0);
}



MarkovParameters & MarkovParameters::operator += (const MarkovParameters & rhs)
{


    if(noMarker!=rhs.noMarker)
    {
        cout<<"\n\n Run Time Error !!! Markov Parameters cannot be copied. Inconsistent number of markers ...\n";
        cout<<" Contact Author for more details ...\n";
        cout<<" Program Exiting ...\n\n";
        abort();
    }

    empiricalCount += rhs.empiricalCount;
    for (int i = 0; i < noMarker - 1; i++)
    {
        empError[i] += rhs.empError[i];
        empRecom[i] += rhs.empRecom[i];
    }

empError[noMarker-1] += rhs.empError[noMarker-1];

    return *this;
}




void MarkovParameters::WriteParameters(vector<string> &markerNames,String prefix, bool gz)
{
    String filename=prefix;
    WriteCrossoverRates(markerNames, filename + ".rec" + (gz ? ".gz" : ""));
    cout<<"\n Recombination Rates : "<< filename + ".rec" + (gz ? ".gz" : "")<<endl;;
    WriteErrorRates(markerNames, filename + ".erate" + (gz ? ".gz" : ""));
    cout<<" Error Rates         : "<< filename + ".erate" + (gz ? ".gz" : "")<<endl;;


}

void MarkovParameters::WriteErrorRates(vector<string> & markerNames, const char * filename)
{

    IFILE output = ifopen(filename, "wb");

    if (output == NULL)
        return;

    ifprintf(output, "MarkerName\tErrorRate\n");
    for (int i = 0; i < noMarker; i++)
        ifprintf(output, "%s\t%.5g\n",  markerNames[i].c_str(), Error[i]);

    ifclose(output);
}



void MarkovParameters::WriteCrossoverRates(vector<string> & markerNames, const char * filename)
{
    IFILE output = ifopen(filename, "wb");

    if (output == NULL)
        return;

    ifprintf(output, "Interval\tSwitchRate\n");

    for (int i = 0; i < noMarker - 1; i++)
        ifprintf(output, "%s-%s\t%.5g\n", markerNames[i].c_str(),
                  markerNames[i+1].c_str(), Recom[i]);

    ifclose(output);
}




void MarkovParameters::UpdateModel()

   {
    double scale = 1.0 /(double)empiricalCount;
    double backgroundR = 0.0;
    double backgroundE = 0.0;
    int backgroundEcount = 0, backgroundRcount = 0;

    for (int i = 0; i < noMarker ; i++)
    {
        if (empError[i] < 1.0)
        {
            backgroundE += empError[i];
            backgroundEcount++;
        }

        if (i < noMarker - 1 && empRecom[i] < 2.0)
        {
            backgroundR += empRecom[i];
            backgroundRcount++;
        }
    }


    backgroundR /= (double)empiricalCount * backgroundRcount + 1e-30;
    backgroundE /= (double)empiricalCount * backgroundEcount + 1e-30;

    for (int i = 0; i < noMarker - 1; i++)
    {
        Recom[i] = empRecom[i] >=  2.0 ? empRecom[i] * scale :backgroundR;
        Error[i] = empError[i] >= 1.0 ? empError[i] * scale : backgroundE;

        empError[i]=0.0;
        empRecom[i]=0.0;

    }

    Error[noMarker-1] = empError[noMarker-1] > 1.0 ? empError[noMarker-1] * scale : backgroundE;
    empError[noMarker-1] = 0;
    empiricalCount = 0;


}


