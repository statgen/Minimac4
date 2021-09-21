#ifndef MY_VARIABLES_INCLUDED
#define MY_VARIABLES_INCLUDED
#include<cmath>
#include<fstream>
#include "StringBasics.h"
#include<vector>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <unordered_set>
#include "assert.h"
using namespace std;




class OutputFormatVariable
{
public:


    int vcfBuffer;
    bool        unphasedOutput;
    bool GT,DS,GP,HDS,SD;
    String OutPrefix;
    string CommandLine;
    char* MyCommandLine;
    string BinaryLocation;
    bool onlyrefmarkers;
    bool gzip,RsId,nobgzip,meta;
//    vector<bool> format;
    String formatString;
    String formatStringForVCF;
    bool verbose;
    int PrintBuffer;
    bool longZero;

    bool memUsage;

    bool vcfOutput,doseOutput,hapOutput,TypedOnly;

    bool CheckValidity()
    {
        string formatPiece,formatTemp=formatString.c_str();
        char *end_str1;


        for(char * pch = strtok_r ((char*)formatTemp.c_str(),",", &end_str1);
            pch!=NULL;
            pch = strtok_r (NULL, ",", &end_str1))
        {


            formatPiece=(string)pch;

            if(formatPiece.compare("GT")==0)
                {
                    GT=true;
                }
            else if(formatPiece.compare("DS")==0)
                {
                    DS=true;
                }
            else if(formatPiece.compare("GP")==0)
                {
                    GP=true;
                }
            else if(formatPiece.compare("HDS")==0)
                {
                    HDS=true;
                }
            else if(formatPiece.compare("SD")==0)
                {
                    SD=true;
                }
            else
            {
                cout << " ERROR !!! \n Cannot identify handle for \"--format\" parameter : "<<formatPiece<<endl;
                cout << " Available handles GT, DS, HDS and GP (for genotype, dosage, haplotype dosage and posterior probability). \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
        }


        if(meta)
        {
            vcfOutput=true;
            HDS=true;
        }


        doseOutput=false;
        hapOutput=false;

        if(PrintBuffer<=100)
        {
            cout << " ERROR !!! \n Invalid input for \"--PrintBuffer\" = "<<PrintBuffer<<"\n";;
            cout << " Buffer for writing output files should be at least 1,000 characters long !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }

//        SD=true;
        bool colonIndex=false;
        if(GT)
        {
            formatStringForVCF+="GT";
            colonIndex=true;
        }
        if(DS)
        {
            formatStringForVCF+= (colonIndex?":DS":"DS");
            colonIndex=true;
        }
        if(HDS)
        {
            formatStringForVCF+= (colonIndex?":HDS":"HDS");
            colonIndex=true;
        }
        if(GP)
        {
            formatStringForVCF+= (colonIndex?":GP":"GP");
            colonIndex=true;
        }
        if(SD)
        {
            formatStringForVCF+= (colonIndex?":SD":"SD");
            colonIndex=true;
        }


        if(vcfBuffer<1)
        {
            cout << " ERROR !!! \n Invalid input for \"--vcfBuffer\" = "<<vcfBuffer<<"\n";;
            cout << " Value must be a positive integer !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }

        if(nobgzip)
            gzip=false;

        return true;

    };


    OutputFormatVariable()
    {
        memUsage=false;
        formatStringForVCF="";
        unphasedOutput=false;
        OutPrefix="Minimac4.Output";
        verbose=false;
        vcfBuffer=200;
        nobgzip=false;
        meta=false;
        PrintBuffer = 100000000;

        formatString = "GT,DS";
        GT=false;
        DS=false;
        GP=false;
        HDS=false;
        SD=false;
        longZero=false;

        hapOutput=false;
        doseOutput=false;
        vcfOutput=true;
        gzip=true;
        RsId=false;
        TypedOnly=false;

    };



void CreateCommandLine(int argc, char ** argv)
{
    int len = 0;

    for (int i=0; i<argc; i++)
        len += strlen(argv[i]) + 1;



    char MyCommandLine[len];
    strcpy(MyCommandLine,argv[0]);

    for (int i=1; i<argc; i++)
    {
        strcat(MyCommandLine, " ");
        strcat(MyCommandLine, argv[i]);
    }
    CommandLine=MyCommandLine;
}


};


class ModelVariable
{
public:

    bool        processReference,  reEstimate, updateModel ;
    double      probThreshold, diffThreshold, topThreshold;
    double constantParam;
    bool lowMemory;
    int rounds, states;
    int transFactor;
    int cisFactor ;
    int cpus;
    int minimac3;
    bool referenceEstimates;
    String intermediate;


    ModelVariable()
    {
        intermediate="";
        referenceEstimates=true;
        constantParam=0.0;
        processReference = false;
        reEstimate=false;
        updateModel = false;
        probThreshold = 0.01;
        diffThreshold = 0.01;
        topThreshold = 0.01;
        lowMemory = false;
        rounds = 5;
        states = 200;
        transFactor = 3;
        cisFactor = 2;
        cpus=1;
        #ifdef _OPENMP
            cpus=5;
        #endif
        minimac3=false;


    };
    bool CheckValidity()
    {

        if(processReference)
        {

            cout<<" NOTE: Since \"--estimate\" is ON, all options under \"Target Haplotypes\" \n";
            cout<<"       will be ignored !!!\n";
            cout<<"       Program will only estimate parameters and create a M3VCF file.\n";
            cout<<"       No imputation will be performed, hence other parameters are unnecessary !!!"<<endl<<endl;


            if(rounds<=0)
            {
                cout << " ERROR !!! \n Invalid input for \"--rounds\" = "<<rounds<<"\n";;
                cout << " Value must be POSITIVE if \"--estimate\" is ON !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
            if(states<=0)
            {
                cout << " ERROR !!! \n Invalid input for \"--states\" = "<<states<<"\n";;
                cout << " Value must be POSITIVE if \"--estimate\" is ON !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }

        }

        if(reEstimate)
        {

            processReference=reEstimate;

            cout<<" NOTE: Since \"--reEstimate\" is ON, all options under \"Target Haplotypes\" \n";
            cout<<"       will be ignored !!!\n";
            cout<<"       Program will only estimate parameters and create a M3VCF file.\n";
            cout<<"       No imputation will be performed, hence other parameters are unnecessary !!!"<<endl<<endl;


            if(rounds<=0)
            {
                cout << " ERROR !!! \n Invalid input for \"--rounds\" = "<<rounds<<"\n";;
                cout << " Value must be POSITIVE if \"--reEstimate\" is ON !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
            if(states<=0)
            {
                cout << " ERROR !!! \n Invalid input for \"--states\" = "<<states<<"\n";;
                cout << " Value must be POSITIVE if \"--reEstimate\" is ON !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }

        }
        if(constantParam>0.0)
            referenceEstimates=true;
        if(constantParam>=0.5)
        {
            cout << " ERROR !!! \n Invalid input for \"--constantParam\" = "<<constantParam<<"\n";;
            cout << " Value must be less than 0.5 !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }

        if(rounds<0)
        {
            cout << " ERROR !!! \n Invalid input for \"--rounds\" = "<<rounds<<"\n";;
            cout << " Value must be non-negative !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }
        if(states<0)
        {
            cout << " ERROR !!! \n Invalid input for \"--states\" = "<<states<<"\n";;
            cout << " Value must be non-negative !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }
        if(probThreshold<0.0 || probThreshold>=1.0)
        {
            cout << " ERROR !!! \n Invalid input for \"--probThreshold\" = "<<probThreshold<<"\n";;
            cout << " Value must be between 0.0 and 1.0 (NOT inclusive) !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }
        if(diffThreshold<0.0 || diffThreshold>=1.0)
        {
            cout << " ERROR !!! \n Invalid input for \"--diffThreshold\" = "<<diffThreshold<<"\n";;
            cout << " Value must be between 0.0 and 1.0 (NOT inclusive) !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }
        if(topThreshold<0.0 || topThreshold>=1.0)
        {
            cout << " ERROR !!! \n Invalid input for \"--topThreshold\" = "<<topThreshold<<"\n";;
            cout << " Value must be between 0.0 and 1.0 (NOT inclusive) !!! \n\n";
            cout<<" Program Exiting ..."<<endl<<endl;
            return false;
        }



        if(lowMemory)
        {
            cout<<" Low Memory Version of Minimac3 initiated  !!! \n"<<endl;
        }



        #ifdef _OPENMP
            omp_set_num_threads(cpus);
        #else
            cpus=1;
        #endif

        return true;
    };


};


class HaplotypeDataVariables
{
    public:

        int ChunkSize, ChunkWindow;
        double ChunkLengthMb, ChunkOverlapMb;

        bool passOnly, ignoreDuplicates, allowRefDuplicates;
        String MyChromosome;

        String chr;
        int start, end, window;
         String CHR;
    int START;
    int END;
    int WINDOW;
    double minRatio;
    String mapFile;
    String build;


        HaplotypeDataVariables()
        {
            ChunkLengthMb=20.0;
            ChunkSize=50000;
            ChunkWindow=1000;
            ChunkOverlapMb=3.0;
            chr = "";
            passOnly=false;
            start = 0;
            end = 0;
            window = 0;
            MyChromosome="";
            ignoreDuplicates=false;
            allowRefDuplicates=false;
            minRatio=0.1;
            mapFile="docs/geneticMapFile.b38.map.txt.gz";
        };

        bool CheckValidity()
        {
            if(ChunkLengthMb<=0.001)
            {
                cout << " ERROR !!! \n Invalid input for \"--ChunkLengthMb\" = "<<ChunkLengthMb<<"\n";;
                cout << " Chunk Length must be at least 0.001 Mb (~1Kb) long !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }

            if(ChunkLengthMb>300)
            {
                cout << " ERROR !!! \n Invalid input for \"--ChunkLengthMb\" = "<<ChunkLengthMb<<"\n";;
                cout << " Chunk Length cannot be longer than 300Mb (greater than length of chromosome 2) !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }

            if(ChunkOverlapMb<=0.001)
            {
                cout << " ERROR !!! \n Invalid input for \"--ChunkOverlapMb\" = "<<ChunkOverlapMb<<"\n";;
                cout << " Chunk Length Overlap must be at least 0.001 Mb (~1Kb) long !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }

            if(ChunkOverlapMb>300)
            {
                cout << " ERROR !!! \n Invalid input for \"--ChunkOverlapMb\" = "<<ChunkOverlapMb<<"\n";;
                cout << " Chunk Length Overlap cannot be longer than 300Mb (greater than length of chromosome 2) !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }


            if(ChunkOverlapMb>(ChunkLengthMb/3.0))
            {
                ChunkOverlapMb=(ChunkLengthMb/3.0);
                cout<<" NOTE: By Default \"--ChunkLengthMb\" must be at least 3 times the \"--ChunkOverlapMb\" \n";
                cout<<"       Value of \"--ChunkOverlapMb\" reduced to = "<< ChunkOverlapMb <<"\n";
            }

            if(minRatio<=0.0 || minRatio >= 1.0)
            {
                cout << " ERROR !!! \n Invalid input for \"--minRatio\" = "<<minRatio<<"\n";;
                cout << " Value must be strictly in between 0 and 1 !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }


            if(window<0)
            {
                cout << " ERROR !!! \n Invalid input for \"--window\" = "<<window<<"\n";;
                cout << " Value must be non-negative !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }

            if(start<0)
            {
                cout << " ERROR !!! \n Invalid input for \"--start\" = "<<start<<"\n";;
                cout << " Value must be non-negative !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
            if(end<0)
            {
                cout << " ERROR !!! \n Invalid input for \"--end\" = "<<end<<"\n";;
                cout << " Value must be non-negative !!! \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
            if(start>0 && end> 0 && start>=end)
            {
                cout << " ERROR !!! \n Invalid Input !!!\n Value of \"--start\" must be less than value of \"--end\"."<<endl;
                cout << " User Input \"--start\" = "<<start<<" and \"--end\" = " <<end<<" \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
            if(chr!="")
            {
                if(start==0)
                {
                    cout << "\n ERROR !!! \n Non-zero value of \"--start\" required parameter if using \"--chr\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
                if(end==0)
                {
                    cout << "\n ERROR !!! \n Non-zero value of \"--end\" required parameter if using \"--chr\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
                if(window==0)
                {
                    window=500000;
                    cout<<" NOTE: Default \"--window\" parameter to be used = 500000\n";
                }
            }
            if(start>0)
            {
                if(chr=="")
                {
                    cout << "\n ERROR !!! \n Missing \"--chr\", a required parameter if using \"--start\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
                if(end==0)
                {
                    cout << "\n ERROR !!! \n Non-zero value of \"--end\" required parameter if using \"--start\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
            }
            if(end>0)
            {
                if(chr=="")
                {
                    cout << "\n ERROR !!! \n Missing \"--chr\", a required parameter if using \"--end\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
                if(start==0)
                {
                    cout << "\n ERROR !!! \n Non-zero value of \"--start\" required parameter if using \"--end\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
            }
            if(window>0)
            {

                if(chr=="")
                {
                    cout << "\n ERROR !!! \n Missing \"--chr\", a required parameter if using \"--end\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
                if(start==0 && end==0)
                {
                    cout << "\n ERROR !!! \n Missing \"--start\" or  \"--end\", a required parameter if using \"--window\" parameter.\n";
                    cout<<" Program Exiting ..."<<endl<<endl;
                    return false;
                }
            }
            else
            {
                if(start>0 || end>0)
                 {
                    cout<<" NOTE: No \"--window\" parameter provided  !!! \n";
                    cout<<"       No buffer region will be used on either side of the chunk"<<endl<<endl;
                }
            }



            CHR=chr;
            START=start;
            END=end;
            WINDOW=window;
            if (CHR!="" && WINDOW > 0)
            {
                if (START-WINDOW < 0)
                    START = 0;
                else
                    START -= WINDOW;

                END += WINDOW;
            }

            return true;

        };



    void GetMapFileLocation(int argc, char ** argv)
    {
        int len = strlen(argv[0]);
        string tempString(argv[0]),tempMap;

        tempMap=tempString.substr(0,len-8)+"../docs/geneticMapFile.b38.map.txt.gz";

        mapFile=tempMap.c_str();
        int k=0;

//    char MyCommandLine;
//    strcpy(MyCommandLine,argv[0]);

    }

};


class AllVariable
{
public:

    OutputFormatVariable myOutFormat;
    ModelVariable myModelVariables;
    HaplotypeDataVariables myHapDataVariables;

};




#endif // MY_VARIABLES_INCLUDED
