#ifndef __MARKOVPARAMETERS_H__
#define __MARKOVPARAMETERS_H__


#include "HaplotypeSet.h"
#include "InputFile.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include "MathVector.h"
#include "StringArray.h"


using namespace std;

class MarkovParameters
   {
   public:
    int                 noMarker;
    vector<double>      Recom,Error;
//    vector<int>      countError,countRecom;

    vector<double>      GWASRecom,GWASError;
    vector<double>      empError,empRecom;
    int                 empiricalCount;

    void CopyParameters(MarkovParameters * rhs);
        MarkovParameters  & operator +=             (const MarkovParameters & rhs);


        void                CopyParameters          (MarkovParameters & rhs);
        void                CopyParametersNew          (MarkovParameters * rhs);


        void                WriteParameters         (vector<string> &markerNames,String prefix, bool gz);
        void                WriteErrorRates         (vector<string> & markerNames, const char * filename);
        void                WriteCrossoverRates     (vector<string> & markerNames, const char * filename);
        void                ReadErrorRates          (String &filename);
        void                ReadCrossoverRates      (String &filename);
        void                UpdateModel             ();
        void ScaffoldParameters(HaplotypeSet &rHap, HaplotypeSet &tHap,MarkovParameters FromParameters);
                            MarkovParameters()
                            {
                                empiricalCount=0;
                                noMarker = 0;
                                Recom.resize(0);
                                Error.resize(0);
                                empError.resize(0);
                                empRecom.resize(0);
                            };
                            MarkovParameters(int Markers)
                            {
                                empiricalCount=0;
                                noMarker = Markers;
                                Recom.resize(Markers-1,0.001);
                                Error.resize(Markers,0.01);
                                empError.resize(Markers,0.0);
                                empRecom.resize(Markers-1,0.0);
                            };
                            ~MarkovParameters()
                           {
                           };



   };

#endif
