#ifndef _GLOBALPARAMETERS
#define _GLOBALPARAMETERS

#include <vector>

enum enumDLDType {
	noDLD, allButOne, halfRandom, fullRandom, all, allCombsMean, allCombsMax, first7 
};

enum enumIDType {
	noID, std01, standardDeviation, relativeOrig, range, P_sensitivity, L_diversity, T_closeness, negL_diversity, negEntropy
};


class globalParameters 
{
public:
	static const int DLDTypesCount;// = 8;

	static enum enumDLDType DLDType;

	static char* DLDnames[];
	static std::vector<std::string> DLDTypeNames;

	static const int IDTypesCount;// = 5;


	static enum enumIDType IDType;

	static char * IDnames[];
	static std::vector<std::string> IDTypeNames;

	static const double mutationRate;          // negative means 1/number of variables
	static const double perturbation;          

	static int DLDhalfRandomSampleCount;// = 50;  // number of samples with half attributes in DLD
	static int DLDfullRandomSampleCount;// = 50;  // number of samples with random attributes in DLD
	static double DLDfullRandomSampleRecordsPercent; // percent of records involved in DLD
	static bool BuildKDTreeOnMaskedData;

	static double sdidRatio;
	static int INC_QUANTIL;

	static double KDEpsilon;  // 
	static bool disclosureAwareAggregation; //
	static bool preferenceBased;
	static double prefDR;
	static double prefIL;

	static const std::string outputFolder;

	static const double defaultSafetyDistance;  // 

	static const double defaultPrecision;  // 10% of defaultSafetyDistance used for best point calculation
	static const double Epsilon;  // 
	static const double c1;// = 1;  // DR importance
	static const double c2;// = 1;  // IL importance
	static int debugLevel;// = 8;
	static const int showProgressPercent;// = 25;
	static double stdDevDLD;

	static const int confidentialValuesCount;  // = 10
	static int  kp; // = 10


	//      double prefDR = 0.2;              // preferred DR for tarragona 5
	//      double prefIL = 30;               // preferred IL for tarragona 5
	//double prefDR = 1.2;              // preferred DR for census 5
	//double prefIL = 38;               // preferred IL for census 5
	//      double prefDR = 0.3;              // preferred DR for eia 5
	//      double prefIL = 34;               // preferred IL for eia 5

};


//#else
//	
//	extern class globalParameters;
#endif