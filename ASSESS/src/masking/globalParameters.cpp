#include "globalParameters.h"

enumDLDType globalParameters::DLDType = enumDLDType::all;
enumIDType globalParameters::IDType = enumIDType::relativeOrig;

const int globalParameters::DLDTypesCount = 8;
char* globalParameters::DLDnames[] = {"noDLD", "allButOne", "halfRandom", "fullRandom", "all", "allCombsMean", "allCombsMax", "first7"};
std::vector<std::string> globalParameters::DLDTypeNames(DLDnames, DLDnames + DLDTypesCount);

const int globalParameters::IDTypesCount = 10;
char* globalParameters::IDnames[] = {"noID", "std01", "standardDeviation", "relativeOrig", "range", "P_sensitivity", "L_diversity", "T_closeness","negL_diversity","negEntropy"};
std::vector<std::string> globalParameters::IDTypeNames(IDnames, IDnames + IDTypesCount);

int globalParameters::DLDhalfRandomSampleCount = 50;  // number of samples with half attributes in DLD
int globalParameters::DLDfullRandomSampleCount = 50;  // number of samples with random attributes in DLD
double globalParameters::DLDfullRandomSampleRecordsPercent = 0.1; // percent of records involved in DLD
bool globalParameters::BuildKDTreeOnMaskedData = false;   // make this true to assess FDM results

double globalParameters::sdidRatio = 5;			//ratio of standard deviation used in sdid calculation (in percent, default=5=5%), zero or negative means the average of 1% to 10% is used
int globalParameters::INC_QUANTIL = 5;

double globalParameters::KDEpsilon = 0;  //   kdtree epsilon for approximate search 
bool globalParameters::disclosureAwareAggregation = true;
bool globalParameters::preferenceBased = false;
double globalParameters::prefDR = 0;
double globalParameters::prefIL = 0;

const std::string globalParameters::outputFolder = "output10/";

const double globalParameters::defaultSafetyDistance = 0.05;  // 
const double globalParameters::defaultPrecision = 0.05;  // 10% of defaultSafetyDistance used for best point calculation
const double globalParameters::Epsilon = 1e-8;  // 

const double globalParameters::c1 = 1;  // DR importance
const double globalParameters::c2 = 1;  // IL importance
int globalParameters::debugLevel = 8;

double globalParameters::stdDevDLD = 0;

const int globalParameters::showProgressPercent = 20;

const int globalParameters::confidentialValuesCount = 10; // number of different values in the last column for p-sensitivity
int globalParameters::kp = 10; // default number of nearest neighbors considered during p-sensitivity computation

//      double prefDR = 0.2;              // preferred DR for tarragona 5
//      double prefIL = 30;               // preferred IL for tarragona 5
//double prefDR = 1.2;              // preferred DR for census 5
//double prefIL = 38;               // preferred IL for census 5
//      double prefDR = 0.3;              // preferred DR for eia 5
//      double prefIL = 34;               // preferred IL for eia 5