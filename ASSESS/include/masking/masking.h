#ifndef _MASKING_HEADER_
#define _MASKING_HEADER_


#include <vector>
#include <string>
#include <fstream>							// file I/O
#include <math.h>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <globalParameters.h>
#include <ANN.h>
#include <random>
#include <chrono>
#include <omp.h>


// define log(2)
#ifndef LOG2
#define LOG2 0.69315
#endif
#ifndef NULL
#define NULL 0
#endif

double get_wall_time();
double get_cpu_time();

bool loadDataFast(double **&allData,std::string datasetName, int NRecords, int NDims);
void injectSyntheticConfidentialValues(double **dataOrig, int NRecords, int NDims, int vCount);

#ifndef _fminfmax
#define fmax(a,b) (((a)>(b))?(a):(b))
#define fmin(a,b) (((a)<(b))?(a):(b))
#endif

class precalculatedStats
{
public:
	double* mean_org;
	double* vec2_org;
	double* vec4_org;
	double* varm01_org;
	double* varm2_org;

	double* correlacions_org;
	double* mat11_org;
	double* mat22_org;
	double* mat40_org;
	double* mat04_org;
	double* mat13_org;
	double* mat31_org;
	double* mat02_org;
	double* mat20_org;
	double* varm11_org;
	double* varr_org;
	double* histo_org;
	double* varq_org;
	int ** comb1;
	int ** comb2;
	int ** comb3;
	int ** comb4;
	int ** comb5;
	int ** comb6;
	int ** comb7;

	int *nchoosekSaved;

	double **origNormal;
	ANNkd_tree *kdtreeOrig;

	double *range;

	precalculatedStats();
	~precalculatedStats();
};

int *getSortedIdxFast(double *A,int n);
double **denormalizeMatrix(double **centers,int *assignment,double *meanOrig,double *stdOrig,int NRecords,int NDims);
double **denormalizeMatrix(double **centers, int *assignment, double *meanOrig, double *stdOrig, int NRecords, int NDims, double **orig);
std::vector<std::vector<double>>doubleAlloc2D_vec(int r,int c,double initialValue=0);
std::vector<std::vector<int>>intAlloc2D_vec(int r,int c,int initialValue=0);
double **normalizeMatrix(double **data,int NObjects,int NDims);
double **normalizeMatrix(double **data,double *&mean,double *&std,int NObjects,int NDims);
double** normalizeMatrix(double **data, double *&Ex, double *&stdDev, int NObjects, int NDims);
double calculateDistance2(double *x,double *y,int NDims);
double calculateDistance(double *x,double *y,int NDims);
double *calculateEx(double **data,int NObjects,int NDims);
double calculateEx(double *data,int NObjects);
double *calculateStdDevN_1(double **allData,int NObjects,int NDims);
double calculateStdDevN_1(double* allData, int NObjects);
double *calculateNEx2(double **allData,int NObjects,int NDims);
double calculateSST(double **allData,int NObjects,int NDims);
double *calculateExOnActiveData(double **data,int *active,int NObjects,int NDims);
void vector_moment_column(double** DATA, int n_records, int n_columns, int r, double* mean);
void vector_moment_central_column(double** DATA, int n_records, int n_columns, int r, double* variancies, double* mean);
void matrix_coef_correl_columns(double** DATA, int n_records, int n_columns, double* matrix, double* mean);
void matrix_moment_central_columns(double** DATA, int n_records, int n_columns, int r1, int r2, double *matrix, double* mean);
void calcular_varm01(double* varm01, double* vec2, int dimensio, int n_records);
void calcular_varm2(double* varm2, double* vec4, double* vec2, int dimensio, int n_records);
void calcular_varm11(double* varm11, double* mat22, double* mat11, int dimensio, int n_records);
void histograma_quantils(double** DATA, int n_columns, int n_records, double* histograma, double* h);
void calcular_varq(double* varq, int n_columns, int n_records, double* histo, double *h);
void crear_vector_orgrdenat(double** DATA, int col, int n_cols, int n_records, double* vector);
double p(double x);
double loss_info_vector(double* v1, double* v2, double* var, int dimensio);
double loss_info_mskatrix(double* mat1, double* mat2, double* var, int dimensio);
void loss_info_mskatrix2(double* mat1, double* mat2, double* var, int dimensio);
double loss_info_quantils(double** DATA_org, double** DATA_msk, int n_columns, int n_records, double* varq_org);
void quicksort(double* vector, int inf, int sup) ;
int *quicksortArray(double *inputArray, int n) ;
void intercanvi(double* arr, int a, int b);
void particio(double* vector, int inf, int sup, double x, int* k) ;
void cerca_dicotomica(double* vector, int mida_vector, double valor_buscat, int* pos);
double dist_records_combination(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles, int* vbles);

double calculateInformationLossOptimized(double** orig, double** mask, precalculatedStats *PS, int NRecords, int NDims) ;
double calculateDisclosureRiskOptimized(double** orig, double** mask, double** origNormal, double** maskNormal, precalculatedStats *PS, int NRecords, int NDims, int *assignment=NULL);

void calculate_dld_combination_optimizedfullRandom(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS);
void calculate_dld_combination_optimizedFirst7(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment=NULL);
void calculate_dld_combination_optimizedAllButOne(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment=NULL);
void calculate_dld_combination_optimizedhalfRandom(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment=NULL);
void calculate_dld_combination_optimizedAllCombs(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment=NULL);

void calculate_dld_combination_optimizedAll(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment) ;
void calculate_dld_combination_optimizedAll(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment, ANNkd_tree *kdtreeOrig) ;

double diversity(double *data,int n);
double entropy(double *data,int n);
double dominantValue(double *data,int n);  // calculates the dominant (mode) value in an array
void calculate_P_Sensitivity(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp, ANNkd_tree *kdtreeOrig = NULL, int *assignment = NULL);  // to implement k'-anonymity p-sensitivity
void calculate_L_Diversity(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp, ANNkd_tree *kdtreeOrig = NULL, int *assignment = NULL);  // to implement k'-anonymity l-diversity
void calculate_T_Closeness(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp, ANNkd_tree *kdtreeOrig = NULL, int *assignment = NULL);  // to implement k'-anonymity t-closeness

void calculate_negL_diversity(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp, ANNkd_tree *kdtreeOrig = NULL, int *assignment = NULL);  // to implement k'-anonymity l-diversity
void calculate_negEntropy(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp, ANNkd_tree *kdtreeOrig = NULL, int *assignment = NULL);  // to implement k'-anonymity l-diversity


int *MDAV_GenericOnNormalizedData(double** allData, int k, int NObjects, int NDims);

void generateCombinations(int n, int k, int **&combs, int &count);
void generateCombinationsRecursive(int n, int k, int **combs, int &curCombIdx, int depth, int *curComb);



std::string int2string(int number);
std::string float2string(float number);
std::string double2string(double number);



void calculate_id(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids);
double **cloneArray(double **data, int nrows, int ncols);
double *cloneArray(double *data, int n);


void print2D(double **info,int r,int c,double offset=0);
void print2D(int **info,int r,int c,int offset=0);
void print2D(FILE *fp,double **info,int r,int c,double offset=0);
void print1D(int *info,int n,int offset=0);
void print1D(double *info,int n,double offset=0);
void print1D(FILE *fp,int *info,int n,int offset=0);
void print1D(FILE *fp, double *info, int n, int offset=0);
void copy2DFrom1D(double **dest,double *src,int r,int c);
void copy2D(double **dest,double **src,int r,int c);
void copy1D(int *dest,int *src,int n);
void copy1D(double *dest,double *src,int n);
bool *boolAlloc1D(int r);
bool *boolAlloc1D(int r,bool initialValue);
double *doubleAlloc1D(int r);
double *doubleAlloc1D(int r,double initialValue);
int *intAlloc1D(int r);
int *intAlloc1D(int r,int initialValue);
void intDelete2D_Fast(int **&pa);
int **intAlloc2D_Fast(int r,int c);
int **intAlloc2D_Fast(int r,int c,int initialValue);
int **intAlloc2D_rawFast(int r,int c);
void doubleDelete2D_Fast(double **&pa);
double **doubleAlloc2D_Fast(int r,int c);
double **doubleAlloc2D_Fast(int r,int c,double initialValue);
double **doubleAlloc2D_rawFast(int r,int c);
float *floatAlloc1D(int r);
float *floatAlloc1D(int r, float initialValue);
double *maxMatrix(double **data,int rowsCount,int colsCount);
double *minMatrix(double **data,int rowsCount,int colsCount);
double **normalize_mM(double **data,long min,long max,int NRecords,int NDims);
int maxArray(int *arr,int n);
double **denormalizeNS_Simple(double **dataNormalizedFull,double *mu1,double *smu2,int NRecords,int NDims);
double **denormalizeNS(double **dataNormalized,double *mu1,double *smu2,int NRecords,int NDims);
double **normalizeNS(double **dataOriginal,double *mu1,double *smu2,int NRecords,int NDims);
long fact(int n);
int nchoosekOld(int n,int k);
unsigned nchoosek( unsigned n, unsigned k );
double* mean(double** data, int rowsCount, int colsCount) ;
double* stDev(double** data, int rowsCount, int colsCount);
double dist_records(double** data1, int i1, double** data2, int i2,int ncols);
void calculate_id_relative(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids);
double look_at_intervals2(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p);
double look_at_intervals2_relative(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p);


void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, int *assignment, int NClusters);
void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, int *assignment, int NClusters, ANNkd_tree *kdtree);
void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, ANNkd_tree *kdtree); 

double calculateDistance(double *x, double *y, int NDims);
double calculateDistance2(double *x, double *y, int NDims);
void calculate_sdid_relative(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids);
double look_at_intervals2_SDIDrelative(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p, double *SDmask);


double **calculateNewCenters(double **allData, int *assignment, int NObjects, int NDims, int NClusters) ;


precalculatedStats *precalculateAll(double **orig,int NRecords,int NDims);

void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims);
void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment);

void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment, precalculatedStats *PS, double* mu1Orig, double* smu2Orig, double **origNormalized, double* mu1Mask, double* smu2Mask, double **maskNormalized);
void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment, precalculatedStats *PS, double* mu1Orig, double* smu2Orig, double **origNormalized, double* mu1Mask, double* smu2Mask, double **maskNormalized, double IL0);
void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment, precalculatedStats *PS, double* mu1Orig, double* smu2Orig, double **origNormalized, double* mu1Mask, double* smu2Mask, double **maskNormalized, double IL0, double IL1);

void myShuffle(int *arr, int n, int shuffleCount=100);


double round(double number);

void blend(double *x, int n, int *initialAssignment, int initialNClusters, double alpha);
int *blendAndAssign(const double *x, int n, int *initialAssignment, int initialNClusters, double alpha);
double getNormalRandorm(double mean0, double std0, unsigned long seed);
int *MHM(double *allData, int k, int NRecords);
double **aggregate(double **dataPtsNorm,int *assignment, double *orgMean, double *orgStdDev, int NRecords, int NDims,int NClusters);
double **aggregateWeighted(double **dataPtsNorm,int *assignment, double *weights, double *orgMean, double *orgStdDev, int NRecords, int NDims,int NClusters);

int* array2Assignment(double* mainAllData, int* path, int k, int NRecords);
int *path2Assignment(double** allData, int* path, int k, int NRecords, int NDims);
int *calculateCenterSizes(int *assignment, int NRecords, int NClusters);
double **calculateNewCenters(double **allData, int *assignment,  int NRecords, int NDims, int NClusters, int *centerSizes);
double **calculateNewCenters(double **allData, int *assignment, int NObjects, int NDims, int NClusters);
double **calculateNewCentersDisclosureAware(double **allData, int *assignment, int NObjects, int NDims, int NClusters);
double calculateSSE(double **centers, double **allData, int *assignment, int NObjects, int NDims, int NClusters);

double *calculateExOnTour(double** allData, int *tour, int start, int count, int NObjects, int NDims);
double calculateExOnTour(double* allData, int *tour, int start, int count, int NObjects);
double calculateEx(double* allData, int start, int count, int NObjects);

double* calculateNEx2OnTour(double** allData, int *tour, int start, int count, int NObjects, int NDims);
double calculateNEx2OnTour(double* allData, int *tour, int start, int count, int NObjects);
double calculateNEx2(double* allData, int start, int count, int NObjects);

double calculateSSTOnTour2(double** allData, int *tour, int start, int count, int NObjects, int NDims, const double *Ex);
double calculateSSTOnTour2(double* allData, int *tour, int start, int count, int NObjects, const double *Ex);
double calculateSST(double* allData, int start, int count, int NObjects, const double *Ex);

int *path2Assignment(double** allData, int* path, int k, int NRecords, int NDims);
int* array2Assignment(double* mainAllData, int* path, int k, int NRecords);
double calculateGCP(double **dataPtsNorm,int *assignment,int NRecords,int NDims, int NClusters);


class indexValueClass ;
struct indexValueCompareStruct;


#ifndef NULL
#define NULL 0
#endif

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef pow2
#define pow2(a) ((a)*(a))
#endif

	//double _tmpDouble_;		//this is a global variable used in pow2C macro
	//#ifndef _pow2C_
	//#define _pow2C_
	//#define pow2C(a) ((_tmpDouble_=(a))*(_tmpDouble_))
	//#endif


#ifndef _pow2C_
#define _pow2C_
#define pow2C(a) ((a)*(a))
#endif


#ifndef DOUBLE_MAX_VALUE
#define DOUBLE_MAX_VALUE DBL_MAX
#endif

//extern enumDLDType DLDType;
//extern enumIDType IDType;
//
//extern char ** IDnames;
//extern std::vector<std::string> IDTypeNames;
//extern enumDLDType DLDType;





#endif

