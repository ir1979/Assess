#include <masking.h>
#include <ANN.h>

using namespace std;
//using namespace globalParameters;


class indexValueClass {
public:
	indexValueClass(){};
	double val;
	int idx;
};

struct indexValueCompareStruct {
	bool operator() (const indexValueClass &A, const indexValueClass &B) {
		return (A.val < B.val);
	}
} indexValueCompare;


//////////////////////// portable fmax/ fmin //////////////////////
#ifndef _PORTABLE
#define _PORTABLE

#ifdef _WIN32
#define fmax(a,b) (((a)>(b))?(a):(b))
#define fmin(a,b) (((a)<(b))?(a):(b))
#endif

#endif
////////////////////////// end of fmax / fmin ////////////////////////


///////////////////////// portable timing ///////////////////////////
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
	LARGE_INTEGER time,freq;
	if (!QueryPerformanceFrequency(&freq)){
		//  Handle error
		return 0;
	}
	if (!QueryPerformanceCounter(&time)){
		//  Handle error
		return 0;
	}
	return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
	FILETIME a,b,c,d;
	if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
		//  Returns total user time.
		//  Can be tweaked to include kernel times as well.
		return
			(double)(d.dwLowDateTime |
			((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
	}else{
		//  Handle error
		return 0;
	}
}
//  Posix/Linux
#else
#include <sys/time.h>
#include <time.h>
double get_wall_time(){
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}
#endif
////////////////////////// end of timing   //////////////////////////////


bool loadDataFast(double **&allData,std::string datasetName, int NRecords, int NDims)
{
	///////////////////////////
	FILE * fp = fopen(&datasetName[0],"r");
	allData = doubleAlloc2D_rawFast(NRecords,NDims);
	int curRecordIdx = 0;
	int curAttribIdx = 0;

	if (fp==NULL)
	{
		std::cout << "Cannot open " << datasetName << " file";
		return false;
	}

	char buffer[2000];
	char temp;
	int i;

	int n_data=0;
	int n_columns=0;

	i=0;
	buffer[i]=fgetc(fp);
	while(!feof(fp))// && buffer[i]!='\n' && buffer[i]!='\r') /* while not end of line */
	{
		/*-- skip initial white --*/
		i=0;
		while(buffer[i]==' ' || buffer[i]=='\t' || buffer[i]==',' || buffer[i]=='\r' || buffer[i]=='\n')
			buffer[i]=fgetc(fp);
		/*-- read the number and array of characters --*/
		while(buffer[i]!=' ' && buffer[i]!='\t' && buffer[i]!=',' && !feof(fp) && buffer[i]!='\n' && buffer[i]!='\r')
		{
			i++;
			buffer[i]=fgetc(fp);
		}
		temp=buffer[i];
		buffer[i]='\0';
		/*-- Converting to integer --*/
		if(i>0)
		{
			allData[curRecordIdx][curAttribIdx++]=atof(buffer);
			if (curAttribIdx==NDims) {
				curAttribIdx=0;
				curRecordIdx++;
			}
			n_columns++;
			n_data++;
		}

		//// I added these lines
		//while(!feof(fp)  && (temp==' ' || temp=='\t' || temp==',' || temp=='\r' || temp=='\n'))
		//	temp=fgetc(fp);

		buffer[0]=temp;
	}


	/* --- Read the rest of values ​​file --- */
	while(!feof(fp))
	{
		fscanf(fp,"%[^\r\n\t, ]%*c",buffer);
		if (curRecordIdx<NRecords) {
			allData[curRecordIdx][curAttribIdx++]=atof(buffer);
			if (curAttribIdx==NDims) {
				curAttribIdx=0;
				curRecordIdx++;
			}
		}
		n_data++;
		// I added these lines
		//do {
		//	temp = fgetc(fp);
		//} while (temp=='\r' || temp=='\n' || temp=='\t' || temp==' ' || temp==',');
		//ungetc (temp,fp);
	}



	std::fclose(fp);
	fp = 0;

	if (!(n_data==NRecords*NDims))
		return false;
	return true;
}





//double** normalizeMatrix(double **data, double *&Ex, double *&stdDev, int NObjects, int NDims) {
//	Ex = calculateEx(data, NObjects, NDims); // new double[NDims];
//	stdDev = calculateStdDevN_1(data, NObjects, NDims);
//	double **allData = doubleAlloc2D_rawFast(NObjects, NDims);
//
//	for (int j = 0; j < NDims; j++) 
//		if (stdDev[j]!=0) 
//			for (int i = 0; i < NObjects; i++)
//				allData[i][j] = (data[i][j] - Ex[j]) / stdDev[j];
//		else 
//			for (int i = 0; i < NObjects; i++)
//				allData[i][j] = (data[i][j] - Ex[j]);
//
//
//	//delete* Ex;
//	//Ex = 0;
//	//delete[] stdDev;
//	//stdDev = 0;
//	return allData;
//}


precalculatedStats::precalculatedStats()
{
	mean_org=NULL;
	vec2_org=NULL;
	vec4_org=NULL;
	varm01_org=NULL;
	varm2_org=NULL;

	correlacions_org=NULL;
	mat11_org=NULL;
	mat22_org=NULL;
	mat40_org=NULL;
	mat04_org=NULL;
	mat13_org=NULL;
	mat31_org=NULL;
	mat02_org=NULL;
	mat20_org=NULL;
	varm11_org=NULL;
	varr_org=NULL;
	histo_org=NULL;
	varq_org=NULL;
	comb1=NULL;
	comb2=NULL;
	comb3=NULL;
	comb4=NULL;
	comb5=NULL;
	comb6=NULL;
	comb7=NULL;

	nchoosekSaved=NULL;

	origNormal = NULL;
	kdtreeOrig = NULL;

	range=NULL;

}

precalculatedStats::~precalculatedStats()
{
	delete []mean_org;
	delete []vec2_org;
	delete []vec4_org;
	delete []varm01_org;
	delete []varm2_org;

	delete []correlacions_org;
	delete []mat11_org;
	delete []mat22_org;
	delete []mat40_org;
	delete []mat04_org;
	delete []mat13_org;
	delete []mat31_org;
	delete []mat02_org;
	delete []mat20_org;
	delete []varm11_org;
	delete []varr_org;
	delete []histo_org;
	delete []varq_org;

	if (comb1)
		intDelete2D_Fast(comb1);
	if (comb2)
		intDelete2D_Fast(comb2);
	if (comb3)
		intDelete2D_Fast(comb3);
	if (comb4)
		intDelete2D_Fast(comb4);
	if (comb5)
		intDelete2D_Fast(comb5);
	if (comb6)
		intDelete2D_Fast(comb6);
	if (comb7)
		intDelete2D_Fast(comb7);

	if (nchoosekSaved) {
		delete []nchoosekSaved;
	}


	if (kdtreeOrig) {
		delete kdtreeOrig;
		annClose();
	}

	if (origNormal) {
		doubleDelete2D_Fast(origNormal);
	}

	delete []range;

}



precalculatedStats *precalculateAll(double** orig, int NRecords, int NDims) 
{
	precalculatedStats *PS = new precalculatedStats();
	int i, j, index;
	double t1, t2, t3;
	double *h = new double[NDims]();//new double[1];

	int n_records_org = NRecords;
	int n_columns_org = NDims;

	PS->mean_org = new double[n_columns_org];

	PS->vec2_org = new double[n_columns_org];
	PS->vec4_org = new double[n_columns_org];

	PS->varm01_org = new double[n_columns_org];
	PS->varm2_org = new double[n_columns_org];

	PS->correlacions_org = new double[n_columns_org * n_columns_org];
	PS->mat11_org = new double[n_columns_org * n_columns_org];
	//mat22_org = new double[n_columns_org*n_columns_org];
	PS->mat40_org = new double[n_columns_org * n_columns_org];
	PS->mat04_org = new double[n_columns_org * n_columns_org];
	PS->mat13_org = new double[n_columns_org * n_columns_org];
	PS->mat31_org = new double[n_columns_org * n_columns_org];
	PS->mat22_org = new double[n_columns_org * n_columns_org];
	PS->mat02_org = new double[n_columns_org * n_columns_org];
	PS->mat20_org = new double[n_columns_org * n_columns_org];


	PS->varm11_org = new double[n_columns_org * n_columns_org];
	PS->varr_org = new double[n_columns_org * n_columns_org];

	PS->histo_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org]();
	PS->varq_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org]();

	vector_moment_column(orig, n_records_org, n_columns_org, 1, PS->mean_org);
	vector_moment_central_column(orig, n_records_org, n_columns_org, 2, PS->vec2_org, PS->mean_org);
	vector_moment_central_column(orig, n_records_org, n_columns_org, 4, PS->vec4_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 1, PS->mat11_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 2, PS->mat22_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 4, 0, PS->mat40_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 4, PS->mat04_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 3, PS->mat13_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 3, 1, PS->mat31_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 0, PS->mat20_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 2, PS->mat02_org, PS->mean_org);
	matrix_coef_correl_columns(orig, n_records_org, n_columns_org, PS->correlacions_org, PS->mean_org);

	//-------------------- Var(m01)
	calcular_varm01(PS->varm01_org, PS->vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m2)
	calcular_varm2(PS->varm2_org, PS->vec4_org, PS->vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m11)
	calcular_varm11(PS->varm11_org, PS->mat22_org, PS->mat11_org, n_columns_org, n_records_org);

	//-------------------- Var(r)

	for (i = 0; i < n_columns_org; i++) {
		for (j = 0; j < n_columns_org; j++) {
			index = i * n_columns_org + j;

			t1 = PS->mat22_org[index] / pow(PS->mat11_org[index], 2);

			t2 = PS->mat40_org[index] / pow(PS->mat20_org[index], 2) + PS->mat04_org[index] / pow(PS->mat02_org[index], 2);
			t2 += 2.0 * PS->mat22_org[index] / (PS->mat20_org[index] * PS->mat02_org[index]);
			t2 = t2 / 4.0;

			t3 = PS->mat31_org[index] / (PS->mat11_org[index] * PS->mat20_org[index]);
			t3 += PS->mat13_org[index] / (PS->mat11_org[index] * PS->mat02_org[index]);

			PS->varr_org[index] = pow(PS->correlacions_org[index], 2) / n_records_org * (t1 + t2 - t3);
		}
	}

	//-- Amount of appearances around each quantil --
	histograma_quantils(orig, n_columns_org, n_records_org, PS->histo_org, h);

	//cout << h[0] << endl;
	//cout << endl;

	//-- var(Q) --
	calcular_varq(PS->varq_org, n_columns_org, n_records_org, PS->histo_org, h);

	/*
	* --- Measures over masked data ---
	*/



	int enc = min(n_columns_org, 7);		//effective number of columns
	PS->nchoosekSaved = new int[8];
	for (int k = 1; k <= 7; k++) {
		PS->nchoosekSaved[k] = nchoosek(enc, k);
	}

	if (globalParameters::DLDType==enumDLDType::first7)
	{
		//producing all PS->combinations
		PS->comb1 = intAlloc2D_Fast(PS->nchoosekSaved[1],1);
		for (i = 0; i < PS->nchoosekSaved[1]; i++) {
			PS->comb1[i][0] = i;
		}

		PS->comb2 = intAlloc2D_Fast(PS->nchoosekSaved[2],2);

		int counter = 0;
		for (int i1 = 0; i1 < enc; i1++) {
			for (int i2 = i1 + 1; i2 < enc; i2++) {
				PS->comb2[counter][0] = i1;
				PS->comb2[counter][1] = i2;
				counter++;
			}
		}

		PS->comb3 = intAlloc2D_Fast(PS->nchoosekSaved[3],3);

		counter = 0;
		for (int i1 = 0; i1 < enc; i1++) {
			for (int i2 = i1 + 1; i2 < enc; i2++) {
				for (int i3 = i2 + 1; i3 < enc; i3++) {
					PS->comb3[counter][0] = i1;
					PS->comb3[counter][1] = i2;
					PS->comb3[counter][2] = i3;
					counter++;
				}
			}
		}

		PS->comb4 = intAlloc2D_Fast(PS->nchoosekSaved[4],4);

		counter = 0;
		for (int i1 = 0; i1 < enc; i1++) {
			for (int i2 = i1 + 1; i2 < enc; i2++) {
				for (int i3 = i2 + 1; i3 < enc; i3++) {
					for (int i4 = i3 + 1; i4 < enc; i4++) {
						PS->comb4[counter][0] = i1;
						PS->comb4[counter][1] = i2;
						PS->comb4[counter][2] = i3;
						PS->comb4[counter][3] = i4;
						counter++;
					}
				}
			}
		}

		PS->comb5 = intAlloc2D_Fast(PS->nchoosekSaved[5],5);

		counter = 0;
		for (int i1 = 0; i1 < enc; i1++) {
			for (int i2 = i1 + 1; i2 < enc; i2++) {
				for (int i3 = i2 + 1; i3 < enc; i3++) {
					for (int i4 = i3 + 1; i4 < enc; i4++) {
						for (int i5 = i4 + 1; i5 < enc; i5++) {
							PS->comb5[counter][0] = i1;
							PS->comb5[counter][1] = i2;
							PS->comb5[counter][2] = i3;
							PS->comb5[counter][3] = i4;
							PS->comb5[counter][4] = i5;
							counter++;
						}
					}
				}
			}
		}

		PS->comb6 = intAlloc2D_Fast(PS->nchoosekSaved[6],6);

		counter = 0;
		for (int i1 = 0; i1 < enc; i1++) {
			for (int i2 = i1 + 1; i2 < enc; i2++) {
				for (int i3 = i2 + 1; i3 < enc; i3++) {
					for (int i4 = i3 + 1; i4 < enc; i4++) {
						for (int i5 = i4 + 1; i5 < enc; i5++) {
							for (int i6 = i5 + 1; i6 < enc; i6++) {
								PS->comb6[counter][0] = i1;
								PS->comb6[counter][1] = i2;
								PS->comb6[counter][2] = i3;
								PS->comb6[counter][3] = i4;
								PS->comb6[counter][4] = i5;
								PS->comb6[counter][5] = i6;
								counter++;
							}
						}
					}
				}
			}
		}

		PS->comb7 = intAlloc2D_Fast(PS->nchoosekSaved[7],7);

		counter = 0;
		for (int i1 = 0; i1 < enc; i1++) {
			for (int i2 = i1 + 1; i2 < enc; i2++) {
				for (int i3 = i2 + 1; i3 < enc; i3++) {
					for (int i4 = i3 + 1; i4 < enc; i4++) {
						for (int i5 = i4 + 1; i5 < enc; i5++) {
							for (int i6 = i5 + 1; i6 < enc; i6++) {
								for (int i7 = i6 + 1; i7 < enc; i7++) {
									PS->comb7[counter][0] = i1;
									PS->comb7[counter][1] = i2;
									PS->comb7[counter][2] = i3;
									PS->comb7[counter][3] = i4;
									PS->comb7[counter][4] = i5;
									PS->comb7[counter][5] = i6;
									PS->comb7[counter][6] = i7;
									counter++;
								}
							}
						}
					}
				}
			}
		}

	} else if (globalParameters::DLDType == enumDLDType::all && globalParameters::BuildKDTreeOnMaskedData==false) {

		double *orgMean = 0, * orgStdDev = 0;

		PS->origNormal = normalizeMatrix(orig, orgMean, orgStdDev, NRecords, NDims);
		PS->kdtreeOrig = new ANNkd_tree(PS->origNormal,NRecords,NDims);

		delete []orgMean;
		delete []orgStdDev;
	} else if (globalParameters::BuildKDTreeOnMaskedData==false && (globalParameters::DLDType == enumDLDType::noDLD && (globalParameters::IDType == enumIDType::P_sensitivity || globalParameters::IDType == enumIDType::L_diversity || globalParameters::IDType == enumIDType::T_closeness ||  globalParameters::IDType == enumIDType::negL_diversity ||   globalParameters::IDType == enumIDType::negEntropy)))
	{
		double *orgMean = 0, * orgStdDev = 0;

		//double **orig_1 = doubleAlloc2D_rawFast(NRecords,NDims-1);
		//for (int i=0;i<NRecords;i++) 
		//	for (int j=0;j<NDims-1;j++)
		//		orig_1[i][j] = orig[i][j];

		PS->origNormal = normalizeMatrix(orig, orgMean, orgStdDev, NRecords, NDims);
		PS->kdtreeOrig = new ANNkd_tree(PS->origNormal,NRecords,NDims-1);

		delete []orgMean;
		delete []orgStdDev;
		//doubleDelete2D_Fast(orig_1);
	}

	PS->range = new double[n_columns_org];
	double* rangeMin = new double[n_columns_org];
	double* rangeMax = new double[n_columns_org];

	for (int iCol = 0; iCol < n_columns_org; iCol++) {
		rangeMax[iCol] = rangeMin[iCol] = orig[0][iCol];
		for (int iRow = 1; iRow < n_records_org; iRow++) {
			if (orig[iRow][iCol] < rangeMin[iCol]) {
				rangeMin[iCol] = orig[iRow][iCol];
			} else if (orig[iRow][iCol] > rangeMax[iCol]) {
				rangeMax[iCol] = orig[iRow][iCol];
			}
		}
		double iRange = rangeMax[iCol] - rangeMin[iCol];
		PS->range[iCol] = iRange > 0 ? iRange : 1;
	}

	delete []rangeMin;
	delete []rangeMax;
	return PS;
}

long fact(int n) {
	if (n <= 1) {
		return 1;
	} else {
		long result;
		result = fact(n - 1) * n;
		return result;
	}
}

int nchoosekOld(int n, int k) {
	if (n < k) 
		return 0;
	else if (k==0 || n==k)
		return 1;
	else if (k==1)
		return n;

	return (int)(fact(n) / fact(k) / fact(n - k));
}

unsigned nchoosek( unsigned n, unsigned k )
{
	if (k > n) return 0;
	if (k * 2 > n) k = n-k;
	if (k == 0) return 1;

	int result = n;
	for( int i = 2; i <= k; ++i ) {
		result *= (n-i+1);
		result /= i;
	}
	return result;
}

double** normalizeNS(double** dataOriginal, double* mu1, double* smu2, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double** dataNormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			dataNormalized[i][j] = (dataOriginal[i][j] - mu1[j]) / (smu2[j] == 0 ? 1 : smu2[j]);
		}
	}

	return dataNormalized;
}

double** denormalizeNS(double** dataNormalized, double* mu1, double* smu2, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double* cmu1 = mean(dataNormalized,NRecords,NDims);
	double* csmu2 = stDev(dataNormalized,NRecords,NDims);

	double** dataDenormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			dataDenormalized[i][j] = (dataNormalized[i][j] - cmu1[j]) * smu2[j] / (csmu2[j] == 0 ? 1 : csmu2[j]) + mu1[j];
		}
	}

	return dataDenormalized;
}

double** denormalizeNS_Simple(double** dataNormalizedFull, double* mu1, double* smu2, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double** dataDenormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			dataDenormalized[i][j] = (dataNormalizedFull[i][j]) * smu2[j] + mu1[j];
		}
	}

	return dataDenormalized;
}

int maxArray(int* arr, int n) {
	int M = arr[0];
	for (int i = 1; i < n; i++) {
		if (arr[i] > M) {
			M = arr[i];
		}
	}
	return M;
}


double maxArray(const double* arr, int n) {
	double M = arr[0];
	for (int i = 1; i < n; i++) {
		if (arr[i] > M) {
			M = arr[i];
		}
	}
	return M;
}

double** normalize_mM(double** data, long min, long max, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double** dataNormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	double* minArray = minMatrix(data, rowsCount, colsCount);
	double* maxArray = maxMatrix(data, rowsCount, colsCount);

	for (int row = 0; row < rowsCount; row++) {
		for (int col = 0; col < colsCount; col++) {
			dataNormalized[row][col] = round((data[row][col] - minArray[col]) / (maxArray[col] - minArray[col]) * (max - min)) + min;
		}
	}

	return dataNormalized;
}

double round(double number)
{
	return (number < 0.0)? ceil(number - 0.5) : floor(number + 0.5);
}


double* minMatrix(double** data, int rowsCount, int colsCount) {
	double* r = new double[colsCount];

	for (int j=0;j<colsCount;j++)
		r[j] = data[0][j];

	for (int row = 1; row < rowsCount; row++) {
		for (int col = 0; col < colsCount; col++) {
			if (data[row][col] < r[col]) {
				r[col] = data[row][col];
			}
		}
	}

	return r;
}

double* maxMatrix(double** data, int rowsCount, int colsCount) {
	double* r = new double[colsCount];

	for (int j=0;j<colsCount;j++)
		r[j] = data[0][j];

	for (int row = 1; row < rowsCount; row++) {
		for (int col = 0; col < colsCount; col++) {
			if (data[row][col] > r[col]) {
				r[col] = data[row][col];
			}
		}
	}

	return r;
}


std::vector<std::vector<double> > doubleAlloc2D_vec(int r, int c, double initialValue)
{
	return std::vector<std::vector<double> >(r,std::vector<double>(c,initialValue));
}

double **doubleAlloc2D_rawFast(int r, int c)		// allocate n pts in dim
{
	double ** pa = new double*[r];			// allocate pointers
	double  *  p = new double [r*c];		// allocate space for the array
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

double **doubleAlloc2D_Fast(int r, int c)		// allocate a 2dimensional array with r x dim zero  entries
{
	double ** pa = new double*[r];			// allocate pointers
	double  *  p = new double [r*c]();		// allocate space for the array initialized to zero
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

double **doubleAlloc2D_Fast(int r, int c, double initialValue)		// allocate a 2dimensional array with r x dim zero  entries
{
	double ** pa = new double*[r];			// allocate pointers
	double  *  p = new double [r*c];		// allocate space for the array initialized to zero
	for (int i=0;i<r*c;i++)
		p[i] = initialValue;
	for (int i = 0; i < r; i++) 
		pa[i] = &(p[i*c]);

	return pa;
}

void doubleDelete2D_Fast(double **&pa)			// deallocate points
{
	delete [] pa[0];							// dealloc reserved storage
	delete [] pa;								// dealloc pointers
	pa = NULL;
}

std::vector<std::vector<int> > intAlloc2D_vec(int r, int c, int initialValue) {
	return std::vector<std::vector<int> > (r,std::vector<int>(c,initialValue));
}

int **intAlloc2D_rawFast(int r, int c) {

	int ** pa = new int*[r];			// allocate pointers
	int  *  p = new int [r*c];			// allocate space for the array
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

int **intAlloc2D_Fast(int r, int c) {
	int ** pa = new int*[r];			// allocate pointers
	int  *  p = new int [r*c]();		// allocate space for the array
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

int **intAlloc2D_Fast(int r, int c, int initialValue) {
	int ** pa = new int*[r];			// allocate pointers
	int  *  p = new int [r*c];		// allocate space for the array
	//std::memset (p,initialValue,r*c);
	for (int i=0;i<r*c;i++)
		p[i] = initialValue;

	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;

}

void intDelete2D_Fast(int **&pa) {
	delete [] pa[0];							// dealloc reserved storage
	delete [] pa;								// dealloc pointers
	pa = NULL;
}

int *intAlloc1D(int r) {
	return new int[r]();
}

int *intAlloc1D(int r, int initialValue) {
	int *t = new int[r];
	std::memset(t,initialValue,r);

	return t;
}

double *doubleAlloc1D(int r) {
	return new double[r]();
}

double *doubleAlloc1D(int r, double initialValue) {
	double *t = new double[r];
	for (int i = 0; i < r; i++)
		t[i] = initialValue;

	return t;
}


float *floatAlloc1D(int r) {
	return new float[r]();
}

float *floatAlloc1D(int r, float initialValue) {
	float *t = new float[r];
	for (int i = 0; i < r; i++)
		t[i] = initialValue;

	return t;
}


bool *boolAlloc1D(int r) {
	return new bool[r]();
}

bool *boolAlloc1D(int r, bool initialValue) {
	bool *t = new bool[r];
	for (int i = 0; i < r; i++)
		t[i] = initialValue;

	return t;
}

void copy1D(int *dest, int *src, int n) {
	const int* s = src;
	int* d = dest;
	const int* const dend = dest + n;
	while ( d != dend )
		*d++ = *s++;	
	//std::memcpy(dest,src,n*sizeof(int));

	//for (int i = 0; i < n; i++)
	//	dest[i] = src[i];
}

void copy1D(double *dest, double *src, int n) {
	const double* s = src;
	double* d = dest;
	const double* const dend = dest + n;
	while (d != dend)
		*d++ = *s++;
	//std::memcpy(dest,src,n*sizeof(double));

	//for (int i = 0; i < n; i++)
	//	dest[i] = src[i];
}

void copy2D(double **dest, double **src, int r, int c) {
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			dest[i][j] = src[i][j];
}

void copy2DFrom1D(double **dest, double *src, int r, int c) {
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			dest[i][j] = src[i * c + j];
}

void print1D(int *info, int n, int offset) {
	for (int i = 0; i < n; i++)
		std::cout << info[i] + offset << " ";
}

void print1D(FILE *fp, int *info, int n, int offset) {
	for (int i = 0; i < n; i++)
		fprintf(fp, "%d ", info[i] + offset);
}

void print1D(FILE *fp, double *info, int n, int offset) {
	for (int i = 0; i < n; i++)
		fprintf(fp, "%lf ", info[i] + offset);
}

void print1D(double *info, int n, double offset) {
	for (int i = 0; i < n; i++)
		std::cout << info[i] + offset << " ";
}

void print2D(double **info, int r, int c, double offset) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++)
			std::cout << info[i][j] + offset << " ";
		std::printf("\n");
	}
}

void print2D(FILE *fp, double **info, int r, int c, double offset) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++)
			if (j<c-1)
				std::fprintf (fp, "%lf, ", info[i][j] + offset);
			else
				std::fprintf (fp, "%lf\n", info[i][j] + offset);

		//std::fprintf(fp, "\n");
	}
}


void print2D(int **info, int r, int c, int offset) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++)
			std::cout << info[i][j] + offset << " ";
		std::printf("\n");
	}
}

double* calculateNEx2(double** allData, int NObjects, int NDims) {
	double *NEx2 = new double[NDims];
	for (int j = 0; j < NDims; j++) {
		double sum = 0;
		for (int i = 0; i < NObjects; i++)
			sum += pow2(allData[i][j]);
		NEx2[j] = sum;
	}
	return NEx2;
}

double calculateEx(double *data, int NObjects) {
	double mean=0;
	int i;
	for (i = 0; i < NObjects; i++)
		mean += data[i];
	mean/=NObjects;
	return mean;
}


double* calculateEx(double **data, int NObjects, int NDims) {
	double *mean = new double[NDims];
	double revNObjects = 1.0 / NObjects;
	int i,j;
	for (j = 0; j < NDims; j++) {
		double sum = 0;
		for (i = 0; i < NObjects; i++)
			sum += data[i][j];
		mean[j] = sum * revNObjects;
	}
	return mean;
}

// consider rewriting
double* calculateExOnActiveData(double **data, int *active, int NObjects, int NDims) {
	double *mean = new double[NDims];

	int counter = 0;
	double sum;
	int i,j;


	sum = 0;
	for (i = 0; i < NObjects; i++){
		if (active[i]) {
			sum += data[i][0];
			counter++;
		}
	}

	double revCounter = 1.0 / counter;

	mean[0] = sum * revCounter;


	for (j = 1; j < NDims; j++) {
		sum = 0;
		for (i = 0; i < NObjects; i++)
			if (active[i])
				sum += data[i][j];
		mean[j] = sum * revCounter;
	}

	return mean;
}

double calculateSST(double** allData, int NObjects, int NDims) {
	double *Ex = calculateEx(allData, NObjects, NDims); 
	double *NEx2 = calculateNEx2(allData, NObjects, NDims); 

	double SST = 0;
	for (int j = 0; j < NDims; j++) {
		SST += NEx2[j] - NObjects * Ex[j] * Ex[j];
	}

	delete[] Ex;
	delete[] NEx2;

	return SST;
}

double* calculateStdDevN_1(double** allData, int NObjects, int NDims) {
	double *Ex = calculateEx(allData, NObjects, NDims); // new double[NDims];
	double *stdDev = new double[NDims];

	for (int j = 0; j < NDims; j++) {
		double sum = 0;
		for (int i = 0; i < NObjects; i++)
			sum += (allData[i][j] - Ex[j])*(allData[i][j] - Ex[j]);
		stdDev[j] = sqrt(sum / (NObjects - 1));
	}

	delete[] Ex;
	Ex = 0;
	return stdDev;
}

double calculateStdDevN_1(double* allData, int NObjects) {
	double Ex = calculateEx(allData, NObjects); // new double[NDims];
	double stdDev;

	double sum = 0;
	for (int i = 0; i < NObjects; i++)
		sum += (allData[i] - Ex)*(allData[i] - Ex);
	stdDev = sqrt(sum / (NObjects - 1));

	return stdDev;
}


double calculateDistance(double *x, double *y, int NDims) {
	double d = 0;
	for (int i = 0; i < NDims; i++)
		d += pow2C(x[i] - y[i]);
	return sqrt(d);
}

double calculateDistance2(double *x, double *y, int NDims) {
	double d = 0;
	for (int i = 0; i < NDims; i++)
		d += pow2C(x[i] - y[i]);
	return d;
}


double** normalizeMatrix(double **data, int NObjects, int NDims) {
	double *Ex = calculateEx(data, NObjects, NDims); // new double[NDims];
	double *stdDev = calculateStdDevN_1(data, NObjects, NDims);
	double **allData = doubleAlloc2D_rawFast(NObjects, NDims);

	for (int j = 0; j < NDims; j++) 
		if (stdDev[j]!=0) 
			for (int i = 0; i < NObjects; i++)
				allData[i][j] = (data[i][j] - Ex[j]) / stdDev[j];
		else 
			for (int i = 0; i < NObjects; i++)
				allData[i][j] = (data[i][j] - Ex[j]);


	delete[] Ex;
	Ex = 0;
	delete[] stdDev;
	stdDev = 0;
	return allData;
}

double** normalizeMatrix(double **data, double *&mean, double *&stDev, int NObjects, int NDims) {
	mean = calculateEx(data, NObjects, NDims); // new double[NDims];
	stDev = calculateStdDevN_1(data, NObjects, NDims);
	double **allData = doubleAlloc2D_rawFast(NObjects, NDims);

	if (globalParameters::DLDType == enumDLDType::noDLD && (globalParameters::IDType == enumIDType::P_sensitivity || globalParameters::IDType == enumIDType::L_diversity || globalParameters::IDType == enumIDType::T_closeness ||  globalParameters::IDType == enumIDType::negL_diversity ||    globalParameters::IDType == enumIDType::negEntropy))
	{
		for (int j = 0; j < NDims-1; j++) 
			if (stDev[j]!=0) 
				for (int i = 0; i < NObjects; i++)
					allData[i][j] = (data[i][j] - mean[j]) / stDev[j];
			else 
				for (int i = 0; i < NObjects; i++)
					allData[i][j] = (data[i][j] - mean[j]);

		for (int i = 0; i < NObjects; i++)
			allData[i][NDims-1] = data[i][NDims-1];
	}
	else
	{
		for (int j = 0; j < NDims; j++) 
			if (stDev[j]!=0) 
				for (int i = 0; i < NObjects; i++)
					allData[i][j] = (data[i][j] - mean[j]) / stDev[j];
			else 
				for (int i = 0; i < NObjects; i++)
					allData[i][j] = (data[i][j] - mean[j]);
	}
	return allData;
}

double **denormalizeMatrix(double **centers, int *assignment, double *meanOrig, double *stdOrig, int NRecords, int NDims)
{
	if (globalParameters::DLDType == enumDLDType::noDLD && (globalParameters::IDType == enumIDType::P_sensitivity || globalParameters::IDType == enumIDType::L_diversity || globalParameters::IDType == enumIDType::T_closeness ||  globalParameters::IDType == enumIDType::negL_diversity ||   globalParameters::IDType == enumIDType::negEntropy))
	{
		std::cout << "Error! cannot denormalize without knowing the last column value." << endl;
		exit(-1);
	}

	double **result = doubleAlloc2D_Fast(NRecords,NDims);

	for (int i=0;i<NRecords;i++)
		copy1D(result[i],centers[assignment[i]],NDims);

	double *Ex = calculateEx(result, NRecords, NDims); // new double[NDims];
	double *stdDev = calculateStdDevN_1(result, NRecords, NDims);

	for (int j=0;j<NDims;j++)
		if (stdDev[j]!=0)
			for (int i=0;i<NRecords;i++)
				result[i][j] = (result[i][j]-Ex[j])/stdDev[j]*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];
		else
			for (int i=0;i<NRecords;i++)
				result[i][j] = (result[i][j]-Ex[j])*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];

	delete []Ex;
	delete []stdDev;

	return result;
}

double **denormalizeMatrix(double **centers, int *assignment, double *meanOrig, double *stdOrig, int NRecords, int NDims,double **orig)
{
	double **result = doubleAlloc2D_Fast(NRecords,NDims);

	for (int i=0;i<NRecords;i++)
		copy1D(result[i],centers[assignment[i]],NDims);

	double *Ex = calculateEx(result, NRecords, NDims); // new double[NDims];
	double *stdDev = calculateStdDevN_1(result, NRecords, NDims);

	for (int j=0;j<NDims-1;j++)
		if (stdDev[j]!=0)
			for (int i=0;i<NRecords;i++)
				result[i][j] = (result[i][j]-Ex[j])/stdDev[j]*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];
		else
			for (int i=0;i<NRecords;i++)
				result[i][j] = (result[i][j]-Ex[j])*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];
	for (int i=0;i<NRecords;i++)
		result[i][NDims-1] = orig[i][NDims-1];

	delete []Ex;
	delete []stdDev;

	return result;
}


double **denormalizeMatrixWithoutNormalization(double **centers, int *assignment, double *meanOrig, double *stdOrig, int NRecords, int NDims)
{
	double **result = doubleAlloc2D_Fast(NRecords,NDims);

	for (int i=0;i<NRecords;i++)
		copy1D(result[i],centers[assignment[i]],NDims);

	for (int j=0;j<NDims;j++)
		for (int i=0;i<NRecords;i++)
			result[i][j] = (result[i][j])*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];

	return result;
}


int *getSortedIdxFast(double *A, int n) {

	//std::vector<indexValueClass> indexValue(n, indexValueClass()); 
	std::vector<indexValueClass> indexValue(n); 

	for (int idx = 0; idx < n; idx++) {
		indexValue[idx].val = A[idx];
		indexValue[idx].idx = idx;
	}

	std::sort (indexValue.begin(),indexValue.end(), indexValueCompare);

	int *result = intAlloc1D(n);

	for (int i=0;i<n;i++) {
		result[i] = indexValue[i].idx;
		A[i] = indexValue[i].val;
	}

	return result;

}



double sgn(double d) {
	if (d >= 0) {
		return 1;
	} else {
		return -1;
	}
}

void normalizeMm(double** org, double** msk, int n, int d) {
	double* minOrg = new double[d];
	double* maxOrg = new double[d];
	double* minMsk = new double[d];
	double* maxMsk = new double[d];
	for (int k = 0; k < d; k++) {
		minOrg[k] = maxOrg[k] = org[0][k];
		minMsk[k] = maxMsk[k] = msk[0][k];
	}

	for (int n1 = 1; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			if (minOrg[d1] > org[n1][d1]) {
				minOrg[d1] = org[n1][d1];
			}
			if (minMsk[d1] > msk[n1][d1]) {
				minMsk[d1] = msk[n1][d1];
			}
			if (maxOrg[d1] < org[n1][d1]) {
				maxOrg[d1] = org[n1][d1];
			}
			if (maxMsk[d1] < msk[n1][d1]) {
				maxMsk[d1] = msk[n1][d1];
			}
		}

	}

	for (int n1 = 0; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			org[n1][d1] = (org[n1][d1] - minOrg[d1]) / (maxOrg[d1] - minOrg[d1]);
			msk[n1][d1] = (msk[n1][d1] - minMsk[d1]) / (maxMsk[d1] - minMsk[d1]);
		}
	}
}

void normalizeMmGetMm(double** org, double** msk, int n, int d, double* maxOrg, double* minOrg) {
	//        double* minOrg = new double[d];
	//        double* maxOrg = new double[d];
	double* minMsk = new double[d];
	double* maxMsk = new double[d];
	for (int k = 0; k < d; k++) {
		minOrg[k] = maxOrg[k] = org[0][k];
		minMsk[k] = maxMsk[k] = msk[0][k];
	}

	for (int n1 = 1; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			if (minOrg[d1] > org[n1][d1]) {
				minOrg[d1] = org[n1][d1];
			}
			if (minMsk[d1] > msk[n1][d1]) {
				minMsk[d1] = msk[n1][d1];
			}
			if (maxOrg[d1] < org[n1][d1]) {
				maxOrg[d1] = org[n1][d1];
			}
			if (maxMsk[d1] < msk[n1][d1]) {
				maxMsk[d1] = msk[n1][d1];
			}
		}

	}

	for (int n1 = 0; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			org[n1][d1] = (org[n1][d1] - minOrg[d1]) / (maxOrg[d1] - minOrg[d1]);
			msk[n1][d1] = (msk[n1][d1] - minOrg[d1]) / (maxOrg[d1] - minOrg[d1]);
		}
	}
}

double** denormalizeMn(double** data01, int Nrecords, int NDims, double* maxMask, double* minMask) {
	double** data = doubleAlloc2D_Fast(Nrecords,NDims);

	for (int n1 = 0; n1 < Nrecords; n1++) {
		for (int d1 = 0; d1 < NDims; d1++) {
			data[n1][d1] = (maxMask[d1] - minMask[d1]) * data01[n1][d1] + minMask[d1];
		}
	}
	return data;
}

double** addArray(double** y, double** yp, int nrows, int ncols) {
	double**add = doubleAlloc2D_Fast(nrows,ncols);

	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			add[i][j] = y[i][j]+yp[i][j];
		}
	}
	return add;
}

double** addArray(double** y, double** yp, double c1, double c2, int nrows, int ncols) {

	double**add = doubleAlloc2D_Fast(nrows,ncols);

	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			add[i][j] = c1*y[i][j]+c2*yp[i][j];
		}
	}
	return add;
}

double** invertArrayPlus(double** y, int nrows, int ncols) {

	double **result =  doubleAlloc2D_Fast(nrows,ncols);
	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			result[i][j] = (y[i][j]);///(1+exp(-abs(y[i][j])));
		}
	}
	return result;
}


int knnsearch(double** data, int idx, int* exclude, int excludeCount, int nrows, int ncols) {
	double minDist = DOUBLE_MAX_VALUE;
	int nearIdx = -1;

	for (int i=0;i<nrows;i++) {
		bool flag = false;
		for (int cntr=0;cntr<excludeCount;cntr++) {
			if (i==exclude[cntr]) {
				flag=true;
				break;
			}
		}
		if (flag)
			continue;

		double tmpDist = dist_records(data, i, data, idx,ncols);

		if (tmpDist<minDist) {
			minDist = tmpDist;
			nearIdx = i;
		}
	}

	return nearIdx;
}


double** calculateCenters(double** allData, int* assignment, int NObjects, int NDims, int NClusters) {
	double** centers = doubleAlloc2D_Fast(NClusters,NDims);

	int* centerSizes = new int[NClusters];

	for (int i = 0; i < NClusters; i++) {
		centerSizes[i] = 0;
	}

	for (int i = 0; i < NObjects; i++) {
		centerSizes[assignment[i]]++;//=centerSizes[assignment[i]]+1;
	}
	for (int i = 0; i < NObjects; i++) {
		for (int j = 0; j < NDims; j++) {
			centers[assignment[i]][j] += allData[i][j];
		}
	}

	for (int i = 0; i < NClusters; i++) {
		for (int j = 0; j < NDims; j++) {
			centers[i][j] /= centerSizes[i];
		}
	}
	return centers;
}


void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims)
{	
	precalculatedStats *PS = precalculateAll(orig,NRecords,NDims);

	//double* mu1Orig = calculateEx(orig,NRecords,NDims); 
	//double* smu2Orig = calculateStdDevN_1(orig,NRecords,NDims);
	//double **origNormalized = normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	//double* mu1Mask = calculateEx(mask,NRecords,NDims);
	//double* smu2Mask = calculateStdDevN_1(mask,NRecords,NDims);
	//double **maskNormalized = normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);


	double* mu1Orig = NULL;
	double *smu2Orig = NULL;

	double **origNormalized =  normalizeMatrix(orig,mu1Orig,smu2Orig,NRecords,NDims); //  normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	double* mu1Mask = NULL; // calculateEx(mask,NRecords,NDims);
	double* smu2Mask = NULL; // calculateStdDevN_1(mask,NRecords,NDims);
	double **maskNormalized = normalizeMatrix(mask,mu1Mask,smu2Mask,NRecords,NDims); // normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);



	//#pragma omp sections
	{
		//#pragma omp section
		{
			DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims); 
		}

		//#pragma omp section
		{
			IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);
		}
	}
	delete []mu1Orig;
	delete []smu2Orig;

	delete []mu1Mask;
	delete []smu2Mask;

	doubleDelete2D_Fast(origNormalized);
	doubleDelete2D_Fast(maskNormalized);

	return;
}

void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment, precalculatedStats *PS, double* mu1Orig, double* smu2Orig, double **origNormalized, double* mu1Mask, double* smu2Mask, double **maskNormalized)
{	
	/*double* mu1Orig = calculateEx(orig,NRecords,NDims); 
	double* smu2Orig = calculateStdDevN_1(orig,NRecords,NDims);
	double **origNormalized = normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	double* mu1Mask = calculateEx(mask,NRecords,NDims);
	double* smu2Mask = calculateStdDevN_1(mask,NRecords,NDims);
	double **maskNormalized = normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);*/

	//DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims, assignment); // last arg is for speedup
	//IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);

#pragma omp parallel sections
	{
#pragma omp section
		{
			IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);
		}
#pragma omp section
		{
			DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims, assignment); 
		}
	}

	/*delete []mu1Orig;
	delete []smu2Orig;

	delete []mu1Mask;
	delete []smu2Mask;

	doubleDelete2D_Fast(origNormalized);
	doubleDelete2D_Fast(maskNormalized);*/

	return;
}

void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment, precalculatedStats *PS, double* mu1Orig, double* smu2Orig, double **origNormalized, double* mu1Mask, double* smu2Mask, double **maskNormalized, double IL0)
{	
	/*double* mu1Orig = calculateEx(orig,NRecords,NDims); 
	double* smu2Orig = calculateStdDevN_1(orig,NRecords,NDims);
	double **origNormalized = normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	double* mu1Mask = calculateEx(mask,NRecords,NDims);
	double* smu2Mask = calculateStdDevN_1(mask,NRecords,NDims);
	double **maskNormalized = normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);*/

	//DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims, assignment); // last arg is for speedup
	//IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);

	{
		{
			IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);
		}

		if (IL!=IL0)
		{
			DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims, assignment); 
		}
		else
		{
			DR = 100;
		}
	}

	/*delete []mu1Orig;
	delete []smu2Orig;

	delete []mu1Mask;
	delete []smu2Mask;

	doubleDelete2D_Fast(origNormalized);
	doubleDelete2D_Fast(maskNormalized);*/

	return;
}


void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment, precalculatedStats *PS, double* mu1Orig, double* smu2Orig, double **origNormalized, double* mu1Mask, double* smu2Mask, double **maskNormalized, double IL0, double IL1)
{	
	/*double* mu1Orig = calculateEx(orig,NRecords,NDims); 
	double* smu2Orig = calculateStdDevN_1(orig,NRecords,NDims);
	double **origNormalized = normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	double* mu1Mask = calculateEx(mask,NRecords,NDims);
	double* smu2Mask = calculateStdDevN_1(mask,NRecords,NDims);
	double **maskNormalized = normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);*/

	//DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims, assignment); // last arg is for speedup
	//IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);

	{
		{
			IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);
		}
		if (IL!=IL0 && IL!=IL1)
		{
			DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims, assignment); 
		}
		else
		{
			//cout << "*";
			DR = 100;
		}
	}

	/*delete []mu1Orig;
	delete []smu2Orig;

	delete []mu1Mask;
	delete []smu2Mask;

	doubleDelete2D_Fast(origNormalized);
	doubleDelete2D_Fast(maskNormalized);*/

	return;
}


void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims, int *assignment)
{	
	precalculatedStats *PS = precalculateAll(orig,NRecords,NDims);

	//double* mu1Orig = calculateEx(orig,NRecords,NDims); 
	//double* smu2Orig = calculateStdDevN_1(orig,NRecords,NDims);

	double* mu1Orig = NULL;
	double *smu2Orig = NULL;

	double **origNormalized =  normalizeMatrix(orig,mu1Orig,smu2Orig,NRecords,NDims); //  normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	double* mu1Mask = NULL; // calculateEx(mask,NRecords,NDims);
	double* smu2Mask = NULL; // calculateStdDevN_1(mask,NRecords,NDims);
	double **maskNormalized = normalizeMatrix(mask,mu1Mask,smu2Mask,NRecords,NDims); // normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);

	DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims,assignment); // last two args are for speedup
	IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);

	delete []mu1Orig;
	delete []smu2Orig;

	delete []mu1Mask;
	delete []smu2Mask;

	doubleDelete2D_Fast(origNormalized);
	doubleDelete2D_Fast(maskNormalized);

	return;
}

double calculateInformationLossOptimized(double** orig, double** mask, precalculatedStats *PS, int NRecords, int NDims) 
{
	int n_columns_org, n_records_org;
	int n_columns_msk, n_records_msk;
	//        double* h = new double[1];

	double* mean_org = PS->mean_org;
	double* mean_msk;
	double* vec2_org = PS->vec2_org;
	double* vec2_msk;
	//        double* vec4_org=PS->vec4_org;
	double* varm01_org = PS->varm01_org;
	double* varm2_org = PS->varm2_org;

	double* correlacions_org = PS->correlacions_org;
	double* correlacions_msk;
	double* mat11_org = PS->mat11_org;
	double* mat11_msk;
	//        double* mat22_org=PS->mat22_org;
	//        double* mat40_org=PS->mat40_org;
	//        double* mat04_org=PS->mat04_org;
	//        double* mat13_org=PS->mat13_org;
	//        double* mat31_org=PS->mat31_org;
	//        double* mat02_org=PS->mat02_org;
	//        double* mat20_org=PS->mat20_org;
	double* varm11_org = PS->varm11_org;
	double* varr_org = PS->varr_org;
	//        double* histo_org=PS->histo_org;
	double* varq_org = PS->varq_org;

	//        int i, j, index;
	//        double t1, t2, t3;

	//        double* aux_msk_org;
	//        double* aux_v_org;
	//        double* aux_msk_msk;
	//        double* aux_v_msk;

	/*
	* Check for proper number of arguments
	*/

	n_records_org = NRecords;
	n_columns_org = NDims;
	//n_DATA_org = n_records_org * n_columns_org;

	n_records_msk = NRecords;
	n_columns_msk = NDims;
	//n_DATA_msk = n_records_msk * n_columns_msk;


	//mean_org = new double[n_columns_org];

	//vec2_org = new double[n_columns_org];
	//vec4_org = new double[n_columns_org];

	//varm01_org = new double[n_columns_org];
	//varm2_org = new double[n_columns_org];

	//correlacions_org = new double[n_columns_org * n_columns_org];
	//mat11_org = new double[n_columns_org * n_columns_org];
	//mat22_org = new double[n_columns_org*n_columns_org];
	//mat40_org = new double[n_columns_org * n_columns_org];
	//mat04_org = new double[n_columns_org * n_columns_org];
	//mat13_org = new double[n_columns_org * n_columns_org];
	//mat31_org = new double[n_columns_org * n_columns_org];
	//mat22_org = new double[n_columns_org * n_columns_org];
	//mat02_org = new double[n_columns_org * n_columns_org];
	//mat20_org = new double[n_columns_org * n_columns_org];


	//varm11_org = new double[n_columns_org * n_columns_org];
	//varr_org = new double[n_columns_org * n_columns_org];

	//histo_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org];
	//varq_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org];

	//        vector_moment_column(orig, n_records_org, n_columns_org, 1, mean_org);
	//        vector_moment_central_column(orig, n_records_org, n_columns_org, 2, vec2_org, mean_org);
	//        vector_moment_central_column(orig, n_records_org, n_columns_org, 4, vec4_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 1, mat11_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 2, mat22_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 4, 0, mat40_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 4, mat04_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 3, mat13_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 3, 1, mat31_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 0, mat20_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 2, mat02_org, mean_org);
	//        matrix_coef_correl_columns(orig, n_records_org, n_columns_org, correlacions_org, mean_org);

	//-------------------- Var(m01)
	//        calcular_varm01(varm01_org, vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m2)
	//        calcular_varm2(varm2_org, vec4_org, vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m11)
	//        calcular_varm11(varm11_org, mat22_org, mat11_org, n_columns_org, n_records_org);

	//-------------------- Var(r)

	//        for (i = 0; i < n_columns_org; i++) {
	//            for (j = 0; j < n_columns_org; j++) {
	//                index = i * n_columns_org + j;
	//
	//                t1 = mat22_org[index] / pow(mat11_org[index], 2);
	//
	//                t2 = mat40_org[index] / pow(mat20_org[index], 2) + mat04_org[index] / pow(mat02_org[index], 2);
	//                t2 += 2.0 * mat22_org[index] / (mat20_org[index] * mat02_org[index]);
	//                t2 = t2 / 4.0;
	//
	//                t3 = mat31_org[index] / (mat11_org[index] * mat20_org[index]);
	//                t3 += mat13_org[index] / (mat11_org[index] * mat02_org[index]);
	//
	//                varr_org[index] = pow(correlacions_org[index], 2) / n_records_org * (t1 + t2 - t3);
	//            }
	//        }

	//-- Amount of appearances around each quantil --
	//        histograma_quantils(orig, n_columns_org, n_records_org, histo_org, h);

	//-- var(Q) --
	//        calcular_varq(varq_org, n_columns_org, n_records_org, histo_org, h[0]);

	/*
	* --- Measures over masked data ---
	*/

	mean_msk = new double[n_columns_msk];
	vec2_msk = new double[n_columns_msk];
	mat11_msk = new double[n_columns_msk * n_columns_msk];
	correlacions_msk = new double[n_columns_msk * n_columns_msk];

	vector_moment_column(mask, n_records_msk, n_columns_msk, 1, mean_msk);
	vector_moment_central_column(mask, n_records_msk, n_columns_msk, 2, vec2_msk, mean_msk);
	matrix_moment_central_columns(mask, n_records_msk, n_columns_msk, 1, 1, mat11_msk, mean_msk);

	matrix_coef_correl_columns(mask, n_records_msk, n_columns_msk, correlacions_msk, mean_msk);

	//mexPrintf("\nProbabilistic comparison of files %s and %s (PIL(Q) PIL(m^0_1) PIL(m_2) PIL(m_{11}) PIL(r))\n",argv[1],argv[2]);

	//-- Quantiles information loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on quantiles PIL(Q): ");
	double PIL1 = loss_info_quantils(orig, mask, n_columns_org, n_records_org, varq_org);
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Means info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on means PIL(m^0_1): ");
	double PIL2 = loss_info_vector(mean_org, mean_msk, varm01_org, n_columns_org);
	//double PIL2 = 0;
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Variancies info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on variances PIL(m_2): ");
	double PIL3 = loss_info_vector(vec2_org, vec2_msk, varm2_org, n_columns_org);
	//double PIL3 = 0;
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Covariances info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on covariances PIL(m_{11}): ");
	double PIL4 = loss_info_mskatrix(mat11_org, mat11_msk, varm11_org, n_columns_org);
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Correlation info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on Pearson's correlations PIL(r): ");
	double PIL5 = loss_info_mskatrix(correlacions_org, correlacions_msk, varr_org, n_columns_org);

	//if (DETAILED_OUTPUT) mexPrintf("\n");

	if (globalParameters::debugLevel>=10) {
		std::printf("\n PIL1: %f\tPIL2: %f\tPIL3: %f\tPIL4: %f\tPIL5: %f\t", PIL1, PIL2, PIL3, PIL4, PIL5);
	}
	//        std::printf (" IL: %f\t",(PIL1 + PIL2 + PIL3 + PIL4 + PIL5) / 5.0);

	delete []mean_msk;
	delete []vec2_msk;
	delete []mat11_msk;
	delete []correlacions_msk;

	return (PIL1 + PIL2 + PIL3 + PIL4 + PIL5) / 5.0;
	//if (DETAILED_OUTPUT) mexPrintf ("PIL: %lf\%", G_IL);
}


double moment_column(double** DATA, int n_records, int col, int n_columns, int r) {
	int i;
	double acum = 0.0;

	if (r == 1) {
		for (i = 0; i < n_records; i++) {
			acum += DATA[i][col];
		}
	} else {
		for (i = 0; i < n_records; i++) {
			acum += pow(DATA[i][col], r);
		}
	}
	acum /= n_records;
	return (acum);
}

//---------------------------------------------------
// Calcula r-essim moment central de la variable col
// Testejat amb Excel
//---------------------------------------------------
double moment_central_column(double** DATA, int n_records, int col, int n_columns, int r, double mean) {
	int i;
	//double mitja;
	double acum = 0.0;

	//mitja=moment_column(DATA,n_records,col,n_columns,1);

	for (i = 0; i < n_records; i++) {
		acum += pow(DATA[i][col] - mean, r);
	}
	acum /= n_records;
	return (acum);
}

//--------------------------------------------------------
// Calcula r1-r2 moment central de les variables col1, col2
// Testejat amb Excel
//--------------------------------------------------------
double moment_central_columns(double** DATA, int n_records, int col1, int col2, int n_columns, int r1, int r2, double* mean) {
	int i;
	/*
	* double mitja1, mitja2;
	*/
	double value1;
	double value2;
	double acum = 0.0;
	double resultat;

	//mitja1=moment_column(DATA,n_records,col1,n_columns,1);
	//mitja2=moment_column(DATA,n_records,col2,n_columns,1);

	for (i = 0; i < n_records; i++) {
		value1 = pow(DATA[i][col1] - mean[col1], r1);
		value2 = pow(DATA[i][col2] - mean[col2], r2);
		acum += (value1 * value2);
	}

	resultat = acum / n_records;
	return (resultat);
}

//-------------------------------------------------------------
// Calcula el coef. de correlacio de les variables col1, col2
// Testejat amb Excel
//-------------------------------------------------------------
double coef_correl_columns(double** DATA, int n_records, int col1, int col2, int n_columns, double* mean) {
	double mu11, mu20, mu02;
	double resultat;

	mu11 = moment_central_columns(DATA, n_records, col1, col2, n_columns, 1, 1, mean);
	mu20 = moment_central_columns(DATA, n_records, col1, col2, n_columns, 2, 0, mean);
	mu02 = moment_central_columns(DATA, n_records, col1, col2, n_columns, 0, 2, mean);

	resultat = mu11 / sqrt(mu20 * mu02);
	return (resultat);
}

//-------------------------------------------------------------
// Omple un vector amb la mitjana de cada column
//-------------------------------------------------------------
void vector_moment_column(double** DATA, int n_records, int n_columns, int r, double* mean) {
	int col;

#pragma omp parallel for
	for (col = 0; col < n_columns; col++) {
		mean[col] = moment_column(DATA, n_records, col, n_columns, r);
	}
}

//-------------------------------------------------------------
// Omple un vector amb la varian\ca de cada column
//-------------------------------------------------------------
void vector_moment_central_column(double** DATA, int n_records, int n_columns, int r, double* variancies, double* mean) {
	int col;
#pragma omp parallel for
	for (col = 0; col < n_columns; col++) {
		variancies[col] = moment_central_column(DATA, n_records, col, n_columns, r, mean[col]);
	}
}

//------------------------------------------------------------------------------
// Omple una matrix amb els r1-r2 moments centrals de cada parella de columns
//------------------------------------------------------------------------------
void matrix_moment_central_columns(double** DATA, int n_records, int n_columns, int r1, int r2, double *matrix, double* mean) {
	int col1, col2;

	for (col1 = 0; col1 < n_columns; col1++) {
		for (col2 = 0; col2 < n_columns; col2++) {
			matrix[col1 * n_columns + col2] = moment_central_columns(DATA, n_records, col1, col2, n_columns, r1, r2, mean);
		}
	}
}


//------------------------------------------------------------------------------
// Omple una matrix amb els coef. de correlacio de cada parella de columns
//------------------------------------------------------------------------------
void matrix_coef_correl_columns(double** DATA, int n_records, int n_columns, double* matrix, double* mean) {
	int col1, col2;

	for (col1 = 0; col1 < n_columns; col1++) {
		for (col2 = 0; col2 < n_columns; col2++) {
			matrix[col1 * n_columns + col2] = coef_correl_columns(DATA, n_records, col1, col2, n_columns, mean);
		}
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varm01(double* varm01, double* vec2, int dimensio, int n_records) {
	int i;

	for (i = 0; i < dimensio; i++) {
		varm01[i] = vec2[i] / (double) n_records;
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varm2(double* varm2, double* vec4, double* vec2, int dimensio, int n_records) {
	int i;

	for (i = 0; i < dimensio; i++) {
		varm2[i] = (vec4[i] - pow(vec2[i], 2)) / (double) n_records;
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varm11(double* varm11, double* mat22, double* mat11, int dimensio, int n_records) {
	int i, j;
	int index;

	for (i = 0; i < dimensio; i++) {
		for (j = 0; j < dimensio; j++) {
			index = i * dimensio + j;
			varm11[index] = (mat22[index] - pow(mat11[index], 2)) / (double) n_records;
		}
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void histograma_quantils(double** DATA, int n_columns, int n_records, double* histograma, double* h) {
	double* vector;
	int i, j, k;
	double min, max;
	double quantil;
	int appearances;

	vector = new double[n_records];

	//--Per each variable
	for (i = 0; i < n_columns; i++) {
		crear_vector_orgrdenat(DATA, i, n_columns, n_records, vector);

		min = vector[0];
		max = vector[n_records - 1];

		//h[0] = abs(max - min) / 1000.0;
		h[i] = abs(max - min) / 1000.0;

		//cout << "max/min" << max << " " << min << " h[0] " << h[0] << endl;


		//-- Quantiles for each (globalParameters::INC_QUANTIL%) of the variable (I do not think the first or last)
		for (j = globalParameters::INC_QUANTIL; j <= (100 - globalParameters::INC_QUANTIL); j += globalParameters::INC_QUANTIL) {
			quantil = vector[(n_records - 1) * j / 100];

			appearances = 0;
			for (k = 0; k < n_records; k++) {
				if (vector[k] > (quantil - (h[i])) && vector[k] < (quantil + (h[i]))) {
					appearances++;
				} else if (appearances>0) {
					break;
					//k = k;
				}
			}

			histograma[n_columns * (j / globalParameters::INC_QUANTIL) + i] = (double) appearances;
			//	  mexPrintf("%d ",j/globalParameters::INC_QUANTIL);
			//      mexPrintf("%g %d\n",quantil, appearances);

		}
		//      mexPrintf("\n\n");

	}
	delete []vector;
	vector = 0;
}

//-------------------------------------------------------------
//-------------------------------------------------------------
//void calcular_varq(double* varq, int n_columns, int n_records, double* histo, double h) {
void calcular_varq(double* varq, int n_columns, int n_records, double* histo, double *h) {
	int j, i;
	int index;

	for (i = 0; i < n_columns; i++) {
		for (j = globalParameters::INC_QUANTIL; j <= (100 - globalParameters::INC_QUANTIL); j += globalParameters::INC_QUANTIL) {
			index = n_columns * (j / globalParameters::INC_QUANTIL) + i;
			//varq[index] = j * (100 - j) * n_records * pow(2 * h, 2) / (pow(histo[index], 2) * 100.0 * 100.0);
			varq[index] = j * (100 - j) * n_records * pow2((2.0 * h[i])/(histo[index]*100.0));
			//varq[index] = pow2(sqrt(j * (100 - j) * n_records) * (2.0 * h[i])/(histo[index]*100.0));

			if (varq[index]<EPSILON1) {
				//cout << "\n\n** j " << j << "  j * (100 - j)  " <<  j * (100 - j) << ",   n_records  " << n_records << ",  h[i] " << h[i] << "  histo[index]" <<  histo[index]  << ", h/hist" << h[i] / histo[index] <<", pow2((2.0 * h[i])/(histo[index]*100.0)) " << pow((2.0 * h[i])/(histo[index]*100.0),2)  <<  "\t  j * (100 - j) * n_records * pow2(2.0 * h[i]/histo[index]*100.0) " << varq[index] << ",   index" << index << endl;
				varq[index] = 0;
			}
		}
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void crear_vector_orgrdenat(double** DATA, int col, int n_cols, int n_records, double* vector) {
	int i;

	/*
	* --- Primer copio la column sobre el vector ---
	*/
	for (i = 0; i < n_records; i++) {
		vector[i] = DATA[i][col];
	}

	/*
	* --- La ordeno ---
	*/

	std::sort(vector,vector+n_records);
	//quicksort(vector, 0, n_records - 1);


	/*
	* --- Verifico la ordenacio (aixo es podra treure) ---
	*/
	/*
	* for(i=0;i<(n_records-1);i++) if(vector[i]>vector[i+1]) {
	* mexPrintf("\nVector mal ordenat. Verificar quicksort!!!\n"); exit(0);
	* }
	*/
}

double p(double x) {
	double b1 = 0.319381530;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double pi = 3.14159265;

	double t, z, value;

	t = 1 / (1 + p * x);

	z = exp(-pow(x, 2) / 2.0) / sqrt(2.0 * pi);

	value = 0.5 - z * (b1 * t + b2 * pow(t, 2) + b3 * pow(t, 3) + b4 * pow(t, 4) + b5 * pow(t, 5));

	return (value);
}

//-­------------------------------------------------------
//--------------------------------------------------------
double loss_info_vector(double* v1, double* v2, double* var, int dimensio) {
	int i;
	double* perdues;
	double acum = 0.0;

	perdues = new double[dimensio];

	for (i = 0; i < dimensio; i++) {
		perdues[i] = 2 * 100.0 * p(abs(v1[i] - v2[i]) / sqrt(var[i]));
		acum += perdues[i];
	}

	delete []perdues;

	return acum / (double) dimensio;
}

//-­------------------------------------------------------
//--------------------------------------------------------
double loss_info_mskatrix(double* mat1, double* mat2, double* var, int dimensio) {
	int i, j, index;
	double* perdues;
	double acum = 0.0;
	int n_DATA = 0;

	perdues = new double[dimensio * dimensio];

	for (i = 0; i < dimensio; i++) {
		for (j = 0; j < dimensio; j++) {
			if (i != j) {
				index = i * dimensio + j;
				perdues[index] = 2 * 100.0 * p(abs(mat1[index] - mat2[index]) / sqrt(var[index]));

				acum += perdues[index];
				n_DATA++;
			}
		}
	}

	delete []perdues;
	return acum / (double) n_DATA;
}

//-­------------------------------------------------------
//--------------------------------------------------------
void loss_info_mskatrix2(double* mat1, double* mat2, double* var, int dimensio) {
	int i, j, index;
	double* perdues;
	double acum = 0.0;
	int n_DATA = 0;

	perdues = new double[dimensio * dimensio];

	for (i = 0; i < dimensio; i++) {
		for (j = 0; j < dimensio; j++) {
			if (i != j) {
				index = i * dimensio + j;
				perdues[index] = 50.0 * abs(mat1[index] - mat2[index]);

				acum += perdues[index];
				n_DATA++;
			}
		}
	}
}

//-­------------------------------------------------------
//--------------------------------------------------------
double loss_info_quantils(double** DATA_org, double** DATA_msk, int n_columns, int n_records, double* varq_org) {
	double* vector_org;
	double* vector_msk;
	double dada_org, dada_msk, var_org;
	double prob;
	int i, j;
	double acum = 0.0;
	double count = 0;

	vector_org = new double[n_records];
	vector_msk = new double[n_records];

	for (i = 0; i < n_columns; i++) {
		crear_vector_orgrdenat(DATA_org, i, n_columns, n_records, vector_org);
		crear_vector_orgrdenat(DATA_msk, i, n_columns, n_records, vector_msk);

		for (j = globalParameters::INC_QUANTIL; j <= (100 - globalParameters::INC_QUANTIL); j += globalParameters::INC_QUANTIL) {

			dada_org = vector_org[(n_records - 1) * j / 100];
			dada_msk = vector_msk[(n_records - 1) * j / 100];
			var_org = varq_org[n_columns * (j / globalParameters::INC_QUANTIL) + i];

			prob = 2 * 100.0 * p(abs(dada_org - dada_msk) / sqrt(var_org));
			if (_isnan(prob) || prob!=prob)   // prob is indeterminant, i.e., divided by zero
				prob = 0;
			else
			{			
				acum += prob;
				count++;
			}
			//	  mexPrintf("%g %g %g %g -- ",dada_org,dada_msk,var_org,abs(dada_org-dada_msk));
			//	  mexPrintf("%g\n",prob);
		}
		//      mexPrintf("\n");
	}

	delete []vector_org;
	delete []vector_msk;
	double result = 0;
	if (count>0)
		result = acum / (double) count;

	return result;
}

/*
* ---------------------------- quicksort -------------------------------
*/
void quicksort(double* vector, int inf, int sup) {
	int* k = new int[1];
	if (inf <= sup) {
		particio(vector, inf + 1, sup, vector[inf], k);
		intercanvi(vector, inf, k[0]);
		quicksort(vector, inf, k[0] - 1);
		quicksort(vector, k[0] + 1, sup);
	}
}
/*
* ----------------------------- intercanvi -----------------------------
*/

void intercanvi(double* arr, int a, int b) {
	double temp;
	temp = arr[a];
	arr[a] = arr[b];
	arr[b] = temp;
}
/*
* ----------------------------- particio --------------------------------
*/

void particio(double* vector, int inf, int sup, double x, int* k) {
	int k2;

	k[0] = inf - 1;
	k2 = sup + 1;
	while (k2 != (k[0] + 1)) {
		if (vector[k[0] + 1] <= x) {
			(k[0])++;
		} else if (vector[k2 - 1] >= x) {
			k2--;
		} else {
			intercanvi(vector, k[0] + 1, k2 - 1);
			(k[0])++;
			k2--;
		}
	}
}

/*
* ------------------------- cerca dicotomica -----------------------
*/
void cerca_dicotomica(double* vector, int mida_vector, double valor_buscat, int* pos) {
	int inf, sup, mig;

	inf = 0;
	sup = mida_vector - 1;

	while (inf != (sup + 1)) {
		mig = (sup + inf) / 2;

		if (vector[mig] <= valor_buscat) {
			inf = mig + 1;
		} else {
			sup = mig - 1;
		}
	}
	pos[0] = sup;
}


double calculateDisclosureRiskOptimized(double** orig, double** mask, double** origNormal, double** maskNormal, precalculatedStats *PS, int NRecords, int NDims, int *assignment) {

	int n_data, n_columns, n_regs;
	//int i,j;
	double* ids = new double[10];
	double dld, id;
	double DLDstdDev_1 = 0;

	n_regs = NRecords;
	n_columns = NDims;
	n_data = n_regs * n_columns;

	int* match_all = new int[n_regs];

	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
#pragma omp parallel for
	for (int kk = 0; kk < n_regs; kk++) {
		match_all[kk] = kk;
	}


#pragma omp parallel sections
	{
		{
			if (globalParameters::DLDType!=enumDLDType::noDLD) 
			{

				if (globalParameters::DLDType==enumDLDType::first7) {
					double* dlds = new double[7];
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (all combination of the first 7 attributes)...";

					calculate_dld_combination_optimizedFirst7(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS,assignment);
					dld = calculateEx(dlds,min(7, n_columns));
					globalParameters::stdDevDLD = calculateStdDevN_1(dlds,min(7, n_columns));
					delete []dlds;

				} else if (globalParameters::DLDType==enumDLDType::all) {
					double* dlds = new double[1];
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (all attributes)...";

					if (globalParameters::BuildKDTreeOnMaskedData || assignment==NULL || PS->kdtreeOrig==NULL) { // if to search orig "on mask", i.e., the kdtree is build on the masked file
						if (!globalParameters::BuildKDTreeOnMaskedData)
							std::cerr << "*** Warning, globalParameters::BuildKDTreeOnMaskedData is false, but we have to set it to true. ***";
						calculate_dld_combination_optimizedAll (origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
					} else {
						calculate_dld_combination_optimizedAll (origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment, PS->kdtreeOrig);
					}

					dld = dlds[0];
					globalParameters::stdDevDLD=0;

					delete []dlds;

				} else if (globalParameters::DLDType==enumDLDType::allButOne) {
					double* dlds = new double[NDims];

					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (all but one attribute)...";

					calculate_dld_combination_optimizedAllButOne(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS,assignment);
					dld = calculateEx(dlds,NDims);
					globalParameters::stdDevDLD = calculateStdDevN_1(dlds,NDims);

					delete []dlds;
				} else if (globalParameters::DLDType==enumDLDType::allCombsMean) {
					double* dlds = new double[(1<<NDims)-1];
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (Mean of all combinations of attributes)...";

					calculate_dld_combination_optimizedAllCombs(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS,assignment);

					dld = calculateEx(dlds,(1<<NDims)-1);
					globalParameters::stdDevDLD = calculateStdDevN_1(dlds,(1<<NDims)-1);

					delete []dlds;

				} else if (globalParameters::DLDType==enumDLDType::allCombsMax) {
					double* dlds = new double[(1<<NDims)-1];
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (Max of all combinations of attributes)...";

					calculate_dld_combination_optimizedAllCombs (origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS,assignment);

					dld = maxArray(dlds,(1<<NDims)-1);
					globalParameters::stdDevDLD = calculateStdDevN_1(dlds,(1<<NDims)-1);

					delete []dlds;

				} else if (globalParameters::DLDType==enumDLDType::halfRandom) {
					double* dlds = new double[globalParameters::DLDhalfRandomSampleCount];
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (half random attributes)...";

					calculate_dld_combination_optimizedhalfRandom(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS,assignment);

					dld = calculateEx (dlds,globalParameters::DLDhalfRandomSampleCount);
					globalParameters::stdDevDLD = calculateStdDevN_1(dlds,globalParameters::DLDhalfRandomSampleCount);

					//dld = 0;
					//for (int j = 0; j < globalParameters::DLDhalfRandomSampleCount; j++) {
					//	/*std::printf ("\n\tDLD[%d] = %lf \t", j, dlds[j]); ///delme*/
					//	dld += dlds[j];
					//}
					//dld/=globalParameters::DLDhalfRandomSampleCount;
					delete []dlds;

				} else if (globalParameters::DLDType==enumDLDType::fullRandom) {
					double* dlds = new double[globalParameters::DLDfullRandomSampleCount];
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing DLD (full random attributes)...";

					calculate_dld_combination_optimizedfullRandom(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS);

					dld = calculateEx (dlds,globalParameters::DLDfullRandomSampleCount);
					globalParameters::stdDevDLD = calculateStdDevN_1(dlds,globalParameters::DLDfullRandomSampleCount);

					//dld = 0;
					//for (int j = 0; j < globalParameters::DLDfullRandomSampleCount; j++) {
					//	/*std::printf ("\n\tDLD[%d] = %lf \t", j, dlds[j]); ///delme*/
					//	dld += dlds[j];
					//}
					//dld/=globalParameters::DLDfullRandomSampleCount;
					delete []dlds;

				}

				if (globalParameters::debugLevel>=10 && DLDstdDev_1!=0)
					if (globalParameters::debugLevel>=10)
						std::cout << "OK (DLD stddev: " << DLDstdDev_1 << ")." << endl;
			}
			else {
				dld = 0;
				globalParameters::stdDevDLD = 0;
			}
		}

#pragma omp section
		{

			if (globalParameters::IDType!=enumIDType::noID) 
			{
				if (globalParameters::IDType==enumIDType::std01) {
					//Normalize using max and min

					double **orig01 = doubleAlloc2D_Fast(n_regs, n_columns);
					double **mask01 = doubleAlloc2D_Fast(n_regs, n_columns);

					copy2D(orig01,orig,n_regs,n_columns);
					copy2D(mask01,mask,n_regs,n_columns);

					normalizeMm(orig01, mask01, n_regs, n_columns);

					//calculate id on 0-1 normalized data
					calculate_id(orig01, mask01, n_data, n_columns, match_all, ids);

					doubleDelete2D_Fast (orig01);
					doubleDelete2D_Fast (mask01);
				} else if (globalParameters::IDType==enumIDType::standardDeviation) {
					calculate_sdid_relative(orig, mask, n_data, n_columns, match_all, ids);
				} else if (globalParameters::IDType==enumIDType::relativeOrig)
				{
					calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);
				} else if (globalParameters::IDType==enumIDType::P_sensitivity) {
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing ID (P-Sensitivity)...";
					if (globalParameters::BuildKDTreeOnMaskedData || assignment==NULL || PS->kdtreeOrig==NULL)  // if to search orig "on mask", i.e., the kdtree is build on the masked file
					{
						if (!globalParameters::BuildKDTreeOnMaskedData)
							std::cerr << "*** Warning, globalParameters::BuildKDTreeOnMaskedData is false, but we have to set it to true. ***";
						calculate_P_Sensitivity (origNormal, maskNormal, n_data, n_columns, match_all, ids , globalParameters::kp); // kdtree is built on masked data
					} else {
						calculate_P_Sensitivity (origNormal, maskNormal, n_data, n_columns, match_all, ids, globalParameters::kp, PS->kdtreeOrig, assignment);
					}
				}
				else if (globalParameters::IDType==enumIDType::L_diversity) {
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing ID (L-Diversity)...";
					if (globalParameters::BuildKDTreeOnMaskedData || assignment==NULL || PS->kdtreeOrig==NULL)  // if to search orig "on mask", i.e., the kdtree is build on the masked file
					{
						if (!globalParameters::BuildKDTreeOnMaskedData)
							std::cerr << "*** Warning, globalParameters::BuildKDTreeOnMaskedData is false, but we have to set it to true. ***";
						calculate_L_Diversity (origNormal, maskNormal, n_data, n_columns, match_all, ids , globalParameters::kp); // kdtree is built on masked data
					} else {
						calculate_L_Diversity (origNormal, maskNormal, n_data, n_columns, match_all, ids, globalParameters::kp, PS->kdtreeOrig, assignment);
					}
				}
				else if (globalParameters::IDType==enumIDType::negL_diversity) {
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing ID (revL-Diversity)...";
					if (globalParameters::BuildKDTreeOnMaskedData || assignment==NULL || PS->kdtreeOrig==NULL)  // if to search orig "on mask", i.e., the kdtree is build on the masked file
					{
						if (!globalParameters::BuildKDTreeOnMaskedData)
							std::cerr << "*** Warning, globalParameters::BuildKDTreeOnMaskedData is false, but we have to set it to true. ***";
						calculate_negL_diversity (origNormal, maskNormal, n_data, n_columns, match_all, ids , globalParameters::kp); // kdtree is built on masked data
					} else {
						calculate_negL_diversity (origNormal, maskNormal, n_data, n_columns, match_all, ids, globalParameters::kp, PS->kdtreeOrig, assignment);
					}
				}
				else if (globalParameters::IDType==enumIDType::negEntropy) {
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing ID (revL-Diversity)...";
					if (globalParameters::BuildKDTreeOnMaskedData || assignment==NULL || PS->kdtreeOrig==NULL)  // if to search orig "on mask", i.e., the kdtree is build on the masked file
					{
						if (!globalParameters::BuildKDTreeOnMaskedData)
							std::cerr << "*** Warning, globalParameters::BuildKDTreeOnMaskedData is false, but we have to set it to true. ***";
						calculate_negEntropy (origNormal, maskNormal, n_data, n_columns, match_all, ids , globalParameters::kp); // kdtree is built on masked data
					} else {
						calculate_negEntropy (origNormal, maskNormal, n_data, n_columns, match_all, ids, globalParameters::kp, PS->kdtreeOrig, assignment);
					}
				}
				else if (globalParameters::IDType==enumIDType::T_closeness) {
					if (globalParameters::debugLevel>=10)
						std::cout << "Computing ID (T-Closeness)...";
					if (globalParameters::BuildKDTreeOnMaskedData || assignment==NULL || PS->kdtreeOrig==NULL)  // if to search orig "on mask", i.e., the kdtree is build on the masked file
					{
						if (!globalParameters::BuildKDTreeOnMaskedData)
							std::cerr << "*** Warning, globalParameters::BuildKDTreeOnMaskedData is false, but we have to set it to true. ***";
						calculate_T_Closeness (origNormal, maskNormal, n_data, n_columns, match_all, ids , globalParameters::kp); // kdtree is built on masked data
					} else {
						calculate_T_Closeness (origNormal, maskNormal, n_data, n_columns, match_all, ids, globalParameters::kp, PS->kdtreeOrig, assignment);
					}
				} else {
					cerr << "No suitable method for interval disclosure is chosen." << endl;
					exit(-1);
				}

				id = 0;
				for (int j = 0; j < 10; j++) {
					if (globalParameters::debugLevel>=10)
						std::printf ("\n\tIDS[%d] = %lf \t", j, ids[j]); 
					id += ids[j];
				}
				id /= 10.0;
				if (globalParameters::debugLevel>=10)
					std::cout << "OK." << endl;
			} else {
				id = 0;
			}
		}
	}
	///////////////////

	if (globalParameters::debugLevel>=10) {
		std::printf ("\tdld: %f ", dld);
		std::printf ("\tid: %f\n", id);
	}

	delete []ids;
	delete []match_all;

	double dRisk;

	if (globalParameters::DLDType==enumDLDType::noDLD && !globalParameters::IDType==enumIDType::noID)
		dRisk = id;
	else if (!globalParameters::DLDType==enumDLDType::noDLD && !globalParameters::IDType==enumIDType::noID)
		dRisk = (dld + id) / 2.0;
	else if (!globalParameters::DLDType==enumDLDType::noDLD && globalParameters::IDType==enumIDType::noID)
		dRisk = dld;
	else
		dRisk = 100;

	return dRisk;

}


//double calculateDisclosureRiskOptimized(double** orig, double** mask, double** origNormal, double** maskNormal, precalculatedStats *PS, int NRecords, int NDims, int *assignment) {
//
//	int n_data, n_columns, n_regs;
//	//int i,j;
//	double* ids = new double[10];
//	double dld, id;
//	double DLDstdDev_1 = 0;
//
//	n_regs = NRecords;
//	n_columns = NDims;
//	n_data = n_regs * n_columns;
//
//	int* match_all = new int[n_regs];
//
//	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
//	for (int kk = 0; kk < n_regs; kk++) {
//		match_all[kk] = kk;
//	}
//
//
//
//	if (DLDType!=noDLD) {
//
//		if (DLDType==first7) {
//			double* dlds = new double[7];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (all combination of the first 7 attributes)...";
//
//			calculate_dld_combination_optimizedFirst7(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
//			dld = 0;
//			for (int j = 0; j < min(7, n_columns); j++) {
//				/*std::printf ("\n\tDLD[%d] = %lf \t", j, dlds[j]); ///delme*/
//				dld += dlds[j];
//			}
//			dld /= min(7, n_columns);
//			delete []dlds;
//
//		} else if (DLDType==all) {
//			double* dlds = new double[1];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (all attributes)...";
//
//			calculate_dld_combination_optimizedAll (origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
//			dld = dlds[0];
//
//			delete []dlds;
//
//		} else if (DLDType==allButOne) {
//			double* dlds = new double[NDims];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (all but one attribute)...";
//
//			calculate_dld_combination_optimizedAllButOne(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
//			dld = 0;
//			for (int j = 0; j < NDims; j++) {
//				/*std::printf ("\n\tDLD[%d] = %lf \t", j, dlds[j]); ///delme*/
//				dld += dlds[j];
//			}
//			dld/=NDims;
//			delete []dlds;
//
//		} else if (DLDType==allCombsMean) {
//			double* dlds = new double[(1<<NDims)-1];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (Mean of all combinations of attributes)...";
//
//			calculate_dld_combination_optimizedAllCombs(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
//			dld = calculateEx(dlds,(1<<NDims)-1);
//
//			delete []dlds;
//
//		} else if (DLDType==allCombsMax) {
//			double* dlds = new double[(1<<NDims)-1];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (Max of all combinations of attributes)...";
//
//			calculate_dld_combination_optimizedAllCombs (origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
//			dld = maxArray(dlds,(1<<NDims)-1);
//
//			delete []dlds;
//		} else if (DLDType==halfRandom) {
//			double* dlds = new double[globalParameters::DLDhalfRandomSampleCount];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (half random attributes)...";
//
//			calculate_dld_combination_optimizedhalfRandom(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS, assignment);
//
//			dld = calculateEx (dlds,globalParameters::DLDhalfRandomSampleCount);
//			DLDstdDev_1 = calculateStdDevN_1(dlds,globalParameters::DLDhalfRandomSampleCount);
//
//			//dld = 0;
//			//for (int j = 0; j < globalParameters::DLDhalfRandomSampleCount; j++) {
//			//	/*std::printf ("\n\tDLD[%d] = %lf \t", j, dlds[j]); ///delme*/
//			//	dld += dlds[j];
//			//}
//			//dld/=globalParameters::DLDhalfRandomSampleCount;
//			delete []dlds;
//
//		} else if (DLDType==fullRandom) {
//			double* dlds = new double[globalParameters::DLDfullRandomSampleCount];
//			if (globalParameters::debugLevel>=10)
//				std::cout << "Computing DLD (full random attributes)...";
//
//			calculate_dld_combination_optimizedfullRandom(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS);
//
//			dld = calculateEx (dlds,globalParameters::DLDfullRandomSampleCount);
//			DLDstdDev_1 = calculateStdDevN_1(dlds,globalParameters::DLDfullRandomSampleCount);
//
//			//dld = 0;
//			//for (int j = 0; j < globalParameters::DLDfullRandomSampleCount; j++) {
//			//	/*std::printf ("\n\tDLD[%d] = %lf \t", j, dlds[j]); ///delme*/
//			//	dld += dlds[j];
//			//}
//			//dld/=globalParameters::DLDfullRandomSampleCount;
//			delete []dlds;
//
//		}
//
//
//		if (globalParameters::debugLevel>=10)
//			std::cout << "OK (" << DLDstdDev_1 << ")." << endl;
//	}
//	else {
//		dld = 0;
//	}
//
//
//	if (IDType!=noID) {
//		if (globalParameters::standard01) {
//			//Normalize using max and min
//
//			double **orig01 = doubleAlloc2D_Fast(n_regs, n_columns);
//			double **mask01 = doubleAlloc2D_Fast(n_regs, n_columns);
//
//			copy2D(orig01,orig,n_regs,n_columns);
//			copy2D(mask01,mask,n_regs,n_columns);
//
//			normalizeMm(orig01, mask01, n_regs, n_columns);
//
//			//calculate id on 0-1 normalized data
//			calculate_id(orig01, mask01, n_data, n_columns, match_all, ids);
//
//			doubleDelete2D_Fast (orig01);
//			doubleDelete2D_Fast (mask01);
//		} else if (globalParameters::sdid) {
//			calculate_sdid_relative(orig, mask, n_data, n_columns, match_all, ids);
//		} else
//		{
//			calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);
//		}
//
//		//        calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);
//
//		id = 0;
//		for (int j = 0; j < 10; j++) {
//			if (globalParameters::debugLevel>=10)
//				std::printf ("\n\tIDS[%d] = %lf \t", j, ids[j]); ///delme
//			id += ids[j];
//		}
//		id /= 10.0;
//		if (globalParameters::debugLevel>=10)
//			std::cout << "OK." << endl;
//	}else
//	{
//		id = 0;
//	}
//	///////////////////
//
//	if (globalParameters::debugLevel>=9) {
//		std::printf ("\tdld: %f  ", dld);
//		std::printf ("\tid: %f  ", id);
//	}
//
//	delete []ids;
//	delete []match_all;
//
//	double dRisk;
//
//	if (globalParameters::ignoreDLD && !globalParameters::ignoreIDR)
//		dRisk = id;
//	else if (!globalParameters::ignoreDLD && !globalParameters::ignoreIDR)
//		dRisk = (dld + id) / 2.0;
//	else if (!globalParameters::ignoreDLD && globalParameters::ignoreIDR)
//		dRisk = dld;
//	else
//		dRisk = 100;
//
//	return dRisk;
//
//}


//double calculateDisclosureRiskOptimizedRelative(double** orig, double** mask, double** origNormal, double** maskNormal, precalculatedStats *PS, int nrows, int ncols) {
//
//	int n_data, n_columns, n_regs;
//	//int i,j;
//	double* dlds = new double[7];
//	double* ids = new double[10];
//	double dld, id;
//
//	n_regs = nrows;
//	n_columns = ncols;
//	n_data = n_regs * n_columns;
//
//	int* match_all = new int[n_regs];
//
//	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
//	for (int kk = 0; kk < n_regs; kk++) {
//		match_all[kk] = kk;
//	}
//
//	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
//	//if(DEBUG)printf("matching done\n");
//
//	//producing all combinations
//	//int enc = min(n_columns, 7);		//effective number of columns
//	//int** comb1 = PS->comb1;
//
//	//int** comb2 = PS->comb2;
//
//	//int** comb3 = PS->comb3;
//
//	//int** comb4 = PS->comb4;
//
//	//int** comb5 = PS->comb5;
//
//	//int** comb6 = PS->comb6;
//
//	//int** comb7 = PS->comb7;
//
//	//-- with standardized data
//	//        for (int ii=0;ii<3;ii++)
//	//        {
//	//            for (int jj=0;jj<5;jj++)
//	//                std::printf("%f, ", origNormal[ii][jj]);
//	//            System.out.println();
//	//        }
//	//
//	//        System.out.println();
//	//        System.out.println();
//	//
//	//        for (int ii=0;ii<3;ii++)
//	//        {
//	//            for (int jj=0;jj<5;jj++)
//	//                std::printf("%f, ", maskNormal[ii][jj]);
//	//            System.out.println();
//	//        }
//
//
//	calculate_dld_combination_optimized(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS);
//	//calculate_dld(data,mask,n_data,n_columns,match_all,dlds);
//
//	/*
//	* //-- I reset the data without normalize for(j=0;j<n_data;j++)
//	* data[j]=orig[j]; for(j=0;j<n_data;j++) mask[j]=mask_copy[j];
//	*
//	* //Normalize using max and min
//	* normalizeMm(data,mask,n_regs,n_columns);
//	* calculate_id(data,mask,n_data,n_columns,match_all,ids);
//	* if(DEBUG)printf("id calculated\n");
//	*/
//
//	//        for (int ii=0;ii<20;ii++)
//	//        {
//	//            for (int jj=0;jj<2;jj++)
//	//                std::printf("%f ", origNormal[ii][jj]);
//	//            std::printf ("\n");
//	//        }
//	//
//	//        std::printf ("\n");
//	//        std::printf ("\n");
//	//                
//	//        for (int ii=0;ii<20;ii++)
//	//        {
//	//            for (int jj=0;jj<2;jj++)
//	//                std::printf("%f ", maskNormal[ii][jj]);
//	//            std::printf ("\n");
//	//        }
//
//	if (globalParameters::standard01) {
//		//Normalize using max and min
//		double** orig01 = cloneArray(orig,n_regs, n_columns);
//		double** mask01 = cloneArray(mask,n_regs, n_columns);
//
//		normalizeMm(orig01, mask01, n_regs, n_columns);
//
//
//		//calculate id on 0-1 normalized data
//		calculate_id(orig01, mask01, n_data, n_columns, match_all, ids);
//	} else {
//		calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);
//		//            calculate_id_relative(origNormal, maskNormal, n_data, n_columns, match_all, ids);
//	}
//
//
//	dld = 0;
//	for (int j = 0; j < min(7, n_columns); j++) {
//		dld += dlds[j];
//	}
//
//	id = 0;
//	for (int j = 0; j < 10; j++) {
//		id += ids[j];
//	}
//	id /= 10.0;
//
//	///////////////////
//
//	if (globalParameters::debugLevel>=10) {
//		std::printf("dld: %f\t", dld);
//		std::printf("id: %f\t", id);
//	}
//
//	return (dld + id) / 2.0;
//
//}

void do_matching_combination(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}

	/*
	* -- For all registration date2 --
	*/
	for (i = 0; i < n_regs1; i++) {
		min_distancia = 999999;
		quin = -1;

		for (j = 0; j < n_regs2; j++) {
			distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
			//printf ("d(%d,%d)=%f ",i,j,distancia);
			if (quin == -1) {
				min_distancia = distancia;
				quin = j;
			} else {
				if (distancia < min_distancia) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
					min_distancia = distancia;
					quin = j;
				}
			}
		}

		if (quin != -1) {
			//                if (i==0)
			//                    std::printf("%d\t%d\t%g\n",quin,i,min_distancia);
			resultat[i] = quin;
		} else {
			resultat[i] = -i;
		}
	}
}

void do_matching_combination_optimized(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}

	/*
	* -- For all registration date2 --
	*/
	for (i = 0; i < n_regs1; i++) {
		min_distancia = 999999;
		quin = -1;

		for (j = 0; j < n_regs2; j++) {
			//distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
			if (quin == -1) {
				min_distancia = distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
				quin = j;
			} else {
				//if ((distancia = dist_records_combination_with_threshold(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles, min_distancia)) < min_distancia) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
				if ((distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles)) < min_distancia) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
					min_distancia = distancia;
					quin = j;
				}
			}
		}

		if (quin != -1) {
			//                if (i==0)
			//                    std::printf("%d\t%d\t%g\n",quin,i,min_distancia);
			resultat[i] = quin;
		} else {
			resultat[i] = -i;
		}
	}
}

void do_matching_combination_optimized2(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}



	/*
	* -- For all registration date2 --
	*/
	for (i = 0; i < n_regs1; i++) {
		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		for (j = 0; j < n_regs2; j++) {
			if (j == i) {
				continue;
			}
			distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
			if ((distancia < min_distancia) || (distancia == min_distancia && j < i)) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
				quin = -1;
				break;
			}
		}

		if (quin != -1) {
			//                if (i==0)
			//                    std::printf("%d\t%d\t%g\n",quin,i,min_distancia);
			resultat[i] = quin;
		} else {
			resultat[i] = -1;
		}
	}
}

void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) 
{
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/

	double **tmpData = NULL;
	double *tmpQueryPoint = NULL;

	if (n_columns2!=n_vbles) {
		tmpData = doubleAlloc2D_rawFast(n_regs1,n_vbles);
		for (i = 0; i < n_regs1; i++) {
			resultat[i] = -1;

			for (int j=0;j<n_vbles;j++) {
				tmpData[i][j] = data2[i][vbles[j]];
			}
		}
		tmpQueryPoint = new double[n_vbles];
	} else {
		tmpData = data2;
		for (i = 0; i < n_regs1; i++) {
			resultat[i] = -1;
		}
	}



	ANNkd_tree *kdtree = new ANNkd_tree(tmpData,n_regs2,n_vbles);

	double minDist[50];
	int minIdx[50];

	/*
	* -- For all registration date2 --
	*/

	for (i = 0; i < n_regs1; i++) {
		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		if (n_columns2!=n_vbles) {
			for (int j=0;j<n_vbles;j++)
				tmpQueryPoint[j] = data1[i][vbles[j]];
		} else {
			tmpQueryPoint = data1[i];
		}


		kdtree->annkSearch(tmpQueryPoint,1,minIdx,minDist, globalParameters::KDEpsilon);


		if (minDist[0]<min_distancia) {
			resultat[i] = -1;
			continue;
		}

		int found = kdtree->annkFRSearch(tmpQueryPoint,min_distancia,50,minIdx,minDist, globalParameters::KDEpsilon);


		for (int j=0;j<found;j++) {
			if (minIdx[j]<i || minDist[j]<min_distancia) {
				quin = -1;
				break;
			}
		}

		resultat[i] = quin;
	}


	if (n_columns2!=n_vbles) {
		delete []tmpQueryPoint;
		doubleDelete2D_Fast(tmpData);
	} else {
		tmpQueryPoint = NULL;
		tmpData = NULL;
	}

	delete kdtree;
	annClose();

}

void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, ANNkd_tree *kdtree) 
{
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/

	double **tmpData = NULL;
	double *tmpQueryPoint = NULL;

	tmpData = data2;
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}


	double minDist[50];
	int minIdx[50];

	/*
	* -- For all registration date2 --
	*/

	for (i = 0; i < n_regs1; i++) {
		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		tmpQueryPoint = data1[i];

		kdtree->annkSearch(tmpQueryPoint,1,minIdx,minDist,globalParameters::KDEpsilon);


		if (minDist[0]<min_distancia) {
			resultat[i] = -1;
			continue;
		}

		int found = kdtree->annkFRSearch(tmpQueryPoint,min_distancia,50,minIdx,minDist,globalParameters::KDEpsilon);


		for (int j=0;j<found;j++) {
			if (minIdx[j]<i || minDist[j]<min_distancia) {
				quin = -1;
				break;
			}
		}

		resultat[i] = quin;
	}


	tmpQueryPoint = NULL;
	tmpData = NULL;

}


void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, int *assignment, int NClusters, ANNkd_tree *kdtreeOrig)
{
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/

	double *tmpQueryPoint = NULL;

	for (i = 0; i < n_regs2; i++) {
		resultat[i] = -1;
	}

	double minDist[50];
	int minIdx[50];

	/*
	* -- For all registration date2 --
	*/

	bool *checked = boolAlloc1D(NClusters,false);

	for (i = 0; i < n_regs2; i++) {

		/*if (assignment[i]==0)
		assignment[i] = assignment[i];*/


		if (checked[assignment[i]])
			continue;

		checked[assignment[i]] = true;

		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		tmpQueryPoint = data2[i];

		kdtreeOrig->annkSearch(tmpQueryPoint,1,minIdx,minDist, globalParameters::KDEpsilon);

		if (minDist[0]==min_distancia && assignment[minIdx[0]]==assignment[i]) {
			resultat[i] = i;
			//cout << i << " is matched" << endl;
			continue;
		} else if (minDist[0]<min_distancia) {
			if (assignment[minIdx[0]]!=assignment[i]) {
				//resultat[i] = -1;
				continue;
			} else {   // matched to another data point with the same assignment (in the same group)
				//cout << minIdx[0] << " is matched" << endl;
				resultat[minIdx[0]] = minIdx[0];
				continue;
			}
		} else { // minDist[0] == min_distancia is true
			quin = -1;
			int found = kdtreeOrig->annkFRSearch(tmpQueryPoint,min_distancia,50,minIdx,minDist, globalParameters::KDEpsilon);
			for (int j=0;j<found;j++) {
				if (minIdx[j] != i && minDist[j]<=min_distancia && assignment[minIdx[j]]==assignment[i]) {  // in fact, minDist[j] == min_distancia
					quin = i;
					//cout << minIdx[j] << " is matched" << endl;
					break;
				}
			}
		}

		resultat[i] = quin;
	}


	tmpQueryPoint = NULL;
	//tmpData = NULL;

	delete []checked;

}

void do_matching_combination_optimized2_kdtree(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, int *assignment, int NClusters)
{
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/

	double **tmpData = NULL;
	double *tmpQueryPoint = NULL;

	if (n_columns2!=n_vbles) {
		tmpData = doubleAlloc2D_rawFast(n_regs1,n_vbles);
		for (i = 0; i < n_regs1; i++) {
			resultat[i] = -1;

			for (int j=0;j<n_vbles;j++) {
				tmpData[i][j] = data2[i][vbles[j]];
			}
		}
		tmpQueryPoint = new double[n_vbles];
	} else {
		for (i = 0; i < n_regs1; i++) {
			resultat[i] = -1;
		}
		tmpData = data2;
	}



	ANNkd_tree *kdtree = new ANNkd_tree(tmpData,n_regs2,n_vbles);

	double minDist[50];
	int minIdx[50];

	/*
	* -- For all registration date2 --
	*/

	bool *checked = boolAlloc1D(NClusters,false);

	for (i = 0; i < n_regs1; i++) {

		if (checked[assignment[i]])
			continue;

		checked[assignment[i]] = true;

		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		if (n_columns2!=n_vbles) {
			for (int j=0;j<n_vbles;j++)
				tmpQueryPoint[j] = data1[i][vbles[j]];
		} else {
			tmpQueryPoint = data1[i];
		}


		kdtree->annkSearch(tmpQueryPoint,1,minIdx,minDist, globalParameters::KDEpsilon);


		if (minDist[0]<min_distancia) {
			resultat[i] = -1;
			continue;
		}

		int found = kdtree->annkFRSearch(tmpQueryPoint,min_distancia,50,minIdx,minDist, globalParameters::KDEpsilon);


		for (int j=0;j<found;j++) {
			if (minIdx[j]<i || minDist[j]<min_distancia) {
				quin = -1;
				break;
			}
		}

		resultat[i] = quin;
	}


	if (n_columns2!=n_vbles) {
		delete []tmpQueryPoint;
		doubleDelete2D_Fast(tmpData);
	} else {
		tmpQueryPoint = NULL;
		tmpData = NULL;
	}

	delete kdtree;
	annClose();

	delete []checked;

}


void do_matching_combination_optimized2_kdtreeRandomRecords(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles, int &count) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/

	double **tmpData = doubleAlloc2D_rawFast(n_regs1,n_vbles);
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;

		for (int j=0;j<n_vbles;j++) {
			tmpData[i][j] = data2[i][vbles[j]];
		}
	}

	ANNkd_tree *kdtree = new ANNkd_tree(tmpData,n_regs2,n_vbles);

	double minDist[50];
	int minIdx[50];

	/*
	* -- For all registration date2 --
	*/
	double *tmpQueryPoint = new double[n_vbles];

	for (i = 0; i < n_regs1; i++) {

		if ((double)rand()/RAND_MAX > globalParameters::DLDfullRandomSampleRecordsPercent)
			continue;

		count++;

		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		for (int j=0;j<n_vbles;j++)
			tmpQueryPoint[j] = data1[i][vbles[j]];

		kdtree->annkSearch(tmpQueryPoint,1,minIdx,minDist, globalParameters::KDEpsilon);

		if (minDist[0]<min_distancia) {
			resultat[i] = -1;
			continue;
		}

		int found = kdtree->annkFRSearch(tmpQueryPoint,min_distancia,50,minIdx,minDist, globalParameters::KDEpsilon);


		for (int j=0;j<found;j++) {
			if (minIdx[j]<i || minDist[j]<min_distancia) {
				quin = -1;
				break;
			}
		}


		//for (j = 0; j < n_regs2; j++) {
		//	if (j == i) {
		//		continue;
		//	}
		//	distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
		//	if ((distancia < min_distancia) || (distancia == min_distancia && j < i)) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
		//		quin = -1;
		//		break;
		//	}
		//}


		resultat[i] = quin;
		/*if (i%1000==0)
		std::cout << i << " ok.\t";*/
	}


	delete []tmpQueryPoint;
	delete kdtree;
	annClose();

	doubleDelete2D_Fast(tmpData);
}

/*----------------------------------------------------------------------------
Returns the Euclidean distance squared between data1 and register of i1
the register of i2 data2
----------------------------------------------------------------------------*/
double dist_records(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles) {
	double distancia = 0;
	int i;

	for (i = 0; i < n_vbles; i++) {
		distancia += (data1[i1][i] - data2[i2][i]) * (data1[i1][i] - data2[i2][i]);
	}
	return (distancia);
}

double dist_records(double** data1, int i1, double** data2, int i2,int ncols) {
	double distancia = 0;
	int i;

	for (i = 0; i < ncols; i++) {
		distancia += (data1[i1][i] - data2[i2][i]) * (data1[i1][i] - data2[i2][i]);
	}
	return sqrt(distancia);
}


/*
* ----------------------------------------------------------------------------
* Returns the Euclidean distance squared between data1 and register of i1
* the register of i2 data2, the combination is stores in vbles array
----------------------------------------------------------------------------
*/
double dist_records_combination(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles, int* vbles) {
	double distancia = 0;
	int i;
	double tmpValue;

	for (i = 0; i < n_vbles; i++) {
		tmpValue = (data1[i1][vbles[i]] - data2[i2][vbles[i]]);
		distancia += tmpValue * tmpValue;
	}
	return distancia;
}

double dist_records_combination_with_threshold(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles, int* vbles, double threshold) {
	double distancia = 0;
	int i;
	double tmpValue;

	for (i = 0; i < n_vbles; i++) {
		tmpValue = (data1[i1][vbles[i]] - data2[i2][vbles[i]]);
		distancia += tmpValue * tmpValue;
		if (distancia >= threshold) {
			return DOUBLE_MAX_VALUE;
		}
	}
	return distancia;
}

void calculate_dld_combination_optimizedfullRandom(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS) 
{
	int n_regs;
	int* match;
	int i, j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	match = new int[n_regs];

	int *attribIdx = new int[n_columns];

	for (int i=0;i<n_columns;i++)
		attribIdx[i] = i;

	srand(time(NULL));
	for (i = 0; i < globalParameters::DLDfullRandomSampleCount; i++) {

		for (int j=0;j<100;j++) {				// shuffle
			int a = rand()%n_columns;
			int b = rand()%n_columns;

			int tmp = attribIdx[a];
			attribIdx[a] = attribIdx[b];
			attribIdx[b] = tmp;
		}

		equal = 0;

		int activeAttribsCount = rand()%n_columns+1;
		int *curAttribs = new int[activeAttribsCount];

		for (int j=0;j<activeAttribsCount;j++) {
			curAttribs[j] = attribIdx[j];
		}

		int n_regsUtilized = 0;
		do_matching_combination_optimized2_kdtreeRandomRecords(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, activeAttribsCount, curAttribs, n_regsUtilized);

		for (j = 0; j < n_regs; j++) {
			if (match[j] == j) {
				equal++;
			}
		}

		//std::printf ("eq=%d\n\n",equal);
		percentatge = (double) equal * 100.0 / (double) n_regsUtilized ;
		dld[i] = percentatge;

		if (globalParameters::debugLevel>=10)
			std::printf ("\n\tDLD[%d] = %lf \t", i, dld[i]); ///delme
		delete []curAttribs;

	}
	delete []match;
	delete []attribIdx;
}

void calculate_dld_combination_optimizedhalfRandom(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment) 
{
	int n_regs;
	int* match;
	int i, j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	int NClusters = maxArray(assignment, n_regs)+1;

	match = new int[n_regs];

	int activeAttribsCount = (int)(n_columns/2.0+0.5);
	int *curAttribs = new int[activeAttribsCount];

	int *attribIdx = new int[n_columns];

	for (int i=0;i<n_columns;i++)
		attribIdx[i] = i;

	srand(time(NULL));
	for (i = 0; i < globalParameters::DLDhalfRandomSampleCount; i++) {

		for (int j=0;j<100;j++) {				// shuffle
			int a = rand()%n_columns;
			int b = rand()%n_columns;

			int tmp = attribIdx[a];
			attribIdx[a] = attribIdx[b];
			attribIdx[b] = tmp;
		}

		equal = 0;

		for (int j=0;j<activeAttribsCount;j++) {
			curAttribs[j] = attribIdx[j];
		}


		int NClusters = maxArray(assignment,n_regs)+1;
		do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, activeAttribsCount, curAttribs);

		for (j = 0; j < n_regs; j++) {
			if (match[j] == j) {
				equal++;
			}
		}

		//std::printf ("eq=%d\n\n",equal);
		percentatge = (double) equal * 100.0 / (double) n_regs ;
		dld[i] = percentatge;

		if (globalParameters::debugLevel>=10)
			std::printf ("\n\tDLD[%d] = %lf \t", i, dld[i]); ///delme

	}
	delete []curAttribs;
	delete []match;
	delete []attribIdx;
}

void calculate_dld_combination_optimizedAll(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment) 
{
	int n_regs;
	int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	match = new int[n_regs];

	int *curAttribs = new int[n_columns];

	for (int j=0;j<n_columns;j++) {	
		curAttribs[j] = j;
	}

	equal = 0;

	if (assignment==NULL)
		do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns, curAttribs);
	else {
		int NClusters = maxArray(assignment,n_regs)+1;
		do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns, curAttribs, assignment, NClusters);
	}

	for (j = 0; j < n_regs; j++) {
		if (match[j] == j) {
			equal++;
		}
	}

	//std::printf ("eq=%d\n\n",equal);
	percentatge = (double) equal * 100.0 / (double) n_regs;
	dld[0] = percentatge;

	if (globalParameters::debugLevel>=10)
		std::printf ("\n\tDLD[] = %lf \t", dld[0]); ///delme

	delete []curAttribs;
}


void calculate_dld_combination_optimizedAll(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment, ANNkd_tree *kdtreeOrig) 
{
	int n_regs;
	int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	match = new int[n_regs];

	int *curAttribs = new int[n_columns];

	for (int j=0;j<n_columns;j++) {	
		curAttribs[j] = j;
	}

	equal = 0;

	int NClusters = maxArray(assignment,n_regs)+1;
	do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns, curAttribs, assignment, NClusters, kdtreeOrig);

	for (j = 0; j < n_regs; j++) {
		if (match[j] == j) {
			//cout << j << "(" << assignment[j] << ") is matched" << endl;
			equal++;
		}
	}

	//std::printf ("eq=%d\n\n",equal);
	percentatge = (double) equal * 100.0 / (double) n_regs;
	dld[0] = percentatge;

	if (globalParameters::debugLevel>=10)
		std::printf ("\n\tDLD[] = %lf \t", dld[0]); ///delme

	delete []curAttribs;
}

void calculate_dld_combination_optimizedAllButOne(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment) 
{
	int n_regs;
	int* match;
	int i, j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	int NClusters = 0;

	if (assignment!=NULL)
		NClusters = maxArray(assignment, n_regs)+1;

	match = new int[n_regs];

	int *curAttribs = new int[n_columns-1];
	for (i = 0; i < n_columns; i++) {
		equal = 0;

		int idx = 0;
		for (int j=0;j<n_columns;j++) {
			if (j==i)
				continue;
			curAttribs[idx++] = j;
		}



		if (assignment==NULL)
			do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns-1, curAttribs);
		else 
			do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns-1, curAttribs, assignment, NClusters);

		for (j = 0; j < n_regs; j++) {
			if (match[j] == j) {
				equal++;
			}
		}

		//std::printf ("eq=%d\n\n",equal);
		percentatge = (double) equal * 100.0 / (double) n_regs;
		dld[i] = percentatge;

		if (globalParameters::debugLevel>=10)
			std::printf ("\n\tDLD[%d] = %lf \t", i, dld[i]); ///delme

	}
	delete []curAttribs;
}

void calculate_dld_combination_optimizedFirst7(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment) 
{
	int n_regs;
	int* match;
	int i, j, k, equal;
	double percentatge;
	int counter;

	n_regs = n_data / n_columns;

	match = new int[n_regs];

	if (assignment!=NULL) {
		counter = 0;
		int NClusters = maxArray(assignment, n_regs)+1;

		for (i = 1; i <= min(7, n_columns); i++) {
			equal = 0;
			for (k = 0; k < PS->nchoosekSaved[i]; k++, counter++) {
				switch (i) {
				case 1:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb1[k], assignment, NClusters);
					break;
				case 2:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb2[k], assignment, NClusters);
					break;
				case 3:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb3[k], assignment, NClusters);
					break;
				case 4:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb4[k], assignment, NClusters);
					break;
				case 5:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb5[k], assignment, NClusters);
					break;
				case 6:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb6[k], assignment, NClusters);
					break;
				case 7:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb7[k], assignment, NClusters);
					break;
				default:
					break;
				}

				for (j = 0; j < n_regs; j++) {
					if (match[j] == j) {
						equal++;
					}
				}
			}

			//std::printf ("eq=%d\n\n",equal);
			percentatge = (double) equal * 100.0 / (double) n_regs / PS->nchoosekSaved[i];
			dld[i - 1] = percentatge;

			if (globalParameters::debugLevel>=10)
				std::printf ("\n\tDLD[%d] = %lf \t", i, dld[i-1]); ///delme
		}
	} else {
		counter = 0;

		for (i = 1; i <= min(7, n_columns); i++) {
			equal = 0;
			for (k = 0; k < PS->nchoosekSaved[i]; k++, counter++) {
				switch (i) {
				case 1:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb1[k]);
					break;
				case 2:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb2[k]);
					break;
				case 3:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb3[k]);
					break;
				case 4:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb4[k]);
					break;
				case 5:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb5[k]);
					break;
				case 6:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb6[k]);
					break;
				case 7:
					do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb7[k]);
					break;
				default:
					break;
				}

				for (j = 0; j < n_regs; j++) {
					if (match[j] == j) {
						equal++;
					}
				}
			}

			//std::printf ("eq=%d\n\n",equal);
			percentatge = (double) equal * 100.0 / (double) n_regs / PS->nchoosekSaved[i];
			dld[i - 1] = percentatge;

			if (globalParameters::debugLevel>=10)
				std::printf ("\n\tDLD[%d] = %lf \t", i, dld[i-1]); ///delme
		}

	}

}

void calculate_dld_combination_optimizedAllCombs(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS, int *assignment)
{
	int n_regs;
	int* match;
	int i, j, k, equal;
	double percentatge;
	int counter;

	n_regs = n_data / n_columns;

	match = new int[n_regs];

	counter = 0;

	if (assignment!=NULL) {
		int NClusters = maxArray(assignment,n_regs) + 1;

		for (i = 1; i <= n_columns; i++) {

			int **combs = NULL; 
			int combsCount = 0;
			generateCombinations(n_columns, i, combs, combsCount);

			for (int idx =0;idx<combsCount; idx++) {
				equal = 0;
				do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, combs[idx], assignment, NClusters);

				for (j = 0; j < n_regs; j++) {
					if (match[j] == j) {
						equal++;
					}
				}
				dld[counter] = (double) equal * 100.0 / (double) n_regs;

				//std::printf ("\n\tDLD[%d] = %lf \t", counter, dld[counter]); 

				counter++;

				if (globalParameters::debugLevel>=10)
					if (counter%((int)((1<<n_columns)/8.0))==0)
						std::cout << (double)counter/((1<<n_columns))*100 << "% \t";

			}

			intDelete2D_Fast(combs);
		}
	} else {

		for (i = 1; i <= n_columns; i++) {

			int **combs = NULL; 
			int combsCount = 0;
			generateCombinations(n_columns, i, combs, combsCount);

			for (int idx =0;idx<combsCount; idx++) {
				equal = 0;
				do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, combs[idx]);

				for (j = 0; j < n_regs; j++) {
					if (match[j] == j) {
						equal++;
					}
				}
				dld[counter] = (double) equal * 100.0 / (double) n_regs;

				//std::printf ("\n\tDLD[%d] = %lf \t", counter, dld[counter]); 

				counter++;

				if (globalParameters::debugLevel>=10)
					if (counter%((int)((1<<n_columns)/8.0))==0)
						std::cout << (double)counter/((1<<n_columns))*100 << "% \t";

			}

			intDelete2D_Fast(combs);
		}
	}
}

void generateCombinations(int n, int k, int **&combs, int &count)
{
	count = nchoosek(n,k);
	combs = intAlloc2D_rawFast(count,k);

	int curCombIdx = 0;
	int *curComb = new int[k];
	generateCombinationsRecursive(n, k, combs, curCombIdx, 0, curComb);
	delete []curComb;
}


void generateCombinationsRecursive(int n, int k, int **combs, int &curCombIdx, int depth, int *curComb) 
{
	if (depth==k) {
		//std:cout << curCombIdx << ": ";
		for (int j=0;j<k;j++) {
			combs[curCombIdx][j] = curComb[j];
			//std::cout << combs[curCombIdx][j] << " ";
		}
		//std::cout << endl;
		curCombIdx++;
		return;
	}

	int start = (depth==0)?0:curComb[depth-1]+1;
	for (int v=start;v<n;v++) {
		curComb[depth] = v;
		generateCombinationsRecursive(n,k,combs,curCombIdx,depth+1,curComb);
	}

}


/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/
void calculate_id(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids) {
	double p;
	double percentatge;
	int i, n_regs;

	n_regs = n_data / n_columns;

	i = 0;
	p = 1.0;
	while (i < 10) {
		//percentatge=look_at_intervals(data,n_regs,n_columns,data_masked,n_regs,n_columns,match_all,p);
		percentatge = look_at_intervals2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, p);

		ids[i] = percentatge;
		//            std::printf("%d -> %f\t", i, percentatge);
		i++;
		p += 1.0;
	}
}


double dominantValue(double *data,int n)
{
	int *differentValuesCount = new int[globalParameters::confidentialValuesCount]();

	for (int i=0;i<n;i++)
		differentValuesCount[(int)data[i]]++;


	//int maxCount = differentValuesCount[0];
	//double result = 0;  // we assume the last column values (i.e., confidential values) start from 0

	//for (int i=1;i<globalParameters::confidentialValuesCount;i++) {
	//	if (differentValuesCount[i]>maxCount) 
	//	{
	//		maxCount = differentValuesCount[i];
	//		result = i;
	//	}
	//}

	int maxCount = differentValuesCount[(int)data[0]];
	double result = data[0];  // we assume the last column values (i.e., confidential values) start from 0

	for (int i=1;i<n;i++) {
		if (differentValuesCount[(int)data[i]]>maxCount) 
		{
			maxCount = differentValuesCount[(int)data[i]];
			result = data[i];
		}
	}



	return result;
}


void calculate_P_Sensitivity(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp,  ANNkd_tree *kdtreeOrig, int *assignment) 
{
	int n_regs;
	//int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	//match = new int[n_regs];

	equal = 0;

	//int NClusters = maxArray(assignment,n_regs)+1;
	//do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns, curAttribs, assignment, NClusters, kdtreeOrig);


	int *matchedIdx = new int[kp];
	double *matchedDist = new double[kp];
	double *queryPoint = new double[n_columns-1];

	double *confidentialValues = new double[kp];

	if (kdtreeOrig)
	{
		for (int i=0;i<n_regs;i++) 
		{

			for (int j=0;j<n_columns-1;j++)
				queryPoint[j] = data_masked[i][j];

			kdtreeOrig->annkSearch(queryPoint,kp,matchedIdx, matchedDist, globalParameters::KDEpsilon);
			for (int kk=0;kk<kp;kk++)
				confidentialValues[kk] = data[matchedIdx[kk]][n_columns-1];

			if (dominantValue(confidentialValues,kp)==data_masked[i][n_columns-1])
				equal++;
		}
	} 
	else 
	{
		ANNkd_tree *kdtreeMask = new ANNkd_tree(data_masked, n_regs, n_columns-1);
		for (int i=0;i<n_regs;i++) 
		{

			for (int j=0;j<n_columns-1;j++)
				queryPoint[j] = data[i][j];

			kdtreeMask->annkSearch(queryPoint,kp,matchedIdx, matchedDist, globalParameters::KDEpsilon);
			for (int kk=0;kk<kp;kk++)
				confidentialValues[kk] = data_masked[matchedIdx[kk]][n_columns-1];

			if (dominantValue(confidentialValues,kp)==data[i][n_columns-1])
				equal++;
		}
		delete kdtreeMask;
	}

	delete []confidentialValues;
	delete []queryPoint;
	delete []matchedDist;
	delete []matchedIdx;


	//std::printf ("eq=%d\n\n",equal);
	percentatge = (double) equal * 100.0 / (double) n_regs;
	for (int i=0;i<10;i++)
		ids[i] = percentatge;
}


void calculate_L_Diversity(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp,  ANNkd_tree *kdtreeOrig, int *assignment) 
{
	int n_regs;
	//int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	//match = new int[n_regs];

	equal = 0;

	//int NClusters = maxArray(assignment,n_regs)+1;
	//do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns, curAttribs, assignment, NClusters, kdtreeOrig);

	int *matchedIdx = new int[kp];
	double *matchedDist = new double[kp];
	double *queryPoint = new double[n_columns-1];

	double *confidentialValues = new double[kp];

	//double minDiversity=myInf;

	double sumDiversity = 0;

	if (kdtreeOrig)
	{
		for (int i=0;i<n_regs;i++) 
		{

			for (int j=0;j<n_columns-1;j++)
				queryPoint[j] = data_masked[i][j];

			kdtreeOrig->annkSearch(queryPoint,kp,matchedIdx, matchedDist, globalParameters::KDEpsilon);
			for (int kk=0;kk<kp;kk++)
				confidentialValues[kk] = data[matchedIdx[kk]][n_columns-1];

			//int tmp;
			//if ((tmp=diversity(confidentialValues,kp))<minDiversity || i==0)
			//	minDiversity = tmp;
			sumDiversity+=diversity(confidentialValues,kp);
		}
	} 
	else 
	{
		ANNkd_tree *kdtreeMask = new ANNkd_tree(data_masked, n_regs, n_columns-1);
		for (int i=0;i<n_regs;i++) 
		{

			for (int j=0;j<n_columns-1;j++)
				queryPoint[j] = data[i][j];

			kdtreeMask->annkSearch(queryPoint,kp,matchedIdx, matchedDist, globalParameters::KDEpsilon);
			for (int kk=0;kk<kp;kk++)
				confidentialValues[kk] = data_masked[matchedIdx[kk]][n_columns-1];

			//int tmp;
			//if ((tmp=diversity(confidentialValues,kp))<minDiversity || i==0)
			//	minDiversity = tmp;
			sumDiversity+=diversity(confidentialValues,kp);
		}
		delete kdtreeMask;
	}

	delete []confidentialValues;
	delete []queryPoint;
	delete []matchedDist;
	delete []matchedIdx;


	//std::printf ("eq=%d\n\n",equal);
	//percentatge = (1-(double)minDiversity)/(kp-1);  // between -1 to 0

	//percentatge = (double)-minDiversity;   // we want to maximize minDiversity, between -kp to -1

	percentatge =  (kp-(double)sumDiversity/n_regs)/(kp-1)*100; //  = 1+(n_regs - (double)sumDiversity)/(n_regs*(kp-1));

	for (int i=0;i<10;i++)
		ids[i] = percentatge;
}


void calculate_negL_diversity(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp,  ANNkd_tree *kdtreeOrig, int *assignment) 
{
	if (assignment==NULL)
	{
		std::cout << "Error, rev_L_Diversity cannot be calculated when no assignment is provided";
		exit(-1);
	}


	int n_regs;
	//int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	//match = new int[n_regs];

	equal = 0;

	int NClusters = maxArray(assignment,n_regs)+1;

	std::vector<int> *parts = new std::vector<int>[NClusters];
	for (int i=0;i<n_regs;i++)
		parts[assignment[i]].push_back(i);

	double sumDiversity = 0;
	for (int i=0;i<NClusters;i++)
	{
		double *confidentialValues = new double[parts[i].size()];
		for (int j=0;j<parts[i].size();j++)
			confidentialValues[j] = data[parts[i][j]][n_columns-1];

		sumDiversity+=diversity(confidentialValues,parts[i].size());
		delete []confidentialValues;
	}

	delete []parts;

	//std::printf ("eq=%d\n\n",equal);
	//percentatge = (1-(double)minDiversity)/(kp-1);  // between -1 to 0

	//percentatge = (double)-minDiversity;   // we want to maximize minDiversity, between -kp to -1

	percentatge =  -(double)sumDiversity/NClusters; // DV=1/avg(div)

	for (int i=0;i<10;i++)
		ids[i] = percentatge;
}


void calculate_negEntropy(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp,  ANNkd_tree *kdtreeOrig, int *assignment) 
{
	if (assignment==NULL)
	{
		std::cout << "Error, rev_L_Diversity cannot be calculated when no assignment is provided";
		exit(-1);
	}


	int n_regs;
	//int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	//match = new int[n_regs];

	equal = 0;

	int NClusters = maxArray(assignment,n_regs)+1;

	std::vector<int> *parts = new std::vector<int>[NClusters];
	for (int i=0;i<n_regs;i++)
		parts[assignment[i]].push_back(i);

	//double sumDiversity = 0;
	double minEntropy = myInf;
	for (int i=0;i<NClusters;i++)
	{
		double *confidentialValues = new double[parts[i].size()];
		for (int j=0;j<parts[i].size();j++)
			confidentialValues[j] = data[parts[i][j]][n_columns-1];

		//sumDiversity+=diversity(confidentialValues,parts[i].size());
		minEntropy=min(minEntropy,entropy(confidentialValues,parts[i].size()));
		delete []confidentialValues;
		if (minEntropy==0)	// cannot be lower than 0
			break;
	}

	delete []parts;

	//std::printf ("eq=%d\n\n",equal);
	//percentatge = (1-(double)minDiversity)/(kp-1);  // between -1 to 0

	//percentatge = (double)-minDiversity;   // we want to maximize minDiversity, between -kp to -1

	percentatge =  -minEntropy;

	for (int i=0;i<10;i++)
		ids[i] = percentatge;
}


void calculate_T_Closeness(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids, int kp,  ANNkd_tree *kdtreeOrig, int *assignment) 
{
	cout << "Error, not implemented.";
	exit(-1);
	return;

	int n_regs;
	//int* match;
	int j, k, equal;
	double percentatge;

	n_regs = n_data / n_columns;

	//match = new int[n_regs];

	equal = 0;

	//int NClusters = maxArray(assignment,n_regs)+1;
	//do_matching_combination_optimized2_kdtree(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, n_columns, curAttribs, assignment, NClusters, kdtreeOrig);

	int *matchedIdx = new int[kp];
	double *matchedDist = new double[kp];
	double *queryPoint = new double[n_columns-1];

	double *confidentialValues = new double[kp];

	double minDiversity=myInf;

	if (kdtreeOrig)
	{
		for (int i=0;i<n_regs;i++) 
		{

			for (int j=0;j<n_columns-1;j++)
				queryPoint[j] = data_masked[i][j];

			kdtreeOrig->annkSearch(queryPoint,kp,matchedIdx, matchedDist, globalParameters::KDEpsilon);
			for (int kk=0;kk<kp;kk++)
				confidentialValues[kk] = data[matchedIdx[kk]][n_columns-1];

			int tmp;
			if ((tmp=diversity(confidentialValues,kp))<minDiversity || i==0)
				minDiversity = tmp;
		}
	} 
	else 
	{
		ANNkd_tree *kdtreeMask = new ANNkd_tree(data_masked, n_regs, n_columns-1);
		for (int i=0;i<n_regs;i++) 
		{

			for (int j=0;j<n_columns-1;j++)
				queryPoint[j] = data[i][j];

			kdtreeMask->annkSearch(queryPoint,kp,matchedIdx, matchedDist, globalParameters::KDEpsilon);
			for (int kk=0;kk<kp;kk++)
				confidentialValues[kk] = data_masked[matchedIdx[kk]][n_columns-1];

			int tmp;
			if ((tmp=diversity(confidentialValues,kp))<minDiversity || i==0)
				minDiversity = tmp;
		}
		delete kdtreeMask;
	}

	delete []confidentialValues;
	delete []queryPoint;
	delete []matchedDist;
	delete []matchedIdx;


	//std::printf ("eq=%d\n\n",equal);
	percentatge = (1-(double)minDiversity)/(kp-1);
	for (int i=0;i<10;i++)
		ids[i] = percentatge;
}


double diversity(double *data,int n)
{
	int *differentValuesCount = new int[globalParameters::confidentialValuesCount]();

	for (int i=0;i<n;i++)
		differentValuesCount[(int)data[i]]++;

	double result = 0;  // we assume the last column values (i.e., confidential values) start from 0

	for (int i=0;i<globalParameters::confidentialValuesCount;i++) 
	{
		if (differentValuesCount[i]>0) 
		{
			result+=1;
		}
	}

	return result;
}

double entropy(double *data,int n)
{
	int *differentValuesCount = new int[globalParameters::confidentialValuesCount]();

	for (int i=0;i<n;i++)
		differentValuesCount[(int)data[i]]++;

	double result = 0;  // we assume the last column values (i.e., confidential values) start from 0

	for (int i=0;i<globalParameters::confidentialValuesCount;i++) 
	{
		if (differentValuesCount[i]>0) 
		{
			double p = (double)differentValuesCount[i]/n;
			result+=(p*(-log(p)/LOG2));
			//result+=((double)differentValuesCount[i]/n*(-log((double)differentValuesCount[i]/n)/log(2.0)));
		}
	}

	return result;
}


void calculate_sdid_relative(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids) {
	double p;
	double percentatge;
	int i, n_regs;

	n_regs = n_data / n_columns;

	double *sd = calculateStdDevN_1(data_masked,n_regs,n_columns);

	if (globalParameters::sdidRatio<=0) {
		i = 0;
		p = 1.0;
		while (i < 10) {
			//percentatge=look_at_intervals(data,n_regs,n_columns,data_masked,n_regs,n_columns,match_all,p);
			percentatge = look_at_intervals2_SDIDrelative(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, p,sd);

			ids[i] = percentatge;
			i++;
			p += 1.0;
		}
	}
	else {
		percentatge = look_at_intervals2_SDIDrelative(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, globalParameters::sdidRatio,sd);
		for (int i=0;i<10;i++)
			ids[i] = percentatge;
	}

	delete []sd;
}


/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/
void calculate_id_relative(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids) {
	double p;
	double percentatge;
	int i, n_regs;

	n_regs = n_data / n_columns;

	i = 0;
	p = 1.0;
	while (i < 10) {
		//percentatge=look_at_intervals(data,n_regs,n_columns,data_masked,n_regs,n_columns,match_all,p);
		percentatge = look_at_intervals2_relative(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, p);

		ids[i] = percentatge;
		i++;
		p += 1.0;
	}

}

double look_at_intervals2(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p) {
	int vble, i, i_ori;
	double valor_masked, valor_original, v_min, v_max;
	int pos;
	int within = 0, out = 0;


	for (vble = 0; vble < n_vbles1; vble++) 
	{
		for (i = 0; i < n_regs2; i++) {
			valor_original = data1[i][vble];
			valor_masked = data2[i][vble];
			v_min = valor_masked - p * valor_masked / 100.0;
			v_max = valor_masked + p * valor_masked / 100.0;
			//                v_min = valor_masked - p / 100.0;
			//                v_max = valor_masked + p / 100.0;

			if (valor_original >= v_min && valor_original <= v_max) {
				within++;
			} else {
				out++;
			}
		}
	}

	return (100.0 * (double) within / ((double) within + (double) out));
}

double look_at_intervals2_SDIDrelative(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p, double *SDmask) {
	int vble, i, i_ori;
	double valor_masked, valor_original, v_min, v_max;
	int pos;
	int within = 0, out = 0;

	for (vble = 0; vble < n_vbles1; vble++)  {
		for (i = 0; i < n_regs2; i++) {
			valor_original = data1[i][vble];
			valor_masked = data2[i][vble];

			////////////////////////////////////////////////

			if (abs(valor_original-valor_masked)/SDmask[vble] <= p/100.0) {   // in microaggregation, SDmask = SDorig
				within++;
			} else {
				out++;
			}

			///////////////////////////////////////////////////

		}
	}

	return (100.0 * (double) within / ((double) within + (double) out));
}


double look_at_intervals2_relative(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p) {
	int vble, i, i_ori;
	double valor_masked, valor_original, v_min, v_max;
	int pos;
	int within = 0, out = 0;

	for (vble = 0; vble < n_vbles1; vble++)  {
		for (i = 0; i < n_regs2; i++) {
			valor_original = data1[i][vble];
			valor_masked = data2[i][vble];

			////////////////////////////////////////////////

			v_min = valor_original - p * abs(valor_original / 100.0);
			v_max = valor_original + p * abs(valor_original / 100.0);
			//                v_min = valor_masked - p / 100.0;
			//                v_max = valor_masked + p / 100.0;

			if (valor_masked >= v_min && valor_masked <= v_max) {
				within++;
			} else {
				out++;
			}

			///////////////////////////////////////////////////

		}
	}

	return (100.0 * (double) within / ((double) within + (double) out));
}


double* mean(double** data, int rowsCount, int colsCount) 
{

	double* result = new double[colsCount];

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			result[j] += data[i][j];
		}
	}

	for (int j = 0; j < colsCount; j++) {
		result[j] /= rowsCount;
	}

	return result;
}

double* stDev(double** data, int rowsCount, int colsCount) 
{

	double* result = new double[colsCount];
	double* result2 = new double[colsCount];

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			result[j] += data[i][j];
			result2[j] += data[i][j] * data[i][j];
		}
	}

	for (int j = 0; j < colsCount; j++) {
		result[j] = sqrt((result2[j] - result[j] * result[j] / rowsCount) / (rowsCount - 1));
	}

	return result;
}

double **cloneArray(double **data, int nrows, int ncols)
{
	double **result = doubleAlloc2D_Fast(nrows,ncols);
	copy2D(result,data,nrows,ncols);
	return result;
}

double *cloneArray(double *data, int n)
{
	double *result = new double [n];
	copy1D(result,data,n);
	return result;
}


//int *getSortedIdxFast(double *A, int n) {
//
//	//std::vector<indexValueClass> indexValue(n, indexValueClass()); 
//	std::vector<indexValueClass> indexValue(n); 
//
//	for (int idx = 0; idx < n; idx++) {
//		indexValue[idx].val = A[idx];
//		indexValue[idx].idx = idx;
//	}
//
//	std::sort (indexValue.begin(),indexValue.end(), indexValueCompare);
//
//	int *result = intAlloc1D(n);
//
//	for (int i=0;i<n;i++) {
//		result[i] = indexValue[i].idx;
//		A[i] = indexValue[i].val;
//	}
//
//	return result;
//
//}
//
//
//


//inline double calculatedistance(double *x, double *y, int ndims) {
//	double d = 0;
//	for (int i = 0; i < ndims; i++)
//		d += pow2c(x[i] - y[i]);
//	return sqrt(d);
//}

//inline double calculateDistance2(double *x, double *y, int NDims) {
//	double d = 0;
//	for (int i = 0; i < NDims; i++)
//		d += pow2C(x[i] - y[i]);
//	return d;
//}

std::string int2string(int number)
{
	ostringstream convert;   // stream used for the conversion
	convert << number;      // insert the textual representation of 'Number' in the characters in the stream
	return convert.str(); // set 'Result' to the contents of the stream
}

std::string float2string(float number)
{
	ostringstream convert;   // stream used for the conversion
	convert << number;      // insert the textual representation of 'Number' in the characters in the stream
	return convert.str(); // set 'Result' to the contents of the stream
}

std::string double2string(double number)
{
	ostringstream convert;   // stream used for the conversion
	convert << number;      // insert the textual representation of 'Number' in the characters in the stream
	return convert.str(); // set 'Result' to the contents of the stream
}


int middle(double *a, int left, int right) {
	// "sort" three elements and take the middle one
	int i = left;
	int j = (left + right) / 2;
	int k = right;
	// order the first two
	if (a[i]>a[j]) {
		int tmp = j;
		j = i;
		i = tmp;
	}
	// bubble the third down
	if (a[j] > a[k]) {
		if (a[i]>a[k]) 
		{
			return i;
		}
		return k;
	}
	return j;
}

void swap(int i, int j, double *a, int *index) 
{
	int tmp = index[i];
	index[i] = index[j];
	index[j] = tmp;

	double tmp2=a[i];
	a[i]=a[j];
	a[j]=tmp2;            
}

int quickSortPartition(double *a, int left, int right, int *index) 
{
	int mid = middle(a, left, right);
	swap(right, mid, a, index);
	int i = left - 1;
	int j = right;

	while (true) {
		while (a[++i]<a[right]);
		while (a[right]<a[--j]) {
			if (j == left) {
				break;
			}
		}
		if (i >= j) {
			break;
		}
		swap(i, j, a, index);
	}
	swap(i, right, a, index);
	return i;
}


void quicksort(double *a, int left, int right,int *index) 
{
	if (right <= left) 
	{
		return;
	}
	int i = quickSortPartition(a, left, right,index);
	quicksort(a, left, i - 1,index);
	quicksort(a, i + 1, right,index);
}

int *quicksortArray(double *inputArray, int n) 
{
	int *index = new int[n];
	for (int i=0;i<n;i++)
		index[i] = i;
	quicksort(inputArray, 0, n - 1,index);
	return index;
}

int calculateEdgesCount(int n, int k) {
	if (n < k) {
		return 0;
	} else if (n <= 2 * k - 1) {
		return (n - k + 1) * (n - k + 2) / 2;
	} else {
		return (n - 2 * k + 1) * k + k * (k + 1) / 2;
	}
}

int *calculateCenterSizes(int *assignment, int NRecords, int NClusters)
{
	int *centerSizes = intAlloc1D(NClusters);
	for (int j=0;j<NRecords;j++)
		centerSizes[assignment[j]]++;

	return centerSizes;
}

double **calculateNewCenters(double **allData, int *assignment,  int NRecords, int NDims, int NClusters, int *centerSizes) {
	double **centers;
	centers = doubleAlloc2D_Fast(NClusters, NDims);

	for (int i = 0; i < NRecords; i++)
		for (int j = 0; j < NDims; j++)
			centers[assignment[i]][j] += allData[i][j];

	for (int i = 0; i < NClusters; i++)
		for (int j = 0; j < NDims; j++)
			centers[i][j] /= centerSizes[i];

	return centers;
}

double **calculateNewCenters(double **allData, int *assignment, int NObjects, int NDims, int NClusters) 
{
	double **centers;
	centers = doubleAlloc2D_Fast(NClusters, NDims);

	int *centerSizes = new int[NClusters]();

	for (int i = 0; i < NObjects; i++)
		centerSizes[assignment[i]]++; 

	for (int i = 0; i < NObjects; i++)
		for (int j = 0; j < NDims; j++)
			centers[assignment[i]][j] += allData[i][j];

	for (int i = 0; i < NClusters; i++)
		for (int j = 0; j < NDims; j++)
			centers[i][j] /= centerSizes[i];
	delete[] centerSizes;
	return centers;
}


double **calculateNewCentersWeighted(double **allData, int *assignment, double *weights, int NObjects, int NDims, int NClusters) 
{
	double **centers;
	centers = doubleAlloc2D_Fast(NClusters, NDims);

	double *centerSizes = new double[NClusters]();

	for (int i = 0; i < NObjects; i++)
		centerSizes[assignment[i]]+=weights[i]; 

	for (int i = 0; i < NObjects; i++)
		for (int j = 0; j < NDims; j++)
			centers[assignment[i]][j] += weights[i] * allData[i][j];

	for (int i = 0; i < NClusters; i++)
		for (int j = 0; j < NDims; j++)
			centers[i][j] /= centerSizes[i];
	delete[] centerSizes;
	return centers;
}



//double **calculateNewCenters(double **allData, int *assignment,  int NRecords, int NDims, int NClusters, int *centerSizes) {
//	double **centers;
//	centers = doubleAlloc2D_Fast(NClusters, NDims);
//
//	for (int i = 0; i < NRecords; i++)
//		for (int j = 0; j < NDims; j++)
//			centers[assignment[i]][j] += allData[i][j];
//
//	for (int i = 0; i < NClusters; i++)
//		for (int j = 0; j < NDims; j++)
//			centers[i][j] /= centerSizes[i];
//
//	return centers;
//}
//
//double **calculateNewCenters(double **allData, int *assignment, int NObjects, int NDims, int NClusters) {
//	double **centers;
//	centers = doubleAlloc2D_Fast(NClusters, NDims);
//
//	int *centerSizes = new int[NClusters]();
//
//	for (int i = 0; i < NObjects; i++)
//		centerSizes[assignment[i]]++; 
//
//	for (int i = 0; i < NObjects; i++)
//		for (int j = 0; j < NDims; j++)
//			centers[assignment[i]][j] += allData[i][j];
//
//	for (int i = 0; i < NClusters; i++)
//		for (int j = 0; j < NDims; j++)
//			centers[i][j] /= centerSizes[i];
//	delete[] centerSizes;
//	return centers;
//}



double **calculateNewCentersDisclosureAware(double **allData, int *assignment, int NObjects, int NDims, int NClusters) 
{
	double **centers;
	double **centers2;

	centers = doubleAlloc2D_Fast(NClusters, NDims);
	centers2 = doubleAlloc2D_Fast(NClusters, NDims);

	int *centerSizes = new int[NClusters]();

	//std::vector<int> *asnIndices = new std::vector<int>[NClusters];


	double **xMinus = doubleAlloc2D_Fast(NClusters,NDims, myInf);
	double **xPlus = doubleAlloc2D_Fast(NClusters,NDims, -myInf);

	for (int i = 0; i < NObjects; i++) {
		centerSizes[assignment[i]]++; 
		//asnIndices[assignment[i]].push_back(i);

	}

	for (int i = 0; i < NObjects; i++)
		for (int j = 0; j < NDims; j++) {
			centers[assignment[i]][j] += allData[i][j];
		}

		for (int i = 0; i < NClusters; i++)
			for (int j = 0; j < NDims; j++)
				centers[i][j] /= centerSizes[i];

		for (int i = 0; i < NObjects; i++) {
			for (int j = 0; j < NDims; j++) {
				if ((xMinus[assignment[i]][j] == myInf || xMinus[assignment[i]][j]<allData[i][j]) && allData[i][j]<centers[assignment[i]][j]) {
					xMinus[assignment[i]][j]=allData[i][j];
				} else if (centers[assignment[i]][j]<=allData[i][j] && (allData[i][j]<=xPlus[assignment[i]][j] || xPlus[assignment[i]][j]==-myInf)) {
					xPlus[assignment[i]][j]=allData[i][j];
				}
			}
		}


		for (int i=0;i<NClusters;i++) {
			for (int j=0;j<NDims;j++) {
				if (xMinus[i][j] == myInf) {
					centers2[i][j] = min(centers[i][j],xPlus[i][j]-abs(xPlus[i][j]*globalParameters::defaultSafetyDistance));
					continue;
				} else if (xPlus[i][j] == -myInf) {
					centers2[i][j] = max(centers[i][j],xMinus[i][j]+abs(xMinus[i][j]*globalParameters::defaultSafetyDistance));
					continue;
				} else {
					double check = xMinus[i][j] * xPlus[i][j];

					double incident;
					if (check > 0) {
						incident = 2 * check / (xMinus[i][j] + xPlus[i][j]);
					} else {
						incident = 0;
					}

					if (incident < xMinus[i][j] + globalParameters::defaultSafetyDistance * abs(xMinus[i][j])
						&& incident > xPlus[i][j] - globalParameters::defaultSafetyDistance * abs(xPlus[i][j])) {
							centers2[i][j] = incident;
					} else {
						if (centers[i][j] < xMinus[i][j] + globalParameters::defaultSafetyDistance * abs(xMinus[i][j])) {
							centers2[i][j] = xMinus[i][j] + globalParameters::defaultSafetyDistance * abs(xMinus[i][j]);
						} else if (centers[i][j] > xPlus[i][j] - globalParameters::defaultSafetyDistance * abs(xPlus[i][j])) {
							centers2[i][j] = xPlus[i][j] - globalParameters::defaultSafetyDistance * abs(xPlus[i][j]);
						} else {
							centers2[i][j] = centers[i][j];
						}
					}
				}
			}
		}



		//FILE *fp = 0;
		//for (int i=0;i<NClusters;i++)
		//{
		//	fp = fopen ("tmp.dat","w+");
		//	if (!fp)
		//	{
		//		std::cerr << "Error! Cannout write data file.";
		//		std::getchar();
		//		std::exit(-1);
		//	}

		//	fprintf (fp,"param n := %d; \nparam d := %d; \n", centerSizes[i], NDims);

		//	fprintf (fp,"param m := ");
		//	for (int k=0;k<NDims;k++) {
		//		fprintf (fp, "%d %lf ", k+1, centers[i][k]);
		//	}
		//	fprintf (fp,";\n");

		//	fprintf (fp,"param x: ");
		//	for (int k=0;k<NDims;k++) {
		//		fprintf (fp, "%d ", k+1);
		//	}

		//	fprintf (fp," := \n");

		//	for (int j=0;j<centerSizes[i];j++) {
		//		fprintf (fp, "%d ", j+1);
		//		for (int k=0;k<NDims;k++) {
		//			fprintf (fp, "%lf ", allData[asnIndices[i][j]][k]);
		//		}
		//		fprintf (fp,"\n");
		//	}
		//	fprintf (fp,";\n");

		//	std::fclose(fp);

		//	//system("cd optimization");
		//	system("ampl.exe < tmp.mod");

		//	FILE *fp2 = fopen("tmp.out","r+");
		//	for (int k=0;k<NDims;k++)
		//		fscanf(fp2,"%lf ", &centers[i][k]);

		//	fclose(fp2);

		//}

		delete[] centerSizes;
		doubleDelete2D_Fast(centers);
		doubleDelete2D_Fast(xMinus);
		doubleDelete2D_Fast(xPlus);

		return centers2;
}

double calculateSSE(double **centers, double **allData, int *assignment, int NObjects, int NDims, int NClusters) {
	double SSE = 0;
	for (int i = 0; i < NObjects; i++)
		for (int j = 0; j < NDims; j++)
			SSE += pow2C(allData[i][j] - centers[assignment[i]][j]);
	return SSE;
}

double *calculateExOnTour(double** allData, int *tour, int start, int count, int NObjects, int NDims) {
	double *mean = new double[NDims];
	for (int j = 0; j < NDims; j++) {
		double sum = 0;
		for (int i = 0; i < count; i++)
			sum += allData[tour[(start + i) % NObjects]][j];
		mean[j] = sum / count;
	}
	return mean;
}

double calculateExOnTour(double* allData, int *tour, int start, int count, int NObjects) {
	double mean;
	double sum = 0;
	for (int i = 0; i < count; i++)
		sum += allData[tour[(start + i) % NObjects]];
	mean = sum / count;
	return mean;
}

double calculateEx(double* allData, int start, int count, int NObjects) {
	double mean;
	double sum = 0;
	for (int i = 0; i < count; i++)
		sum += allData[(start + i) % NObjects];
	mean = sum / count;
	return mean;
}



double* calculateNEx2OnTour(double** allData, int *tour, int start, int count, int NObjects, int NDims) {
	double *NEx2 = new double[NDims];
	for (int j = 0; j < NDims; j++) {
		double sum = 0;
		for (int i = 0; i < count; i++)
			sum += pow2C(allData[tour[(start + i) % NObjects]][j]);
		NEx2[j] = sum;
	}
	return NEx2;
}

double calculateNEx2OnTour(double* allData, int *tour, int start, int count, int NObjects) {
	double NEx2;
	double sum = 0;
	for (int i = 0; i < count; i++)
		sum += pow2C(allData[tour[(start + i) % NObjects]]);
	NEx2 = sum;
	return NEx2;
}


double calculateNEx2(double* allData, int start, int count, int NObjects) {
	double NEx2;
	double sum = 0;
	for (int i = 0; i < count; i++)
		sum += pow2C(allData[(start + i) % NObjects]);
	NEx2 = sum;
	return NEx2;
}


double calculateSSTOnTour2(double** allData, int *tour, int start, int count, int NObjects, int NDims, const double *Ex) {
	double *NEx2 = calculateNEx2OnTour(allData, tour, start, count, NObjects, NDims); // new double[NDims];

	double SST = 0;
	for (int j = 0; j < NDims; j++) 
		SST += NEx2[j] - count * pow2(Ex[j]);

	delete[] NEx2;

	return SST;
}

double calculateSSTOnTour2(double* allData, int *tour, int start, int count, int NObjects, const double Ex) 
{
	double NEx2 = calculateNEx2OnTour(allData, tour, start, count, NObjects); // new double[NDims];
	double SST = NEx2 - count * pow2(Ex);
	return SST;
}

double calculateSST(double* allData, int start, int count, int NObjects, const double Ex) 
{
	double NEx2 = calculateNEx2(allData, start, count, NObjects); 
	double SST = NEx2 - count * pow2(Ex);
	return SST;
}



int *path2Assignment(double** allData, int* path, int k, int NRecords, int NDims) 
{
	/////////////////  generate graph    ////////////////
	int edgesCount = calculateEdgesCount(NRecords, k);

	int* rows = new int[edgesCount];
	int* cols = new int[edgesCount];
	double* vals = new double[edgesCount];
	double* g = new double[edgesCount];
	double** groupCenters = doubleAlloc2D_Fast(edgesCount,NDims);

	int counter = 0;
	double* mainMM = new double[NDims];

	double* mm = NULL;

	double* delta = new double[NDims];

	double mainSM = 0;
	double sm = 0;

	for (int i = 0; i < NRecords - k + 1; i++) {
		for (int j = i + k - 1; j < min(NRecords, i + 2 * k - 1); j++) {
			if (i == 0 && j == k - 1) {
				mm = calculateExOnTour(allData, path, i, k, NRecords, NDims);
				sm = calculateSSTOnTour2(allData, path, i, k, NRecords, NDims, mm);
				mainMM = cloneArray(mm,NDims);
				mainSM = sm;

				for (int idx = 0; idx < NDims; idx++) {
					groupCenters[counter][idx] = mm[idx];
				}

			} else if (j == i + k - 1) {
				for (int idx = 0; idx < NDims; idx++) {
					delta[idx] = allData[path[j]][idx] - allData[path[i - 1]][idx];
					mm[idx] = mainMM[idx] + (double) delta[idx] / k;
				}

				for (int idx = 0; idx < NDims; idx++) {
					groupCenters[counter][idx] = mm[idx];
				}

				sm = mainSM;
				for (int idx = 0; idx < NDims; idx++) {
					sm += (delta[idx] * (delta[idx] * (1 - k) + 2 * (allData[path[j]][idx] - mm[idx]) * k)) / k;
				}
				mainMM = cloneArray(mm,NDims);
				mainSM = sm;
			} else {
				for (int idx = 0; idx < NDims; idx++) {
					sm += (pow(mm[idx] - allData[path[j]][idx], 2)) * (j - i) / (j - i + 1);
					mm[idx] = mm[idx] + (allData[path[j]][idx] - mm[idx]) / (j - i + 1);
				}

				for (int idx = 0; idx < NDims; idx++) {
					groupCenters[counter][idx] = mm[idx];
				}
			}

			rows[counter] = i;
			cols[counter] = j + 1;

			vals[counter] = sm;
			counter++;
		}
	}

	/*for (int idx = 0; idx < edgesCount; idx++)
	{
	CGroupEdge curGroupEdge = new CGroupEdge();
	curGroupEdge.weight = vals[idx];
	curGroupEdge.center = new double[NDims];
	for (int idx2 = 0; idx2 < NDims; idx2++)
	curGroupEdge.center[idx2] = groupCenters[idx][idx2];
	g[rows[idx]*k+cols[idx]] = curGroupEdge;
	}
	*/

	for (int idx = 0; idx < edgesCount; idx++) {
		g[idx] = vals[idx];
	}

	/////////////////////////////////////////////////////
	double* costs = new double[NRecords + 1];
	int* parents = new int[NRecords + 1];

	for (int i = 0; i <= NRecords; i++) {
		costs[i] = myInf;
		parents[i] = -1;
	}

	costs[0] = 0;
	int iCounter = 0;

	for (int i = 0; i <= NRecords - k; i++) {
		for (int j = i + k; j <= min(i + 2 * k - 1, NRecords); j++, iCounter++) {
			//                    System.out.printf ("%d ", iCounter);
			//                    System.out.flush();
			if (costs[j] > costs[i] + g[iCounter]) {
				costs[j] = costs[i] + g[iCounter];
				parents[j] = i;
			}
		}
	}

	int groupsCount = 0;
	int curNodeIdx = NRecords;
	while (curNodeIdx != 0) {
		groupsCount++;
		curNodeIdx = parents[curNodeIdx];
	}

	int* shortestPath = new int[groupsCount + 1];

	curNodeIdx = NRecords;
	int groupsCountCopy = groupsCount;
	while (curNodeIdx != 0) {
		shortestPath[groupsCountCopy--] = curNodeIdx;
		curNodeIdx = parents[curNodeIdx];
	}
	shortestPath[0] = 0;

	int* assignment = new int[NRecords];

	for (int i = 0; i < groupsCount; i++) //number of nodes in the path are equal to length + 1
	{
		for (int j = shortestPath[i]; j < shortestPath[i + 1]; j++) {
			assignment[path[j]] = i;
		}
	}

	return assignment;
}

int* array2Assignment(double* mainAllData, int* path, int k, int NRecords) 
{

	int i, j;

	double** allData2d = doubleAlloc2D_rawFast(NRecords,1);

	for (i = 0, j = 0; i < NRecords; i++) {
		allData2d[i][j] = mainAllData[i];
	}

	int *assignment = path2Assignment(allData2d, path, k, NRecords, 1);
	doubleDelete2D_Fast(allData2d);
	return assignment;
}

//path dagShortestPathParallel(double* g, int NObjects, int k) 
//{
//	int idx, idx1, idx2;
//	path finalPath = new path();
//
//	double** costs = new double[2 * k - 1][NObjects + 2 * k - 1];
//	int** parents = new int[2 * k - 1][NObjects + 2 * k - 1];
//
//	for (idx = 0; idx < 2 * k - 1; idx++) {
//		costs[idx][idx] = 0;
//		for (idx2 = idx + 1; idx2 <= NObjects + idx; idx2++) {
//			costs[idx][idx2] = Double.MAX_VALUE;
//			parents[idx][idx2] = -1;
//		}
//	}
//
//	int minIdx = -1;
//	finalPath.weight = Double.MAX_VALUE;
//
//	for (idx = 0; idx < 2 * k - 1; idx++) {
//		for (idx1 = idx; idx1 < NObjects + 0 * idx; idx1++) {
//			for (idx2 = idx1 + k; idx2 <= Math.min(idx1 + 2 * k - 1, idx + NObjects); idx2++) {
//				if (g[idx1 * k + idx2 - idx1 - k] > 0 && (costs[idx][idx2] == Double.MAX_VALUE || (costs[idx][idx2] > costs[idx][idx1] + g[idx1 * k + idx2 - idx1 - k]))) {
//					costs[idx][idx2] = costs[idx][idx1] + g[idx1 * k + idx2 - idx1 - k];
//					parents[idx][idx2] = idx1;
//				}
//			}
//
//		}
//		if (finalPath.weight > costs[idx][NObjects + idx]) {
//			minIdx = idx;
//			finalPath.weight = costs[idx][NObjects + idx];
//		}
//	}
//
//	finalPath.length = 0;
//	int curNodeIdx = minIdx + NObjects;
//
//	while (curNodeIdx != minIdx) {
//		finalPath.length++;
//		//if (DETAILS) std::printf ("%d\t",curNodeIdx);
//		curNodeIdx = parents[minIdx][curNodeIdx];
//	}
//
//	finalPath.nodesInPath = new int[finalPath.length + 1];
//	curNodeIdx = minIdx + NObjects;
//
//	for (idx = finalPath.length; idx >= 0; idx--) {
//		finalPath.nodesInPath[idx] = curNodeIdx;
//		curNodeIdx = parents[minIdx][curNodeIdx];
//	}
//
//	return finalPath;
//}

double calculateGCP(double **dataPtsNorm,int *assignment,int NRecords,int NDims, int NClusters)
{
	double **minG = doubleAlloc2D_Fast(NClusters,NDims, myInf);
	double **maxG = doubleAlloc2D_Fast(NClusters,NDims, -myInf);
	double *MM = doubleAlloc1D(NDims,-myInf);
	double *mm = doubleAlloc1D(NDims, myInf);

	for (int i=0;i<NRecords;i++)
	{
		for (int j=0;j<NDims;j++)
		{
			if (dataPtsNorm[i][j] < minG[assignment[i]][j]) {
				minG[assignment[i]][j] = dataPtsNorm[i][j];
				if (dataPtsNorm[i][j] < mm[j])
					mm[j] = dataPtsNorm[i][j];
			}
			if (dataPtsNorm[i][j] > maxG[assignment[i]][j]) {
				maxG[assignment[i]][j] = dataPtsNorm[i][j];
				if (dataPtsNorm[i][j] > MM[j])
					MM[j] = dataPtsNorm[i][j];
			}
		}
	}


	double penalty = 0;
	for (int i=0;i<NRecords;i++)
	{
		for (int j=0;j<NDims;j++)
		{
			penalty+=(maxG[assignment[i]][j] - minG[assignment[i]][j])/(MM[j]-mm[j]);
		}
	}

	doubleDelete2D_Fast(minG);
	doubleDelete2D_Fast(maxG);
	delete []MM;
	delete []mm;

	return penalty/(NRecords*NDims)*100;
}


int *MDAV_GenericOnNormalizedData(double** allData, int k, int NObjects, int NDims) 
{
	int* assignment = new int[NObjects];

	int* unused = new int[NObjects];
	for (int i = 0; i < NObjects; i++) {
		unused[i] = 1;
	}

	double* distances = new double[NObjects];
	double* Ex = NULL;
	int lastClusterNumber = 0;
	int remainedObjects = NObjects;


	int* neighborsIdxList = new int[k];
	double* neighborsDistList = new double[k];

	while (remainedObjects >= 3 * k) {
		int curFarIdx = -1;
		double curFarDist = -myInf;
		if (Ex!=NULL)
			delete []Ex;

		Ex = calculateExOnActiveData(allData, unused, NObjects, NDims);

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(Ex, allData[i], NDims);
			} else {
				distances[i] = 0;
			}
			if (distances[i] > curFarDist) {
				curFarDist = distances[i];
				curFarIdx = i;
			}
		}

		int counter = 0;		//number of assigned records

		neighborsIdxList[counter] = curFarIdx;
		neighborsDistList[counter++] = 0;
		unused[curFarIdx] = 0;


		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(allData[curFarIdx], allData[i], NDims);
				if (counter < k) {
					int insertIdx = counter - 1;
					while (insertIdx >= 0 && distances[i] < neighborsDistList[insertIdx]) {
						neighborsDistList[insertIdx + 1] = neighborsDistList[insertIdx];
						neighborsIdxList[insertIdx + 1] = neighborsIdxList[insertIdx];
						insertIdx--;
					}

					neighborsIdxList[insertIdx + 1] = i;
					neighborsDistList[insertIdx + 1] = distances[i];
					counter++;

					//neighborsIdxList[counter] = i;
					//neighborsDistList[counter++] = distances[i];
				} else if (distances[i] < neighborsDistList[counter - 1]) {
					int insertIdx = k - 1;
					while (insertIdx > 0 && distances[i] < neighborsDistList[insertIdx - 1]) {
						neighborsDistList[insertIdx] = neighborsDistList[insertIdx - 1];
						neighborsIdxList[insertIdx] = neighborsIdxList[insertIdx - 1];
						insertIdx--;
					}

					neighborsIdxList[insertIdx] = i;
					neighborsDistList[insertIdx] = distances[i];
				}
			} else {
				distances[i] = myInf; // Double.MAX_VALUE;
			}
		}

		for (int idx = 0; idx < k; idx++) {
			assignment[(neighborsIdxList[idx])] = lastClusterNumber;
			unused[neighborsIdxList[idx]] = 0;

			//std::cout << neighborsIdxList[idx] << " " << assignment[(neighborsIdxList[idx])] << std::endl;
			//getchar();

		}

		lastClusterNumber++;

		//print1D(assignment,NObjects);
		//getchar();

		//////////////////////////
		double nextFarDist = 0;
		int nextFarIdx = -1;
		double tmpDist;
		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				if ((tmpDist = calculateDistance(allData[i], allData[curFarIdx], NDims)) > nextFarDist) {
					nextFarDist = tmpDist;
					nextFarIdx = i;
				}
			}
		}

		counter = 0;		//number of assigned records

		//neighborsIdxList = new int[k];
		//neighborsDistList = new double[k];

		neighborsIdxList[counter] = nextFarIdx;
		neighborsDistList[counter++] = 0;
		unused[nextFarIdx] = 0;

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(allData[nextFarIdx], allData[i], NDims);
				if (counter < k) {
					int insertIdx = counter - 1;
					while (insertIdx >= 0 && distances[i] < neighborsDistList[insertIdx]) {
						neighborsDistList[insertIdx + 1] = neighborsDistList[insertIdx];
						neighborsIdxList[insertIdx + 1] = neighborsIdxList[insertIdx];
						insertIdx--;
					}

					neighborsIdxList[insertIdx + 1] = i;
					neighborsDistList[insertIdx + 1] = distances[i];
					counter++;
				} else if (distances[i] < neighborsDistList[counter - 1]) {
					int insertIdx = k - 1;
					while (insertIdx > 0 && distances[i] < neighborsDistList[insertIdx - 1]) {
						neighborsDistList[insertIdx] = neighborsDistList[insertIdx - 1];
						neighborsIdxList[insertIdx] = neighborsIdxList[insertIdx - 1];
						insertIdx--;
					}

					neighborsIdxList[insertIdx] = i;
					neighborsDistList[insertIdx] = distances[i];
				}
			} else {
				distances[i] = myInf; // Double.MAX_VALUE;
			}
		}

		for (int idx = 0; idx < k; idx++) {
			assignment[neighborsIdxList[idx]] = lastClusterNumber;
			unused[neighborsIdxList[idx]] = 0;
		}

		lastClusterNumber++;
		remainedObjects -= 2 * k;
	}

	if (remainedObjects >= 2 * k && remainedObjects < 3 * k) {
		int curFarIdx = -1;
		double curFarDist = -myInf; // Double.NEGATIVE_INFINITY;
		Ex = calculateExOnActiveData(allData, unused, NObjects, NDims);

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(Ex, allData[i], NDims);
			} else {
				distances[i] = 0;
			}
			if (distances[i] > curFarDist) {
				curFarDist = distances[i];
				curFarIdx = i;
			}
		}

		int counter = 0;		//number of assigned records

		//int *neighborsIdxList = new int[k];
		//double *neighborsDistList = new double[k];

		neighborsIdxList[counter] = curFarIdx;
		neighborsDistList[counter++] = 0;
		unused[curFarIdx] = 0;

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(allData[curFarIdx], allData[i], NDims);
				if (counter < k) {
					int insertIdx = counter - 1;
					while (insertIdx >= 0 && distances[i] < neighborsDistList[insertIdx]) {
						neighborsDistList[insertIdx + 1] = neighborsDistList[insertIdx];
						neighborsIdxList[insertIdx + 1] = neighborsIdxList[insertIdx];
						insertIdx--;
					}

					neighborsIdxList[insertIdx + 1] = i;
					neighborsDistList[insertIdx + 1] = distances[i];
					counter++;
				} else if (distances[i] < neighborsDistList[counter - 1]) {
					int insertIdx = k - 1;
					while (insertIdx > 0 && distances[i] < neighborsDistList[insertIdx - 1]) {
						neighborsDistList[insertIdx] = neighborsDistList[insertIdx - 1];
						neighborsIdxList[insertIdx] = neighborsIdxList[insertIdx - 1];
						insertIdx--;
					}

					neighborsIdxList[insertIdx] = i;
					neighborsDistList[insertIdx] = distances[i];
				}
			} else {
				distances[i] = myInf; // Double.MAX_VALUE;
			}
		}

		for (int idx = 0; idx < k; idx++) {
			assignment[neighborsIdxList[idx]] = lastClusterNumber;
			unused[neighborsIdxList[idx]] = 0;
		}

		lastClusterNumber++;
	}

	for (int i = 0; i < NObjects; i++) {
		if (unused[i] != 0) {
			assignment[i] = lastClusterNumber;
			unused[i] = 0;
		}
	}

	/*
	for (int i=0;i<NObjects;i++)
	if (assignment[i]>=NClusters)
	assignment[i] = NClusters-1;
	*/
	lastClusterNumber++;

	//        for (int i=0;i<NObjects;i++)
	//            System.out.printf ("%d ",assignment[i]);

	delete []unused;
	delete []distances;
	delete []neighborsDistList;
	delete []neighborsIdxList;
	if (Ex!=NULL)
		delete []Ex;
	return assignment;
}

int *MDAV_GenericOnNormalizedDataKdtree(double** allData, int k, int NObjects, int NDims) 
{
	int* assignment = new int[NObjects];

	int* unused = new int[NObjects];
	for (int i = 0; i < NObjects; i++) {
		unused[i] = 1;
	}

	double* distances = new double[NObjects];
	double* Ex = NULL;
	int lastClusterNumber = 0;
	int remainedObjects = NObjects;


	int* neighborsIdxList = new int[k];
	double* neighborsDistList = new double[k];

	while (remainedObjects >= 3 * k) {
		int curFarIdx = -1;
		double curFarDist = -myInf;
		Ex = calculateExOnActiveData(allData, unused, NObjects, NDims);

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(Ex, allData[i], NDims);
			} else {
				distances[i] = 0;
			}
			if (distances[i] > curFarDist) {
				curFarDist = distances[i];
				curFarIdx = i;
			}
		}

		int counter = 0;		//number of assigned records

		neighborsIdxList[counter] = curFarIdx;
		neighborsDistList[counter++] = 0;
		unused[curFarIdx] = 0;


		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(allData[curFarIdx], allData[i], NDims);
				if (counter < k) {
					int insertIdx = counter - 1;
					while (insertIdx >= 0 && distances[i] < neighborsDistList[insertIdx]) {
						neighborsDistList[insertIdx + 1] = neighborsDistList[insertIdx];
						neighborsIdxList[insertIdx + 1] = neighborsIdxList[insertIdx];
						insertIdx--;
					}

					neighborsIdxList[insertIdx + 1] = i;
					neighborsDistList[insertIdx + 1] = distances[i];
					counter++;

					//neighborsIdxList[counter] = i;
					//neighborsDistList[counter++] = distances[i];
				} else if (distances[i] < neighborsDistList[counter - 1]) {
					int insertIdx = k - 1;
					while (insertIdx > 0 && distances[i] < neighborsDistList[insertIdx - 1]) {
						neighborsDistList[insertIdx] = neighborsDistList[insertIdx - 1];
						neighborsIdxList[insertIdx] = neighborsIdxList[insertIdx - 1];
						insertIdx--;
					}

					neighborsIdxList[insertIdx] = i;
					neighborsDistList[insertIdx] = distances[i];
				}
			} else {
				distances[i] = myInf; // Double.MAX_VALUE;
			}
		}

		for (int idx = 0; idx < k; idx++) {
			assignment[(neighborsIdxList[idx])] = lastClusterNumber;
			unused[neighborsIdxList[idx]] = 0;

			//std::cout << neighborsIdxList[idx] << " " << assignment[(neighborsIdxList[idx])] << std::endl;
			//getchar();

		}

		lastClusterNumber++;

		//print1D(assignment,NObjects);
		//getchar();

		//////////////////////////
		double nextFarDist = 0;
		int nextFarIdx = -1;
		double tmpDist;
		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				if ((tmpDist = calculateDistance(allData[i], allData[curFarIdx], NDims)) > nextFarDist) {
					nextFarDist = tmpDist;
					nextFarIdx = i;
				}
			}
		}

		counter = 0;		//number of assigned records

		//neighborsIdxList = new int[k];
		//neighborsDistList = new double[k];

		neighborsIdxList[counter] = nextFarIdx;
		neighborsDistList[counter++] = 0;
		unused[nextFarIdx] = 0;

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(allData[nextFarIdx], allData[i], NDims);
				if (counter < k) {
					int insertIdx = counter - 1;
					while (insertIdx >= 0 && distances[i] < neighborsDistList[insertIdx]) {
						neighborsDistList[insertIdx + 1] = neighborsDistList[insertIdx];
						neighborsIdxList[insertIdx + 1] = neighborsIdxList[insertIdx];
						insertIdx--;
					}

					neighborsIdxList[insertIdx + 1] = i;
					neighborsDistList[insertIdx + 1] = distances[i];
					counter++;
				} else if (distances[i] < neighborsDistList[counter - 1]) {
					int insertIdx = k - 1;
					while (insertIdx > 0 && distances[i] < neighborsDistList[insertIdx - 1]) {
						neighborsDistList[insertIdx] = neighborsDistList[insertIdx - 1];
						neighborsIdxList[insertIdx] = neighborsIdxList[insertIdx - 1];
						insertIdx--;
					}

					neighborsIdxList[insertIdx] = i;
					neighborsDistList[insertIdx] = distances[i];
				}
			} else {
				distances[i] = myInf; // Double.MAX_VALUE;
			}
		}

		for (int idx = 0; idx < k; idx++) {
			assignment[neighborsIdxList[idx]] = lastClusterNumber;
			unused[neighborsIdxList[idx]] = 0;
		}

		lastClusterNumber++;
		remainedObjects -= 2 * k;
	}

	if (remainedObjects >= 2 * k && remainedObjects < 3 * k) {
		int curFarIdx = -1;
		double curFarDist = -myInf; // Double.NEGATIVE_INFINITY;
		Ex = calculateExOnActiveData(allData, unused, NObjects, NDims);

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(Ex, allData[i], NDims);
			} else {
				distances[i] = 0;
			}
			if (distances[i] > curFarDist) {
				curFarDist = distances[i];
				curFarIdx = i;
			}
		}

		int counter = 0;		//number of assigned records

		//int *neighborsIdxList = new int[k];
		//double *neighborsDistList = new double[k];

		neighborsIdxList[counter] = curFarIdx;
		neighborsDistList[counter++] = 0;
		unused[curFarIdx] = 0;

		for (int i = 0; i < NObjects; i++) {
			if (unused[i] != 0) {
				distances[i] = calculateDistance(allData[curFarIdx], allData[i], NDims);
				if (counter < k) {
					int insertIdx = counter - 1;
					while (insertIdx >= 0 && distances[i] < neighborsDistList[insertIdx]) {
						neighborsDistList[insertIdx + 1] = neighborsDistList[insertIdx];
						neighborsIdxList[insertIdx + 1] = neighborsIdxList[insertIdx];
						insertIdx--;
					}

					neighborsIdxList[insertIdx + 1] = i;
					neighborsDistList[insertIdx + 1] = distances[i];
					counter++;
				} else if (distances[i] < neighborsDistList[counter - 1]) {
					int insertIdx = k - 1;
					while (insertIdx > 0 && distances[i] < neighborsDistList[insertIdx - 1]) {
						neighborsDistList[insertIdx] = neighborsDistList[insertIdx - 1];
						neighborsIdxList[insertIdx] = neighborsIdxList[insertIdx - 1];
						insertIdx--;
					}

					neighborsIdxList[insertIdx] = i;
					neighborsDistList[insertIdx] = distances[i];
				}
			} else {
				distances[i] = myInf; // Double.MAX_VALUE;
			}
		}

		for (int idx = 0; idx < k; idx++) {
			assignment[neighborsIdxList[idx]] = lastClusterNumber;
			unused[neighborsIdxList[idx]] = 0;
		}

		lastClusterNumber++;
	}

	for (int i = 0; i < NObjects; i++) {
		if (unused[i] != 0) {
			assignment[i] = lastClusterNumber;
			unused[i] = 0;
		}
	}

	/*
	for (int i=0;i<NObjects;i++)
	if (assignment[i]>=NClusters)
	assignment[i] = NClusters-1;
	*/
	lastClusterNumber++;

	//        for (int i=0;i<NObjects;i++)
	//            System.out.printf ("%d ",assignment[i]);

	return assignment;
}


void myShuffle(int *arr, int n, int shuffleCount)
{
	for (int j=0;j<shuffleCount;j++) 
	{
		int a = rand()%n;
		int b = rand()%n;

		int tmp = arr[a];
		arr[a] = arr[b];
		arr[b] = tmp;
	}
}

void blend(double *x, int n, int *initialAssignment, int initialNClusters, double alpha)
{
	if (alpha<1 && alpha>=0 && initialNClusters>1)
		for (int i = 0; i < n; i++) {
			if (x[i] >= alpha) {
				double mean0 = initialAssignment[i]/(initialNClusters-1.0);
				double std0 = 1.0/initialNClusters;
				unsigned long seed0 = (unsigned long)(myInf * ((alpha - x[i])/(alpha - 1))); 
				double randNormal = getNormalRandorm(mean0,std0,seed0);

				/*if(randNormal>=0) {
				x[i] = randNormal-(int)randNormal;
				} else {
				x[i] = -floor(randNormal)+randNormal;
				}*/
				x[i] = randNormal;
			}
			else {
				x[i] = initialAssignment[i]/(initialNClusters-1.0);
			}

		}
	else if (alpha==1 && initialNClusters>1)
	{
		for (int i = 0; i < n; i++) {
			x[i] = initialAssignment[i]/(initialNClusters-1.0);
		}
	} else if (alpha<1 && alpha>=0 && initialNClusters==1) 
	{
		for (int i = 0; i < n; i++) 
			x[i] = (alpha - x[i])/(alpha - 1);
	} else if (alpha==1 && initialNClusters==1) {
		double max = maxArray(x,n);
		for (int i = 0; i < n; i++) 
			x[i] = x[i]/max;
	} else {
		cerr << "Values of alpha or NClusters of initial assignment are not usable." << endl;
		exit(-1);
	}
}


int * blendAndAssign(const double *x, int n, int *initialAssignment, int initialNClusters, double alpha)
{
	return 0;  // not correctly implemented
	int *assignment = new int [n];

	if (alpha<1 && alpha>=0 && initialNClusters>1)
		for (int i = 0; i < n; i++) {
			if (x[i] >= alpha) {
				double mean0 = initialAssignment[i]/(initialNClusters-1.0);
				double std0 = 1.0/initialNClusters;
				unsigned long seed0 = (unsigned long)(myInf * ((alpha - x[i])/(alpha - 1))); 
				double randNormal = getNormalRandorm(mean0,std0,seed0);

				/*if(randNormal>=0) {
				x[i] = randNormal-(int)randNormal;
				} else {
				x[i] = -floor(randNormal)+randNormal;
				}*/

				//x[i] = randNormal;
				if (randNormal>=1) {
					assignment[i] = initialAssignment[n-1];
				} else if (randNormal<=0) {
					assignment[i] = initialAssignment[0];
				} else {
					assignment[i] = initialAssignment[(int)round(randNormal * (initialNClusters-1.0))];
				}
				//if (initialAssignment[i] != assignment[i])
				//	cout << "* " << abs(initialAssignment[i] != assignment[i]) << "\t";
			}
			else {
				//x[i] = initialAssignment[i]/(initialNClusters-1.0);
				assignment[i] = initialAssignment[i];
				//cout << "***";
			}

		}
	else if (alpha==1 && initialNClusters>1)
	{
		for (int i = 0; i < n; i++) {
			//x[i] = initialAssignment[i]/(initialNClusters-1.0);
			assignment[i] = initialAssignment[i];
		}
	} else if (alpha<1 && alpha>=0 && initialNClusters==1) 
	{
		for (int i = 0; i < n; i++) {
			//x[i] = (alpha - x[i])/(alpha - 1);
			assignment[i] = initialAssignment[(int)(x[i] * (initialNClusters))];  // check later
		}
	} else if (alpha==1 && initialNClusters==1) {
		double max = maxArray(x,n);
		for (int i = 0; i < n; i++) {
			//x[i] = x[i]/max;
			assignment[i] = initialAssignment[i];
		}
	} else {
		cerr << "Values of alpha or NClusters of initial assignment are not usable." << endl;
		exit(-1);
	}

	return assignment;
}


double getNormalRandorm(double mean0, double std0, unsigned long seed)
{
	std::normal_distribution<double> dist(mean0,std0);
	std::mt19937_64 eng;
	eng.seed(seed);

	std::normal_distribution<double>::result_type tmpNormal;
	tmpNormal = dist(eng);

	//double *arr = new double[1000];
	//std::mt19937_64 eng2;
	//for (int ii=0;ii<1000;ii++) {
	//	eng2.seed(100);
	//	arr[ii] = dist(eng2);
	//}
	//FILE *fp  =fopen("testRandom.txt", "wt");
	//print1D(fp,arr,1000);
	//delete[] arr;
	//fclose(fp);

	return tmpNormal;
}

int *MHM(double *allData, int k, int NRecords)
{
	int *sortIdx = getSortedIdxFast(allData,NRecords);

	//sort(allData,allDataCopy+NRecords);


	/////////////////  generate graph    ////////////////
	int edgesCount = calculateEdgesCount(NRecords, k);

	//int* rows = new int[edgesCount];
	//int* cols = new int[edgesCount];
	double* vals = new double[edgesCount];
	//double* g = new double[edgesCount];
	double* groupCenters = new double[edgesCount];

	int counter = 0;
	double mainMM;

	double mm;

	double delta;

	double mainSM = 0;
	double sm = 0;

	for (int i = 0; i < NRecords - k + 1; i++) {
		for (int j = i + k - 1; j < min(NRecords, i + 2 * k - 1); j++) {
			if (i == 0 && j == k - 1) {
				mm = calculateEx(allData, i, k, NRecords);
				sm = calculateSST(allData, i, k, NRecords, mm);
				mainMM = mm;
				mainSM = sm;

				groupCenters[counter] = mm;

			} else if (j == i + k - 1) {
				delta = allData[j] - allData[i - 1];
				mm = mainMM + (double) delta / k;



				groupCenters[counter] = mm;


				sm = mainSM;

				sm += (delta * (delta * (1 - k) + 2 * (allData[j] - mm) * k)) / k;

				mainMM = mm;
				mainSM = sm;
			} else {

				sm += (pow(mm - allData[j], 2)) * (j - i) / (j - i + 1);
				mm = mm + (allData[j] - mm) / (j - i + 1);



				groupCenters[counter] = mm;

			}


			vals[counter] = sm;
			counter++;
		}
	}


	/////////////////////////////////////////////////////
	double* costs = new double[NRecords + 1];
	int* parents = new int[NRecords + 1];

	for (int i = 0; i <= NRecords; i++) {
		costs[i] = myInf;
		parents[i] = -1;
	}

	costs[0] = 0;
	int iCounter = 0;

	for (int i = 0; i <= NRecords - k; i++) {
		for (int j = i + k; j <= min(i + 2 * k - 1, NRecords); j++, iCounter++) {
			//                    System.out.printf ("%d ", iCounter);
			//                    System.out.flush();
			if (costs[j] > costs[i] + vals[iCounter]) {
				costs[j] = costs[i] + vals[iCounter];
				parents[j] = i;
			}
		}
	}

	int groupsCount = 0;
	int curNodeIdx = NRecords;
	while (curNodeIdx != 0) {
		groupsCount++;
		curNodeIdx = parents[curNodeIdx];
	}

	int* shortestPath = new int[groupsCount + 1];

	curNodeIdx = NRecords;
	int groupsCountCopy = groupsCount;
	while (curNodeIdx != 0) {
		shortestPath[groupsCountCopy--] = curNodeIdx;
		curNodeIdx = parents[curNodeIdx];
	}
	shortestPath[0] = 0;

	int* assignment = new int[NRecords];

	for (int i = 0; i < groupsCount; i++) //number of nodes in the path are equal to length + 1
	{
		for (int j = shortestPath[i]; j < shortestPath[i + 1]; j++) {
			assignment[sortIdx[j]] = i;
		}
	}

	delete []vals;
	delete []groupCenters;
	delete []shortestPath;
	delete []parents;
	delete []costs;
	delete []sortIdx;

	return assignment;
}

double **aggregate(double **dataPtsNorm,int *assignment, double *orgMean, double *orgStdDev, int NRecords, int NDims,int NClusters)
{
	double **tc = NULL;

	if (globalParameters::disclosureAwareAggregation)
		tc = calculateNewCentersDisclosureAware(dataPtsNorm,assignment,NRecords,NDims,NClusters);
	else
		tc = calculateNewCenters(dataPtsNorm,assignment,NRecords,NDims,NClusters);

	double **dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims);
	//double **dataPtsMask = denormalizeMatrixWithoutNormalization(tc,assignment,orgMean,orgStdDev,NRecords,NDims);

	doubleDelete2D_Fast(tc);
	return dataPtsMask;
}

double **aggregateWeighted(double **dataPtsNorm,int *assignment, double *weights, double *orgMean, double *orgStdDev, int NRecords, int NDims,int NClusters)
{
	double **tc = NULL;

	tc = calculateNewCentersWeighted(dataPtsNorm,assignment, weights, NRecords,NDims,NClusters);

	double **dataPtsMask = 0;

	if (globalParameters::DLDType == enumDLDType::noDLD && (globalParameters::IDType == enumIDType::P_sensitivity || globalParameters::IDType == enumIDType::L_diversity || globalParameters::IDType == enumIDType::T_closeness ||  globalParameters::IDType == enumIDType::negL_diversity ||   globalParameters::IDType == enumIDType::negEntropy))
	{
		dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims, dataPtsNorm);   // this will save the value of the last column
	}
	else 
	{
		dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims);
	}
	//double **dataPtsMask = denormalizeMatrixWithoutNormalization(tc,assignment,orgMean,orgStdDev,NRecords,NDims);



	doubleDelete2D_Fast(tc);
	return dataPtsMask;
}


float *calculateExOnTourFloat(double** allData, int *tour, int start, int count, int NObjects, int NDims) {
	float *mean = new float[NDims];
	for (int j = 0; j < NDims; j++) {
		float sum = 0;
		for (int i = 0; i < count; i++)
			sum += allData[tour[(start + i) % NObjects]][j];
		mean[j] = sum / count;
	}
	return mean;
}


void injectSyntheticConfidentialValues(double **dataOrig, int NRecords, int NDims, int vCount)
{
	if (vCount<=1) {
		cout << "Error in number of synthetic confidential values count. It must be > 1.\n";
		exit(-1);
	}

	double *vOrig = new double[NRecords];
	for (int i=0;i<NRecords;i++)
		vOrig[i] = dataOrig[i][NDims-1];

	int *sIdx = getSortedIdxFast(vOrig,NRecords);

	int lastIdx = 0;

	for (int i=0;i<NRecords;i++)
	{
		dataOrig[sIdx[i]][NDims-1] = lastIdx;
		if (i%(NRecords/vCount)==((NRecords/vCount)-1) && lastIdx<vCount-1)
			lastIdx++;
	}


	delete []sIdx;
	delete []vOrig;
}
