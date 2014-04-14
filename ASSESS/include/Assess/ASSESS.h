#ifndef _ASSESS_
#define _ASSESS_

#include <masking.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>							// file I/O
#include <ANN.h>
#include <checkActiveClass2.h>
#include <checkActiveClass3.h>
#include <DisjointSets.h>
#include <boost/algorithm/string.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <time.h> 


//#include <cstdlib>						// C standard library
//#include <cstdio>							// C I/O (for sscanf)
//#include <cstring>						// string manipulation
//#include <sparseMatrix.h>
//#include <windows.h>
//#include <math.h>
//#include <time.h>
//#include <string>
//#include <vector>
//#include <array>
//#include <boost/lexical_cast.hpp>
//#include <boost/numeric/ublas/io.hpp>


using namespace std;
using namespace boost::numeric::ublas;


#ifndef remSimple
#define remSimple(a,n) (((a)<(n))?(a):(a-n))
//#define remSimple(a,n) ((a)<(n)?a:a-n)
#endif

int NNeighbors = 1; // number of nearest neighbors
long maxPts = 10000000; // maximum number of data points
int NDims = 0; // dimension
int NRecords = 0; // number of records
int verbosePeriod = 0;
bool generateOutput = false;
int fullAssess = 0;					// if full assess of information loss and disclosure risk


extern enum enumDLDType DLDType;

#ifdef __ADD_LOG__
	int totalPushCount = 0;
	int totalUselessCount = 0;
	int *updateSavingsCallCount = 0;
	double totalKdtreeSearchTime = 0;
	double totalKdtreeDeleteTime = 0;
	double kdtreeConstructionTime = 0;
	double sortTime = 0;
	double aggregationTime = 0;


	double wTmp01;
	double wTmp02;
	double wTmp03;
	
#endif

	
typedef struct {
	double val;
	unsigned int index;
} t_savings;

typedef struct {
	double val;
	unsigned int i;
	unsigned int j;
} t_savings2;


int *nnIdxAll = 0;
double *distsAll = 0;

string datasetName = "";
string datasetPartName = "";
string datasetMaskName = "";
string datasetAsnName = "";

double wall1, wall2;
double wall11, wall22;


istream* dataIn = NULL; // input for data points
istream* dataInAsn = NULL; // input for data points
istream* dataInPart = NULL; // input for data points
istream* dataInMask = NULL; // input for data points


typedef struct {
	double saveVal;
	int start;
	int end;
} pqElementType;

struct compareSaving
{
	bool operator() (const pqElementType *lhs, const pqElementType *rhs) const
	{
		return lhs->saveVal < rhs->saveVal;
	}
};

struct compareHeap
{
	bool operator() (const  boost::heap::fibonacci_heap<pqElementType*, boost::heap::compare<compareSaving> > *lhs, const boost::heap::fibonacci_heap<pqElementType*, boost::heap::compare<compareSaving> > *rhs) const
	{
		return lhs->top()->saveVal  < rhs->top()->saveVal;
	}
};

typedef boost::heap::fibonacci_heap<pqElementType*, boost::heap::compare<compareSaving> > innerHeapType;
typedef boost::heap::fibonacci_heap<innerHeapType *, boost::heap::compare<compareHeap> >  outerHeapType;

#endif