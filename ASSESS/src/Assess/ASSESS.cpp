/////////////////////////////////////////////////
/*
* File:   ASSESS.cpp
* Author: Reza Mortazavi
*
* Created on July 23, 2013, 12:20 AM
*/
////////////////////////////////////////////////

//#include <vld.h>    visual leak detector
#include <ASSESS.h>

using namespace std;
using namespace boost::numeric::ublas;

//int **calculateCenterSizes_RangeK(int **assignment, int *NClusters, int kMin, int kMax)
//{
//	int **centerSizes= new int*[kMax-kMin+1];
//	for (int curK=kMin;curK<=kMax;curK++) {
//		centerSizes[curK-kMin] = intAlloc1D(NClusters[curK-kMin]);
//		for (int j=0;j<NRecords;j++)
//			centerSizes[curK-kMin][assignment[curK-kMin][j]]++;
//	}
//
//	return centerSizes;
//}
//
//int *calculateCenterSizes(int *assignment, int NClusters)
//{
//	int *centerSizes = intAlloc1D(NClusters);
//	for (int j=0;j<NRecords;j++)
//		centerSizes[assignment[j]]++;
//
//	return centerSizes;
//}
//
//double **calculateNewCenters(double **allData, int *assignment,  int NClusters, int *centerSizes) {
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
//
//double **calculateNewCentersDisclosureAware(double **allData, int *assignment, int NObjects, int NDims, int NClusters) 
//{
//	double **centers;
//	centers = doubleAlloc2D_Fast(NClusters, NDims);
//
//	int *centerSizes = new int[NClusters]();
//
//	std::vector<int> *asnIndices = new std::vector<int>[NClusters];
//
//	for (int i = 0; i < NObjects; i++) {
//		centerSizes[assignment[i]]++; 
//		asnIndices[assignment[i]].push_back(i);
//	}
//
//	for (int i = 0; i < NObjects; i++)
//		for (int j = 0; j < NDims; j++)
//			centers[assignment[i]][j] += allData[i][j];
//
//	for (int i = 0; i < NClusters; i++)
//		for (int j = 0; j < NDims; j++)
//			centers[i][j] /= centerSizes[i];
//
//	FILE *fp = 0;
//	for (int i=0;i<NClusters;i++)
//	{
//		fp = fopen ("tmp.dat","w+");
//		if (!fp)
//		{
//			std::cerr << "Error! Cannout write data file.";
//			std::getchar();
//			std::exit(-1);
//		}
//
//		fprintf (fp,"param n := %d; \nparam d := %d; \n", centerSizes[i], NDims);
//
//		fprintf (fp,"param m := ");
//		for (int k=0;k<NDims;k++) {
//			fprintf (fp, "%d %lf ", k+1, centers[i][k]);
//		}
//		fprintf (fp,";\n");
//
//		fprintf (fp,"param x: ");
//		for (int k=0;k<NDims;k++) {
//			fprintf (fp, "%d ", k+1);
//		}
//
//		fprintf (fp," := \n");
//
//		for (int j=0;j<centerSizes[i];j++) {
//			fprintf (fp, "%d ", j+1);
//			for (int k=0;k<NDims;k++) {
//				fprintf (fp, "%lf ", allData[asnIndices[i][j]][k]);
//			}
//			fprintf (fp,"\n");
//		}
//		fprintf (fp,";\n");
//
//		std::fclose(fp);
//
//		//system("cd optimization");
//		system("ampl.exe < tmp.mod");
//
//		FILE *fp2 = fopen("tmp.out","r+");
//		for (int k=0;k<NDims;k++)
//			fscanf(fp2,"%lf ", &centers[i][k]);
//
//		fclose(fp2);
//
//	}
//
//	delete[] centerSizes;
//	return centers;
//}


//double calculateSSE(double **centers, double **allData, int *assignment, int NObjects, int NDims, int NClusters) {
//	double SSE = 0;
//	for (int i = 0; i < NObjects; i++)
//		for (int j = 0; j < NDims; j++)
//			SSE += pow2C(allData[i][j] - centers[assignment[i]][j]);
//	return SSE;
//}

//double *calculateExOnTour(double** allData, int *tour, int start, int count, int NObjects, int NDims) {
//	double *mean = new double[NDims];
//	for (int j = 0; j < NDims; j++) {
//		double sum = 0;
//		for (int i = 0; i < count; i++)
//			sum += allData[tour[(start + i) % NObjects]][j];
//		mean[j] = sum / count;
//	}
//	return mean;
//}
//
//double* calculateNEx2OnTour(double** allData, int *tour, int start, int count, int NObjects, int NDims) {
//	double *NEx2 = new double[NDims];
//	for (int j = 0; j < NDims; j++) {
//		double sum = 0;
//		for (int i = 0; i < count; i++)
//			sum += pow2C(allData[tour[(start + i) % NObjects]][j]);
//		NEx2[j] = sum;
//	}
//	return NEx2;
//}
//
//double calculateSSTOnTour2(double** allData, int *tour, int start, int count, int NObjects, int NDims, const double *Ex) {
//	double *NEx2 = calculateNEx2OnTour(allData, tour, start, count, NObjects, NDims); // new double[NDims];
//
//	double SST = 0;
//	for (int j = 0; j < NDims; j++) 
//		SST += NEx2[j] - count * pow2(Ex[j]);
//
//	delete[] NEx2;
//
//	return SST;
//}

//bool loadDataFast(double **&allData,std::string datasetName)
//{
//	///////////////////////////
//	FILE * fp = fopen(&datasetName[0],"r");
//	allData = doubleAlloc2D_rawFast(NRecords,NDims);
//	int curRecordIdx = 0;
//	int curAttribIdx = 0;
//
//	if (fp==NULL)
//	{
//		std::cout << "Cannot open " << datasetName << " file";
//		return false;
//	}
//
//	char buffer[2000];
//	char temp;
//	int i;
//
//	int n_data=0;
//	int n_columns=0;
//
//	i=0;
//	buffer[i]=fgetc(fp);
//	while(!feof(fp))// && buffer[i]!='\n' && buffer[i]!='\r') /* while not end of line */
//	{
//		/*-- skip initial white --*/
//		i=0;
//		while(buffer[i]==' ' || buffer[i]=='\t' || buffer[i]==',' || buffer[i]=='\r' || buffer[i]=='\n')
//			buffer[i]=fgetc(fp);
//		/*-- read the number and array of characters --*/
//		while(buffer[i]!=' ' && buffer[i]!='\t' && buffer[i]!=',' && !feof(fp) && buffer[i]!='\n' && buffer[i]!='\r')
//		{
//			i++;
//			buffer[i]=fgetc(fp);
//		}
//		temp=buffer[i];
//		buffer[i]='\0';
//		/*-- Converting to integer --*/
//		if(i>0)
//		{
//			allData[curRecordIdx][curAttribIdx++]=atof(buffer);
//			if (curAttribIdx==NDims) {
//				curAttribIdx=0;
//				curRecordIdx++;
//			}
//			n_columns++;
//			n_data++;
//		}
//
//		//// I added these lines
//		//while(!feof(fp)  && (temp==' ' || temp=='\t' || temp==',' || temp=='\r' || temp=='\n'))
//		//	temp=fgetc(fp);
//
//		buffer[0]=temp;
//	}
//
//
//	/* --- Read the rest of values ​​file --- */
//	while(!feof(fp))
//	{
//		fscanf(fp,"%[^\r\n\t, ]%*c",buffer);
//		if (curRecordIdx<NRecords) {
//			allData[curRecordIdx][curAttribIdx++]=atof(buffer);
//			if (curAttribIdx==NDims) {
//				curAttribIdx=0;
//				curRecordIdx++;
//			}
//		}
//		n_data++;
//		// I added these lines
//		//do {
//		//	temp = fgetc(fp);
//		//} while (temp=='\r' || temp=='\n' || temp=='\t' || temp==' ' || temp==',');
//		//ungetc (temp,fp);
//	}
//
//
//
//	std::fclose(fp);
//	fp = 0;
//
//	if (!(n_data==NRecords*NDims))
//		return false;
//	return true;
//}


bool loadAssignmentFast(int *&allData,std::string datasetName)
{
	///////////////////////////
	FILE * fp = fopen(&datasetName[0],"r");
	allData = intAlloc1D(NRecords);
	int curRecordIdx = 0;

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
	while(!feof(fp) && buffer[i]>=0) // && buffer[i]!='\n' && buffer[i]!='\r') /* while not end of line */
	{
		//if (n_data==149640)
		//	n_data = n_data;

		/*-- skip initial white --*/
		i=0;
		while(buffer[i]==' ' || buffer[i]=='\t' || buffer[i]==',') {// || buffer[i]=='\r' || buffer[i]=='\n')
			buffer[i]=fgetc(fp);
			if (buffer[i]<0 || buffer[i]==EOF) 
				break;
		}
		/*-- read the number and array of characters --*/
		while(buffer[i]!=' ' && buffer[i]!='\t' && buffer[i]!=',' && !feof(fp) && buffer[i]!='\n' && buffer[i]!='\r' && buffer[i]>=0 && buffer[i]!=EOF)
		{
			i++;
			buffer[i]=fgetc(fp);
		}
		temp=buffer[i];
		buffer[i]='\0';
		/*-- Converting to integer --*/
		if(i>0)
		{
			allData[curRecordIdx]=atoi(buffer);
			curRecordIdx++;
			n_columns++;
			n_data++;
		}
		else
			break;
		//// I added these lines
		//while(!feof(fp)  && (temp==' ' || temp=='\t' || temp==',' || temp=='\r' || temp=='\n'))
		//	temp=fgetc(fp);
		i=0;
		buffer[0]=temp;
	}


	/* --- Read the rest of values ​​file --- */
	while(!feof(fp) && fscanf(fp,"%[^\r\n\t, ]%*c",buffer) == 1)
	{
		if (curRecordIdx<NRecords) {
			allData[curRecordIdx]=atoi(buffer);
			curRecordIdx++;
		}
		n_data++;
		// I added these lines
		//do {
		//	temp = fgetc(fp);
		//} while (temp=='\r' || temp=='\n' || temp=='\t' || temp==' ' || temp==',');
		//ungetc (temp,fp);
	}


	//if((n_data)%(n_columns) != 0 || (n_columns)==1)
	//	n_data--;


	int curNRecords = n_data;


	std::fclose(fp);
	fp = 0;

#ifdef DEBUG 
	std::cout << curNRecords << " records with " << curNDims << " attributes are read.\n";
#endif

	if (curNRecords!=NRecords)
		return false;
	return true;
}

bool loadPartFast(int *&allData,std::string datasetName)
{
	///////////////////////////
	FILE * fp = fopen(&datasetName[0],"r");
	allData = intAlloc1D(NRecords);
	int curNRecords = 0;

	if (fp==NULL)
	{
		std::cout << "Cannot open " << datasetName << " file";
		return false;
	}


	int clusterCount;

	fscanf (fp, "%d", &clusterCount);

	int dataPointsCount;
	int rcrd;
	for (int counter=0;counter<clusterCount;counter++) {
		fscanf (fp, "%d", &dataPointsCount);
		for (int j=0;j<dataPointsCount;j++) {
			fscanf (fp, "%d", &rcrd);
			allData[rcrd] = counter;
			curNRecords++;
		}
	}

	fclose(fp);
	fp = 0;

	if (curNRecords!=NRecords)
		return false;
	return true;
}

void getArgs(int argc, char **argv) 
{
	static ifstream dataStream; // data file stream
	static ifstream dataStreamMask; // data file stream
	static ifstream dataStreamAsn; // data file stream
	static ifstream dataStreamPart; // data file stream

	globalParameters::DLDType = enumDLDType::all;  // default
	globalParameters::IDType = enumIDType::relativeOrig; // default

	if (argc < 9) { // no arguments
		cerr << "Usage:\n\n"
			<< "ASSESS -f dataFile -n NRecords -d dim [-p partFile] [-a asnFile] [-m maskFile] [-tdld DLDtype=4 (all)] [-tid IDtype=3 (relativeOrig)] [-kp kpValue=10] [-v DebugLevel=8]\n\n"
			<< "where:\n"
			<< " dataFile    Name of file containing data points\n"
			<< " partFile    Name of file containing partition names\n"
			<< " asnFile     Name of file containing assignments\n"
			<< " maskFile    Name of file containing masked data\n"
			<< " NRecords    Number of records\n"
			<< " dim         Dimension of the space\n"
			<< " kpValue     the value of kp used in kp-anonymity P-sensitivity\n"
			<< " v           Debug Level\n"
			<< " DLDtype     Type of distance based disclosure risk computation\n" 
			<< "               available types are (0=noDLD, 1=allButOne, 2=halfRandom, 3=fullRandom, 4=all, 5=allCombsMean, 6=allCombsMax, 7=first7)\n"
			<< " IDtype      Type of interval disclosure computation\n" 
			<< "               available types are (0=noID, 1=std01, 2=standardDeviation, 3=relativeOrig, 4=range, 5=P-Sensitivity (DLDType must be set to noDLD))\n"
			<< "\nResults are sent to the standard output.\n";
		getchar();
		exit(0);
	}
	int i = 1;
	while (i < argc) { // read arguments
		//std::cout << argv[i] << ' ' << argv[i+1] << endl;
		if (!strcmp(argv[i], "-d")) { // -d option
			NDims = atoi(argv[++i]); // get dimension to dump
		} else if (!strcmp(argv[i], "-n")) { // -n option
			NRecords = atoi(argv[++i]); // get NRecords to dump
		} else if (!strcmp(argv[i], "-v")) { // -v option
			globalParameters::debugLevel = atoi(argv[++i]); // get Debug Level
		} else if (!strcmp(argv[i], "-tdld")) { // -tdld option
			switch (atoi(argv[++i]))
			{
			case 0:
				globalParameters::DLDType = enumDLDType::noDLD;
				break;
			case 1:
				globalParameters::DLDType = enumDLDType::allButOne;
				break;
			case 2:
				globalParameters::DLDType = enumDLDType::halfRandom;
				break;
			case 3:
				globalParameters::DLDType = enumDLDType::fullRandom;
				break;
			case 4:
				globalParameters::DLDType = enumDLDType::all;
				break;
			case 5:
				globalParameters::DLDType = enumDLDType::allCombsMean;
				break;
			case 6:
				globalParameters::DLDType = enumDLDType::allCombsMax;
				break;
			case 7:
				globalParameters::DLDType = enumDLDType::first7;
				break;
			default:
				std::cout << "DLDType is set to all." << endl;
				globalParameters::DLDType = enumDLDType::all;
				break;
			}
		} else if (!strcmp(argv[i], "-tid")) { // -tid option
			switch (atoi(argv[++i]))
			{
			case 0:
				globalParameters::IDType = enumIDType::noID;
				break;
			case 1:
				globalParameters::IDType = enumIDType::std01;
				break;
			case 2:
				globalParameters::IDType = enumIDType::standardDeviation;
				break;
			case 3:
				globalParameters::IDType = enumIDType::relativeOrig;
				break;
			case 4:
				globalParameters::IDType = enumIDType::range;
				break;
			case 5:
				globalParameters::IDType = enumIDType::P_sensitivity;
				break;
			case 6:
				globalParameters::IDType = enumIDType::L_diversity;
				break;
			case 7:
				globalParameters::IDType = enumIDType::T_closeness;
				break;
			case 8:
				globalParameters::IDType = enumIDType::negL_diversity;
				break;
			case 9:
				globalParameters::IDType = enumIDType::negEntropy;
				break;
			default:
				std::cout << "ID is set to relativeOrig." << endl;
				globalParameters::IDType = enumIDType::relativeOrig;
				break;
			}
		} else if (!strcmp(argv[i], "-m")) { // -lowMem option
			datasetMaskName = argv[++i];
			dataStreamMask.open(argv[i], ios::in); // open data file
			if (!dataStreamMask) {
				cerr << "Cannot open " << datasetMaskName << " file" << endl;
				cerr << "You are here: " << argv[0] << endl;

				getchar();
				exit(1);
			}
			dataInMask = &dataStreamMask; // make this the data stream
		} else if (!strcmp(argv[i], "-a")) { // -fullAssess option
			datasetAsnName = argv[++i];
			dataStreamAsn.open(argv[i], ios::in); // open data file
			if (!dataStreamAsn) {
				cerr << "Cannot open " << datasetAsnName << " file" << endl;
				cerr << "You are here: " << argv[0] << endl;

				getchar();
				exit(1);
			}
			dataInAsn = &dataStreamAsn; // make this the data stream
		} else if (!strcmp(argv[i], "-p")) { // -fullAssess option
			datasetPartName = argv[++i];
			dataStreamPart.open(argv[i], ios::in); // open data file
			if (!dataStreamPart) {
				cerr << "Cannot open " << datasetPartName << " file" << endl;
				cerr << "You are here: " << argv[0] << endl;

				getchar();
				exit(1);
			}
			dataInPart = &dataStreamPart; // make this the data stream
		} else if (!strcmp(argv[i], "-f")) { // -f option
			datasetName = argv[++i];
			dataStream.open(argv[i], ios::in); // open data file
			if (!dataStream) {
				cerr << "Cannot open " << datasetName << " file" << endl;
				cerr << "You are here: " << argv[0] << endl;

				getchar();
				exit(1);
			}
			dataIn = &dataStream; // make this the data stream
		} else { // illegal syntax
			cerr << "Unrecognized option.\n";
			cerr << argv[i] << endl;
			getchar();
			exit(1);
		}
		i++;
	}

	if (dataIn == NULL) {
		cerr << "Cannot open data file.\n";
		exit(1);
	}

	if ((globalParameters::IDType==enumIDType::P_sensitivity || globalParameters::IDType==enumIDType::L_diversity || globalParameters::IDType==enumIDType::T_closeness || globalParameters::IDType==enumIDType::negL_diversity || globalParameters::IDType==enumIDType::negEntropy) && globalParameters::DLDType!=enumDLDType::noDLD)
	{
		std::cout << "Warning: dldType is not compatible with idType. The dldType is set to noDLD.\n";
		globalParameters::DLDType = enumDLDType::noDLD;
	}


	//switch

}

int main(int argc, char** argv) 
{
	double **dataPtsOrig=0;	// data points
	double **dataPtsNorm=0;	// data points

	double **dataPtsMask=0;	// data points
	double **dataPtsMaskNorm=0;	// data points


	/////////////////////////////////////////////
	getArgs(argc, argv);		// read command-line arguments


	if (globalParameters::debugLevel>=10)
		std::cout << "Reading Data Points...";

	wall1 = get_wall_time();

	////////////////////////////////////////////
	if (!loadDataFast(dataPtsOrig,datasetName, NRecords, NDims)) {
		std::cout << "Error loading data points. Please try again.\n";
		getchar();
		exit(-1);
	}

	wall2 = get_wall_time();

	if (globalParameters::debugLevel>=10)
		std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

	if (globalParameters::DLDType==enumDLDType::noDLD && (globalParameters::IDType==enumIDType::P_sensitivity || globalParameters::IDType==enumIDType::L_diversity || globalParameters::IDType==enumIDType::T_closeness || globalParameters::IDType==enumIDType::negL_diversity || globalParameters::IDType==enumIDType::negEntropy))
	{
		injectSyntheticConfidentialValues(dataPtsOrig,NRecords,NDims,globalParameters::confidentialValuesCount);
	}



	////////////////////////////////////////////
	if (globalParameters::debugLevel>=10)
		std::cout << "Normalizing original file...";
	wall1 = get_wall_time();

	double *orgMean = 0, *orgStdDev = 0;

	dataPtsNorm = normalizeMatrix(dataPtsOrig, orgMean, orgStdDev, NRecords, NDims);

	wall2 = get_wall_time();
	if (globalParameters::debugLevel>=10)
		std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

	////////////////////////////////////////////

	double *mskMean = 0, *mskStdDev = 0;

	int *assignment = NULL;

	if (dataInMask!=NULL) {

		if (globalParameters::debugLevel>=10)
			std::cout << "Reading Data Points...";

		wall1 = get_wall_time();

		////////////////////////////////////////////
		if (!loadDataFast(dataPtsMask,datasetMaskName, NRecords, NDims)) {
			std::cout << "Error loading data points. Please try again.\n";
			getchar();
			exit(-1);
		}

		wall2 = get_wall_time();

		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

		if (globalParameters::debugLevel>=10)
			std::cout << "Normalizing masked file...";
		wall1 = get_wall_time();

		dataPtsMaskNorm = normalizeMatrix(dataPtsMask, mskMean, mskStdDev, NRecords, NDims);

		wall2 = get_wall_time();
		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;
	} else if (dataInAsn!=NULL) {

		if (globalParameters::debugLevel>=10)
			std::cout << "Reading Assignment Points...";

		wall1 = get_wall_time();


		////////////////////////////////////////////
		if (!loadAssignmentFast(assignment,datasetAsnName)) {
			std::cout << "Error loading assignments. Please try again.\n";
			getchar();
			exit(-1);
		}

		wall2 = get_wall_time();

		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

		if (globalParameters::debugLevel>=10)
			std::cout << "Generating masked Points...";

		wall1 = get_wall_time();

		double **tc = calculateNewCenters(dataPtsNorm,assignment,NRecords,NDims,maxArray(assignment,NRecords)+1);

		//double **tmpRawMask = doubleAlloc2D_rawFast(NRecords,NDims);
		//for (int i=0;i<NRecords;i++) {
		//	for (int j=0;j<NDims;j++) {
		//		tmpRawMask[i][j] = tc[assignment[i]][j];
		//	}
		//}

		if (globalParameters::DLDType==enumDLDType::noDLD && (globalParameters::IDType==enumIDType::P_sensitivity || globalParameters::IDType==enumIDType::L_diversity || globalParameters::IDType==enumIDType::T_closeness || globalParameters::IDType==enumIDType::negL_diversity || globalParameters::IDType==enumIDType::negEntropy)) 
		{
			dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims, dataPtsNorm);   // this will save the value of the last column
		}
		else 
		{
			dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims);
		}
		doubleDelete2D_Fast(tc);


		wall2 = get_wall_time();

		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;


		if (globalParameters::debugLevel>=10)
			std::cout << "Normalizing masked file...";
		wall1 = get_wall_time();

		dataPtsMaskNorm = normalizeMatrix(dataPtsMask, mskMean, mskStdDev, NRecords, NDims);

		wall2 = get_wall_time();
		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

	} else if (dataInPart!=NULL) {

		if (globalParameters::debugLevel>=10)
			std::cout << "Reading Partitioning Data...";

		wall1 = get_wall_time();


		////////////////////////////////////////////
		if (!loadPartFast(assignment,datasetPartName)) {
			std::cout << "Error loading partitioning data. Please try again.\n";
			getchar();
			exit(-1);
		}

		wall2 = get_wall_time();

		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

		if (globalParameters::debugLevel>=10)
			std::cout << "Generating masked Points...";

		wall1 = get_wall_time();

		double **tc = calculateNewCenters(dataPtsNorm,assignment,NRecords,NDims,maxArray(assignment,NRecords)+1);

		if (globalParameters::DLDType==enumDLDType::noDLD && (globalParameters::IDType==enumIDType::P_sensitivity || globalParameters::IDType==enumIDType::L_diversity || globalParameters::IDType==enumIDType::T_closeness || globalParameters::IDType==enumIDType::negL_diversity || globalParameters::IDType==enumIDType::negEntropy)) 
		{
			dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims, dataPtsNorm);   // this will save the value of the last column
		}
		else 
		{
			dataPtsMask = denormalizeMatrix(tc,assignment,orgMean,orgStdDev,NRecords,NDims);
		}

		doubleDelete2D_Fast(tc);


		wall2 = get_wall_time();

		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;


		if (globalParameters::debugLevel>=10)
			std::cout << "Normalizing masked file...";
		wall1 = get_wall_time();

		dataPtsMaskNorm = normalizeMatrix(dataPtsMask, mskMean, mskStdDev, NRecords, NDims);

		wall2 = get_wall_time();
		if (globalParameters::debugLevel>=10)
			std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;
	}



	//double SST = calculateSST(dataPtsNorm, NRecords, NDims);
	if (globalParameters::debugLevel>=10)
		std::cout << "Dataset: " << datasetName << " (" << NRecords << "x" << NDims << ")" << endl;

	/////////////////////////////////////////

	if (globalParameters::debugLevel>=10) {
		std::cout << "DLD Type: " << globalParameters::DLDTypeNames[(int)globalParameters::DLDType] << endl;
		std::cout << " ID Type: " << globalParameters::IDTypeNames[(int)globalParameters::IDType] << endl;
	}


	if (globalParameters::debugLevel>=10)
		std::cout << "Assessing anonymized dataset...";

	double DR=0,IL=0;

	wall1 = get_wall_time();

	if (assignment==NULL) {
		calculateDRIL (dataPtsOrig,dataPtsMask, DR, IL, NRecords, NDims);
	} else {
		calculateDRIL (dataPtsOrig,dataPtsMask, DR, IL, NRecords, NDims, assignment);
		delete []assignment;
	}

	wall2 = get_wall_time();
	if (globalParameters::debugLevel>=10)
		std::cout << "Done.  (" << (wall2 - wall1) << " seconds)" << endl ;

	std::string processedFilename="";
	if (datasetAsnName!="") {
		size_t found = datasetAsnName.find_last_of("/\\");
		processedFilename = datasetAsnName.substr(found+1) ;
	} else if (datasetPartName!="") {
		size_t found = datasetPartName.find_last_of("/\\");
		processedFilename = datasetPartName.substr(found+1) ;
	} else if (datasetMaskName!="") {
		size_t found = datasetMaskName.find_last_of("/\\");
		processedFilename = datasetMaskName.substr(found+1) ;
	} else {
		processedFilename = "Unknown";
	}


	if (globalParameters::debugLevel>=10)
		std::cout << "*************************************************************************" << endl;

	if (globalParameters::debugLevel>=8)
		if (globalParameters::DLDType==enumDLDType::allButOne ||
			globalParameters::DLDType==enumDLDType::allCombsMax ||
			globalParameters::DLDType==enumDLDType::allCombsMean ||
			globalParameters::DLDType==enumDLDType::first7 ||
			globalParameters::DLDType==enumDLDType::fullRandom ||
			globalParameters::DLDType==enumDLDType::halfRandom)
		{
			//std::cout << processedFilename << "\tDisclosure Risk (+/-stdDev) = " << DR << " (" << globalParameters::stdDevDLD << ")" << ",\tInformation Loss = " << IL << ",\t(DR+IL)/2 = " << (DR+IL)/2.0 << endl;
			std::cout << processedFilename << "\tDisclosure Risk (" << globalParameters::DLDTypeNames[globalParameters::DLDType] << ", " << globalParameters::IDTypeNames[globalParameters::IDType]  << ") , Information Loss (PIL) , SI=(DR+IL)/2 = \t" << DR << "±" << globalParameters::stdDevDLD << " , \t" << IL << ", \t" << (DR+IL)/2.0 << endl;
		}
		else
		{
			std::cout << processedFilename << "\tDisclosure Risk (" << globalParameters::DLDTypeNames[globalParameters::DLDType] << ", " << globalParameters::IDTypeNames[globalParameters::IDType]  << ") , Information Loss (PIL) , SI=(DR+IL)/2 = \t" << DR << " , \t" << IL << ", \t" << (DR+IL)/2.0 << endl;
		}
	if (globalParameters::debugLevel>=10)
		std::cout << "*************************************************************************" << endl;

	if (globalParameters::debugLevel>=10)
		std::cout << "Finished.\a\n";
	//std::getchar();

	doubleDelete2D_Fast(dataPtsOrig);
	doubleDelete2D_Fast(dataPtsNorm);
	doubleDelete2D_Fast(dataPtsMask);
	doubleDelete2D_Fast(dataPtsMaskNorm);

	delete []orgMean;
	orgMean = 0;

	delete []orgStdDev;
	orgStdDev = 0;

	delete []mskMean;
	mskMean = 0;

	delete []mskStdDev;
	mskStdDev = 0;

	return EXIT_SUCCESS;
}
