
#ifndef DATA
#define DATA

#include "lem.h"
#include "MapInfo.h"
//#include "production.h"
// #include "memory.h"


/*! This is the main data structure used to pass values across the sub-models */



typedef struct DataS {
    double* dem;
    double* prevdem;
	char* demfile;
	char* clim_file;
	const char* maskfile;
	const char* logfile;
	const char* outfilename;
	const char* matrixDIR;

	char lines[39][128]; // used for readin of par file

	char* modelcode;
	char* dummystring;

	MapInfo mapInfo;

	int restart;
	int start_iter;
    int runiniter;
	int max_iterations;
	int flowdirectiontype;

	double rescale;

	int* fd;
	int* prevfd;
	int* SFD;
	int* fdmod;
	double* fa;
	int* catchmentarea;

	double cellarea;
	double* runoffweight ;
	double* SlopePtr; // slope grid for SFD (x,y)
	double* Slopes; // slope-grid for MFD (x,y,z)
	int cpucontribA_max;

	double* finesPtr /*! fines_content% */ ;
	double fineval;
	double* stonePtr /*! stone content% */ ;
	double stoneval;
	double* soilTPtr /*! soil thickness [mm] */ ;
	double soilTval;
	double* soilMPtr /*! soil moisture % */ ;
	double soilMval;
	double* soilBPtr  /*! soil biomass % */ ;
	double soilBval;
	double* oldsoilBPtr /*! old soil biomass % */ ;
	double* nutPtr   /*! nutrient content % */ ;
	double nutval;
	double* TotBPtr  /*! total biomass calculated load in kg per ??*/ ;

	double* eroPtr /*! diffuse erosion grid */ ;
	double* geliPtr /*! gelifluction grid */ ;
	double* inciPtr /*! concentrated erosion grid */ ;
	double* depoPtr /*! deposition grid */ ;

	double totE /*! total annual diffuse erosion */ ;
	double totI /*! total annual concentrated erosion */ ;
	double totD /*! total annual deposition */ ;
	double totG /*! total annual gelifluction */ ;
	double totbio /*! total annual biomass */ ;
	double totweatherP;
	double totweatherC;

	double* summary /*! used to perform reduction stats using thurst */;

	int soil_struct  /*! soil structure unused at present */ ;
	int profile_permeability  /*! soil permeability unused at present */ ;

	double* weatherP /*! physical weathering */ ;
	double* weatherC /*! chemical weathering */ ;
	//double* evapPtr;

	double  adfGeoTransform[6]; // needed for gdal write operations to contain projection parameters

	int* dx;
	int* dy;
	double* dz;
	double* diffdem;
	float* shortest_paths;
    int* watershed_id;
    double* lowHeight;

	int* contribA;
	double Qthreshold;
	double Gslope;

	int outletcellidx;

	double rain;
	double temp;
	double* rainmat;
	double* tempmat;
	double last_rain;
	double last_temp;
	double* last_rainmat;
	double* last_tempmat;
	double* rainchange;
	double* tempchange;
	double* slopetotal;

	int* year;
	int* start_year;

	int* mask;
	int* flatmask;
	int* rivermask;

	int sinkcount;

	const char* logfileOut;
	FILE* outlog;

	double* prop;
	double FA_max;

	double temp_inflex;
	double temp_sens;

	int activecells;

	double kappa; //needed for Gelifluction
	double Gflowlim; // limit for channels where G is turned off
	double dummy3;
	double dummy4;

	double Flow_A; // conc erosion "soft" sediment erodibility
	double Flow_B; // conc erosion discharge exponent
	double* Flow_C; // conc erosion erodibility parameter i.e.bedrock erodibility

	double max_incision;
	double ddk;
	double dck;
	double dgk;

	//matrix output file names
	char* heightfile;
	char* dzfile;
	char* diff_file;
	char* FDfile ;
    char* FAfile ;

    char* Precipfile;
    char* Tempfile;

	char* erofile;
	char* incifile ;
	char* gelifile ;
	char* depofile;
	char* slopefile;
	char* finesfile ;
	char* stonesfile;
	char* totbiofile ;
	char* soilTfile ;
	char* nutfile;
	char* soilMfile;
	char* soilBfile;
	char* wCfile;
	char* wPfile;
	char* Burnfile;

	char* catchmap;
	char* catchmask;
	char* contrib ;
	char* rivermaskfile;
	char* flatfile;
	char* outputfilefile;

	char* bedrockfile; // for Flow_C

} Data, DataCopy;


typedef struct CatchmentS {

	int* watershed_id;
	int* catchmentarea;
	int* toptencatchmentArea;
	int* toptencatchmentID;
	int* toptenrootcell;
	int* mask;

} Catchment;


int checkparread(Data* data);

int SetClimateMatrices(Data* data, int iteration);

int setClimate(Data* data, int iterations);

int zerogrids(Data* data);

//int readgdalDEMfromFile(Data* data);

//int readBurnFDfromFile(Data* data);

int init_setoutletcell(Data* data);

//int readgdalBedrockfromFile(Data* data, Data* device);

void createSoilTfromformula(Data* data);

int createmask(Data* data);

//int readgdalmaskfromFile(Data* data) ;

//void savegrids(Data* data, Catchment* catchments,  int iteration);

void budget_calc(Data* data);

void copyprevht(Data* data);

//int writegdalGRIDtoFile(Data* data, Catchment* catchments, char* fileName, int whichtype, int what) ;

//int writeGRIDtoFile(Data* data, char* fileName, int whichtype, int what);

int writeSummaryDataToFile(Data* data, int iteration);

int setnodatavalues(Data* data);

int setProcessMatrices(Data* data);

#endif
