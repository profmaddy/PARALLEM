#include "io.h"
#include <errno.h>
#include "Data.h"

/******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full path AND if path leaf is a dir.
 *
 * @return  >0 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          <0 if an error occurred (errno is also set)
 *****************************************************************************/
int dirExists(const char* const path)
{
	struct stat info;

	int statRC = stat(path, &info);
	if (statRC != 0)
	{
		if (errno == ENOENT) { return 0; } // something along the path does not exist
		if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
		return -1;
	}

	return (info.st_mode & S_IFDIR) ? 1 : 0;
}

int readgrid(Data* data)
{
	struct header dem_head;
    double *next;
    int x,y,z;
    int cell;

	//printf("%s", data->demfile);
    FILE *in = fopen(data->demfile, "r");
	//FILE* in = fopen("30mfam.asc", "r");
	if(in == NULL) {
		perror("Cannot open file: ");
		exit(0);
	}

/* typical header looks like this
	ncols         1542
	nrows         2032
	xllcorner     707862.85644531
	yllcorner     4467835.0488281
	cellsize      30
	NODATA_value  -9999
*/
	fscanf(in, "%*s%d", &dem_head.ncols);
    fscanf(in, "%*s%d", &dem_head.nrows);
    fscanf(in, "%*s%f", &dem_head.xllcorner);
    fscanf(in, "%*s%f", &dem_head.yllcorner);
    fscanf(in, "%*s%f", &dem_head.cellsize);
    fscanf(in, "%*s%f", &dem_head.nodata);
    
    data->mapInfo.width = dem_head.ncols;
    data->mapInfo.height = dem_head.nrows;
    data->mapInfo.xllcorner =  dem_head.xllcorner;
    data->mapInfo.yllcorner = dem_head.yllcorner;
    data->mapInfo.cellsize = dem_head.cellsize;
    data->mapInfo.nodata = dem_head.nodata;

    fprintf(data->outlog,"ncols %d \n", data->mapInfo.width);
    fprintf(data->outlog,"nrows %d \n", data->mapInfo.height);
    fprintf(data->outlog,"xllcorner %lf \n", data->mapInfo.xllcorner);
    fprintf(data->outlog,"yllcorner %lf \n", data->mapInfo.yllcorner);
    fprintf(data->outlog,"cellsize %lf \n", data->mapInfo.cellsize);
    fprintf(data->outlog,"nodata %lf \n", data->mapInfo.nodata);

    int fullsize;
    fullsize = data->mapInfo.width * data->mapInfo.height;
	data->dem = (double *) malloc (fullsize * sizeof(double));

    for (x = 0; x < (data->mapInfo.height); ++x)
    {
    	for (y= 0; y < (data->mapInfo.width); ++y)
			{
    		cell = 	x*data->mapInfo.width + y;
            fscanf(in, "%lf", &data->dem[cell]);
            //if (x<10) printf ("%lf ,",data->dem[cell]); // note this is reading an integer into a double
			}
     }


    fclose(in);
	printf("\nFinished reading %s...\n", data->demfile);
	fprintf(data->outlog, "Finished reading %s...\n", data->demfile);
	fflush(data->outlog);
    return 0;
}

int write_float(Data *data, float *grid, char *filename)
{
    int elements;
    double *next;
    int x, y;
    int cell;
    FILE *out = fopen(filename, "w+");

    fprintf(out, "ncols %d\n", data->mapInfo.width);
    fprintf(out, "nrows %d\n", data->mapInfo.height);
    fprintf(out, "xllcorner %lf\n", data->mapInfo.xllcorner);
    fprintf(out, "yllcorner %lf\n", data->mapInfo.yllcorner);
    fprintf(out, "cellsize %lf\n", data->mapInfo.cellsize);
    fprintf(out, "NODATA_value %lf\n", data->mapInfo.nodata);

    for (x = 0; x < (data->mapInfo.height); ++x)
    {
		for (y = 0; y < data->mapInfo.width; ++y)
		{
			cell = x * data->mapInfo.width + y;
			if (data->mask[cell] == 1) 
				{ 
					fprintf(out, "%f ", grid[cell]); 
				}
			else 
				{
					fprintf(out, "%d ", -9999);
				}
		
		}
    	fprintf(out,"\n");
    }
    fclose(out);
    return 0;
}

int write_int(Data *data, int *grid, char *filename)
{
    int elements;
    double *next;
    int x, y;
    int cell;
    FILE *out = fopen(filename, "w+");

    fprintf(out, "ncols %d\n", data->mapInfo.width);
    fprintf(out, "nrows %d\n", data->mapInfo.height);
    fprintf(out, "xllcorner %lf\n", data->mapInfo.xllcorner);
    fprintf(out, "yllcorner %lf\n", data->mapInfo.yllcorner);
    fprintf(out, "cellsize %lf\n", data->mapInfo.cellsize);
    fprintf(out, "NODATA_value %lf\n", data->mapInfo.nodata);

    for (x = 0; x < (data->mapInfo.height); ++x)
    {
    	for (y= 0; y < data->mapInfo.width; ++y)
			{
    		cell = 	x*data->mapInfo.width + y;
			if (data->mask[cell] == 1)
				{
					fprintf(out, "%d ", grid[cell]);
				}
			else
				{
					fprintf(out, "%d ", -9999);
				}
			}
    	fprintf(out,"\n");
    }
    fclose(out);
    return 0;
}

int write_double(Data *data, double *grid, char *filename)
{
    int elements;
    double *next;
    int x, y;
    int cell;
    FILE *out = fopen(filename, "w+");

    fprintf(out, "ncols %d\n", data->mapInfo.width);
    fprintf(out, "nrows %d\n", data->mapInfo.height);
    fprintf(out, "xllcorner %lf\n", data->mapInfo.xllcorner);
    fprintf(out, "yllcorner %lf\n", data->mapInfo.yllcorner);
    fprintf(out, "cellsize %lf\n", data->mapInfo.cellsize);
    fprintf(out, "NODATA_value %lf\n", data->mapInfo.nodata);

    for (x = 0; x < (data->mapInfo.height); ++x)
    {
    	for (y= 0; y < data->mapInfo.width; ++y)
			{
    		cell = 	x*data->mapInfo.width + y;
			if (data->mask[cell] == 1)
				{
					fprintf(out, "%lf ", grid[cell]);
				}
			else
				{
					fprintf(out, "%d ", -9999);
				}
			}
    	fprintf(out,"\n");
    }
    
    fclose(out);
    return 0;
}

void writegrids(Data* data, int iteration)
{
	// Write the recalculated DEM
	sprintf(data->heightfile, "%s/%d_height.asc", data->matrixDIR, iteration);
	sprintf(data->FDfile, "%s/%d_FD.asc", data->matrixDIR, iteration);
	sprintf(data->FAfile, "%s/%d_FA.asc", data->matrixDIR, iteration);
	sprintf(data->erofile, "%s/%d_ero.asc", data->matrixDIR, iteration);
	sprintf(data->incifile, "%s/%d_inci.asc", data->matrixDIR, iteration);
	sprintf(data->gelifile, "%s/%d_geli.asc", data->matrixDIR, iteration);
	sprintf(data->depofile, "%s/%d_depo.asc", data->matrixDIR, iteration);
	sprintf(data->slopefile, "%s/%d_slope.asc", data->matrixDIR, iteration);
	sprintf(data->Precipfile, "%s/%d_precip.asc", data->matrixDIR, iteration);
	sprintf(data->Tempfile, "%s/%d_temp.asc", data->matrixDIR, iteration);
	sprintf(data->soilTfile, "%s/%d_soilT.asc", data->matrixDIR, iteration);
	
	if ((iteration % 50) == 0)
	{
		write_int(data, data->fd, data->FDfile);
		write_double(data, data->fa, data->FAfile);
		write_double(data, data->dem, data->heightfile);
	}

	if ((iteration % 100) == 0)
	{
		write_double(data, data->eroPtr, data->erofile);
		write_double(data, data->inciPtr, data->incifile);
		write_double(data, data->geliPtr, data->gelifile);
		write_double(data, data->depoPtr, data->depofile);
	}
}


int writeSummaryDataToFile(Data* data, int iteration)
{
FILE* out1;

			if (iteration == data->start_iter)
			{
				out1 = fopen(data->outfilename, "w");
				fprintf(out1, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n",
						      "iteration","year", "temp", "rain", "diffuse", "concentrated", "gelifluction", "chemicalW", "physicalW", "deposition", "totbio", "FA_max", "sedoutflow", "waterout", "sinkcount");
				fclose(out1);
			}


	int thisyear;
	thisyear =  iteration;
	// if (iteration > 0) thisyear = data->year[iteration];


	out1=fopen(data->outfilename,"a");
	double sedoutflow = ( (data->totE + data->totI + data->totG - data->totD) / 1300000 )  ; // outputs average per cell in mm, erosion already expressed in mm?

	fprintf(data->outlog, "Ave erosion rate mm/a : %10.8lf \n", sedoutflow);
	double totalrunoff = ((data->FA_max)* 1000); // 3982559 ; // water depth is in m [see runoffweights] therefore * 1000 for mm div no of cells to get rainfall equivalent
	fprintf(out1, "%d, %d, %f, %f, %10.8lf, %10.81lf, %10.8lf, %10.8lf, %10.8lf, %10.8lf, %10.8lf, %10.3lf, %10.8lf,  %10.5lf, %d \n",
			        iteration, thisyear, data->temp, data->rain, data->totE, data->totI, data->totG, data->totweatherC, data->totweatherP, data->totD, data->totbio, data->FA_max,  sedoutflow, totalrunoff, data->sinkcount);
	fflush(out1);
	fclose(out1);

	//free(outfilename);
	return 0;
}

int ReadClimateDataFromFile(Data* data)
{
	char a[20];
	char b[20];
	char c[20];
	int year;
	float changetemp;
	float changerain;

	// create sufficient space for climate data
	data->rainchange = (double *) calloc (data->max_iterations, sizeof(double));
	data->tempchange = (double *) calloc (data->max_iterations, sizeof(double));
	data->year =       (int *)    calloc (data->max_iterations, sizeof(int));

	// now read in the data
	printf("Climate file identified: %s \n", data->clim_file);
	FILE* in = fopen(data->clim_file, "r");
	//FILE* in = fopen("climateF.txt", "r");
	if (in == NULL) {
		perror("Cannot open file: ");
		return 1;
	}
	/*
    if (!data->clim_file)
        {
			printf("Climate file does not exist \n ");
	        exit (0);
        }*/


	// read the header
	fscanf(in, "%s", a); // year
	fscanf(in, "%s", b); // ppt
	fscanf(in, "%s", c); // temp

	printf("headings are %s %s %s \n", a, b, c);



	for(int i = 0; i < data->max_iterations; i++)
	{

		fscanf(in, "%d", &year);
		data->year[i] = (int) year;
		fscanf(in, "%f", &changerain);
		data->rainchange[i] = (double) changerain;
		fscanf(in, "%f", &changetemp);
		data->tempchange[i] = (double) changetemp;
	}

	fclose(in);
	fprintf(data->outlog,"climate file read \n");
	printf("Read Climate Data\n");
	return 0;
}

int checkparread(Data* data)
{
	printf("This is model : %s\n", data->modelcode);
	printf("Landscape file is: %s\n", data->demfile);
	printf("Climate file is: %s\n", data->clim_file);
	printf("Log file is: %s\n", data->logfileOut);
	printf("Summary csv file is: %s\n\n", data->outfilename);
	printf("restart is : %d\n", data->restart);
	printf("start is     :  %d\n", data->start_iter);
    printf("finish is    : %d\n", data->max_iterations );
    printf("no. of runin iterations is   :  %d\n", data->runiniter );
    printf("flow directions type  is  : %d\n", data->flowdirectiontype);
	printf("start rain is  : %f\n", data->rain);
    printf("start temp is: %lf\n", data->temp );
    printf("stone   : %lf\n", data->stoneval);
    printf("fines     : %lf\n", data->fineval );
    printf("struct   : %d\n", data->soil_struct );
    printf("perm    : %d\n", data->profile_permeability);
    printf("soilT     : %lf\n", data->soilTval );
    printf("soilM    : %lf\n", data->soilMval );
    printf("soilB     : %lf\n", data->soilBval );
    printf("nut        : %lf\n", data->nutval );
    printf("Qthres  : %lf\n", data->Qthreshold );
    printf("Gslope  : %lf\n", data->Gslope );
    printf("kappa     : %lf\n", data->kappa );
    printf("Gflowlim  : %lf\n", data->Gflowlim );
    printf("Flow_A  : %lf\n", data->Flow_A );
    printf("Flow_B  : %lf\n", data->Flow_B );
    printf("Max_in   : %lf\n", data->max_incision );
    printf("dum3     : %lf\n", data->dummy3 );
    printf("dum4     : %lf\n", data->dummy4 );
    printf("ddk        : %lf\n", data->ddk );
    printf("dck        : %lf\n", data->dck );
    printf("dgk        : %lf\n", data->dgk );
	printf("t_inf       : %lf\n", data->temp_inflex );
    printf("t_sens   : %lf\n\n", data->temp_sens);

    return(1);
}

int stripblank(char* q, char* w)
{
    do {
      if (*q != ' ')
        *(w++) = *q;
    } while (*(q++));

    return 0;
}

void choppy(char* s)
{
	s[strcspn(s, "\n")] = '\0';
}

void readinpar(Data* data, const char* file)
{

	  FILE *parfile = fopen ( file, "r" );
	  char s[39];
	  char* fname1= (char*) malloc(sizeof(char) *50);
	  char* fname2= (char*) malloc(sizeof(char) *50);
	  char* fname3= (char*) malloc(sizeof(char) *50);

	  if (parfile != NULL)
	  {
		  int i=0;
		  //if you use **lines allocate size for all lines with the loop or else you can allocate size inside loop and then read.
		  while((fgets(data->lines[i], 128, parfile)!=NULL) && (i<39))
		      {		     // printf("%s",lines[i++]);
			  i++;
		      }
		  fclose(parfile);

		  int parcount = 0;

          //#What is this file
		  strncpy_s(s, 39, data->lines[parcount], 38);
		  //choppy(s); // removbe the \n from the string
	      printf("\n\n%s\n",s );	      parcount++ ;

	      //#Date of file
	      strncpy_s(s, 39, data->lines[parcount], 38);
		  choppy(s);
	      printf("%s\n",s );	      parcount++ ;

	      //#model identifier
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 38);	      parcount++ ;
		  stripblank(data->dummystring, data->modelcode);
		  choppy(data->modelcode);

		  data->matrixDIR = data->modelcode;

		  /*
			  DIR* dir = opendir(data->matrixDIR);
			  if (dir)
			  {
				  // Directory exists. //
				  closedir(dir);
			  }
			  else if (ENOENT == errno)
			  {
				  printf("\n\nOutput directory %s does not exist\n", data->matrixDIR);
				  exit(0);
			  }
		  */

		  printf("%d\n", dirExists(data->matrixDIR));

	      // use model code to setup output files and directory
	      sprintf(fname1, "%s/%s_logfile.txt", data->matrixDIR, data->modelcode);
		  data->logfileOut = fname1;
	      
		  sprintf(fname2, "%s/%s_budgets.csv", data->matrixDIR, data->modelcode);
		  data->outfilename= fname2;
          
 	      //#Initial DEM file
		  char temp[39];
		  strncpy_s(data->dummystring, 39, data->lines[parcount], 38);	      parcount++ ;
          stripblank(data->dummystring, data->demfile);
		  choppy(data->demfile);

	      //#climate filename
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 38);	      parcount++ ;
          stripblank(data->dummystring, data->clim_file);
		  choppy(data->clim_file);

	      //#bedrock filename
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 38);	      parcount++ ;
          stripblank(data->dummystring, data->bedrockfile);
		  choppy(data->bedrockfile);

	      //restart yes(1) or no(0)
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 10);          parcount++ ;
		  data->restart = (int)strtol(data->dummystring, (char**)NULL, 10);

		  //start iteration
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 10); parcount++;
	      data->start_iter = (int)strtol(data->dummystring, (char**)NULL, 10);

	       //#finish iteration (Max iteration)
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 10); parcount++;
	      data->max_iterations = (int)strtol(data->dummystring, (char**)NULL, 10);

	      //#number of runin iterations
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 10); parcount++;
	      data->runiniter = (int)strtol(data->dummystring, (char**)NULL, 10);

	      //flow direction type
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 10); parcount++;
	      data->flowdirectiontype= (int)strtol(data->dummystring, (char**)NULL, 10);
		  printf("flowdirectiontype %d \n", data->flowdirectiontype);

          //#rain for runin
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->rain = strtod(data->dummystring, NULL);

	      //#temperature for runin
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->temp = strtod(data->dummystring, NULL);

	  	  // stonePtr
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->stoneval = strtod(data->dummystring, NULL);

	  	  // finePtr
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->fineval = strtod(data->dummystring, NULL);

	  	 // soil structure
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->soil_struct= (int)strtol(data->dummystring, (char**)NULL, 10);

	  	 // soil permeability
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->profile_permeability= (int)strtol(data->dummystring, (char**)NULL, 10);

	  	 // soil thickness
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20); parcount++;
	      data->soilTval = strtod(data->dummystring, NULL);

	  	 // soil moisture
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->soilMval = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // soil organics (bio)
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->soilBval = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // soil nutrients
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->nutval = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // Q threshold
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->Qthreshold = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // gelifluction slope threshold
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->Gslope = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // gelifluction kappa
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->kappa = strtod(data->dummystring, NULL);	      parcount++ ;

	  	  // geliflucution Gflowlim
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->Gflowlim = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // Flow A
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->Flow_A = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // Flow B
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->Flow_B = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // maximum incision
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->max_incision = strtod(data->dummystring, NULL);	      parcount++ ;

	  	// dummy3
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->dummy3 = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // dummy4
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->dummy4 = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // Diffused ddk
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->ddk = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // Concentrated dck
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->dck = strtod(data->dummystring, NULL);	      parcount++ ;

	  	 // Gelifluction dgk
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
	      data->dgk = strtod(data->dummystring, NULL);	      parcount++ ;

		 // temp_inflex for vegetation growth
	      strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
		  data->temp_inflex = strtod(data->dummystring, NULL);		  parcount++ ;

	    // temp_sesitivity for vegetation growth
		  strncpy_s(data->dummystring, 39, data->lines[parcount], 20);
		  data->temp_sens = strtod(data->dummystring, NULL);		  parcount++ ;

	  printf("All %d parameters read from file\n", parcount);
	  }
	  else {
		  printf("parameter file %s does not exist\n", file);
		  exit(0);
	  }
		data->last_rain = data->rain; // set up defaults
		data->last_temp = data->temp;
 return;

}


