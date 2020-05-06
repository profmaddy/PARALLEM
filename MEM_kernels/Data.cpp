
#include "Data.h"
#include <sys/types.h>
#include "io.h"


#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include "helper_cuda.h"

//#include <dirent.h>
//#include <libgen.h>


int setClimate(Data* data, int iteration)
{	

	//double mean_rain = 850. ;
	//double mean_temp = 10.0 ;

	data->last_rain = data->rain; // set up defaults
	data->last_temp = data->temp; 
	//data->rain = mean_rain + mean_rain * data->rainchange[iteration-1];
	//data->temp = mean_temp + data->tempchange[iteration-1];

	data->rain = data->rainchange[iteration]; // from -1
	data->temp = data->tempchange[iteration]; // from 1

	fprintf(data->outlog,"rain: %f. temp: %f, last_rain: %f, last temp: %f \n", data->rain, data->temp, data->last_rain, data->last_temp);

	return 0;
}



int init_setoutletcell(Data* data) 
{
 double lowest ;
  lowest = 1000.;
  int cell;
  int index;
  index = 0;

  // what if the lowest cell is currently in the grid!
  for (int i = 0; i < data->mapInfo.height; i++)
  {
	  for (int j = 0; j < data->mapInfo.width; j++)
	  {
		  cell = i * data->mapInfo.width + j;
		  if (data->mask[cell] == 1) // I am in the catchment
		  {
			  if ((data->dem[cell]) < lowest) {
				  lowest = data->dem[cell];
				  index = cell;
			  }
		  }
	  }
  }
  data->outletcellidx = index;
  printf("identified outlet cell has index %d and value %f \n", index, data->dem[data->outletcellidx]);
  fprintf(data->outlog, "identified outlet cell has index %d and value %f \n", index, data->dem[data->outletcellidx]);

  return(0);
}


void budget_calc(Data* data){

	 double fullsize =  data->mapInfo.width * data->mapInfo.height;

	 int x;
	  for (x = 0; x < fullsize; x++)
	  {
	    data->diffdem[x] = data->dem[x] - data->prevdem[x];
	  }
}

void copyprevht(Data* data){

	int x;
	 double fullsize =  data->mapInfo.width * data->mapInfo.height;

	  for (x = 0; x < fullsize; x++)
	  {
	    data->prevdem[x] = data->dem[x];
	  }
}


int setnodatavalues(Data* data)
{
	  int fullsize;
	  int x;

	  fullsize =  data->mapInfo.width * data->mapInfo.height;

	  for (x = 0; x < fullsize; x++)
	  {
		  if (data->mask[x] == 0) // not inside the catchment(s)
		  {
			  data->fa[x]= -9999;
			  data->fd[x] = -9999;

		  }
	  }
	  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Process Data Matrices
//////////////////////////////////////////////////////////////////////////////

int setProcessMatrices(Data* data)
{
  int fullsize;
  int x; int y;

  //double calinitBio;
  double tempeffect;
  
  fullsize =  data->mapInfo.width * data->mapInfo.height;
  data->rescale = (10000 / ( data->mapInfo.cellsize * data->mapInfo.cellsize ));
  printf("rescale =  div %f \n",  data->rescale);
  
  //calinitBio = ( 3000.0 + (5.0 * (data->rain - 100.0) ) ) / data->rescale;  // scalar added 15/01/16

  for (x = 0; x < fullsize; x++)
  {
	  if (data->mask[x] == 1) // inside the catchment(s)
	  {
		  data->fd[x] = 0; // now set for MFD notation
		  data->SFD[x] = 0;
		  data->fa[x] = 0.0;
		  data->SlopePtr[x] = 0.0;
		  data->dz[x] =0.0;
		  data->prevdem[x] = 0.0;
		  data->diffdem[x] = 0.0;

		  data->runoffweight[x] = 	  1.0 ;
		  //data->soilTPtr[x] = 		 data->soilTval ;
		  data->stonePtr[x] = 		 data->stoneval ;
		  data->finesPtr[x] = 		 data->fineval ;
		  data->soilMPtr[x] = 		 data->soilMval ;
		  data->nutPtr[x] = 		 data->nutval ;
		  data->soilBPtr[x] = 		 data->soilBval ;

		  tempeffect = 100 + (100 / ( data->tempmat[x] + 15)) ; // make the second parameter 100 a variable
		  //data->TotBPtr[x] = 		( 3000.0 + (5.0 * (data->rainmat[x] - 100.0) ) ) / data->rescale;  // scalar added 15/01/16
		  data->TotBPtr[x] = 		( 3000.0 + (5.0 * (data->rainmat[x] - tempeffect) ) ) / data->rescale;  // added 11/11/17

	  } else // if not in catchment(s) set values to zero
	  {
		  data->fd[x] = 0;// altered for mfd 8
		  data->SFD[x] = 0;
		  data->fa[x] = 0.0 ;
		  data->dz[x]=0.0;
		  data->prevdem[x] = 0.0;
		  data->SlopePtr[x] = 0.0;
		  data->diffdem[x] = 0.0;

		  data->runoffweight[x] = 	  0.0 ;
		  data->soilTPtr[x] = 		  0.0 ;
		  data->stonePtr[x] = 		  0.0 ;
		  data->finesPtr[x] = 		  0.0 ;
		  data->soilMPtr[x] = 		  0.0 ;
		  data->nutPtr[x] =   		  0.0 ;
		  data->soilBPtr[x] = 		  0.0 ;
		  data->TotBPtr[x] =  		  0.0 ;
	  }
		  data->eroPtr[x] = 		  0.0 ;
		  data->geliPtr[x] = 		  0.0 ;
		  data->inciPtr[x] = 		  0.0 ;
		  data->depoPtr[x] = 		  0.0 ;
		  data->weatherC[x] = 		  0.0 ;
		  data->weatherP[x] = 		  0.0 ;
  }

  // set mfd slopes and props to zero
  for (y = 0; y<fullsize*8; y++)
    {
	  data->Slopes[y]= 0.0;
	  data->prop[y]= 0.0;
    }


  // data->outletcellidx = 0;// already set

  printf("All process grids now allocated \n\n");

  fprintf(data->outlog, "All process grids now allocated and default values set \n\n");

  int fdp;  //firstdatapoint
  for (int i=0; i<1000000; i++)
  {
	  if ((data->mask[i]) == 1) {
		  fdp = i;
		  printf("first data cell is %d \n\n", fdp);
		  break;
	  }
  }

  fprintf(data->outlog,"runoffweight: %f \nsoilThickness: %f \nstone%: %f, fines%: %f \nsoilMoisture: %f, nutrients: %f, soilBio: %f, TotalBio: %f \n",
		               data->runoffweight[fdp], data->soilTPtr[fdp], data->stonePtr[fdp], data->finesPtr[fdp], data->soilMPtr[fdp], data->nutPtr[fdp], data->soilBPtr[fdp], data->TotBPtr[fdp]);

  //write_double(data, data->TotBPtr, "bio.txt");


  return 0;
}

int SetClimateMatrices(Data* data, int iteration)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;	  int fullsize;
	int dataSize;
	double therain;
	double thetemp;

	printf("Setting climate gradients\n");
	fullsize = data->mapInfo.width * data->mapInfo.height;
	dataSize = fullsize * sizeof(double);

	if (iteration == -10000) // do this once only!
	{
		checkCudaErrors(cudaMallocHost((void**)&data->rainmat, dataSize));
		checkCudaErrors(cudaMallocHost((void**)&data->tempmat, dataSize));
		checkCudaErrors(cudaMallocHost((void**)&data->last_rainmat, dataSize));
		checkCudaErrors(cudaMallocHost((void**)&data->last_tempmat, dataSize));
	}

	data->last_rain = data->rain; // set up defaults
	data->last_temp = data->temp;

	if (iteration > -10000) // need to do this only to setup initial matrices
	{
		data->rain = data->rainchange[iteration]; // from -1
		data->temp = data->tempchange[iteration]; // from 1
	}

	printf("rain %f\n", data->rain);
	printf("temp %f\n", data->temp);
	printf("last rain %f\n", data->last_rain);
	printf("last temp %f\n", data->last_temp);


	double rainmax = data->rain * 2;// this multiplier needs to be a parameter
	double rainmin = data->rain;
	thetemp = data->temp;

	double increment = (rainmax - rainmin) / ncell_y;

	double minheight = 467.228; // this needs to be a parameter


	for (int i = 0; i < ncell_y; i++) {
		therain = rainmin + (i * increment);

		for (int j = 0; j < ncell_x; j++) {

			if (data->mask[i * ncell_x + j] == 1) // inside the catchment(s)
			{
				data->rainmat[i * ncell_x + j] = therain;
				data->tempmat[i * ncell_x + j] = thetemp - ((data->dem[i * ncell_x + j] - minheight) * 0.0065);
			}
			else // if not in catchment(s) set values to non value
			{
				data->rainmat[i * ncell_x + j] = -9999.;
				data->tempmat[i * ncell_x + j] = -9999.;
			}
		}
	}

	if (iteration == -10000) // need to do this only to setup initial matrices
	{


		double Lrainmax = data->last_rain * 2;// this multiplier needs to be a parameter
		double Lrainmin = data->last_rain;
		thetemp = data->last_temp;
		increment = (Lrainmax - Lrainmin) / ncell_y;

		for (int i = 0; i < ncell_y; i++) {
			therain = Lrainmin + (i * increment);

			for (int j = 0; j < ncell_x; j++) {

				if (data->mask[i * ncell_x + j] == 1) // inside the catchment(s)
				{
					data->last_rainmat[i * ncell_x + j] = therain;
					data->last_tempmat[i * ncell_x + j] = thetemp - ((data->dem[i * ncell_x + j] - minheight) * 0.0065);
				}
				else // if not in catchment(s) set values to non value
				{
					data->last_rainmat[i * ncell_x + j] = -9999.;
					data->last_tempmat[i * ncell_x + j] = -9999.;
				}
			}
		}

	}

	return 0;
}
