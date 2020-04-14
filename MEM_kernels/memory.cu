#include "memory.h"


int createcontribAspace(Data* data)
{
	int fullsize;
	int dataSize;
	fullsize =  data->mapInfo.width * data->mapInfo.height;
	dataSize = fullsize * sizeof(int);
	data->contribA = (int *) malloc(dataSize);
	fprintf(data->outlog,"Host memory allocation for contribA  \n");
	return 0;
}


int clearcontribAspace(Data* data)
{
	free(data->contribA);
	//free(data->watershed_id); // need to clear this?

	return 0;
}

int createfilenamespace(Data* data)
{
	    data->heightfile = (char*) malloc(sizeof(char) *100);
	    data->diff_file= (char*) malloc(sizeof(char) *100);
	    data->FDfile = (char*) malloc(sizeof(char) *100);
	    data->FAfile = (char*) malloc(sizeof(char) *100);

	    data->Precipfile = (char*) malloc(sizeof(char) *100);
	    data->Tempfile = (char*) malloc(sizeof(char) *100);

	    data->erofile = (char*) malloc(sizeof(char) *100);
	    data-> incifile = (char*) malloc(sizeof(char) *100);
	    data->gelifile = (char*) malloc(sizeof(char) *100);
	    data->depofile = (char*) malloc(sizeof(char) *100);
	    data->slopefile = (char*) malloc(sizeof(char) *100);

	    data->finesfile = (char*) malloc(sizeof(char) *100);
	    data->stonesfile = (char*) malloc(sizeof(char) *100);
	    data->totbiofile = (char*) malloc(sizeof(char) *100);
	    data->soilTfile = (char*) malloc(sizeof(char) *100);
	    data->nutfile = (char*) malloc(sizeof(char) *100);
	    data->soilMfile = (char*) malloc(sizeof(char) *100);
	    data->soilBfile = (char*) malloc(sizeof(char) *100);

	    data->wCfile = (char*) malloc(sizeof(char) *100);
	    data->wPfile = (char*) malloc(sizeof(char) *100);

	    data->catchmap = (char*) malloc(sizeof(char) *100);
	    data->catchmask = (char*) malloc(sizeof(char) *100);
	    data->contrib = (char*) malloc(sizeof(char) *100);
	    data->rivermaskfile = (char*) malloc(sizeof(char) *100);

	    data->flatfile = (char*) malloc(sizeof(char) *100);

	    data-> logfile = (char*) malloc(sizeof(char) *100);
	    data->outfilename = (char*) malloc(sizeof(char) *100);
	    data->matrixDIR = (char*) malloc(sizeof(char) *100);
	    data->modelcode = (char*) malloc(sizeof(char) *100);
	    data->outputfilefile = (char*) malloc(sizeof(char) *100);

	    data->bedrockfile = (char*) malloc(sizeof(char) *100);

	    data->demfile = (char*) malloc(sizeof(char) *100);
	    data->clim_file = (char*) malloc(sizeof(char) *100);
	    data->dummystring = (char*) malloc(sizeof(char) *100);

	    data->Burnfile = (char*) malloc(sizeof(char) *100);



	    return(1);

}
int createProcessMatrices(Data* data)
{
  int fullsize;
  int dataSize;
  int dataSizeInt;

  fullsize =  data->mapInfo.width * data->mapInfo.height;
  dataSize = fullsize * sizeof(double);
  dataSizeInt = fullsize * sizeof(int);

// these are the static grids in which data is stored from one iteration to the next ie. these are ONLY freed at the end of the simulation

  checkCudaErrors(cudaMallocHost((void **)&data->prevdem, dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->diffdem, dataSize));

  checkCudaErrors(cudaMallocHost((void **)&data->fd, dataSizeInt ));
  checkCudaErrors(cudaMallocHost((void **)&data->SFD, dataSizeInt ));
  fprintf(data->outlog, "Flow direction (fd and SFD) space on host allocated %s\n", cudaGetErrorString(cudaGetLastError()));

  checkCudaErrors(cudaMallocHost((void **)&data->fa,           dataSize));
  fprintf(data->outlog, "Flow accumulation space on host allocated %s\n", cudaGetErrorString(cudaGetLastError()));

  checkCudaErrors(cudaMallocHost((void **)&data->SlopePtr,     dataSize));

  checkCudaErrors(cudaMallocHost((void **)&data->flatmask, dataSizeInt  ));

// room to store the slopes and proportions in all directions
  checkCudaErrors(cudaMallocHost((void **)&data->Slopes,  dataSize * 8 ));
  checkCudaErrors(cudaMallocHost((void **)&data->prop,  dataSize * 8  ));


  checkCudaErrors(cudaMallocHost((void **)&data->runoffweight, dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->stonePtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->finesPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->soilMPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->soilBPtr,     dataSize));
  //checkCudaErrors(cudaMallocHost((void **)&data->soilTPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->nutPtr,       dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->TotBPtr,      dataSize));

  checkCudaErrors(cudaMallocHost((void **)&data->eroPtr,       dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->geliPtr,      dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->inciPtr,      dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->depoPtr,      dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->dz,           dataSize));

  checkCudaErrors(cudaMallocHost((void **)&data->weatherC,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->weatherP,     dataSize));

  fprintf(data->outlog, "All hosts matrices memory allocated %s\n", cudaGetErrorString(cudaGetLastError()));

  return 0;

}


int deleteProcessMatrices(Data* data)
{

	  //checkCudaErrors(cudaFreeHost(data->dem));
	  checkCudaErrors(cudaFreeHost(data->prevdem));
	  checkCudaErrors(cudaFreeHost(data->diffdem));
	  checkCudaErrors(cudaFreeHost(data->dz));
	  checkCudaErrors(cudaFreeHost(data->fd));
	  checkCudaErrors(cudaFreeHost(data->fa));
	  checkCudaErrors(cudaFreeHost(data->SlopePtr));

	  checkCudaErrors(cudaFreeHost(data->rainmat));
	  checkCudaErrors(cudaFreeHost(data->tempmat));
	  checkCudaErrors(cudaFreeHost(data->last_rainmat));
	  checkCudaErrors(cudaFreeHost(data->last_tempmat));

	  checkCudaErrors(cudaFreeHost(data->mask));
	  checkCudaErrors(cudaFreeHost(data->flatmask));
	  //checkCudaErrors(cudaFreeHost(data->Flow_C));
	  // free the slope and prop matrices
	  checkCudaErrors(cudaFreeHost(data->Slopes));
	  checkCudaErrors(cudaFreeHost(data->prop));

	  checkCudaErrors(cudaFreeHost(data->runoffweight));
	  checkCudaErrors(cudaFreeHost(data->stonePtr));
	  checkCudaErrors(cudaFreeHost(data->finesPtr));
	  checkCudaErrors(cudaFreeHost(data->soilMPtr));
	  checkCudaErrors(cudaFreeHost(data->soilBPtr));
	  checkCudaErrors(cudaFreeHost(data->soilTPtr));
	  checkCudaErrors(cudaFreeHost(data->nutPtr));
	  checkCudaErrors(cudaFreeHost(data->TotBPtr));

	  checkCudaErrors(cudaFreeHost(data->eroPtr)); 
	  checkCudaErrors(cudaFreeHost(data->geliPtr));
	  checkCudaErrors(cudaFreeHost(data->inciPtr));
	  checkCudaErrors(cudaFreeHost(data->depoPtr));
	  // checkCudaErrors(cudaFreeHost(data->dz));
	  checkCudaErrors(cudaFreeHost(data->weatherC));
	  checkCudaErrors(cudaFreeHost(data->weatherP));

	  fprintf(data->outlog, "All hosts matrices memory freed \n");

	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Setup store for catchment data ( needed for summary outputs etc)
//////////////////////////////////////////////////////////////////////////////

int createCatchmentSpace(Data* data, Catchment* Catchments) {
	//allocate space for catchment data and selective list and set values to zero
	Catchments->watershed_id = (int *) calloc(sizeof(int) , data->mapInfo.height * data->mapInfo.width);
	Catchments->mask = (int *) calloc(sizeof(int),  data->mapInfo.height * data->mapInfo.width); // all mask values set to zero

	fprintf(data->outlog, "Catchment space allocated \n");
	return 0;
}


void createSoilTfromformula(Data* data){

	int cell;

	checkCudaErrors(cudaMallocHost((void **)&data->soilTPtr, data->mapInfo.width * data->mapInfo.height * sizeof(double)));

	  // what if the lowest cell is currently in the grid!
	  for (int i = 0; i < data->mapInfo.height; i++){
			for (int j = 0; j < data->mapInfo.width; j++) {
				cell = 	i*data->mapInfo.width + j;
					if ( (data->dem[cell]) > 0 ) { data->soilTPtr[cell] = ( (data->dem[cell] - 400) / 1400 ) * 5;}
						else { data->soilTPtr[cell] = 0.0;}

				  }
	  }
	  printf( "Soil Thickness Data Created \n");
	  return;
}


int createmask(Data* data)
{

int width = data->mapInfo.width;
int height = data->mapInfo.height ;
int fullsize = width * height;
double nodataold;

nodataold = data->mapInfo.nodata;
printf("DEM old no data value = %.6f will be reset to -9999\n", data->mapInfo.nodata);

checkCudaErrors(cudaMallocHost((void **)&data->mask,     fullsize*sizeof(int)  ));
	for (int x = 0; x < fullsize; x++) {
		data->mask[x] = 1;

		if (data->dem[x] == -9999)
		{
			data->mask[x] = 0;
			data->dem[x] = -9999; //reset the no data value to -9999
		}

	}
	data->mapInfo.nodata = -9999;

 return 1;
}
