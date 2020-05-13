
#include "runoffweight.h"
#include "io.h"

struct is_not_zero
{
	__host__ __device__ bool operator() (double x)
	{
		return (x != 0.0);
	}
};

struct is_not_negative
{
	__host__ __device__ bool operator() (double x)
	{
		return (x >= 0.0);
	}
};

__global__ void calcrunoffweight(int ncell_x, int ncell_y, double* ppn, int* mask, double* stonePtr, double* TotBPtr, double* soilMPtr, double* runoffweight, double rescale)
{
  double runoff;
  double runoff_coeff;
  double myppt;

  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  // If outside of DEM return (nothing to do)
  if(icol >= ncell_x || irow >= ncell_y)
    return;
  
  int self = irow * ncell_x + icol;

	  if (mask[self] == 0)
		  {
		  runoffweight[self] = -9999;
		  return; // don't calculate if not in catchment(s) of interest
		  }
  
    runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	* tanh (0.6 + (soilMPtr[self]) / 0.2);
 // runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	; // newer version is simplified?
 
	myppt = ppn[self] / 1000;// ppn/1000 is depth of rain in m, runoff is therefore a depth in m
		
	runoffweight[self] = myppt * runoff_coeff;
	
	soilMPtr[self] += myppt * (1. - runoff_coeff);

	if ((soilMPtr[self]) > 0.5)
				{
					soilMPtr[self] = 0.5;
				}
	if (runoffweight[self] == 0) printf("weight = 0 IN CATCHMENT");
	
}


int computeRunOffWeights(Data* data, Data* device)
{
	//double currentRain;
	//double ppn;

	data->cellarea = data->mapInfo.cellsize * data->mapInfo.cellsize ;

	double rescale = 1 / (10000 / (data->mapInfo.cellsize * data->mapInfo.cellsize));


	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int fullsize = ncell_x * ncell_y;

	dim3  threads( BLOCKCOLS, BLOCKROWS);
	dim3  grid((ncell_x/BLOCKCOLS) +1, (ncell_y/ BLOCKROWS) +1);
	

	calcrunoffweight<<<grid, threads>>>(ncell_x, ncell_y, device->rainmat, device->mask, device->stonePtr, device->TotBPtr, device->soilMPtr, device->runoffweight, rescale);
	fprintf(data->outlog, "FA: calcrunoffweight:%s\n", cudaGetErrorString(cudaGetLastError()));

	thrust::device_ptr<double> runoff_d = thrust::device_pointer_cast(device->runoffweight);
	thrust::device_ptr<double> summary_d = thrust::device_pointer_cast(device->summary);

	//predicate function is_not_zero defined in header
	thrust::copy_if(runoff_d, runoff_d + fullsize, summary_d, is_not_negative());

	double maxRO;
	double minRO;
	double sumRO;
	double aveRO;

	//data->activecells = 1337310;

	maxRO = thrust::reduce(summary_d, summary_d + data->activecells, (double) 0, thrust::maximum<double>());
	minRO = thrust::reduce(summary_d, summary_d + data->activecells, (double) 10, thrust::minimum<double>());
	sumRO = thrust::reduce(summary_d, summary_d + data->activecells);
	aveRO = sumRO / data->activecells;

	fprintf(data->outlog, "FA: Max RunOff: %lf,  Min RunOff: %lf, Ave Runoff: %lf \n", maxRO, minRO, aveRO);
	printf("FA: Max RunOff: %lf,  Min RunOff: %lf, Ave Runoff: %lf \n", maxRO, minRO, aveRO);

	cudaMemcpy(data->runoffweight, device->runoffweight, sizeof(double)* ncell_x* ncell_y, cudaMemcpyDeviceToHost); // here just for checking
	fprintf(data->outlog, "FA: runoff weights calculated :%s\n", cudaGetErrorString(cudaGetLastError()));

	fflush(data->outlog);

	//write_double(data, data->runoffweight, "row.txt");

	thrust::fill(summary_d, summary_d + data->activecells, 0.0); // reset the summary grid

	return 1;
}

void calcprops(Data* data) {

	double slopetot;
	int width = data->mapInfo.width;
	int height = data->mapInfo.height;
	int fullsize = width * height;
	double propsum;
	int counter;
	int bigcounter;
	int zeroslope;
	counter = 0;
	bigcounter = 0;
	zeroslope = 0;
	double propdiff;
	int count;
	int slopeidx;
	int propdircount;

	slopeidx = 0;
	count = 0;

	for (int i = 0; i < fullsize; i++)
	{
		slopetot = 0;
		propsum = 0;
		propdiff = 0;
		
		if (data->mask[i] == 1) {

			slopeidx = i * 8;
			count++;
			propdircount = 0;

			//if (count < 30) printf("fd: %d \n", data->fd[i]);
			//if (count < 30) printf("slopes [0][1][2][3][4][5][6][7] :");
			for (int k = 0; k < 8; k++)
			{
				slopetot += data->Slopes[slopeidx + k];
				if ((data->Slopes[slopeidx + k]) < 0) printf("slope is negative %lf \n", data->Slopes[slopeidx + k]);

			}
			//if (count < 30) printf("total = %lf\n", slopetot);

			if (slopetot == 0) {
				zeroslope++;  // this is counting the flats with no slope?
				//printf("zeroslope at %d with fd[%d]\n", i, data->fd);
			}

			//if (count <30) printf("props @ %d [0][1][2][3][4][5][6][7] :", i);

			

			for (int k = 0; k < 8; k++)
			{
				data->prop[slopeidx + k] = data->Slopes[slopeidx + k] / slopetot;
				if (data->prop[slopeidx + k] > 0) propdircount++;
				propsum += data->prop[slopeidx + k];
			}

			if ((data->fd[i] == 1) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 2) 
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else 
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}
			
			if ((data->fd[i] == 2) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 3)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}

			if ((data->fd[i] ==4) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 4)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}

			if ((data->fd[i] == 8) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 5)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}
			
			if ((data->fd[i] == 16) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 6)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}

			if ((data->fd[i] == 32) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 7)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}

			if ((data->fd[i] == 64) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 0)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}
			
			if ((data->fd[i] == 128) && (propdircount > 1))
			{
				counter++;
				for (int z = 0; z < 8; z++) {
					if (z == 1)
					{
						data->Slopes[slopeidx + z] = 0.0001;
						data->prop[slopeidx + z] = 1;
					}
					else
					{
						data->Slopes[slopeidx + z] = 0;
						data->prop[slopeidx + z] = 0;
					}
				}
			}


			if (propsum == 0)printf("cell %d, fd[%d],[%lf] \n]", i, data->fd[i], propsum);
			//if (count < 30 || slopetot == 0) printf("total prop %lf\n", propsum);

			propdiff = abs( 1 - propsum) ;

			if ( propsum !=1.0 ) 
			{
				// printf("illegal propsum for cell %d diff is [%15.12f] \n", i, propdiff);
				if (propdiff > 0.01) bigcounter++;
			}
		}
	}

	printf("zero slope totals  %d \n", zeroslope);
	printf("total number of non zero diff proportion entries = %d \n", counter);
	printf("total number where propdiff > 0.01  = %d \n", bigcounter);

}

void calcwater(Data* data)
{
	double sumwater;
	sumwater = 0;
	int fullsize;
	fullsize = data->mapInfo.height * data->mapInfo.width;
	double sum;
	sum = 0;
	
	for (int i = 0; i < fullsize; i++)
	{
		if (data->mask[i] == 1)
		{
			sum = data->runoffweight[i];
			sumwater += sum;
		}
	}

	printf("Total amount of water in grid for runoff = %lf \n", sumwater);

}

