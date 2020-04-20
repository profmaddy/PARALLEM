
#include "runoffweight.h"
#include "io.h"
/*!

  Calculates the actual runoff from individual cell in m.

  OLD VERSION:: runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	* tanh (0.6 + (soilMPtr[self]) / 0.2);

  runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	; // newer version is simplified?

  runoff = ppn * runoff_coeff;  ppn is depth of rain in m, runoff is therefore a depth in m

*/

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
   
    //runoff = (ppn[self] * 0.001) * runoff_coeff; // ppn/1000 is depth of rain in m, runoff is therefore a depth in m
    runoff = ( 600 * 0.001) * runoff_coeff; // ppn/1000 is depth of rain in m, runoff is therefore a depth in m
    runoffweight[self] = runoff; // runoff; //* rescale; //actual runoff depth per cell in m altered 14/01/16


  soilMPtr[self] += (ppn[self]/1000) * (1. - runoff_coeff);

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

	data->activecells = 1337310;

	maxRO = thrust::reduce(summary_d, summary_d + data->activecells, (double) 0, thrust::maximum<double>());
	minRO = thrust::reduce(summary_d, summary_d + data->activecells, (double) 0, thrust::minimum<double>());
	sumRO = thrust::reduce(summary_d, summary_d + data->activecells);
	aveRO = sumRO / data->activecells;

	fprintf(data->outlog, "FA: Max RunOff: %lf,  Min RunOff: %lf, Ave Runoff: %lf \n", maxRO, minRO, aveRO);
	printf("FA: Max RunOff: %lf,  Min RunOff: %lf, Ave Runoff: %lf \n", maxRO, minRO, aveRO);

	cudaMemcpy(data->runoffweight, device->runoffweight, sizeof(double)* ncell_x* ncell_y, cudaMemcpyDeviceToHost); // here just for checking
	fprintf(data->outlog, "FA: runoff weights calculated :%s\n", cudaGetErrorString(cudaGetLastError()));

	write_double(data, data->runoffweight, "row.txt");

	thrust::fill(summary_d, summary_d + data->activecells, 0.0); // reset the summary grid

	fflush(data->outlog);

	return 1;
}
