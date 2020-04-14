#include "eroincidep.h"
#include "updates.h"


__global__ void Update_Surface (int *mask, double base_level, double *weatherC_d, double *weatherP_d, double* zs_d, double* dz_d, double cell_size, int rows, int cols)
{

//  double surface_flow;
//  double evap_trans;
  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

  int self = irow * cols + icol;
  if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

 // update_surface
    //double dz = depoPtr_d[self] - eroPtr_d[self] - inciPtr_d[self];

	zs_d[self] = zs_d[self] + dz_d[self] - weatherC_d[self] - weatherP_d[self];


	if (zs_d[self] <= base_level) zs_d[self] = base_level + 0.01 ;//base_level;


}

__global__ void update_nutrient_content (int *mask, double* temperature, double* nutPtr_d, double* soilTPtr_d, double* soilBPtr_d, double* TotBPtr_d, int rows, int cols)
{

	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;

    if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

    int self = irow * cols + icol;
	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	//double old_soil_biomass = soil_biomass; // removed - variable not used anywhere
	if (soilTPtr_d[self] <= 0.0)
		nutPtr_d[self] = 0.0;
	double leaf_fall = TotBPtr_d[self] / (72.8 * log (TotBPtr_d[self] + 1.4));
	double decomposition = soilBPtr_d[self] * exp (0.12 * temperature[self] - 5.75);
	soilBPtr_d[self] = soilBPtr_d[self] + leaf_fall - decomposition;
	if (soilTPtr_d[self] == 0.0 && soilBPtr_d[self] > 0.0)
		nutPtr_d[self] = 100.0;
	else if (soilTPtr_d[self] == 0.0)
		nutPtr_d[self] = 0.0;
	else
		nutPtr_d[self] = (soilBPtr_d[self] / (1300. * soilTPtr_d[self])) * 100.;

	if (nutPtr_d[self] > 100) nutPtr_d[self] = 100.;
	//return old_soil_biomass;  //removed - does nothing// required so vegetation growth does calculations on correct values

}

__global__  void grow_vegetation (int *mask, double* temperature, double* last_temperature, double* rain, double* last_rain, double* TotBPtr, int rows, int cols,
		                                                            double temp_inflex, double temp_sens, double rescale)
{
	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;

    if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

    int self = irow * cols + icol;

	double delta_temp;
	delta_temp = temperature[self] - last_temperature[self];
	double delta_rain;
	delta_rain = rain[self] - last_rain[self];

	double total_cc;
	double tempeffect;
	tempeffect = 100 + (100 / ( temperature[self] + 15)) ;
	total_cc =  ( 3000.0 + (5.0 * (rain[self] - tempeffect) ) ) / rescale;

	double r_plus;
	double r_minus;

	if (mask[self] == 0) {
		TotBPtr[self]=0;
		return; // don't calculate if not in catchment(s) of interest
	}


	//  Temperature values centred around 6 degrees not 12 as in previous incarnation
	if (delta_temp > 0)
		{
            // increase
            r_plus = delta_temp * (1.0 + tanh ((temperature[self] - temp_inflex) / temp_sens));
            TotBPtr[self] += ((r_plus * TotBPtr[self]) * ( 1. - (TotBPtr[self] / total_cc))) ;
        }

	if (delta_temp < 0)
		{
            // decrease
             r_minus = delta_temp * (1.0 - tanh ((temperature[self] - temp_inflex) / temp_sens));
            TotBPtr[self] += ((r_minus * TotBPtr[self]) * (TotBPtr[self] / total_cc)) ;
        }

	if (delta_temp == 0)
	{
        // No change
        //double r_static = delta_rain * (1.0 - tanh ((temperature - 12.) / 2.5));
        //total_biomass += ((delta_rain * total_biomass) * ((total_biomass / total_cc)));
        //total_biomass += ((0.000001 * total_biomass) * ( 1. - (total_biomass / total_cc)));
		TotBPtr[self] =  total_cc ; //( 3000.0 + (5.0 * (rain - 100.0) ) ) / rescale ;
		//(((3000. * (1.0 - exp (-0.000664 * delta_rain)))) * (TotBPtr[self] / total_cc)) ;
	}

	if (TotBPtr[self] < 1)
	{
		TotBPtr[self] = 1.0;
	}

	//TotBPtr[self] = 22;// delta_temp; // used for checking kernel execution
}


void update_newSurface(Data* data, Data* device, int iter)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	dim3  threads( BLOCKCOLS, BLOCKROWS);
	dim3  grid((ncell_x/BLOCKCOLS)+1, (ncell_y/ BLOCKROWS)+1);

	double base_level = data->dem[data->outletcellidx] - 0.00012; // do not allow any point to go lower than the outlet cell but let it incise in line with uplift

	Update_Surface<<<grid, threads>>>(device->mask, base_level, device->weatherC, device->weatherP, device->dem, device->dz, data->mapInfo.cellsize,data->mapInfo.height, data->mapInfo.width);
	fprintf(data->outlog, "MOD: update surface :%s\n", cudaGetErrorString(cudaGetLastError()));

	checkCudaErrors( cudaMemcpy ( data->dem,  device->dem,  full_size * sizeof(double), cudaMemcpyDeviceToHost) ); //UPDATED SURFACE
	fprintf(data->outlog, "MOD: update surface memcopy back :%s\n", cudaGetErrorString(cudaGetLastError()));

	double max_dem;
	thrust::device_ptr<double> max_dem_d = thrust::device_pointer_cast(device->dem);
	max_dem = thrust::reduce(max_dem_d, max_dem_d + full_size, (double) 0, thrust::maximum<double>() );
	printf("update : Maximum dem height = %f \n", max_dem);

	// reset the outlet to lowest once the data is back from the GPU!
	data->dem[data->outletcellidx] = base_level ;

	thrust::device_ptr<double> demlow_d = thrust::device_pointer_cast(device->dem);
	thrust::device_ptr<double> summary_d = thrust::device_pointer_cast(device->summary);

	//predicate function is_not_zero defined in header
	thrust::copy_if(demlow_d, demlow_d + full_size, summary_d, is_not_negative());

	if (data->dem[data->outletcellidx] >= 467.22)  data->dem[data->outletcellidx] = 467.22 ; // do not allow outlet cell to rise
	data->dem[data->outletcellidx] = data->dem[data->outletcellidx] - 0.00012 ;

	// outlet cell is fixed in position now
	fprintf(data->outlog, "MOD: Min DEM ht : %f \n", data->dem[data->outletcellidx]);
	//printf("MOD: Min DEM ht : %f \n", data->dem[data->outletcellidx]);

}




void update_nutrients(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;

	dim3  threads( BLOCKCOLS, BLOCKROWS);
	dim3  grid((ncell_x/BLOCKCOLS)+1, (ncell_y/BLOCKROWS)+1);

	update_nutrient_content<<<grid, threads>>>(device->mask, device->tempmat, device->nutPtr, device->soilTPtr, device->soilBPtr, device->TotBPtr, data->mapInfo.height, data->mapInfo.width);
	fprintf(data->outlog, "MOD: update_nutrient_content :%s\n", cudaGetErrorString(cudaGetLastError()));
}




void update_vegetation(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	dim3  threads( BLOCKCOLS, BLOCKROWS);
	dim3  grid((ncell_x/BLOCKCOLS)+1, (ncell_y/ BLOCKROWS)+1);

	//printf("rain: %f ,  last_rain: %f, temp: %f, last_temp: %f, tempinflex: %f, tempsens:%f, rescale %f",
			 //  data->rain, data->last_rain, data->temp, data->last_temp, data->temp_inflex, data->temp_sens, data->rescale);

	grow_vegetation<<<grid, threads>>>(device->mask, device->tempmat, device->last_tempmat, device->rainmat, device->last_rainmat, device->TotBPtr, data->mapInfo.height, data->mapInfo.width,
			                                                                 data->temp_inflex, data->temp_sens, data->rescale);

	fprintf(data->outlog, "MOD: grow vegetation :%s\n", cudaGetErrorString(cudaGetLastError()));

	thrust::device_ptr<double> biotot_d = thrust::device_pointer_cast(device->TotBPtr);
	//cudaSetDevice(0);
	data->totbio = thrust::reduce(biotot_d, biotot_d + full_size, (double) 0);
	fprintf(data->outlog, "total TotBio from thrust is %10.8lf \n", data->totbio);



}

