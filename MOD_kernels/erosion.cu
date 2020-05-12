
//#include "lem.h"
#include "erosion.h"
#include <math.h>
#include "device_constants.cuh"
#include "config.h"

#include "updates.h"
#include "mfd_accum.h"
#include <thrust/reduce.h>
#include <iostream>
#include "io.h"

void fillthesinks(Data* data) 
{
	int ncell_y = data->mapInfo.height;
	int ncell_x = data->mapInfo.width;
	int self;
	int cellx, celly, dcell;
	int upslope, downslope,flatcount, sinkcounter;
	double lowest;
	double thiscellht, targetcellht;
	int xmove[9] = { 0,1,1,1,0,-1,-1,-1,0 };
	int ymove[9] = { -1,-1,0,1,1,1,0,-1,0 };
	int* dx;
	int* dy;
	dx = &xmove[0];
	dy = &ymove[0];
	data->dx = dx;
	data->dy = dy;
	sinkcounter = 0;
	
	for (int irow = 0; irow < ncell_y; irow++)  //irow loop
	{
		for (int icol = 0; icol < ncell_x; icol++) //icol loop
		{
			self = irow * ncell_x + icol;
			if (data->mask[self] == 1) // in catchment
			{
				lowest = 2000;
				upslope = 0;
				downslope = 0;
				flatcount = 0;
				for (dcell = 0; dcell < 8; dcell++)
				{
					cellx = icol + data->dx[dcell];
					celly = irow + data->dy[dcell];
					targetcellht = data->dem[celly * ncell_x + cellx];
					if (targetcellht != -9999) // dont look outside the grid
					{
						if ((data->dem[self]) < targetcellht) upslope++;
						if ((data->dem[self]) > targetcellht) downslope++;
						if ((data->dem[self]) == targetcellht) flatcount++;
						if (targetcellht < lowest) lowest = targetcellht; // find the lowest neighbour
					}
				}
				if  (upslope > 7) {
					data->dem[self] = lowest; // +0.00000001; // fill the sink to level of lowest neighbour and create a slope
					sinkcounter++;
				}
			}
		}
	}
	data->dem[data->outletcellidx] = data->dem[data->outletcellidx] - 0.00012;
	printf("number of sinks filled %d \n", sinkcounter);
	write_double(data, data->dem, "filleddem.asc");
}

void erosionGPU(Data* data, Data* device, int iter)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	cudaEvent_t start, stop;
	float time;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

   if (cudaSuccess != cudaSetDevice(CUDA_DEVICE)){
    printf("Unable to access CUDA card\n");
    exit(0);
  }

	size_t freenow, total; 
		
	fprintf(data->outlog, "MOD: Starting Model Process Routines \n");

	//calc averageslope SlopePtr and transfer to device
	aveslope(data, device);

	calc_diff_erosion(data, device);
		
		thrust::device_ptr<double> difftot_d = thrust::device_pointer_cast(device->eroPtr);
		data->totE = thrust::reduce(difftot_d, difftot_d + full_size);
		fprintf(data->outlog, "total concentrated  from thrust is %10.8lf \n", data->totE);
		printf("total concentrated  from thrust is %f \n", data->totE);
		if (iter == -1) write_double(data, data->eroPtr, "differo.asc");

	calc_conc_erosion(data, device);
		
		thrust::device_ptr<double> incitot_d = thrust::device_pointer_cast(device->inciPtr);
		cudaSetDevice(0);
		data->totI = thrust::reduce(incitot_d, incitot_d + full_size, (double)0);
		fprintf(data->outlog, "total Incision from thrust is %10.8lf \n", data->totI);
		printf("total Incision from thrust is %10.8lf \n", data->totI);
		if (iter == -1) write_double(data, data->inciPtr, "concero.asc");

	calc_gelifluction(data, device);
		
		thrust::device_ptr<double> gelitot_d = thrust::device_pointer_cast(device->geliPtr);
		cudaSetDevice(0);
		data->totG = thrust::reduce(gelitot_d, gelitot_d + full_size, (double)0);
		fprintf(data->outlog, "total gelifluction from thrust is %10.8lf \n", data->totG);
		printf("total gelifluction from thrust is %10.8lf \n", data->totG);
		if (iter == -1) write_double(data, data->geliPtr, "geliero.asc");
	
	fflush(data->outlog);

	checkCudaErrors( cudaMemcpy ( device->mask, data->mask, full_size * sizeof(int), cudaMemcpyHostToDevice) );

	sedmfdaccum(data, device);
	
	if (iter == -1) write_double(data, data->depoPtr, "depo.asc");
	
	fprintf(data->outlog, "MOD: returned from sedmfdaccum :%s\n", cudaGetErrorString(cudaGetLastError()));

	checkCudaErrors( cudaMemcpy ( data->geliPtr,  device->geliPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	checkCudaErrors( cudaMemcpy ( data->depoPtr,  device->depoPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors(cudaMemcpy(data->SlopePtr, device->SlopePtr, full_size * sizeof(double), cudaMemcpyDeviceToHost));
	fprintf(data->outlog, "MOD: ero/inc/dep/slope memcopy :%s\n", cudaGetErrorString(cudaGetLastError()));

	calc_dz(data,device); // now includes gelifluction erosion

	checkCudaErrors( cudaMemcpy ( data->dz, device->dz, full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	// Now add in weathering products and update cell calibre and cell moisture data

	calc_weathering(data, device);

	// now copy back all updated matrices
	checkCudaErrors( cudaMemcpy ( data->finesPtr, device->finesPtr, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->stonePtr, device->stonePtr, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->soilTPtr, device->soilTPtr, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->soilMPtr, device->soilMPtr, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->weatherC, device->weatherC, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->weatherP, device->weatherP, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	fprintf(data->outlog, "MOD: fines/stone/soilT/soilM/weatherC and P memcopy back :%s\n", cudaGetErrorString(cudaGetLastError()));
	
	// Now update the surface height
	update_newSurface(data, device, iter);

	// Now update the nutrients on surface and in soil profile
	update_nutrients(data, device);

	checkCudaErrors( cudaMemcpy ( data->soilBPtr, device->soilBPtr, full_size * sizeof(double), cudaMemcpyDeviceToHost)) ;
	checkCudaErrors( cudaMemcpy ( data->nutPtr,   device->nutPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	fprintf(data->outlog, "MOD: conc_soilB/nutB copyback :%s\n", cudaGetErrorString(cudaGetLastError()));

	fflush(data->outlog);
	// Now grow the vegetation
	update_vegetation(data,device);

	checkCudaErrors( cudaMemcpy( data->TotBPtr,  device->TotBPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	fprintf(data->outlog, "MOD: mem copyback TotBn :%s\n", cudaGetErrorString(cudaGetLastError()));

	//checkCudaErrors(cudaMemcpy(data->dem, device->dem, full_size * sizeof(double), cudaMemcpyDeviceToHost));
	//write_double(data, data->dem, "dem2.asc");

	//fillthesinks(data); // use this until flooding is working
	//checkCudaErrors(cudaMemcpy(device->dem, data->dem, full_size * sizeof(double), cudaMemcpyHostToDevice));

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

#ifndef PRODUCTION_RUN
	printf("Time to complete model calculations %.6f s\n\n", time / 1000.0);
#endif

	fprintf(data->outlog, "MOD: time to complete flow accumulation %.6f s\n", time / 1000.0);

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog, "MOD: Memory on CUDA card free at end of erosion: %zd total: %zd\n\n",freenow/1024,total/1024);

}


