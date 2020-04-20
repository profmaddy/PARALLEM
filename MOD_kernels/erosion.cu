
//#include "lem.h"
#include "erosion.h"
#include <math.h>
#include "device_constants.cuh"
#include "config.h"

#include "updates.h"
#include "mfd_accum.h"

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

	calc_diff_erosion(data, device);

	calc_conc_erosion(data, device);
	calc_gelifluction(data, device);

	fflush(data->outlog);

	//calc_sedflux(data, device);

	checkCudaErrors( cudaMemcpy ( device->mask, data->mask, full_size * sizeof(int), cudaMemcpyHostToDevice) );

	// sedmfdaccum(data, device);
	fprintf(data->outlog, "MOD: returned from sedmfdaccum :%s\n", cudaGetErrorString(cudaGetLastError()));

	checkCudaErrors( cudaMemcpy ( data->eroPtr,   device->eroPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->geliPtr,  device->geliPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->inciPtr,  device->inciPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->depoPtr,  device->depoPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy ( data->SlopePtr, device->SlopePtr, full_size * sizeof(double), cudaMemcpyDeviceToHost) );
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

	// Now grow the vegetation
	update_vegetation(data,device);

	checkCudaErrors( cudaMemcpy( data->TotBPtr,  device->TotBPtr,  full_size * sizeof(double), cudaMemcpyDeviceToHost) );
	fprintf(data->outlog, "MOD: mem copyback TotBn :%s\n", cudaGetErrorString(cudaGetLastError()));

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


