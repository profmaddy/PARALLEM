
#include "memory_dev.h"
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include "helper_cuda.h"

//#include "helper_timer.h"
//#include "helper_functions.h"
//#include "cuda_runtime_api.h"
//#include "cuda_runtime.h"
//#include <device_launch_parameters.h>



void setslopesspace(Data* data, Data* device)
{
	 size_t freenow, total;
	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y * 9;

	 //checkCudaErrors( cudaMalloc((void**)&(device->Slopes), fullsize * sizeof(double)) );


	  cudaMemGetInfo(&freenow, &total);
	  fprintf(data->outlog, "CUDA card free after Slopes space allocated: %zd total: %zd \n",freenow/1024,total/1024);
	  fprintf(data->outlog, "slopes space allocated:on host and device  %s\n", cudaGetErrorString(cudaGetLastError()));
}

void setdevicespace_FD(Data* data, Data* device)
{
	 size_t freenow, total;
	 int fullsize;
	 int doublefull;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;
	 doublefull = fullsize* sizeof(double) *8;

	 checkCudaErrors( cudaMalloc((void**)&(device->fd), fullsize * sizeof(int)) );
	 checkCudaErrors( cudaMalloc((void**)&(device->SFD), fullsize * sizeof(int)) );

	 checkCudaErrors( cudaMalloc((void**)&(device->dx), 9 * sizeof(int)) );
	 checkCudaErrors( cudaMalloc((void**)&(device->dy), 9 * sizeof(int)) );
	 checkCudaErrors( cudaMalloc((void**)&(device->shortest_paths), fullsize * sizeof(float))  );
	 checkCudaErrors( cudaMalloc((void**)&(device->lowHeight),      fullsize * sizeof(double)) );
	 checkCudaErrors( cudaMalloc((void**)&(device->watershed_id),     fullsize * sizeof(int))    );
	 checkCudaErrors( cudaMalloc((void**)&(device->flatmask),     fullsize * sizeof(int))    );
	 fprintf(data->outlog, "FD: setdevicespace0:%s\n", cudaGetErrorString(cudaGetLastError()));

	 checkCudaErrors( cudaMalloc((void**)&(device->Slopes), doublefull) );
	 checkCudaErrors( cudaMalloc((void**)&(device->prop),   doublefull) );
	 fprintf(data->outlog, "FD: setdevicespace1:%s\n", cudaGetErrorString(cudaGetLastError()));

	  cudaMemGetInfo(&freenow, &total);
	  fprintf(data->outlog, "Memory on CUDA card free after FD space allocated: %zd total: %zd \n",freenow/1024,total/1024);
	  fprintf(data->outlog, "FD: setdevicespace2:%s\n", cudaGetErrorString(cudaGetLastError()));

}


void cleardevicespace_FD(Data* data, Data* device)
{
	size_t freenow, total;

		cudaFree(device->fd);
		cudaFree(device->SFD);
		// cudaFree(device->Slopes); //keep on device
		// cudaFree(device->prop); //keep on device

		cudaFree(device->dx);
		cudaFree(device->dy);
		cudaFree(device->shortest_paths);
		cudaFree(device->lowHeight);
		cudaFree(device->watershed_id);
		cudaFree(device->flatmask);
		fprintf(data->outlog, "FD: error after FD clear :%s\n", cudaGetErrorString(cudaGetLastError()));

		cudaMemGetInfo(&freenow, &total);
		fprintf(data->outlog, "FD: Memory on CUDA card free after FD space freed: %zd total: %zd \n\n",freenow/1024,total/1024);

}

void setdevicespace_FA(Data* data, Data* device)
{
	int full_size;
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	full_size= ncell_x * ncell_y;

	checkCudaErrors( cudaMalloc( (void**) &device->runoffweight, full_size * sizeof(double)) );

	checkCudaErrors( cudaMalloc( (void**) &device->rainmat, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->tempmat, full_size * sizeof(double)) );

	checkCudaErrors( cudaMalloc( (void**) &device->fa, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->fd, full_size * sizeof(int)) );

	//checkCudaErrors( cudaMalloc( (void**) &device->Slopes, full_size * 8 * sizeof(double)) );

	checkCudaErrors( cudaMalloc( (void**) &device->stonePtr, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->TotBPtr, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->soilMPtr, full_size * sizeof(double) ));

	// now copy the necessary data - these will not overlap becasue they are all on the same stream

	//checkCudaErrors(cudaSetDevice(0));
	checkCudaErrors( cudaMemcpy( device->dem, data->dem, full_size * sizeof(double), cudaMemcpyHostToDevice)) ; // copy the non-raised DEM back to GPU
	checkCudaErrors( cudaMemcpy( device->fd, data->fd, full_size * sizeof(int), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->runoffweight, data->runoffweight, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;

	checkCudaErrors( cudaMemcpy( device->rainmat, data->rainmat, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->tempmat, data->tempmat, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;

	checkCudaErrors( cudaMemcpy( device->stonePtr, data->stonePtr, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->TotBPtr, data->TotBPtr, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->soilMPtr, data->soilMPtr, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;

#ifndef PRODUCTION_RUN
	printf("FA: setdevicespace_FA:%s\n", cudaGetErrorString(cudaGetLastError()));
#endif

	fprintf(data->outlog, "FA: setdevicespace_FA:%s\n", cudaGetErrorString(cudaGetLastError()));
	fflush(data->outlog);
}

void cleardevicespace_FA(Data* data, Data* device)
{
	int full_size;
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	full_size = ncell_x * ncell_y;

	size_t freenow, total;
	checkCudaErrors(cudaMemcpy(data->prevfd, device->fd, full_size * sizeof(int), cudaMemcpyDeviceToHost));
	cudaFree(device->fd);
	cudaFree(device->runoffweight);
	cudaFree(device->fa);

	//cudaFree(device->contribA); // free it here as it is no longer needed)

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog, "FA: Memory on CUDA card free after FA space freed: %zd total: %zd \n\n",freenow/1024,total/1024);

	fprintf(data->outlog, "FA: cleardevicespace_FA:%s\n", cudaGetErrorString(cudaGetLastError()));
}


void setdevicespace_Process(Data* data, Data* device)
{
	size_t freenow, total;
	int full_size;
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	full_size= ncell_x * ncell_y;



	checkCudaErrors( cudaMalloc( (void**) &device->fa,       full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->fd,       full_size * sizeof(int))    );
	checkCudaErrors( cudaMalloc( (void**) &device->SFD,      full_size * sizeof(int))    );

	checkCudaErrors( cudaMalloc( (void**) &device->dz,       full_size * sizeof(double)) ); // create room for product dz
	checkCudaErrors( cudaMalloc( (void**) &device->finesPtr, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->soilTPtr, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->nutPtr,   full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->soilBPtr, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->eroPtr,   full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->geliPtr,  full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->inciPtr,  full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->depoPtr,  full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->weatherC, full_size * sizeof(double)) );
	checkCudaErrors( cudaMalloc( (void**) &device->weatherP, full_size * sizeof(double)) );

		fprintf(data->outlog, "MOD: setdevicespace_Process :%s\n", cudaGetErrorString(cudaGetLastError()));

		// stones, TotBio, soilM plus dem, slopes and mask still on device
		checkCudaErrors( cudaMemcpy ( device->fa,       data->fa,         full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->fd,       data->fd,         full_size * sizeof(int),    cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->SFD,      data->fd,         full_size * sizeof(int),    cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->SlopePtr,  data->SlopePtr,  full_size * sizeof(double), cudaMemcpyHostToDevice) );

		checkCudaErrors( cudaMemcpy ( device->finesPtr, data->finesPtr,   full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->soilTPtr, data->soilTPtr,   full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->nutPtr,   data->nutPtr,     full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->soilBPtr, data->soilBPtr,   full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->eroPtr,   data->eroPtr,     full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->geliPtr,  data->geliPtr,    full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->inciPtr,  data->inciPtr,    full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->depoPtr,  data->depoPtr,    full_size * sizeof(double), cudaMemcpyHostToDevice) );

		checkCudaErrors( cudaMemcpy ( device->weatherC,  data->weatherC,    full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->weatherP,  data->weatherP,    full_size * sizeof(double), cudaMemcpyHostToDevice) );

		fprintf(data->outlog, "MOD: Matrix memcopy operations :%s\n", cudaGetErrorString(cudaGetLastError()));

		cudaMemGetInfo(&freenow, &total);
		fprintf(data->outlog, "MOD: Memory on CUDA card free after model matrix space allocated: %zd total: %zd \n",freenow/1024,total/1024);
}

void cleardevicespace_Process(Data* data, Data* device)
{
	size_t freenow, total;

	cudaFree(device->fa);
	cudaFree(device->fd);
	cudaFree(device->SFD);
	cudaFree(device->fdmod);

	cudaFree(device->Slopes);
	//cudaFree(device->SlopePtr); // do not free here as it will not be redeclared.
	cudaFree(device->prop);

	cudaFree(device->rainmat);
	cudaFree(device->tempmat);

	cudaFree(device->dz);

	cudaFree(device->finesPtr);
	cudaFree(device->soilTPtr);
	cudaFree(device->nutPtr);
	cudaFree(device->soilBPtr);
	cudaFree(device->eroPtr);
	cudaFree(device->geliPtr);
	cudaFree(device->inciPtr);
	cudaFree(device->depoPtr);
	cudaFree(device->weatherC);
	cudaFree(device->weatherP);

	// free after being left at end of FA routines.
	cudaFree(device->stonePtr);
	cudaFree(device->TotBPtr);
	cudaFree(device->soilMPtr);

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog, "MOD: Memory on CUDA card free after model space freed: %zd total: %zd \n",freenow/1024,total/1024);
	fprintf(data->outlog, "MOD: Clear matrix operations :%s\n", cudaGetErrorString(cudaGetLastError()));

}

int copyMask(Data* data, Data* device)
{


	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;

	 checkCudaErrors( cudaMalloc( (void**) &device->mask, fullsize * sizeof(int)) ); // create space for the mask

	 checkCudaErrors( cudaMemcpy(device->mask, data->mask, fullsize * sizeof(int), cudaMemcpyHostToDevice) );  // copy back flag
	 fprintf(data->outlog, "Mask data sent to device %s\n", cudaGetErrorString(cudaGetLastError()));

	thrust::device_ptr<int> activecells = thrust::device_pointer_cast(device->mask);
	data->activecells  = thrust::count(activecells, activecells + fullsize, 1);

#ifndef PRODUCTION_RUN
	printf("No of active cells = %d \n", data->activecells);
#endif

	return 0;
}

int copylastclimate(Data* data, Data* device)
{

	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;

		checkCudaErrors( cudaMalloc( (void**) &device->last_rainmat, fullsize * sizeof(double)) );
		checkCudaErrors( cudaMalloc( (void**) &device->last_tempmat, fullsize * sizeof(double)) );

		checkCudaErrors( cudaMemcpy( device->last_rainmat, data->last_rainmat, fullsize * sizeof(double), cudaMemcpyHostToDevice)) ;
		checkCudaErrors( cudaMemcpy( device->last_tempmat, data->last_tempmat, fullsize * sizeof(double), cudaMemcpyHostToDevice)) ;

		fprintf(data->outlog, "Last climate data sent to device %s\n", cudaGetErrorString(cudaGetLastError()));

		//initial host matrices no longer needed. Will direct copy on only on GPU subsequently
		cudaFree(data->last_rainmat);
		cudaFree(data->last_tempmat);

		fprintf(data->outlog, "Error here? %s\n", cudaGetErrorString(cudaGetLastError()));

	return 0;
}



void createDeviceSpace(Data* data, Data* device)
{
	size_t freenow, total;

	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;

	 //cudaDeviceReset();

	  checkCudaErrors( cudaMalloc((void **)&(device->dem), fullsize * sizeof(double)) );
	  checkCudaErrors( cudaMalloc((void **)&(device->SlopePtr), fullsize * sizeof(double)) );
	  checkCudaErrors( cudaMalloc((void **)&(device->summary), fullsize * sizeof(double)) );

	fprintf(data->outlog,"Allocated DEM and slope matrices on device :%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog,"Memory on CUDA card free after device DEM and slope grids allocated: %zd total: %zd \n",freenow/1024,total/1024);

#ifndef PRODUCTION_RUN
	printf("Device space created \n");
#endif

}


int clearDeviceSpace(Data* data, Data* device)
{
	size_t freenow, total;

		//cudaFree(device->dem);
		cudaFree(device->SlopePtr);
		cudaFree(device->summary);
		cudaFree(device->Slopes);


	cudaMemGetInfo(&freenow, &total);
#ifndef PRODUCTION_RUN
	printf("Memory on CUDA card free after DEM and slope device grids space freed: %d total: %d \n",freenow/1024,total/1024);
#endif

	fprintf(data->outlog,"Memory on CUDA card free after DEM and slope device grids space freed: %zd total: %zd \n",freenow/1024,total/1024);





	return 0;
}

int zerogrids(Data* data)
{

	memset(data->eroPtr, 0.0, sizeof(data->eroPtr));
	memset(data->geliPtr, 0.0, sizeof(data->eroPtr));
	memset(data->inciPtr, 0.0, sizeof(data->inciPtr));
	memset(data->depoPtr, 0.0, sizeof(data->depoPtr));

	memset(data->fa, 0.0, sizeof(data->fa));

	return 0;
}

int copytolastclimate(Data* data, Data* device)
{
	int fullsize;
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	fullsize= ncell_x * ncell_y;

	checkCudaErrors( cudaMemcpy( device->last_rainmat, device->rainmat, fullsize * sizeof(double), cudaMemcpyDeviceToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->last_tempmat, device->tempmat, fullsize * sizeof(double), cudaMemcpyDeviceToDevice)) ;

	return 0;
}
