#include "lem.h"
#include "Data.h"
#include "io.h"
#include "Directions.h"
#include "mfd_simple.h"


__device__ int nghbrindex(int client, int nghbr, int gridCols)
{
    switch (nghbr)
    {
        case EAST:
            return client + 1;
        case SOUTHEAST:
            return client + gridCols + 1;
        case SOUTH:
            return client + gridCols;
        case SOUTHWEST:
            return client + gridCols - 1;
        case WEST:
            return client - 1;
        case NORTHWEST:
            return client - gridCols - 1;
        case NORTH:
            return client - gridCols;
        case NORTHEAST:
            return client - gridCols + 1;
        default:
            return -1;
    }
}

__global__ void CalcMFD(int *mask, double *hv, int *fd, double *fa, double* props, int *ok, unsigned int *progressd,  int* localOK, int gridCols, int gridRows, double *weights)
{
    int self;
    int proploc;
    double accum;
    double addhere;
	unsigned int theval = 99999999;
    int   nie, nise, nis, nisw, niw, ninw, nin, nine;
    int cnfd;
    int irow, icol;

    irow = blockIdx.y * blockDim.y + threadIdx.y;
    icol = blockIdx.x * blockDim.x + threadIdx.x;
    self = irow * gridCols + icol;

    // Are we outside the DEM?
    if(icol >= gridCols || irow >= gridRows) return; // if outside of DEM nothing to do
    //if (icol == 0) || irow == 0 ) return; // if outside of DEM nothing to do

    if (mask[self] != 1 ) return;


    if (ok[self] != 0) return;
	if (localOK[self] == 1) {
      ok[self] = 1;
      return;
    }

    __syncthreads();
    nie  = nghbrindex(self, EAST,      gridCols);
    nise = nghbrindex(self, SOUTHEAST, gridCols);
    nis  = nghbrindex(self, SOUTH,     gridCols);
    nisw = nghbrindex(self, SOUTHWEST, gridCols);
    niw  = nghbrindex(self, WEST,      gridCols);
    ninw = nghbrindex(self, NORTHWEST, gridCols);
    nin  = nghbrindex(self, NORTH,     gridCols);
    nine = nghbrindex(self, NORTHEAST, gridCols);

    accum = 0.;
    fa[self] = weights[self];

	  // ORIGINAL CODE directions
	  // NW  N NE    7  0  1
	  //  W  *  E    6  8  2
	  // SW  S SE    5  4  3

    // revised code directions
    // NW  N NE    32  64   128
    //  W  *  E    16   0   1
    // SW  S SE     8   4   2

    /*
    if (dcell == 0)  newaspect = 64;
    if (dcell == 1)  newaspect = 128;
    if (dcell == 2)  newaspect = 1;
    if (dcell == 3)  newaspect = 2;
    if (dcell == 4)  newaspect = 4;
    if (dcell == 5)  newaspect = 8;
    if (dcell == 6)  newaspect = 16;
    if (dcell == 7)  newaspect = 32;
    */

    cnfd = fd[nie];
    if (cnfd & WEST)
    {
        if (!ok[nie]) return;
        proploc = (nie * 8) + 6;// flowing west
        addhere = fa[nie] * props[proploc];
        accum += addhere; 
    }
    
    cnfd = fd[nise];
    if (cnfd & NORTHWEST)
    {
        if (!ok[nise]) return;
        proploc = (nise * 8) + 7;// flowing northwest
        addhere = fa[nise] * props[proploc];
        accum += addhere;
    }
    
    cnfd = fd[nis];
    if (cnfd & NORTH)
    {
        if (!ok[nis]) return;
        proploc = (nis * 8) + 0;//flowing north
        addhere = fa[nis] * props[proploc];
        accum += addhere;
    }
    
    cnfd = fd[nisw];
    if (cnfd & NORTHEAST)
    {
        if (!ok[nisw]) return;
        proploc = (nisw * 8) + 1;//flowing northeast
        addhere = fa[nisw] * props[proploc];
        accum += addhere;
    }
    
    cnfd = fd[niw];
    if (cnfd & EAST)
    {
        if (!ok[niw]) return;
        proploc = (niw * 8) + 2; // flowing east
        addhere = fa[niw] * props[proploc];
        accum += addhere;
    }
    
    cnfd = fd[ninw];
    if (cnfd & SOUTHEAST)
    {
        if (!ok[ninw]) return;
        proploc = (ninw * 8) + 3;//flowing southeast
        addhere = fa[ninw] * props[proploc];
        accum += addhere;
    }
    
    cnfd = fd[nin];
    if (cnfd & SOUTH)
    {
        if (!ok[nin]) return;
        proploc = (nin * 8) + 4;//flowing south
        addhere = fa[nin] * props[proploc];
        accum += addhere;
    }

    cnfd = fd[nine];
    if (cnfd & SOUTHWEST)
    {
        if (!ok[nine]) return;
        proploc = (nine * 8) + 5;//flowing southwest
        addhere = fa[nine] * props[proploc];
        accum += addhere;
    }

    if (accum < 0.0) accum = 0.0;
    fa[self] += accum;
    localOK[self] = 1;
    atomicInc(progressd, theval);
}


int processtheGrid(Data* data, Data* device, int loopMax, int percent, int gridRows, int gridColumns, int* okGrid,  int* localOK, int blockRows, int blockColumns, int dimBlock3, int* doneP)
{
  int loopForever = (loopMax < 0) ? 1 : 0;
  int *ok;
  int *localOK_d;
  unsigned int progressh, *progressd;
  int gridProgress;
  


  //allocate GPU memory
    checkCudaErrors(cudaMalloc((void **) &ok,                   gridRows * gridColumns * sizeof(int))    );
    checkCudaErrors(cudaMalloc((void **) &localOK_d,     gridRows * gridColumns * sizeof(int))   );
    checkCudaErrors(cudaMalloc((void **) &progressd,      sizeof(progressh)) );
    fprintf(data->outlog, "FA: FA memory allocation :%s\n", cudaGetErrorString(cudaGetLastError()));
      

  gridProgress = 0;
  int loop = 0;
  int grid1 = (data->mapInfo.width  / blockColumns ) +1 ;
  int grid2 = (data->mapInfo.height / blockRows )  +1 ;
  do {
	  	  dim3 dimGrid(grid1, grid2, 1);
	  	  dim3 dimBlock(blockColumns, blockRows, dimBlock3);
	  	  printf("Grid is %d by %d by 1\n", grid1, grid2);
	  	  int oneZero = 0;
	      int oneZero2 = 0;
	      //copy grids to GPU
	        checkCudaErrors(cudaMemcpy(ok,                    okGrid,    gridRows * gridColumns * sizeof(int),    cudaMemcpyHostToDevice) );
	        checkCudaErrors(cudaMemcpy(localOK_d,      localOK,  gridRows * gridColumns * sizeof(int),    cudaMemcpyHostToDevice) );
	        do {
	        		progressh = 0;
	        		oneZero2 = oneZero;

	        		checkCudaErrors(cudaMemcpy(progressd, &progressh, sizeof(progressh), cudaMemcpyHostToDevice) );

	        		//CalcMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SFD, device->fa, device->Slopes, ok, progressd, localOK_d, gridColumns, gridRows, device->runoffweight);
	        		CalcMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->fd, device->fa, device->prop, ok, progressd, localOK_d, gridColumns, gridRows, device->runoffweight);

        			checkCudaErrors(cudaMemcpy(&progressh, progressd, sizeof(progressh), cudaMemcpyDeviceToHost) );
        			//printf("FA: progressed copy back :%s\n", cudaGetErrorString(cudaGetLastError()));
        			gridProgress += progressh;
	        			if (progressh == 0)
	        					oneZero = 1;
	        					else {
	        							oneZero = 0;
	        							oneZero2 = 0;
	        					}
            } while ((progressh > 10 ) || (!oneZero2) ); // && ((double) *doneP * 100 / (gridRows * gridColumns)) < 99.0);

	        checkCudaErrors(cudaMemcpy(progressd, &progressh, sizeof(progressh), cudaMemcpyHostToDevice) );
	        loop++;
  } while (loop < loopMax || (loopForever && gridProgress) );

    checkCudaErrors(cudaMemcpy(okGrid, ok, gridRows * gridColumns * sizeof(int),    cudaMemcpyDeviceToHost) ) ;
	fprintf(data->outlog, "FA: okGrid copy back :%s\n", cudaGetErrorString(cudaGetLastError()));
  /* Free the GPU copies */
  cudaFree(ok);
  cudaFree(localOK_d);
  cudaFree(progressd);
  return gridProgress; // this does not seem to be copied back in this code??????
}


void correctmfdflow(Data* data, Data* device, int iter)
{
  int ncols = data->mapInfo.width;
  int nrows = data->mapInfo.height;
  int fullsize = nrows * ncols;
  int x;
  float time;
  int percent;
  int gridprogress;
  int doneP;
  int *okgrid;
  int *localOK;

  cudaEvent_t start, stop;
  percent = 0;
  gridprogress = 0;
  doneP = 0;
	// Set all values to 0.0
	for (int x = 0; x < fullsize; ++x) {
		data->fa[x] = 0.0;
	}

   int doublefull;
   doublefull = fullsize * sizeof(double) * 8;

   cudaMemcpy(device->fa, data->fa, fullsize * sizeof(double), cudaMemcpyHostToDevice);
   checkCudaErrors(cudaMemcpy(device->prop, data->prop, doublefull, cudaMemcpyHostToDevice));


  fflush(data->outlog);
 //* start the timer for correctflow;
 cudaEventCreate(&start);
 cudaEventCreate(&stop);
 cudaEventRecord(start, 0);
  okgrid      = (int *)    malloc(fullsize * sizeof(int));
  localOK     = (int *)    malloc(fullsize * sizeof(int));
  if (okgrid == NULL || localOK == NULL ) {
    printf("Not enough memory to allocate grids 1\n");
    return ;
  }
  //set all host machine okgrid to zero in the catchment and 1 outside the catchment
  for (x = 0; x < fullsize; ++x) {
		  okgrid[x] = 0;
		  localOK[x] = 0;
		  if (data->dem[x] == -9999){
			  okgrid[x] = 1;
			  localOK[x] = 1;
		  }



		  //}
  }
  int loop = 0;
  int totalGP = 0;
  do {
	  gridprogress = processtheGrid(data, device, 1, percent, nrows, ncols, okgrid,  localOK, 16, 16, 1, &doneP) ;

	  if (gridprogress != 0)
         gridprogress += processtheGrid(data, device, 1, percent, nrows, ncols, okgrid,  localOK, 16, 16, 1, &doneP) ;

	//printf("gridProgress = %d\n", gridprogress);
	//printf("cuda error NOW2:%s\n", cudaGetErrorString(cudaGetLastError()));
    totalGP += gridprogress;
    printf("%d / %d (%f %%)\n", totalGP, data->activecells, (double) totalGP / data->activecells * 100);
    /* ---- end of second run ---- */
    loop ++;
  } while ((gridprogress > 0));// && (((double)totalGP / (nrows * ncols)) < 0.99));

#ifndef PRODUCTION_RUN
  printf("About to finish FA_MFD\n");
#endif

  //#DEBUG-CHECK
  int count = 0;
  for (int r = 0; r < nrows; r++) {
	  for (int c = 0; c < ncols; c++) {
		  data->flatmask[r * ncols +c] = okgrid[r * ncols + c]; // copy okgrid for export here
		  if (okgrid[r * ncols + c] != 1) {
			  //printf("Cell at [%d,%d] has not been computed!\n", r, c);
			  if (data->mask[r * ncols + c] !=0) count ++; // count not computed in grid
		  }
	  }
  }
  printf("Number of actual cells not computed = %d\n", count);
  fprintf(data->outlog, "Number of actual cells not computed = %d\n", count);
  fprintf(data->outlog, "FA: FA memcopy operation 0 :%s\n", cudaGetErrorString(cudaGetLastError()));
  checkCudaErrors( cudaMemcpy( data->fa, device->fa, fullsize * sizeof(double), cudaMemcpyDeviceToHost)) ;
  fprintf(data->outlog, "FA: FA memcopy operation 1 :%s\n", cudaGetErrorString(cudaGetLastError()));

  //data->flatfile = "okgrid.txt";
  //write_int(data, data->flatmask, data->flatfile);
  //writeGRIDtoFile(data, "okgrid.tif", 1, 4); // this is the flatmask data

  double FA_max;
  	int FAindex = 0;
  	double cpuFA_max = 0.0;

  	int firstiter = data->start_iter - data->runiniter;
  	if (iter == firstiter) // cpu calculation otherwise we cannot locate the outletcell index
  	{
  		for (int i = 0; i < nrows; i++) {
  			for (int j = 0; j < ncols; j++) {
  				if (data->fa[i * ncols + j] > cpuFA_max)
  					{
  					cpuFA_max = data->fa[i * ncols + j];
  					FAindex = i * ncols + j;
  					}
  			}
  		}
  		data->FA_max = cpuFA_max;
  		data->outletcellidx = FAindex; // this is the outlet cell which will be maintained throughout the simulation
  	} else // do it faster using GPU in all subsequent iterations

  	   {
  			int fullsize = data->mapInfo.width * data->mapInfo.height;
  			thrust::device_ptr<double> max_FA = thrust::device_pointer_cast(device->fa);
  			FA_max = thrust::reduce(max_FA, max_FA + fullsize, (double) 0, thrust::maximum<double>());
  			data->FA_max = FA_max;
  		}
  	fprintf(data->outlog, "FA: Maximum FA is  %.6f \n\n", data->FA_max);
  	fprintf(data->outlog, "FA: Outletcell index is  %d \n\n", data->outletcellidx);
#ifndef PRODUCTION_RUN
  	printf("Maximum FA is  %.6f s\n\n", data->FA_max);
#endif
  	//#DEBUG-CHECK
  	count = 0;
  	for (int i = 0; i < nrows; i++) {
  		for (int j = 0; j < ncols; j++) {
  				if (data->fa[i * ncols + j] < 0) {
  				count++;
  				}
  		}
  	}
  	fprintf(data->outlog, "FA: Bad value count (i.e. not in catchment(s) = %d\n", count);
  free(okgrid);
  free(localOK);
  cudaFree(device->slopetotal); // remove slopes grid from GPU
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  fprintf(data->outlog, "time to complete correctflowMFD algorithm %.6f s\n\n", time / 1000.0);
  fprintf(data->outlog, "FA: End of FA :%s\n", cudaGetErrorString(cudaGetLastError()));

  return;
}
