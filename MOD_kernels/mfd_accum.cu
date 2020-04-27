#include "mfd_accum.h"
//#include "Directions.h"
//#include "Data.h"
//#include "config.h"
//#include "lem.h"

__device__ int inghbr(int client, int nghbr, int gridCols)
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

__global__ void MFDsedbud(int *mask, double *hv, int *fd, double* props, double *diffuse, double *conc, double *geli,  double *depo, int *ok,
						  unsigned int *progressd2,  int* localOK, double ddk, double dck, double dgk,
						  double * finesPtr, double *stonePtr, double *nutPtr, double *soilTPtr , int gridCols, int gridRows)
{
	  int irow, icol;
	  irow = blockIdx.y * blockDim.y + threadIdx.y;
	  icol = blockIdx.x * blockDim.x + threadIdx.x;
	  // Are we outside the DEM?
	  if(icol >= gridCols || irow >= gridRows) // if outside of DEM nothing to do
	    return;
      if(icol == 0 || irow == 0 ) return; // if outside of DEM nothing to do

	int self;
    self = irow * gridCols + icol;
	if (mask[self] == 0) return;

    double accum_d, accum_c, accum_g;
	unsigned int theval = 99999999;
    int   nie, nise, nis, nisw, niw, ninw, nin, nine;
    double nce_from, nce_nutfrom;
    double ncse_from, ncse_nutfrom;
    double ncs_from, ncs_nutfrom;
    double ncsw_from, ncsw_nutfrom;
    double ncw_from, ncw_nutfrom;
    double ncnw_from, ncnw_nutfrom;
    double ncn_from, ncn_nutfrom;
    double ncne_from, ncne_nutfrom;
    double  fromEprop,  fromSEprop,  fromSprop,  fromSWprop,  fromWprop , fromNWprop ,  fromNprop ,  fromNEprop;
    int cnfd;
    int proploc;
    
    accum_d = 0.0;
    accum_c = 0.0;
    accum_g = 0.0;
        
    nce_from = nce_nutfrom = 0.0;
    ncse_from = ncse_nutfrom =  0.0;
    ncs_from = ncs_nutfrom = 0.0;
    ncsw_from =ncsw_nutfrom = 0.0;
    ncw_from = ncw_nutfrom = 0.0;
    ncnw_from = ncnw_nutfrom = 0.0;
    ncn_from = ncn_nutfrom = 0.0;
    ncne_from = ncne_nutfrom = 0.0;

    if (ok[self] != 0) return;

	if (localOK[self] == 1) {
      ok[self] = 1;
      return;
    }

    nie  = inghbr(self, EAST,      gridCols);
    nise = inghbr(self, SOUTHEAST, gridCols);
    nis  = inghbr(self, SOUTH,     gridCols);
    nisw = inghbr(self, SOUTHWEST, gridCols);
    niw  = inghbr(self, WEST,      gridCols);
    ninw = inghbr(self, NORTHWEST, gridCols);
    nin  = inghbr(self, NORTH,     gridCols);
    nine = inghbr(self, NORTHEAST, gridCols);

    /* ---- */
    cnfd = fd[nie];  // cnfd is current flow direction of contributing cell to east
    
    if (cnfd & WEST)
    {
        if (!ok[nie]) return;
            /* if (!ok[nie]) return;
            proploc = (nie * 8) + 6;// flowing west
            addhere = fa[nie] * props[proploc];
            accum += addhere;*/
                proploc = (nie * 8) + 6;// flowing west
                accum_d  += props[proploc]  * (diffuse[nie] * (1-ddk));
		        accum_c  += props[proploc] * (conc[nie]    * (1-dck)) ;
		        accum_g  += props[proploc] * (geli[nie]    * (1-dgk));
		        nce_from       = props[proploc] * stonePtr[nie] ;
		        nce_nutfrom    = props[proploc] * nutPtr[nie] ;
    }

    cnfd = fd[nise];
    
    if (cnfd & NORTHWEST)
    {
        if (!ok[nise]) return;
            /*if (!ok[nise]) return;
            proploc = (nise * 8) + 7;// flowing northwest
            addhere = fa[nise] * props[proploc];
            accum += addhere;*/
                proploc = (nise * 8) + 7;// flowing northwest
        	    accum_d  += props[proploc] * (diffuse[nise] * (1-ddk));
        	    accum_c  += props[proploc] * (conc[nise]     * (1-dck)) ;
        	    accum_g  += props[proploc] * (geli[nise]       * (1-dgk));
        	    ncse_from       = props[proploc] * stonePtr[nise] ;
        	    ncse_nutfrom    = props[proploc] * nutPtr[nise] ;
    }

    cnfd = fd[nis];
    if (cnfd & NORTH)
    {
        if (!ok[nis]) return;
        /*if (!ok[nis]) return;
        proploc = (nis * 8) + 0;//flowing north
        addhere = fa[nis] * props[proploc];
        accum += addhere;*/
                proploc = (nis * 8) + 0;//flowing north
        	    accum_d  += props[proploc] * (diffuse[nis] * (1-ddk));
        	    accum_c  += props[proploc] * (conc[nis]    * (1-dck)) ;
        	    accum_g  += props[proploc] * (geli[nis]      * (1-dgk));
        	    ncs_from       = props[proploc] * stonePtr[nis] ;
        	    ncs_nutfrom    = props[proploc] * nutPtr[nis] ;
    }

    cnfd = fd[nisw];
    if (cnfd & NORTHEAST)
    {
        if (!ok[nisw]) return;
        /*if (!ok[nisw]) return;
        proploc = (nisw * 8) + 1;//flowing northeast
        addhere = fa[nisw] * props[proploc];
        accum += addhere;*/
                proploc = (nisw * 8) + 1;//flowing northeast
                accum_d  += props[proploc] *  (diffuse[nisw] * (1-ddk));
                accum_c  += props[proploc] * (conc[nisw]     * (1-dck)) ;
                accum_g  += props[proploc] * (geli[nisw]       * (1-dgk));
	            ncsw_from       = props[proploc] * stonePtr[nisw] ;
	            ncsw_nutfrom    = props[proploc] * nutPtr[nisw] ;
    }

    cnfd = fd[niw];
    if (cnfd & EAST)
    {
        if (!ok[niw]) return;
        /*if (!ok[niw]) return;
        proploc = (niw * 8) + 2; // flowing east
        addhere = fa[niw] * props[proploc];
        accum += addhere;*/
                proploc = (niw * 8) + 2; // flowing east
                accum_d  += props[proploc] * (diffuse[niw] * (1-ddk));
                accum_c  += props[proploc] * (conc[niw]     * (1-dck)) ;
                accum_g  += props[proploc] * (geli[niw]       * (1-dgk));
                ncw_from       = props[proploc] * stonePtr[niw] ;
                ncw_nutfrom    = props[proploc] * nutPtr[niw] ;
    }

    cnfd = fd[ninw];
    if (cnfd & SOUTHEAST)
    {
        if (!ok[ninw]) return;
        /*if (!ok[ninw]) return;
        proploc = (ninw * 8) + 3;//flowing southeast
        addhere = fa[ninw] * props[proploc];
        accum += addhere;*/
                proploc = (ninw * 8) + 3;//flowing southeast
        	    accum_d  += props[proploc] * (diffuse[ninw] * (1-ddk));
        	    accum_c  += props[proploc] * (conc[ninw]     * (1-dck)) ;
        	    accum_g  += props[proploc] * (geli[ninw]       * (1-dgk));
        	    ncnw_from       = props[proploc] * stonePtr[ninw] ;
        	    ncnw_nutfrom    = props[proploc] * nutPtr[ninw] ;
    }

    cnfd = fd[nin];
    if (cnfd & SOUTH)
    {
        if (!ok[nin]) return;
        /*if (!ok[nin]) return;
        proploc = (nin * 8) + 4;//flowing south
        addhere = fa[nin] * props[proploc];
        accum += addhere;*/
                proploc = (nin * 8) + 4;//flowing south
        	    accum_d  += props[proploc] * (diffuse[nin] * (1-ddk));
        	    accum_c  += props[proploc] * (conc[nin]     * (1-dck)) ;
        	    accum_g  += props[proploc] * (geli[nin]       * (1-dgk));
        	    ncn_from       = props[proploc] * stonePtr[nin] ;
        	    ncn_nutfrom    = props[proploc] * nutPtr[nin] ;
    }

    cnfd = fd[nine];
    if (cnfd & SOUTHWEST)
    {
        if (!ok[nine]) return;
       /* if (!ok[nine]) return;
        proploc = (nine * 8) + 5;//flowing southwest
        addhere = fa[nine] * props[proploc];
        accum += addhere;*/
                proploc = (nine * 8) + 5;//flowing southwest
        	    accum_d  += props[proploc] * (diffuse[nine] * (1-ddk));
        	    accum_c  += props[proploc] * (conc[nine]     * (1-dck)) ;
        	    accum_g  += props[proploc] * (geli[nine]       * (1-dgk));
        	    ncne_from       = props[proploc] * stonePtr[nine] ;
        	    ncne_nutfrom    = props[proploc] * nutPtr[nine] ;
    }

    depo[self] =  ((accum_d + diffuse[self])* ddk) + ((accum_c + conc[self]) * dck) + ((accum_g + geli[self]) * dgk);

    // replace original deposit routine to calculate new proportions of fines, stones and nutrients

    double sedthick = (soilTPtr[self]) + (depo[self]) ;
	if (sedthick != 0.0)
	{
		double old = (soilTPtr[self]) / sedthick ;

		double input = 1-old;
		//double input = (depo[self]) /weight;

		double propStonefrom = ( nce_from + ncse_from + ncs_from + ncsw_from + ncw_from + ncnw_from + ncn_from + ncne_from) ;/// ndirections ;
        double propNutfrom = ( nce_nutfrom + ncse_nutfrom + ncs_nutfrom + ncsw_nutfrom + ncw_nutfrom + ncnw_nutfrom + ncn_nutfrom + ncne_nutfrom) ;

        stonePtr[self] =   (old * (stonePtr[self]))   + (input * (propStonefrom)) ;

		if (stonePtr[self] < 0.01) stonePtr[self] = 0.01;
		else if (stonePtr[self] > 99.9) stonePtr[self] = 99.9;

		finesPtr[self] = 100 - stonePtr[self];  // now same as calculation in weathering
		//finesPtr[self] = old * (finesPtr[self]) + input * (finesPtr[from]);

		//nutPtr[self] =   (old * (nutPtr[self]))   + (input * propNutfrom) ;
	}
	double toterosion = diffuse[self]+conc[self]+geli[self];
	if (toterosion >= sedthick) soilTPtr[self]= 0.0;
	if (sedthick > toterosion) soilTPtr[self] = sedthick - toterosion;

    localOK[self] = 1;
    atomicInc(progressd2, theval);
}

int processtheGrid(Data* data, Data* device, int loopMax, int percent, int gridRows, int gridColumns, int* okGrid,  int* localOK, int blockRows, int blockColumns, int dimBlock3, int* doneP,
													double ddk, double dck, double  dgk)

{
  int loopForever = (loopMax < 0) ? 1 : 0;

  int *ok;
  int *localOK_d;


  unsigned int progressh2, *progressd2;
  int gridProgress;

  //allocate GPU memory
    checkCudaErrors(cudaMalloc((void **) &ok,                   gridRows * gridColumns * sizeof(int))    );
    checkCudaErrors(cudaMalloc((void **) &localOK_d,     gridRows * gridColumns * sizeof(int))   );
    checkCudaErrors(cudaMalloc((void **) &progressd2,      sizeof(progressh2)) );
    fprintf(data->outlog, "MOD: sedmfd memory allocation :%s\n", cudaGetErrorString(cudaGetLastError()));

  gridProgress = 0;
  int loop = 0;

  int grid1 = data->mapInfo.width  / (blockColumns ) + 1;
  int grid2 = data->mapInfo.height / (blockRows ) + 1;

  //printf("blockColumns = %d  blockRows = %d  grid1 = %d  grid2 = %d \n", blockColumns, blockRows, grid1, grid2);

  do {
	  	  dim3 dimGrid(grid1, grid2, 1);
	  	  dim3 dimBlock(blockColumns, blockRows, dimBlock3);
  //printf("Grid is %d by %d by %d\n", dimGrid1, dimGrid2, dimGrid3);
	      int oneZero = 0;
	      int oneZero2 = 0;

	      //copy grids to GPU
	        checkCudaErrors(cudaMemcpy(ok,                    okGrid,    gridRows * gridColumns * sizeof(int),    cudaMemcpyHostToDevice) );
	        checkCudaErrors(cudaMemcpy(localOK_d,      localOK,  gridRows * gridColumns * sizeof(int),    cudaMemcpyHostToDevice) );

	        do {
	        		progressh2 = 0;
	        		oneZero2 = oneZero;

	        		checkCudaErrors(cudaMemcpy(progressd2, &progressh2, sizeof(progressh2), cudaMemcpyHostToDevice) );

//	        		MFDsedbud(double *hv, int *fd, double *diffuse, double *conc, double *geli,  int *ok, unsigned int *progressd,  int* localOK,
// 						double ddk, double dck, double dgk, double * finesPtr, double *stonePtr, double *nutPtr, double *soilTPtr , int gridCols)

	        		MFDsedbud<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->fd, device->prop, device->eroPtr, device->inciPtr, device->geliPtr, device->depoPtr,  ok, progressd2, localOK_d,
	        				                                                    ddk, dck, dgk, device->finesPtr, device->stonePtr, device->nutPtr, device->soilTPtr, gridColumns, gridRows);

	        		//fprintf(data->outlog, "MOD: MFDsedbud rtn :%s\n", cudaGetErrorString(cudaGetLastError()));
	        		//printf("MOD: MFDsedbud rtn :%s\n", cudaGetErrorString(cudaGetLastError()));
	        		checkCudaErrors(cudaMemcpy(&progressh2, progressd2, sizeof(progressh2), cudaMemcpyDeviceToHost) );
        			//fprintf(data->outlog, "MOD: MFDsedbud loop :%s\n", cudaGetErrorString(cudaGetLastError()));
	        		//printf("MOD: MFDsedbud loop :%s\n", cudaGetErrorString(cudaGetLastError()));

        			gridProgress += progressh2;

	        			if (progressh2 == 0)
	        					oneZero = 1;
	        					else {
	        							oneZero = 0;
	        							oneZero2 = 0;
	        					}

	        }  while ((progressh2 > 10 || !oneZero2)); // && ((double) *doneP * 100 / (gridRows * gridColumns)) < 99.0);

	        checkCudaErrors(cudaMemcpy(progressd2, &progressh2, sizeof(progressh2), cudaMemcpyHostToDevice) );

	        fprintf(data->outlog, "MOD:  loop :%s\n", cudaGetErrorString(cudaGetLastError()));

	        loop++;

} while (loop < loopMax || (loopForever && gridProgress) );

    checkCudaErrors(cudaMemcpy(okGrid, ok, gridRows * gridColumns * sizeof(int),    cudaMemcpyDeviceToHost) ) ;
    //checkCudaErrors(cudaMemcpy(data->fa, device->fa, gridRows * gridColumns * sizeof(double), cudaMemcpyDeviceToHost) ); //chnged from int !!!!!!!!!

  /* Free the GPU copies */
  cudaFree(ok);
  cudaFree(localOK_d);
  cudaFree(progressd2);
  //free(&progressh);

  //printf("sum of gridProgress = %d\n", gridProgress);
  return gridProgress; // this does not seem to be copied back in this code??????
}


void sedmfdaccum(Data* data, Data* device)
{
  int ncols = data->mapInfo.width;
  int nrows = data->mapInfo.height;
  int fullsize = nrows * ncols;

  if (cudaSuccess != cudaSetDevice(CUDA_DEVICE))
  {
    printf("Unable to access CUDA card\n");
    return ;
  }

#ifndef PRODUCTION_RUN
  printf("\nGPU Card set for correctmfdflow\n\n");
#endif

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

 // start the timer for correctflow;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  okgrid      = (int *)    malloc(fullsize * sizeof(int));
  localOK     = (int *)    malloc(fullsize * sizeof(int));

  if (okgrid == NULL || localOK == NULL ) {
    printf("Not enough memory to allocate grids 1\n");
    return ;
  }

  //set all host machine grids to zero
  for (x = 0; x < fullsize; ++x) {
    okgrid[x] = 0;
    localOK[x] = 0;
    if (data->mask[x] == 1) data->depoPtr[x] = 0.00001; //?????
  }

  cudaMemcpy(device->depoPtr, data->depoPtr, fullsize * sizeof(double),   cudaMemcpyHostToDevice);
#ifndef PRODUCTION_RUN
   printf("cuda error NOW1:%s\n", cudaGetErrorString(cudaGetLastError()));
#endif
   //  transport distance exponents
   //  double ddk = 0.95 * (23. / cell_size2); // travel distance parameter for diffuse flow -- normalized relative to 23 m plot
   //  double dck = 0.05 * (100. / cell_size2); // travel distance parameter for concentrated flow -- normalized relative to 100 m
   //  is the aim is to export ~5% of diffuse and 95% of concentrated after 100m?

    //probably not needed in this context??
  	if (data->ddk > 0.95)	data->ddk = 0.95;
  	if (data->dck > 0.9)	data->dck = 0.9;


  int loop = 0;

  int totalGP = 0;
  do {
    /* ---- beginning of first run ---- */
	//gridprogress = processWithPartitions(1, percent, 1568, 2704, GRDROWS, GRDCOLS, okgrid, fagrid, fdgrid, hvgrid, localOK, 169, 98, 1, BLOCKROWS, BLOCKCOLS, 1, &doneP);
	//printf(" arguments for gridprocess; Abr: %d, Abc: %d, GRDROWS: %d, GRDCOLS: %d, okgrid: %d, nrows: %d, ncols: %d \n",Abr, Abc, GRDROWS, GRDCOLS, okgrid, nrows, ncols);

	  //int processtheGrid(Data* data, Data* device, int loopMax, int percent, int gridRows, int gridColumns, int* okGrid,  int* localOK, int blockRows, int blockColumns, int dimBlock3, int* doneP)

	  gridprogress = processtheGrid(data, device, 1, percent, nrows, ncols, okgrid,  localOK, 16, 16, 1, &doneP, data->ddk, data->dck, data->dgk ) ;

	  /*printf("gridProgress = %d\n", gridprogress);*/
    /* ---- gap between runs ---- */
	//printf("cuda error NOW1:%s\n", cudaGetErrorString(cudaGetLastError()));
	//printf("half way there \n");

	if (gridprogress != 0)
	// gridprogress += processWithPartitions(1, percent,  672, 5408, GRDROWS, GRDCOLS, okgrid, fagrid, fdgrid, localOK, 338, 42, 1, BLOCKROWS, BLOCKCOLS, 1, &doneP);
	gridprogress += processtheGrid(data, device, 1, percent, nrows, ncols, okgrid,  localOK, 16, 16, 1, &doneP, data->ddk, data->dck, data->dgk) ;
   // printf("gridProgress = %d\n", gridprogress);
	//printf("cuda error NOW2:%s\n", cudaGetErrorString(cudaGetLastError()));
    totalGP += gridprogress;
    //printf("%d / %d (%f %%)\n", totalGP, nrows * ncols, (double) totalGP / (nrows * ncols) * 100);
    /* ---- end of second run ---- */
    loop ++;
  } while (gridprogress > 0); // && (double) totalGP / (GRDCOLS * GRDROWS)  < 0.99);


  cudaMemcpy(data->depoPtr, device->depoPtr, fullsize * sizeof(double),   cudaMemcpyDeviceToHost);
#ifndef PRODUCTION_RUN
  printf("cuda error NOW3:%s\n", cudaGetErrorString(cudaGetLastError()));
 // printf("Left to get = %d\n", nrows * ncols - totalGP);
  printf("About to finish sediment accumulation \n");
  //fflush(stdout);
#endif

  int count = 0;
  for (int r = 0; r < nrows; r++) {
	  for (int c = 0; c < ncols; c++) {
		  if (okgrid[r * ncols + c] != 1) {
			 // printf("Cell at [%d,%d] has not been computed!\n", r, c);
			  count ++;
		  }
	  }
  }
  fprintf(data->outlog, "Number of actual cells not computed = %d\n", count);


	thrust::device_ptr<double> deptot_d = thrust::device_pointer_cast(device->depoPtr);
	cudaSetDevice(0);
	data->totD = thrust::reduce(deptot_d, deptot_d + fullsize, (double) 0);
	fprintf(data->outlog, "total Dep from thrust is %10.8lf \n", data->totD);

	if (data->totD == NAN)
	{
		printf("data from totD is NaN \n");
		exit(0);
	}

	  free(okgrid);
	  free(localOK);


  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  fprintf(data->outlog, "time to complete sediment accumulation algorithm %.6f s\n\n", time / 1000.0);

  return;
}
