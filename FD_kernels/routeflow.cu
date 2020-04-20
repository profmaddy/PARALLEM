#include "lem.h"
#include "routeflow.h"
#include "io.h"

double calcprops(Data* data) {

  double slopetot;
  int width = data->mapInfo.width;
  int height = data->mapInfo.height ;
  int fullsize = width * height;
  double propsum;
  int counter;
  int bigcounter;
  int zeroslope;
  counter = 0;
  bigcounter = 0;
  zeroslope = 0;
  double propdiff;

	  for (int i = 0; i < fullsize; i++)
	  {
			  slopetot = 0;
			  propsum = 0;
			  propdiff = 0;

			  if (data->mask[i] == 1){

					for(int k = 0; k < 8; k++)
					{
						slopetot += data->Slopes[i+k];
					}

                    if (slopetot == 0) zeroslope++;

					for(int k = 0; k < 8; k++)
					{
						data->prop[i+k] = data->Slopes[i+k] / slopetot;
						propsum += data->prop[i+k];
					}

                    if ((propsum > 1.0) || (propsum < 0.0)){
						propdiff = 1.0 - propsum;
					}

					if (propdiff != 0.0)
						{
							// printf("illegal propsum for cell %d diff is [%15.12f] \n", i, propdiff);
							counter++;
							if (propdiff > 0.01) bigcounter ++;
						}
			  }
	  }

  printf("zero slope totals  %d \n", zeroslope);
  printf("total number of non zero diff proportion entries = %d \n", counter);
  printf("total number where propdiff > 0.01  = %d \n", bigcounter);

  return(0);
}


double testoutletislowest(Data* data) {
 double lowest ;
  lowest = 1000.;
  //int cell;
  int width = data->mapInfo.width;
  int height = data->mapInfo.height ;
  int fullsize = width * height;

	for (int x = 0; x < fullsize; x++) {
		if (data->mask[x] == 1) {
			if (data->dem[x] < lowest) lowest = data->dem[x];}
	}
  return(lowest);
}

int createflatmask(Data* data)
{

int width = data->mapInfo.width;
int height = data->mapInfo.height ;
int fullsize = width * height;
#ifndef PRODUCTION_RUN
printf("no data value %f\n", data->mapInfo.nodata);
#endif
	for (int x = 0; x < fullsize; x++) {
		data->flatmask[x] = 1;
		if (data->dem[x] == data->mapInfo.nodata) data->flatmask[x] = 0;
	}
 return 1;
}

void floodingDriver(dim3 dimGrid, dim3 dimBlock, Data* data, Data* device, int ncell_x, int ncell_y, int cell_size, int iter) {

//#ifndef PRODUCTION_RUN
  printf("********************\n");
  printf("Computing watershed\n");
  printf("********************\n");
//#endif
  int width = data->mapInfo.width;
  int height = data->mapInfo.height ;
  int fullsize = width * height;

  fprintf(data->outlog,"FD: ************************\n");
  fprintf(data->outlog,"FD: * Computing watersheds *\n");
  fprintf(data->outlog,"FD: ************************\n");

  int block_ncell_y = 16;
  int block_ncell_x = 16;
  int *counter;
  int *counter_d;
  int count = 0;
  int *change_flag_d;
  int *change_flag_h;

  int *change_flag_d2;
  int *change_flag_h2;

  // ** Flats (0) can still exist - but these are now parts of sinks...so deal with these
  counter = (int*) malloc(sizeof(int)); //cleared
  *counter = 0;
  checkCudaErrors( cudaMalloc((void**)&counter_d, sizeof(int)) );
  checkCudaErrors( cudaMemcpy(counter_d, counter, sizeof(int), cudaMemcpyHostToDevice) );

  // Identify each sink cell and assign it a watershed value
  // This kernel identifies each cell which doesn't have a flow direction in SFD and assigns it a unique watershed value.
  // Hence the value of counting is actually the number of cells without a flow direction
  init_watershedMFD<<<dimGrid, dimBlock>>>(device->mask, device->SFD, device->watershed_id, device->dem, device->shortest_paths, device->dx, device->dy, counter_d, ncell_x, ncell_y);
  fprintf(data->outlog,"FD: init_watershed :%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("FD: init_watershed :%s\n", cudaGetErrorString(cudaGetLastError()));

#ifndef PRODUCTION_RUN
  // Copy the watershed_id matrix back to see if there are any mistakes in it
  // create memory to store this:
  int size = data->mapInfo.height * data->mapInfo.width;
  int* local_ws = (int*) malloc(size * sizeof(int));
  checkCudaErrors( cudaMemcpy(local_ws, device->watershed_id, size * sizeof(int), cudaMemcpyDeviceToHost) );
  int irow, icol;
  printf("Before run\n");
  for (irow = 0; irow < data->mapInfo.height; irow++) {
	  for (icol = 0; icol < data->mapInfo.width; icol++) {
		  int this_idx = irow * ncell_x + icol;
		  if (local_ws[this_idx] > 10000) { // || local_ws[this_idx] < -10000) {
			  printf("watershed1 created with invalid value [%d][%d]%d\n", icol, irow, local_ws[this_idx]);
		  }
		  //printf("%d, ", local_ws[this_idx]);
	  }
  }
  free(local_ws);
#endif

//#ifndef PRODUCTION_RUN
  checkCudaErrors( cudaMemcpy(counter, counter_d, sizeof(int), cudaMemcpyDeviceToHost));
  fprintf(data->outlog, "FD: Total number of cells which are being processed as potential sinks = %d\n", *counter);
  printf("FD: Total number of cells which are being processed as potential sinks = %d\n", *counter);
//#endif


  change_flag_h = (int*) malloc(sizeof(int));//cleared
  change_flag_h2 = (int*) malloc(sizeof(int));//cleared

  checkCudaErrors( cudaMalloc((void**) &change_flag_d, sizeof(int)) );
  checkCudaErrors( cudaMalloc((void**) &change_flag_d2, sizeof(int)) );

  // We now go through all the cells marked as watersheds and merge them together. This is done by saying any two
  // cells which are neighbours and are both marked as watersheds are in fact the same watershed. They both get allocated
  // the same watershed id - the larger of the two watershed values.
  // This is repeated until in a given iteration no cell has its watershed value changed - thus we have identified the
  // minimum number of watersheds
  do {
    *change_flag_h = 0; // clear flag
    checkCudaErrors( cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice) );

    //** Process Watersheds
    //**********************************************************************
    // Reduce number of watersheds - neighbouring cells in the same plateau
    // adopt the higher watershed values...

    init_watershed_sinkMFD<<<dimGrid, dimBlock>>>(device->mask, device->SFD, device->watershed_id, device->dem, device->dx, device->dy, ncell_x, ncell_y, change_flag_d);
    //if(count <3) fprintf(data->outlog,"FD: init_watershed_sink :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors(cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost));

	count ++;

	//printf("FD: Merging watersheds iterations %d\n", count);

  } while(*change_flag_h == 1);  // while we're still doing things...


  fprintf(data->outlog,"FD: Merged watersheds in %d iterations\n", count);
  printf("FD: Merged watersheds in %d iterations\n", count);

  fflush(data->outlog);

  checkCudaErrors( cudaMemcpy((void *)data->watershed_id, device->watershed_id, sizeof(int) * fullsize, cudaMemcpyDeviceToHost) );
  //writeGRIDtoFile(data, "ws_merged.tif", 1, 1);


  count = 0;
  // We have now defined all unique sinks, we now need to identify the full watershed for each sink - all cells which flow
  // into that sink. To do this we ???
  do {
    // Reset change flag
    *change_flag_h2 = 0;

    checkCudaErrors(cudaMemcpy(change_flag_d2, change_flag_h2, sizeof(int), cudaMemcpyHostToDevice));
    // *********************************************************************
    // ** identify watershed
    // *********************************************************************
    // Assign watershed values to cells that don't currently have values
    // assume same watershed value as cell I flow into...
    // Q: Does this allocate a watershed value to all cells in the DEM?

    identify_watershedMFD<<<dimGrid, dimBlock>>>(device->mask, device->SFD, device->watershed_id, device->dx, device->dy, ncell_x, ncell_y, change_flag_d2);
    //if(count <3) fprintf(data->outlog,"FD: identify_watersheds :%s\n", cudaGetErrorString(cudaGetLastError()));
    //printf("FD: identify_watersheds :%s\n", cudaGetErrorString(cudaGetLastError()));

    checkCudaErrors( cudaMemcpy(change_flag_h2, change_flag_d2, sizeof(int), cudaMemcpyDeviceToHost) );

    //printf("Reported last error = %s\n", cudaGetErrorString(cudaGetLastError()));

    count ++;

    //printf("FD: identify_watershedMFD loop iteration %d\n", count);

  } while(*change_flag_h2 == 1);

#ifndef PRODUCTION_RUN
  // Copy the watershed_id matrix back to see if there are any mistakes in it
  // create memory to store this:
  int size2 = data->mapInfo.height * data->mapInfo.width;
  int* local_ws2 = (int*) malloc(size2 * sizeof(int));
  checkCudaErrors( cudaMemcpy(local_ws2, device->watershed_id, size2 * sizeof(int), cudaMemcpyDeviceToHost) );
  int irow2, icol2;
  for (irow2 = 0; irow2 < data->mapInfo.height; irow2++)
	  for (icol2 = 0; icol2 < data->mapInfo.width; icol2++) {
		  int this_idx = irow2 * ncell_x + icol2;
		  if (local_ws[this_idx] > 10000 || local_ws[this_idx] < -10000) {
			  printf("watershed created with invalid value %d\n", local_ws[this_idx]);
		  }
	  }
  free(local_ws2);
#endif

  fprintf(data->outlog,"FD: Identified watersheds in %d iterations\n", count);
  printf("FD: Identified watersheds in %d iterations\n", count);


  // From here on in we're doing a hash table solution...
  // Is this a good idea...?

  fflush(data->outlog);



#ifndef PRODUCTION_RUN
  printf("Working on hash\n");
#endif

  // The values here need to be adjusted for each DEM size
  int table_size = 30000; //increased from 10000 14/11/17 was 32767;
  int list_size =800; //800
  //HASH_TYPE *hash_table;
  //+HASH_TYPE *hash_table_d;
  HASH_TYPE *hash_table;
  HASH_TYPE *hash_table_d;

  int *is_succeeded;
  int *is_succeeded_d2;
  int *table_counter;
  int *table_counter_d;

  fprintf(data->outlog, "FD: table size %d, list size %d", table_size, list_size);
  fprintf(data->outlog, "FD: Creating hash table on GPU of size %zd (= %f Mb)\n", table_size * list_size *3 * sizeof(HASH_TYPE), (table_size * list_size * sizeof(HASH_TYPE) * 3.0 / 1024 / 1024));
  //printf("HASH_TYPE = %lu   double = %lu\n", sizeof(HASH_TYPE), sizeof(double));

  hash_table = (HASH_TYPE*) calloc(table_size*list_size*3, sizeof(HASH_TYPE)); //cleared
  table_counter = (int*) malloc(sizeof(int) * table_size); //cleared

  checkCudaErrors(cudaMalloc((void **)&hash_table_d, sizeof(HASH_TYPE) * table_size * list_size * 3)); //checks added 30/01/16
  //checkCudaErrors(cudaMalloc((void **)&hash_table_d, 3)); //checks added 30/01/16
  checkCudaErrors(cudaMalloc((void **)&table_counter_d, sizeof(int) * table_size));//checks added 30/01/16

  dim3 dimTableGrid(list_size/block_ncell_x + 1, table_size/block_ncell_y + 1);

//#ifndef PRODUCTION_RUN
  printf("Working on hash 2\n");
//#endif

  // Set all values in hash table to -1
  //** Initialise hash table...
  // Set all values in hash table to -1
  init_hash_tableMFD<<<dimTableGrid, dimBlock>>>(hash_table_d, table_counter_d, table_size, list_size);
  fprintf(data->outlog,"FD: init_hash_table :%s\n", cudaGetErrorString(cudaGetLastError()));

#ifndef PRODUCTION_RUN
  printf("Working on hash 3\n");
#endif

  is_succeeded = (int*) malloc(sizeof(int)); //cleared
  *is_succeeded = 1;
  checkCudaErrors(cudaMalloc((void **)&is_succeeded_d2, sizeof(int)));//checks added 30/01/16
  checkCudaErrors(cudaMemcpy((void*) is_succeeded_d2, (void*) is_succeeded, sizeof(int), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpy((void*) is_succeeded, (void*) is_succeeded_d2, sizeof(int), cudaMemcpyDeviceToHost) );
  fprintf(data->outlog, "FD: Allocate memory for is_succeeded  = %s\n", cudaGetErrorString(cudaGetLastError()));
//  checkCudaErrors( cudaMemcpy((void*) is_succeeded_d, (void*) is_succeeded, sizeof(int), cudaMemcpyHostToDevice) );

  // Place boundary watershed indexes into hash table along with heights
  identify_watershed_boundary_write_backMFD<<<dimGrid, dimBlock>>>(device->mask, device->watershed_id, device->dem, hash_table_d, table_counter_d, table_size, list_size, device->dx, device->dy, ncell_x, ncell_y, is_succeeded_d2, counter_d);

#ifndef PRODUCTION_RUN
  printf("Working on hash 4\n");
#endif

  fprintf(data->outlog,"FD: identify_watershed_boundary_write_back  :%s\n", cudaGetErrorString(cudaGetLastError()));
  // Copy hash table back to main memory...
//  printf("Reported last error B= %s\n", cudaGetErrorString(cudaGetLastError()));

#ifndef PRODUCTION_RUN
  printf("FD: identify_watershed_boundary_write_back  :%s\n", cudaGetErrorString(cudaGetLastError()));
#endif

  checkCudaErrors( cudaMemcpy((void*) is_succeeded, (void*) is_succeeded_d2, sizeof(int), cudaMemcpyDeviceToHost) );

#ifndef PRODUCTION_RUN
  printf("Working on hash 4a\n");
#endif

  checkCudaErrors( cudaMemcpy((void*) counter, (void*) counter_d, sizeof(int), cudaMemcpyDeviceToHost) );

#ifndef PRODUCTION_RUN
  printf("Working on hash 4b\n");
#endif

  checkCudaErrors( cudaMemcpy((void*) table_counter, (void*) table_counter_d, sizeof(int) * table_size, cudaMemcpyDeviceToHost));

#ifndef PRODUCTION_RUN
  printf("Working on hash 4c\n");
#endif

  checkCudaErrors( cudaMemcpy((void *) hash_table, hash_table_d, sizeof(HASH_TYPE) * table_size * list_size * 3, cudaMemcpyDeviceToHost));
  fprintf(data->outlog,"FD: idevice hash memcopies :%s\n", cudaGetErrorString(cudaGetLastError()));

  //int edgeNum = 0;
  //
  //// go through each table entry... to count up how many pairs we have.
  //for(int i = 0; i < table_size; i++) {
  //  if(table_counter[i] > 0)
  //    // go through each entry within the table
  //    for(int j = 0; j < (table_counter[i] + 1) * 3; j+=3) {
  //      edgeNum++;
  //    }
  //}

#ifndef PRODUCTION_RUN
  printf("Working on hash 5\n");
#endif

  int sorted_list_size = 0;
  int max_list = -1;
  int max_watershed = -1;
  for(int i = 0; i < table_size; i++) {
	  if (max_list < table_counter[i])
		  max_list = table_counter[i];
	  //if (table_counter[i] >= 0) {
	  //printf("table_counter[%d]=%d\n", i, table_counter[i]);
	  //}
    sorted_list_size += table_counter[i] + 1;
  }
  // create new array just large enough to store the edges
  fprintf(data->outlog,"FD: amount of data held in hash table = %d [%d,%d]\n", sorted_list_size, table_size, max_list);

#ifndef PRODUCTION_RUN
  printf("Working on hash 6\n");
#endif

  watershed_edge *sorted_list = new watershed_edge[sorted_list_size];  // cleared

  sorted_list_size = init_and_sort_watershed_bound_labelsMFD(data, hash_table, table_counter, table_size, list_size, sorted_list, sorted_list_size, &max_watershed, *counter);

  fprintf(data->outlog,"FD: ##max_watershed %d \n", max_watershed);
  //int cccount = 0;
  //int nncount = 0;
  //for(int i = 0; i < sorted_list_size; i++) {
  //	  double height1 = sorted_list[i].get_weight();
  //	  int    label1a = sorted_list[i].get_node1();
  //    int    label2a = sorted_list[i].get_node2();
  //	  double raise_d1 =  height1;
  //	  if (raise_d1 == sorted_list[i].get_weight()) {  // we were an integer to start with
  //		  cccount ++;
  //	  }
  //	  else
  //		nncount ++;
  //}

  union_set *set = new union_set(max_watershed + 1); //cleared

  fprintf(data->outlog,"FD: max_watershed %d \n", max_watershed);
#ifndef PRODUCTION_RUN
  printf("FD: max_watershed %d \n", max_watershed);
#endif
  int *done;
  done = (int*) calloc(max_watershed + 1, sizeof(int)); //cleared

#ifndef PRODUCTION_RUN
  printf("Working on hash 7\n");
#endif

  HASH_TYPE *raise, *raise_d;
  raise = (HASH_TYPE*) calloc(max_watershed + 1, sizeof(HASH_TYPE)); //cleared

  cudaMalloc((void**)&raise_d, (max_watershed + 1) * sizeof(HASH_TYPE));
  cudaDeviceSynchronize();

  done[0] = 1;

  for(int i = 0; i < sorted_list_size; i++) {
    int label1, label2, root1, root2;
	HASH_TYPE height;

    label1 = sorted_list[i].get_node1();
    label2 = sorted_list[i].get_node2();
    height = (HASH_TYPE) sorted_list[i].get_weight();
    if(!(set->in_set(label1)))
      set->make_set(label1);
    if(!(set->in_set(label2)))
      set->make_set(label2);
    root1 = set->find_root(label1);
    root2 = set->find_root(label2);
    // If both labels are done, ignore the edge
    if(done[root1] && done[root2]) {
      continue;
    }
    // If only one label is done, assigne the other one to be done, and the other watershed is to be flooded
    if(done[root1] || done[root2]) {
      if(done[root1]) {
	done[root2] = 1;
	raise[root2] = height;
      }
      else {
	done[root1] = 1;
	raise[root1] = height;
      }

      continue;
    }
    //* If both labels are not done, merge them
    if(root1 != root2) {
      set->merge(root1, root2);
      raise[root2] = raise[root1] = height;
    }
  }
  for(int i = 0; i < max_watershed + 1; i++) {
    if(raise[i] != 0)
      {
	raise[i] = raise[set->find_root(i)];
      }
  }

#ifndef PRODUCTION_RUN
  printf("Working on hash 8\n");
#endif

  //int ccount = 0;
  //int ncount = 0;
  //for(int ii = 0; ii < max_watershed; ii++) {
  //	  int raise_i = (int) raise[ii];
  //	  double raise_d = (double) raise_i;
  //	  if (raise_d == raise[ii]) {  // we were an integer to start with
  //		  if (raise_i > 0)
  //		  ccount ++;
  //	  }
  //	  else
  //		ncount ++;
  //}

  cudaMemcpy((void*)raise_d, (void*) raise, (max_watershed + 1) * sizeof(HASH_TYPE), cudaMemcpyHostToDevice);

  raise_watershedMFD<<<dimGrid, dimBlock>>>(device->mask, device->watershed_id, raise_d, device->dem, ncell_x, ncell_y);
  fprintf(data->outlog,"FD: raise_watershed  :%s\n", cudaGetErrorString(cudaGetLastError()));

  cudaDeviceSynchronize();

  // THIS NEEDS CHECKING FOR MFD - add in flow type to call set here to 1 MFD


  // FlowDirs(int *mask, int *flatmask, double *zs, double *slopes, double *prop, int *aspects,int flowtype, int cell_size, int ncell_x, int ncell_y, int *dx, int *dy, int last);
  //FlowDirs<<<dimGrid, dimBlock>>>(device->mask, device->flatmask, device->dem, device->Slopes, device->SFD, device->prop, device->fd, 1, csize, ncell_x, ncell_y, device->dx, device->dy, 0, sinkcounter_d);
  fprintf(data->outlog,"FD: second single flow direction  :%s\n", cudaGetErrorString(cudaGetLastError()));

 // flow_boundary<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, ncell_x, ncell_y, device->dx, device->dy);
 // fprintf(data->outlog,"FD: second flow boundary  :%s\n", cudaGetErrorString(cudaGetLastError()));

  //printf("after flow boundary:%s\n", cudaGetErrorString(cudaGetLastError()));

  shortest_paths_plateaus_initMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->fd, device->SFD, device->shortest_paths, cell_size, ncell_x, ncell_y, device->lowHeight);
  fprintf(data->outlog,"FD: second shortest_paths_plateaus  :%s\n", cudaGetErrorString(cudaGetLastError()));


  do {

    *change_flag_h = 0;
    cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice);
    route_plateausMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->Slopes, device->prop, device->fd, device->SFD, device->shortest_paths,  ncell_x, ncell_y, device->dx, device->dy, change_flag_d, device->lowHeight);
    cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost);
  } while(*change_flag_h == 1);
  fprintf(data->outlog,"FD: after route_plateaus:%s\n", cudaGetErrorString(cudaGetLastError()));

  // Now work out the slopes on cells in a plateaux
  slope_plateausMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->shortest_paths, ncell_x, ncell_y, device->lowHeight, device->SlopePtr, cell_size);
  fprintf(data->outlog,"FD: slope pleateaus  :%s\n", cudaGetErrorString(cudaGetLastError()));


//  Clean up memory usage for flooding driver

  free(counter);
  free(change_flag_h);
  free(change_flag_h2);
  free(is_succeeded);
  free(hash_table);
  free(table_counter);
  free(raise);
  free(done);

  delete [] sorted_list;
  delete  set;

  cudaFree(counter_d);
  cudaFree(change_flag_d);
  cudaFree(change_flag_d2);
  cudaFree(raise_d);
  cudaFree(hash_table_d);
  cudaFree(table_counter_d);
  cudaFree(is_succeeded_d2);

  fprintf(data->outlog,"FD: *********************************\n");
  fprintf(data->outlog,"FD: * Finished Computing watersheds *\n");
  fprintf(data->outlog,"FD: *********************************\n");
}

//******************************************************
//**
//** Main Code for SFD calculation
//**
//****
//** This code sets up and calls all kernels
//**
//******************************************************
void cuFlowDirection(Data* data, Data* device, int iter)

{
#ifndef PRODUCTION_RUN
	printf("cuFlowDirection\n");
#endif

  fprintf(data->outlog, "FD: starting cuFlowDirection :%s\n", cudaGetErrorString(cudaGetLastError()));

  cudaEvent_t start, stop;
  float time;
  size_t freenow, total;

 // start the timer for SFD;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  int xmove [9] = {0,1,1,1,0,-1,-1,-1,0};
  int ymove [9] = {-1,-1,0,1,1,1,0,-1,0};
  int *dx;
  int *dy;
  dx = &xmove[0];
  dy = &ymove[0];

  data->dx = dx;
  data->dy = dy;

  int block_ncell_y = 16;
  int block_ncell_x = 16;
  int fullsize;
  int ncell_x = data->mapInfo.width;
  int ncell_y = data->mapInfo.height;
  int csize = (int) data->mapInfo.cellsize;

  dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
  dim3 dimBlock(block_ncell_x, block_ncell_y);
  fullsize= ncell_x * ncell_y;

#ifndef PRODUCTION_RUN
  printf("#########################\n");
  printf("Max x value = %d\n", (ncell_x/block_ncell_x+1)*block_ncell_x);
  printf("Max y value = %d\n", (ncell_y/block_ncell_y+1)*block_ncell_y);
  printf("#########################\n");
#endif

  // declare transient memory to be freed on exit
  data->shortest_paths = (float*) malloc(fullsize * sizeof(float));
  data->fdmod = (int *) malloc(sizeof(int) * data->mapInfo.height * data->mapInfo.width);
  // *** Set up memory space to hold watershed id's
  data->watershed_id = (int*) malloc(ncell_x * ncell_y * sizeof(int));
  fprintf(data->outlog, "FD: dev mem allocate shortest_path, fdmod, watershed_id :%s\n", cudaGetErrorString(cudaGetLastError()));


  checkCudaErrors(cudaMemcpy((void *) device->dx, data->dx, 9 * sizeof(int), cudaMemcpyHostToDevice) );
  checkCudaErrors(cudaMemcpy((void *) device->dy, data->dy, 9 * sizeof(int), cudaMemcpyHostToDevice) );
  checkCudaErrors(cudaMemcpy((void *) device->dem, data->dem, fullsize * sizeof(double), cudaMemcpyHostToDevice) );
  checkCudaErrors(cudaMemcpy((void *) device->SlopePtr, data->SlopePtr, fullsize * sizeof(double), cudaMemcpyHostToDevice) ); //new 30/01/16
  checkCudaErrors(cudaMemcpy((void *) device->Slopes, data->Slopes, 8* fullsize * sizeof(double), cudaMemcpyHostToDevice) );
  checkCudaErrors(cudaMemcpy((void *) device->prop, data->prop, 8* fullsize * sizeof(double), cudaMemcpyHostToDevice) );

  fprintf(data->outlog, "FD: memcopy dx, dy, dem, slopeptr, Slopes and prop :%s\n", cudaGetErrorString(cudaGetLastError()));

  //check lowest point in DEM is the outlet cell and not in the grid
  double lowest ;
  lowest = testoutletislowest(data);
  if (lowest < data->dem[data->outletcellidx]) data->dem[data->outletcellidx] = lowest - 0.0001; //make sure outlet cell is the lowest for flooding routine

//#ifndef PRODUCTION_RUN
   printf("Lowest DEM value = %f, outletcell index %d, height at outletcell = %f\n",  lowest, data->outletcellidx, data->dem[data->outletcellidx] );
//#endif


// copy the zeroed fd matrix to the GPU (all set for MFD i.e. 0 in setprocessmatices )
   checkCudaErrors(cudaMemcpy((void *)device->fd, data->fd, sizeof(int) * fullsize, cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy((void *)device->SFD, data->SFD, sizeof(int) * fullsize, cudaMemcpyHostToDevice));
   fprintf(data->outlog, "FD: fd and SFD to GPU :%s\n", cudaGetErrorString(cudaGetLastError()));

   int *sinkcounter_h, *sinkcounter_d;
   int *flatcounter_h, *flatcounter_d;

   sinkcounter_h = (int *)malloc(sizeof(int));
   checkCudaErrors( cudaMalloc((void **)&sinkcounter_d, sizeof(int)) );

   flatcounter_h = (int *)malloc(sizeof(int));
   checkCudaErrors( cudaMalloc((void **)&flatcounter_d, sizeof(int)) );

/*   for (int i= 0; i<2; i++)
   { */
	   *sinkcounter_h = 0;
	   *flatcounter_h = 0;
	   checkCudaErrors( cudaMemcpy(sinkcounter_d, sinkcounter_h, sizeof(int), cudaMemcpyHostToDevice) );
	   checkCudaErrors( cudaMemcpy(flatcounter_d, flatcounter_h, sizeof(int), cudaMemcpyHostToDevice) );

	   //calculate the flow directions for both SFD and MFD using MFD coding - n.b. new type parameter set here to 1 for MFD
	   FlowDirs<<<dimGrid, dimBlock>>>(device->mask, device->flatmask, device->dem, device->Slopes, device->SFD, device->prop, device->fd, 1, csize, ncell_x, ncell_y, device->dx, device->dy, 0, sinkcounter_d, flatcounter_d);

	   fprintf(data->outlog,"FD: first flow routing :%s\n", cudaGetErrorString(cudaGetLastError()));
       fflush(data->outlog);
	   checkCudaErrors( cudaMemcpy(sinkcounter_h, sinkcounter_d, sizeof(int) ,cudaMemcpyDeviceToHost) );
	   checkCudaErrors( cudaMemcpy(flatcounter_h, flatcounter_d, sizeof(int) ,cudaMemcpyDeviceToHost) );
	   printf("FD: initial sinks :%d \n", *sinkcounter_h);
	   printf("FD: initial flats :%d \n", *flatcounter_h);
	   fprintf(data->outlog,"FD: initial sinks :%d \n", *sinkcounter_h);
	   fprintf(data->outlog,"FD: initial flats :%d \n", *flatcounter_h);

   //now take care of the boundary cells with unallocated flow direction

  flow_boundaryMFD<<<dimGrid, dimBlock>>>( device->mask, device->dem, device->Slopes, device->SFD, device->fd, ncell_x, ncell_y);

  // printf("After flow_boundary :%s\n", cudaGetErrorString(cudaGetLastError()));

  //copy back new matrices
  //checkCudaErrors( cudaMemcpy((void *)data->flatmask, device->flatmask, sizeof(int) * fullsize, cudaMemcpyDeviceToHost) );
  checkCudaErrors( cudaMemcpy((void *)data->SFD, device->SFD, sizeof(int) * fullsize, cudaMemcpyDeviceToHost) );
  checkCudaErrors(cudaMemcpy((void *) data->Slopes, device->Slopes, 8 * fullsize * sizeof(double), cudaMemcpyDeviceToHost) );
  checkCudaErrors( cudaMemcpy((void *)data->fd, device->fd, sizeof(int) * fullsize, cudaMemcpyDeviceToHost) );

   // checkCudaErrors(cudaMemcpy((void *) data->prop, device->prop, 8* fullsize * sizeof(double), cudaMemcpyDeviceToHost) );

// #ifndef PRODUCTION_RUN
  printf("Checking FD  %d\n", data->mapInfo.width * data->mapInfo.height);
  // check the values coming back
   int flatandsinkcount ;
   int negingrid;
   int cellsincatchment  ;
   flatandsinkcount = 0;
   cellsincatchment = 0;
   negingrid = 0;
	// Validate Flow Directions and mapout cells which need to be ignored for later MFD

   for (int j = 0; j < fullsize; j++){
		if ((data->mask[j] == 1))
			  {
				 cellsincatchment ++ ;
					if (data->SFD[j] == 0) flatandsinkcount ++;
					if (data->SFD[j] == -1) negingrid ++; // no longer used
			  }
	}

	printf("Edge Flat and Sink count in catchment: %d of %d cells in catchment\n", flatandsinkcount, cellsincatchment );
	printf("Negative FD in catchment: %d \n", negingrid);
	fprintf(data->outlog,"Edge Flat and Sink count in catchment: %d of %d cells in catchment\n", flatandsinkcount, cellsincatchment );
	fprintf(data->outlog,"Negative FD in catchment: %d \n", negingrid);
    fflush(data->outlog);

// #endif

	//data->FDfile = "fdwithzeros.txt";
	//write_int(data, data->fd, data->FDfile);
	//writeGRIDtoFile(data, "init_fd.tif", 1, 0);

	shortest_paths_plateaus_initMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->fd, device->SFD, device->shortest_paths, csize, data->mapInfo.width, data->mapInfo.height, device->lowHeight);

	checkCudaErrors( cudaMemcpy((void *)data->shortest_paths, device->shortest_paths, sizeof(float) * fullsize, cudaMemcpyDeviceToHost) );

    fprintf(data->outlog, "FD: return from shortest_paths_plateaus :%s\n", cudaGetErrorString(cudaGetLastError()));


	int numsfdpath;
	numsfdpath= 0;

	for (int j = 0; j < fullsize; j++){
        if ((data->mask[j] == 1))
        {
            if (data->shortest_paths[j] > 0) {
                numsfdpath++;
                //printf("%f, ", data->shortest_paths[j]);
            }
        }
	}
	printf("Number of shortest_paths: %d \n", numsfdpath );
	printf("FD: first shortest_path_plateaus_init :%s\n", cudaGetErrorString(cudaGetLastError()));
	fprintf(data->outlog,"FD: first shortest_path_plateaus_init :%s\n", cudaGetErrorString(cudaGetLastError()));

	//write_float(data, data->shortest_paths, "sp.txt");

  // Set up a flag to indicate if we have had any change since the last iteration?
  int *change_flag_h, *change_flag_d;
  change_flag_h = (int *)malloc(sizeof(int));
  checkCudaErrors( cudaMalloc((void **)&change_flag_d, sizeof(int)) );
  int count = 0;

  do
  {
    // Clear flag
    *change_flag_h = 0;
    checkCudaErrors( cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice) );

    // Each cell looks at the cells around itself finds the one with the same
    // height as itself and shortest path to exit and flows towards that cell
    // When we adopt the shortest route we adopt its lowHeight value too
    route_plateausMFD<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->Slopes, device->prop, device->fd, device->SFD, device->shortest_paths,  data->mapInfo.width, data->mapInfo.height, device->dx, device->dy, change_flag_d, device->lowHeight);
    if (count <3) fprintf(data->outlog,"FD: first route_plateaus_init :%s\n", cudaGetErrorString(cudaGetLastError()));

    checkCudaErrors( cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost) );  // copy back flag

    //printf("Any change in plateaus? %d\n", *change_flag_h);
    count ++;

    //if(*change_flag_h == 1)
  } while(*change_flag_h == 1); // while at least some cells have been modified, repeat

  fprintf(data->outlog,"FD: route plateaus :%s\n", cudaGetErrorString(cudaGetLastError()));
  fprintf(data->outlog,"FD: Routed plateaus in %d iterations. Now starting flooding driver... \n", count);

  printf("FD: route plateaus in %d iterations :%s\n", count, cudaGetErrorString(cudaGetLastError()));

  //copy back new matrices
  checkCudaErrors( cudaMemcpy((void *)data->fd, device->fd, sizeof(int) * fullsize, cudaMemcpyDeviceToHost) );
  checkCudaErrors( cudaMemcpy((void *)data->SFD, device->SFD, sizeof(int) * fullsize, cudaMemcpyDeviceToHost) );
  checkCudaErrors( cudaMemcpy((void *)data->shortest_paths, device->shortest_paths, sizeof(float) * fullsize, cudaMemcpyDeviceToHost) );


  //write_float(data, data->shortest_paths, "sp2.txt");


  //checkCudaErrors(cudaMemcpy((void *) data->Slopes, device->Slopes, 8* fullsize * sizeof(double), cudaMemcpyDeviceToHost) );

// #ifndef PRODUCTION_RUN
    printf("Checking FD after route plateau  %d\n", data->mapInfo.width * data->mapInfo.height);
    // check the values coming back
     int sinkcount ;
     int mfdsinkcount;
     sinkcount = 0;
     mfdsinkcount = 0;
  	// Validate Flow Directions and mapout cells which need to be ignored for later MFD

     for (int j = 0; j < fullsize; j++){
  		if ((data->mask[j] == 1))
  			  {
  					if (data->SFD[j] == 0) sinkcount ++;
  					if (data->fd[j] == 0) mfdsinkcount ++;
  			  }
  	}
     printf("SFD sinkcount : %d , MFD sinkcount : %d \n", sinkcount, mfdsinkcount);

  	if (sinkcount > *sinkcounter_h) {
  		printf ("Something went wrong, sinks remaining = %d but sinks in initial DEM are %d\n", sinkcount, *sinkcounter_h);
  		//writeGRIDtoFile(data, "fd_sinksleft.tif", 1, 0);
  		//exit(0);
  	}
 // #endif

  	calcprops(data); // not needed until later?

    fflush(data->outlog);

	//data->FDfile = "fdwithsinks.txt";
	//write_int(data, data->fd, data->FDfile);

    if (sinkcount > 1){
  	 //floodingDriver(dimGrid, dimBlock, data, device, data->mapInfo.width, data->mapInfo.height, csize /*data->mapInfo.cellsize*/, iter);
    	fprintf(data->outlog,"FD: errors after return from flooding driver :%s\n", cudaGetErrorString(cudaGetLastError())); }

  	//checkCudaErrors(cudaMemcpy(progress_d, temp, sizeof(unsigned int), cudaMemcpyHostToDevice) );cudaMemcpy(progress_d, temp, sizeof(unsigned int), cudaMemcpyHostToDevice) );
  	checkCudaErrors(cudaMemcpy(data->SlopePtr, device->SlopePtr, data->mapInfo.height * data->mapInfo.width * sizeof(double), cudaMemcpyDeviceToHost));
#ifndef PRODUCTION_RUN
  	printf(":sfd copy error %s\n", cudaGetErrorString(cudaGetLastError()));
#endif


 // bring back altered data
	checkCudaErrors( cudaMemcpy(data->fd, device->fd, data->mapInfo.height * data->mapInfo.width * sizeof(int), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy(data->SlopePtr, device->SlopePtr, data->mapInfo.height * data->mapInfo.width * sizeof(int), cudaMemcpyDeviceToHost) );

	 printf("outlet cell direction = %d\n", data->fd[data->outletcellidx])  ; // check the outlet cell direction was set by flow_boundary routine

#ifndef PRODUCTION_RUN
    int nonvaluecount = 0;
	int gridColumns = data->mapInfo.width;
	 for (int i = 0; i < fullsize; i++) {

		 if ( (data->fd[i] <= 0)  && (data->mask[i] == 1)){
			 	 nonvaluecount++ ;
			 //	printf("index of cell with zero fd = %d\n", i);
			}
	 }
	 printf("Total number of FD = 0: %d \n", nonvaluecount);
#endif

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(change_flag_d);

  //cudaFree(device->flatmask); // clear the original flatmap

  free(data->shortest_paths);
  free(data->fdmod);
  free(change_flag_h);

#ifndef PRODUCTION_RUN
  printf("Time to complete flow routing routine algorithm %.6f s\n", time / 1000.0);
#endif
  fprintf(data->outlog,"FD: time to complete flow routing routine algorithm %.6f s\n", time / 1000.0);

  cudaMemGetInfo(&freenow, &total);
  fprintf(data->outlog,"FD: Memory on CUDA card free at end of routing: %zd total: %zd\n",freenow/1024,total/1024);
}
