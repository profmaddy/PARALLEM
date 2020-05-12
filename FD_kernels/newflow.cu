#include "newflow.h"


__device__ int look(int client, int nghbr, int gridCols)
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


//  Before MFD can be called sinks and flats must have been routed. Should they be carved for a number of iterations using SFD?
//  Slopes in all directions need to be stored for use in accumulation partitioning
//  In MFD we must abandon the 0-7 directions in favour of the byte notation

// NW  N NE    32  64   128
//  W  *  E    16   0   1
// SW  S SE     8   4   2


//**************************************************
// ORIGINAL CODE directions

// NW  N NE    7  0  1
//  W  *  E    6  8  2
// SW  S SE    5  4  3


//called FlowDirs<<<dimGrid, dimBlock>>>(device->mask, device->flatmask, device->dem, device->Slopes, device->SFD, device->prop, device->fd, 1, csize, ncell_x, ncell_y, device->dx, device->dy, 0);

__global__ void FlowDirs(int *mask, int *flatmask, double *zs, double *slopes, int *SFD, double *prop, int *mfd, int flowtype, int cell_size, int ncell_x, int ncell_y, int *dx, int *dy, int last, int *sinkcounter, int *flatcounter)
{
  int cellx, celly;
  int icol, irow, dcell;
  int aspect, aspectsfd, newaspect;
  int aspectF, aspectsfdF;

  double thiscellht, targetcellht;
  int upslopecount;
  int flatcount;
  int flatneighbour;
  int slopeidx;

  double smax, stemp, slopetot;
  float dc;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM - nothing to do
    return;

  int self = irow * ncell_x + icol;
  slopeidx = 8 * self;

  // force all aspects to 8 i.e. no direction
  mfd[self] = -9999; // set to nodataValue
  SFD[self] = -9999; // set to nodataValue

  if (self >= ncell_x * ncell_y) { return;}
  if (mask[self] != 1) { return; } // return no data value to cells outside grid

  
  mfd[self] = 0; // set to zero default
  SFD[self] = 0; // set to zero default
     
  smax = 0.0;// needed for SFD aspect
  slopetot = 0.0; // for down slope directions and mfd only
  aspect = 0; // we will now code using MFD aspect coding -default aspect is none i.e. zero, this will signify sink or flat for routing
  aspectsfd = 0;
  upslopecount = 0;
  flatcount = 0;
  stemp = 0;

	  for (dcell = 0; dcell < 8; dcell ++)
	  {
		cellx = icol + dx [dcell];
		celly = irow + dy [dcell];
		if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y) // for each of my neighbours
		{
			//newaspect = 2^((dcell + 6) % 8); // this does not work on the GPU
			  if (dcell == 0)  newaspect = 64;
			  if (dcell == 1)  newaspect = 128;
			  if (dcell == 2)  newaspect = 1;
			  if (dcell == 3)  newaspect = 2;
			  if (dcell == 4)  newaspect = 4;
			  if (dcell == 5)  newaspect = 8;
			  if (dcell == 6)  newaspect = 16;
			  if (dcell == 7)  newaspect = 32;

		  if (dcell % 2 == 0)
			dc = 1.0;   //even directions are cardinal
		  else
			dc = 1.41;  //odd directions are diagonals

		  //calculate the slope
		  thiscellht = zs[self];
		  targetcellht = zs[celly * ncell_x + cellx];

		  stemp = (thiscellht - targetcellht) / (cell_size * dc);

		  if (targetcellht == -9999) {
			  stemp = -0.0000001 ;
		  	  //aspectsfd = newaspect;
		  	  //aspect = newaspect;
			  } // needed to identify and route sinks on boundary

		  //store the slope and aspect value if it flows here i.e. down hill
		  if ( (thiscellht > targetcellht) && (targetcellht != -9999) )
		  {
			  slopes[slopeidx+dcell] = stemp;
			  slopetot += slopes[slopeidx+dcell];
				  if (stemp > smax) // default aspect set to aspect of max slope for SFD
				  {
					  smax = stemp;
					  aspectsfd = newaspect;
					  //if (aspectsfd == 0) return;
				  }

			  aspect = aspect + newaspect; // for MFD
			  //if (aspect == 0) return;

		  }
		  if (thiscellht < targetcellht) upslopecount += 1 ;

		  if (thiscellht == targetcellht)  {
			  flatcount += 1;
	  	  	  aspectF = 0;
	  	  	  aspectsfdF = 0;
		  }

		  }
	  }

	  if (upslopecount >7) {  
              atomicAdd(sinkcounter,1);
	  	  	  aspect = 0;
	  	  	  aspectsfd = 0;
	  } // I got nowhere to go because I am a sink

	  if ((flatcount + upslopecount) > 7)
		{ if(upslopecount <=7){ // do not double count the sinks
								atomicAdd(flatcounter,1);
								aspect = 0;
								aspectsfd = 0;
							  }
	  } // I got nowhere to go because I am a flat next to either adjacent higher or equal height cells i.e. no lower cells

	  if ((flatcount < 8) && (flatcount + upslopecount > 7)){
          aspect = 0; // aspectF;
          aspectsfd = 0; // aspectsfdF;
	  }

	  mfd[self] = aspect;
  	  SFD[self] = aspectsfd;// set the SFD direction - this is used for all processes i.e. they are not MFD

}

// __global__ void flow_boundaryMFD( int *mask, double *zs, double *slopes, int *SFD, int *mfd, int ncell_x, int ncell_y);
__global__ void flow_boundaryMFD( int *mask, double *dem, double *slopes, int *SFD, int *mfd, int ncell_x, int ncell_y) //directions altered to correspond better with slim
{
    int irow, icol;
    irow = blockIdx.y * blockDim.y + threadIdx.y;
    icol = blockIdx.x * blockDim.x + threadIdx.x;

    if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM - nothing to do
      return;

    int self ;
    self = irow * ncell_x + icol;
    int slopeidx;
    slopeidx = 0;

    // if (mask[self] == 0) return; // go back if not in grid
     if (mfd[self] != 0) return ;  // only looking at non defined edge cells

     if (mask[self] != 1) return;

    // __syncthreads();

       //slopes[self] = 0.000001; // give the direction a slope

        int north, south, east, west, northeast, northwest, southeast, southwest;
        int gridCols;
        int lookhere;
        gridCols = ncell_x;

        east  = look(self, EAST,      gridCols);
        south  = look(self, SOUTH,     gridCols);
        west  = look(self, WEST,      gridCols);
        north  = look(self, NORTH,     gridCols);
        northeast = look(self, NORTHEAST, gridCols);
        northwest = look(self, NORTHWEST, gridCols);
        southeast = look(self, SOUTHEAST, gridCols);
        southwest = look(self, SOUTHWEST, gridCols);

        // NW  N NE    32  64   128
        //  W  *  E    16   0   1
        // SW  S SE     8   4   2

       SFD[self] = 0; mfd[self] = 0;
       slopeidx = self * 8;
       if (dem[east] == -9999.) { SFD[self] = 1; mfd[self] = 1; }; // I can go east and don't include me in future calculations
       if (dem[southeast] == -9999.) { SFD[self] = 2;mfd[self] = 2;  }; // I can go southeast
       if (dem[south] == -9999.) { SFD[self] = 4;mfd[self] = 4; };; // I can go south
       if (dem[southwest] == -9999.) { SFD[self] = 8;mfd[self] = 8; };; // go southwest
       if (dem[west] ==-9999.) { SFD[self] = 16;mfd[self] = 16; };; // go west
       if (dem[northwest] ==-9999.) { SFD[self] = 32;mfd[self] = 32; };; // go northwest
       if (dem[north] == -9999.)  { SFD[self] = 64;mfd[self] = 64; };; // go north
       if (dem[northeast] == -9999.) { SFD[self] = 128;mfd[self] = 128; }; // go northeast
}

//***************************************************
// Set up for working out shortest paths on plateaus
// If cell is not a plateau (aspect = 0) then set shortest path to 0 else set shortest path to large number

__global__ void shortest_paths_plateaus_initMFD(int *mask, double *zs, int *fd, int* SFD, float* shortest_paths, int cell_size,
				int ncell_x, int ncell_y, double *lowHeight)
{
	int irow, icol;
	irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;
	if(icol >= ncell_x || irow >= ncell_y)
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest


	// Set the lowHeight to 0 for each cell
	lowHeight[self] = 0.0;

	//lowHeight[self] = zs[self];
    
    shortest_paths[self] = 20000000;
    if (fd[self] > 0) shortest_paths[self] = 0;
    if (SFD[self] > 0) shortest_paths[self] = 0;
        
}

//**************************************************************
// Route water over plateaus
//**************************************************************
// find the cell around self which has same height and a lower distance to the exit
// route to that cell

__global__ void route_plateausMFD(int *mask, double *zs, double *slopes, double *props, int *fd, int* SFD, float* shortest_paths,
				int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag, double *lowHeight)
{
	int irow, icol, dcell, cellx, celly;
	int newaspect;
	float dc, dis;
    int index;

	irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;

  // If outside of DEM return (nothing to do)
	if(icol >= ncell_x || irow >= ncell_y)
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] !=1) return; // don't calculate if not in catchment(s) of interest

	//if (SFD[self] < 0) return;

	__syncthreads();

  // Get the current shortest path for this cell
	float min_distance = shortest_paths[self];

	// If it's zero then we're not part of a plateau so finish
	if(min_distance == 0)
		return;

	// double stemp; // use the same slope nomenclature;
    index = self * 8;

 
	// double this_ele = zs[irow * ncell_x + icol];
	for (dcell = 0; dcell < 8; dcell ++)
	{
        //slopes[index + dcell] = 0.0;
        //newaspect = 2^((dcell + 6) % 8);
		  if (dcell == 0)  newaspect = 64;
		  if (dcell == 1)  newaspect = 128;
		  if (dcell == 2)  newaspect = 1;
		  if (dcell == 3)  newaspect = 2;
		  if (dcell == 4)  newaspect = 4;
		  if (dcell == 5)  newaspect = 8;
		  if (dcell == 6)  newaspect = 16;
		  if (dcell == 7)  newaspect = 32;

		cellx = icol + dx [dcell];
		celly = irow + dy [dcell];
		//stemp = 0.0;
		if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y)
		{
			dc = 1.0;
			if (dcell % 2 == 0)
			{
				dc = 1.0;   //even directions
										//are cardinal
			}
			else
			{
			  dc = 1.41;  //odd directions are diagonals
			}
			dis = shortest_paths[celly * ncell_x + cellx] + dc;

			if (zs[self] == zs[celly * ncell_x + cellx] && min_distance > dis) // if I am a flat

              
            // the following is taken from Shuns dissertation
            //if ((((zs[self] == zs[celly * ncell_x + cellx]) ||
                //zs[self] - zs[celly * ncell_x + cellx] < 0.000001 && zs[self] > zs[celly * ncell_x + cellx]) || (zs[self] - zs[celly * ncell_x + cellx] > -0.000001 &&
                  //  zs[self] < zs[celly * ncell_x + cellx])) && min_distance > dis)



		    {
				min_distance = dis;
				*change_flag = 1;  // don't care how many cells have changed - just that cells have changed
				fd[self] = newaspect;
				SFD[self] = newaspect;
                slopes[index+dcell] = 0.0001;
                //props[self + dcell] = 1;
				shortest_paths[self] = min_distance;
				lowHeight[self] = lowHeight[celly * ncell_x + cellx];
			}
		}
	}
}

__global__ void init_watershedMFD(int *mask, int *SFD, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  //if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM nothing to do
      if (irow >= ncell_y) // if outside of DEM nothing to do
    return;

  int self = irow * ncell_x + icol;
  if (mask[self] != 1)
  {
    watershed_id[self] = 0;
    return;
  }

  //watershed_id[this_idx] = -1;
  //if((icol == 0) || (icol == ncell_x - 1) || (irow == 0) || (irow == ncell_y - 1)) // if we're on the outside of the DEM use a special watershed (ncell_x*ncell_y)
  if ((icol == ncell_x - 1) || (irow == 0) || (irow == ncell_y - 1)) // if we're on the outside of the DEM use a special watershed (ncell_x*ncell_y)
  {
    watershed_id[self] = 0;
	return;
  }

  if(SFD[self] == 0) // this is a sink as all plateaus are now routed (coded here for MFD)
    { watershed_id[self] = atomicAdd(counter,1) + 1;}
  else {
	// aspects in flooding MUST be SFD in old money to work!
	// find cell we flow into
	  //**************************************************
	  // ORIGINAL CODE directions
	  // NW  N NE    7  0  1
	  //  W  *  E    6  8  2
	  // SW  S SE    5  4  3

	  // NW  N NE    32  64   128
	  //  W  *  E    16   0   1
	  // SW  S SE     8   4   2

	 // int dir = ((int) log2(SFD[self]) + 2) % 8;

	  int dir = SFD[self];
      int newdir;
	  if (dir == 64)  newdir = 0;
	  if (dir == 128) newdir = 1;
	  if (dir == 1)   newdir = 2;
	  if (dir == 2)   newdir = 3;
	  if (dir == 4)   newdir = 4;
	  if (dir == 8)   newdir = 5;
	  if (dir == 16)  newdir = 6;
	  if (dir == 32)  newdir = 7;

	int cellx = icol + dx[newdir];
    int celly = irow + dy[newdir];
    // check we are still in the DEM
#ifndef PRODUCTION_RUN
    if (cellx < 0 || cellx >= ncell_x || celly < 0 || celly >= ncell_y) {
    	printf("HOW THE HECK DID WE GET HERE? \n\n");
      watershed_id[self] = ncell_x * ncell_y; // use special watershed
    }
    else
#endif
      // save this as a negative number so we don't confuse it with a watershed
      // each value is shifted down by 1 so that index 0 isn't used (it could be a
      // watershed in the DEM)
      watershed_id[self] = -(celly * ncell_x + cellx + 1);
  }
}

//**********************************************************
// Process Watersheds
//**********************************************************
// If any cell around me part of the same plateau has a higher
// watershed index use that one instead of mine.
// Q: Is this thread safe? What if two neigbouring cells update
// at the same time?

__global__ void init_watershed_sinkMFD(int *mask, int *SFD, int *watershed_id, double *zs, int *dx, int *dy, int ncell_x, int ncell_y, int *changed) {
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  //if(icol >= ncell_x || irow >= ncell_y) // Outside DEM - return as no work to do...
      if (irow >= ncell_y) // Outside DEM - return as no work to do...
    return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

    //if (watershed_id[self] == 0) return;

  if(SFD[self] == 0) { // if we're part of a flat - now coded for MFD
    int w_id = -1;
    for (int dcell = 0; dcell < 8; dcell ++) { // for each of the cells around me...
      int cellx = icol + dx [dcell];
      int celly = irow + dy [dcell];
      if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y) { // check we're still inside of the DEM
        int idx = celly * ncell_x + cellx;
        if(zs[idx] == zs[self] && watershed_id[idx] > watershed_id[self] ) {
        // if there's a watershed next to me (within the same plateau) with a higher index then use the index of that watershed
          w_id = watershed_id[idx];
          *changed = 1; // something has changed...
        }
      }
    }
    if(w_id != -1)
      watershed_id[self] = w_id; // I've got a new watershed value - so update...
  }
}

//*****************************************************************
//** Identify Watersheds
//*****************************************************************
// for each cell which doesn't have a watershed yet...
// assign watershed value for cell I flow into (if assigned)
// if we don't have a watershed value from this point at the
// cell that the cell I point at points to...
__global__ void identify_watershedMFD(int *mask, int *SFD, int *watershed_id, int *dx, int *dy, int ncell_x, int ncell_y, int *changed)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y)   // If outside of DEM nothing to do...
    return;

	int self = irow * ncell_x + icol;
    if (mask[self] != 1)  return; // don't calculate if not in catchment(s) of interest
    
    if (watershed_id[self] < 0)
    {
        //watershed_id[self] = watershed_id[-watershed_id[self] - 1];
        watershed_id[self] = watershed_id[-watershed_id[self] - 1];;
        *changed = 1;
    }

}


//**************************** EDITED TO HERE 23/02/20      *******************************************







//**************************************************************
// Compute the slope for cells in the plateau
//**************************************************************
__global__ void slope_plateausMFD(int *mask, double *zs, float* shortest_paths, int ncell_x, int ncell_y, double *lowHeight, double *slopes, int cell_size)
{
  int irow, icol; // dcell, cellx, celly
  //float dis; // dc
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  // If outside of DEM return (nothing to do)
  if(icol >= ncell_x || irow >= ncell_y)
    return;

  int self = irow * ncell_x + icol;
  if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

  // Get the current shortest path for this cell
  float min_distance = shortest_paths[irow * ncell_x + icol];

  // if we're not part of the plateau then nothing to do
  if (min_distance == 0)
    return;

  // Get the height at which this cell flows to (when it leaves the plateau)
  double low = lowHeight[irow * ncell_x + icol];

  // Get the height of this cell
  double height = zs[irow * ncell_x + icol];

  // The slope is (height - low) / (min_distance * cell_size)
  slopes[irow * ncell_x + icol] = (height - low) / ((min_distance * cell_size)*10000000); // make the slope small
}

//**************************************************************
//** Compute shortest paths in sinks
//**************************************************************
// Same as route water over plateaus except we don't change the aspect

__global__ void comp_shortest_paths_sinkMFD(int *mask, double *zs, int *aspects, float* shortest_paths,
				int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag)
{
	int irow, icol, dcell, cellx, celly;
	float dc, dis;
	irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;

	if(icol >= ncell_x || irow >= ncell_y) // if we're outside of the DEM - nothing to do so return
		return;

	float min_distance = shortest_paths[irow * ncell_x + icol];
	if(min_distance == 0)  // if we're not part of a flat
	  return;

	int flow_nei_idx = -1;
	int this_ele = zs[irow * ncell_x + icol];
	for (dcell = 0; dcell < 8; dcell ++)
	{
		cellx = icol + dx [dcell];
		celly = irow + dy [dcell];
		if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y)
		{
			dc = 1.0;
			if (dcell % 2 == 0)
			{
			  dc = 1.0;   //even directions
							      //are cardinal
			}
			else
			{
				dc = 1.41;  //odd directions
                    //are diagonals
			}
			dis = shortest_paths[celly * ncell_x + cellx] + dc;

            if (zs[this_ele] == zs[celly * ncell_x + cellx] && min_distance > dis) // if I am a flat
			//if((((this_ele == zs[celly * ncell_x + cellx]) ||
			  //this_ele - zs[celly * ncell_x + cellx] < 0.000001 && this_ele > zs[celly * ncell_x + cellx]) || (this_ele - zs[celly * ncell_x + cellx] > -0.000001 &&
				//this_ele < zs[celly * ncell_x + cellx])) && min_distance > dis)
				// Q: again - can this be replaced with:
				// if (abs(this_ele - zs[celly * ncell_x + cellx) < 0.000001 && min_distance > dis)
			{
				flow_nei_idx = dcell;
				min_distance = dis;
				*change_flag = 1;  // don't care which one has changed - just that there's been change
				shortest_paths[irow * ncell_x + icol] = min_distance;
			}
		}
	}
}


//********************************************
// Set what happens on the boundary of the DEM
//********************************************


// *****************************************
//            From floodingDriver
// *****************************************

using namespace std;

watershed_edge::watershed_edge()
{
	node1 = 0;
	node2 = 0;
	weight = 0.0;
}

watershed_edge::~watershed_edge()
{
}

watershed_edge::watershed_edge(int n1, int n2, HASH_TYPE w)
{
    node1 = n1;
    node2 = n2;
    weight = w;
}

int watershed_edge::get_node1()
{
    return node1;
}

int watershed_edge::get_node2()
{
    return node2;
}

HASH_TYPE watershed_edge::get_weight() const
{
    return weight;
}

void watershed_edge::set_node1(int n)
{
    node1 = n;
}

void watershed_edge::set_node2(int n)
{
    node2 = n;
}

void watershed_edge::set_weight(HASH_TYPE w)
{
    weight = w;
}

bool watershed_edge::operator<(const watershed_edge &other)
{
    return ((*this).weight < other.weight);
}

union_set::union_set(int s)
{
    size = s;
    parent = (int*) calloc(size, sizeof(int));
    rank = (int*) calloc(size, sizeof(int));
}

union_set::~union_set()
{
    free(parent);
    free(rank);
}

bool union_set::in_set(int x)
{
    if(x == 0)
        return true;
    if(x > 0 && x < size)
        if(parent[x] > 0)
            return true;
    return false;
}

int union_set::find_root(int x)
{
    if(parent[x] != x)
        return find_root(parent[x]);
    else return x;
}

void union_set::make_set(int x)
{
    if(x > 0 && x < size)
    {
        parent[x] = x;
        rank[x] = 0;
    }
}

void union_set::merge(int a, int b)
{
    int aroot = find_root(a);
    int broot = find_root(b);
    if(aroot == broot)
        return;
    if(rank[aroot] > rank[broot])
        parent[broot] = aroot;
    else
    {
        parent[aroot] = broot;
        if(rank[aroot] == rank[broot])
            rank[broot]++;
    }
}


__device__ double atomicMin(double* address, double val) {
	// HACK!!
	// As there is currently no atomic min for doubles use the following...
	double min = *address;
	double assumed;
	double old = *address;
    //double val = *ptr;
    do {       // If we have a value lower than that stored at the moment
		assumed = old;
		min = fmin(old, val);
		if (min == old)
			break;
        old = __longlong_as_double(atomicCAS((unsigned long long int*)address, __double_as_longlong(assumed), __double_as_longlong(min)));
    //    // When finished if it worked val will contain the old value and ptr will hold the new low value
    //	// if it fails - something else got in there first! - we get the new value held in prt back
    } while (assumed != old); // keep going until val <= min
    return old;
}



//****************************************************************
//** Initialise the hash table
//****************************************************************
//
__global__ void init_hash_tableMFD(HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol > list_size || irow > table_size) // If I'm not in the space for the hash table return...
    return;

  hash_table[(irow * list_size + icol ) * 3] = -1.0;     // set all values to -1
  hash_table[(irow * list_size + icol ) * 3 + 1] = -1.0;
  hash_table[(irow * list_size + icol ) * 3 + 2] = -1.0;
  table_counter[irow] = -1;
}

///#include "floodingdriver.h"
/// using namespace std;

void outputEdges(const char* filename, int table_size, int* table_counter, int* hash_table, int list_size) {
  FILE* out = fopen(filename, "w");
  for(int i = 0; i < table_size; i++)
  {
    if(table_counter[i] > 0)
      // go through each entry within the table
      for(int j = 0; j < (table_counter[i] + 1) * 3; j+=3)
      {
        fprintf(out, "Watershed 1: %d, Watershed 2: %d, Height: %d\n", hash_table[i * list_size * 3 + j], hash_table[i * list_size * 3 + j + 1], hash_table[i * list_size * 3 + j + 2]);
      }
  }
  fclose(out);
}

void outputEdges(const char* filename, int table_size, int* table_counter, double* hash_table, int list_size) {
  FILE* out = fopen(filename, "w");
  for(int i = 0; i < table_size; i++)
  {
    if(table_counter[i] > 0)
      // go through each entry within the table
      for(int j = 0; j < (table_counter[i] + 1) * 3; j+=3)
      {
        fprintf(out, "Watershed 1: %d, Watershed 2: %d, Height: %f\n", (int) hash_table[i * list_size * 3 + j], (int) hash_table[i * list_size * 3 + j + 1], hash_table[i * list_size * 3 + j + 2]);
      }
  }
  fclose(out);
}

//****************************************************************
//** Identify watershed boundaries
//****************************************************************
// If we're part of a boundary between two watersheds:
// Try to find boundary indexes in the hash table - if present ensure height for pair is minimum of old value and this value
// If not in hash table, add to hash table

__global__ void identify_watershed_boundary_write_backMFD(int *mask, int *watershed_id, double *zs, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, int *dx, int *dy, int ncell_x, int ncell_y, int *is_succeeded, int* counter)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y) // If I'm outside of the DEM nothing to do so return...
    return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return;

  int dcell, cellx, celly;
  int this_idx = irow * ncell_x + icol;
  int this_id = watershed_id[this_idx]; // my watershed index
  double height = -1.0;

#ifndef PRODUCTION_RUN
      // is this cell id too high?
      if (this_id > 10000) {
    	  printf("This cell index too high %d mask[%d]=%d, cellx=%d celly=%d\n", this_id, this_idx, mask[this_idx], icol, irow);
    	  //nei_id = nei_id;
      }
#endif

  //*is_succeeded = 0;

  for (dcell = 0; dcell < 8; dcell ++) // for each cell around me...
  {
    cellx = icol + dx [dcell];
    celly = irow + dy [dcell];
    if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y)
    {
      //atomicAdd(counter, 1);
      int idx = celly * ncell_x + cellx;
      int nei_id = watershed_id[idx];       // watershed index of neighbour
#ifndef PRODUCTION_RUN
      // is the neigbour cell part of the masked out area?
      if (mask[idx] != 1) {
    	  nei_id = -nei_id;
      }
      // is the neighbour id too high?
      if (nei_id > 10000) {
    	  printf("Neighbour cell index too high %d mask[%d]=%d, cellx=%d celly=%d\n", nei_id, idx, mask[idx], cellx, celly);
    	  //nei_id = nei_id;
      }
#endif
      if(this_id != nei_id) // if I'm part of a boundary between watersheds...
      {
        if(zs[this_idx] >= zs[idx]) // if I'm the higher cell of the boundary
          height = zs[this_idx];// store my height
        else
	  continue; // Ignore here - when neighbour is central cell it will be processed there...
        atomicAdd(counter, 1);
        int hash_val = (this_id + nei_id) % table_size;  // compute hash value for this boundary
        int found = 0;
        // try to find this hash value in those stored already...
        for(int i = 0; i < table_counter[hash_val] * 3; i += 3)
        {
          int label1 = (int) hash_table[hash_val * (list_size * 3) + i];
          int label2 = (int) hash_table[hash_val * (list_size * 3) + i + 1];
          //HASH_TYPE edge_height = (HASH_TYPE) hash_table[hash_val * (list_size * 3) + i + 2];
          if((label1 == this_id && label2 == nei_id) || (label1 == nei_id && label2 == this_id))
          {
            // if this boundary pair exists already...
            // make sure we have the lower height of the previous one and this one.
            atomicMin(hash_table + (hash_val * (list_size * 3) + i + 2), height);
            found = 1; // we've found this pair
          }
        }
        if(found == 0) // did we find this?
        {  // no
          // ask for a unique storage space...
          int list_idx = (atomicAdd(table_counter + hash_val, 1) + 1) * 3;
          if(list_idx < list_size * 3) // is there enough space to store another item in the hash?
          { // yes
            //list[list_idx] = this_id;
            //list[list_idx + 1] = nei_id;
            //list[list_idx + 2] = height;
            // Store the three peices of information about this watershed boundary
            hash_table[hash_val * (list_size * 3) + list_idx] = (HASH_TYPE) this_id;
            hash_table[hash_val * (list_size * 3) + list_idx + 1] = (HASH_TYPE) nei_id;
            hash_table[hash_val * (list_size * 3) + list_idx + 2] = (HASH_TYPE) height;
          }
          else

           *is_succeeded = 0; // mark that we couldn't store this data...
        }
      }
    }
  }
}

__global__ void raise_watershedMFD(int *mask, int *watershed_id, HASH_TYPE *raise, double *zs, int ncell_x, int ncell_y)
{
    int irow, icol;
    irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;
	if(icol >= ncell_x || irow >= ncell_y)
        return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

    int this_idx = irow * ncell_x + icol;
    int label = watershed_id[this_idx];
    if(zs[this_idx] < raise[label])
        zs[this_idx] = raise[label];
}

bool operator<(const watershed_edge& a, const watershed_edge &b)
{
    return (a.get_weight()) < (b.get_weight());
}

//*******************************************************************
//** Initialise and sort watershed boundary labels
//*******************************************************************
//
int init_and_sort_watershed_bound_labelsMFD(Data* data, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, watershed_edge *sorted_list, int sorted_list_size, int *max_watershed_p, int max_sink_id)
{
#ifndef PRODUCTION_RUN
    for (int ii = 0; ii < table_size; ii++) {
    	//printf("==========\n");
        for (int jj = 0; jj < table_counter[ii]; jj++) {
        	if ((int)hash_table[(ii*list_size+jj)*3] > max_sink_id || (int)hash_table[(ii*list_size+jj)*3+1] > max_sink_id)
            printf("[%d, %d, %f]  \n", (int)hash_table[(ii*list_size+jj)*3], (int)hash_table[(ii*list_size+jj)*3+1], hash_table[(ii*list_size+jj)*3+2]);
        	if ((int)hash_table[(ii*list_size+jj)*3] < 0 || (int)hash_table[(ii*list_size+jj)*3+1] < 0) {
        		printf("Mask watershed value used\n");
        		exit(0);
        	}
        }
       // if (table_counter[ii] >=0)
       // 	printf("\n\n");
    }
#endif

  //printf("Allocate memory for edge list.\n");
  int max_watershed = 0;
  int max_ws_index = -1;
  int current_idx = 0;
  //watershed_edge* memory_index_of_last;
  // store the labels of all watersheds that have been found

  //int* watershed_labels = (int*) malloc(sizeof(int) * sorted_list_size);
    int* watershed_labels = (int*) calloc (sorted_list_size, sizeof(int));

    fprintf(data->outlog,"FD: table_size = %d\n", table_size);
  //printf("Memory allocated for watersheds = %d\n", sorted_list_size);
  for(int i = 0; i < table_size; i++)
  {
    int last_start = current_idx;
    int edge_num = table_counter[i] + 1;
    //printf("i: %d\n", i);
    //printf("i = %d, Boundary to be exceeded? %d\n", i, current_idx + counter > sorted_list_size);
    for(int j = 0; j < edge_num * 3; j += 3)
    {
      //printf("i = %d j = %d\n", i, j);

	  int n1 =  (int) hash_table[i * list_size * 3 + j];  // added casting
      int n2 =  (int) hash_table[i * list_size * 3 + j + 1];

      HASH_TYPE w = hash_table[i * list_size * 3 + j + 2];

	  	  if (n1 < 0 )
	  {
		  printf("n1 < 0");
	  }



      // perform a merge sort to see if this value already exists in the list...
      int placed = 0;
      for (int pos = last_start; pos < current_idx; pos++) {
        if (n1 == sorted_list[pos].get_node1() && n2 == sorted_list[pos].get_node2() ||
            n2 == sorted_list[pos].get_node1() && n1 == sorted_list[pos].get_node2() )
		{
          // This pair already exists in the list so just check to see if we have a
          // lower value for weight
	      //printf("found (%d,%d,%d) at %d\n", n1, n2, w, pos);
          if (w < sorted_list[pos].get_weight()) {
            sorted_list[pos].set_weight(w);
          } // otherwise the value is larger. In eather case we're finished looking
          placed = 1;
          break;
		}
      }

      // if we didn't find the watershed pair in there already...
      if (placed == 0)
	  {
        sorted_list[current_idx].set_node1(n1);
        sorted_list[current_idx].set_node2(n2);
        sorted_list[current_idx].set_weight(w);
        //memory_index_of_last = &sorted_list[current_idx];

        current_idx++;

        // if placed == 1 then we've already stored these
        int found_n1 = 0;
        int found_n2 = 0;
        for (int pp = 0; pp < max_watershed; pp++)
		{
          if (watershed_labels[pp] == n1)
            found_n1 = 1;
          if (watershed_labels[pp] == n2)
            found_n2 = 1;
          if (found_n1 && found_n2)
            break;
		}
        if (found_n1 == 0)
		{
			watershed_labels[max_watershed] = n1;
			max_watershed++;
				if (max_ws_index < n1)
				max_ws_index = n1;
		}
		if (found_n2 == 0)
		{
			watershed_labels[max_watershed] = n2;
			max_watershed++;
			if (max_ws_index < n2)
			 max_ws_index = n2;
		}
      }
    }
  }

  *max_watershed_p = max_ws_index;
  fprintf(data->outlog,"FD: Sort edge list.\n");
  fprintf(data->outlog,"FD: Sorted_list_size: %d, Current_idx: %d\n", sorted_list_size, current_idx);
  //qsort(sorted_list, sorted_list_size, sizeof(watershed_edge*), compare_watershed_edge);
 // printf("Memory_index_of_last = %u\n", memory_index_of_last);
  fprintf(data->outlog,"FD: sorted_list + sorted_list_size = %u\n", sorted_list + sorted_list_size);
  fprintf(data->outlog,"FD: current idx = : %d \n", current_idx);

  sort(sorted_list, sorted_list + current_idx);

#ifndef PRODUCTION_RUN
  // check the sorted list
   for(int i = 0; i < sorted_list_size; i++)
   {
	  double height1 = sorted_list[i].get_weight();
	  int    label1a = sorted_list[i].get_node1();
      int    label2a = sorted_list[i].get_node2();

	  if (label1a < 0 )
	  {
		  fprintf(data->outlog,"FD: label1 < 0");
//		  _getch();
	  }
   }
#endif

   fprintf(data->outlog,"FD: List sorted.\n");
//  for (int i = 0; i < current_idx; i++) {
  //  printf("[%d] W1: %d W2: %d H: %f\n", i, sorted_list[i].get_node1(), sorted_list[i].get_node2(), sorted_list[i].get_weight());
// }
   fprintf(data->outlog,"FD: Number of watersheds = %d\n", max_watershed);
   data->sinkcount = max_watershed;

  //_getch();
  free(watershed_labels);
  return current_idx;
}



