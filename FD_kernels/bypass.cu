#include "bypass.h"



void transferFAtoFAgrid(Data* data)
{
	  int cell;
	  int counter;
	  counter = 0;

	  // what if the lowest cell is currently in the grid!
	  for (int i = 0; i < data->mapInfo.height; i++){
			for (int j = 0; j < data->mapInfo.width; j++) {
				cell = 	i*data->mapInfo.width + j;

				if (data->mask[cell] !=0)
					{	counter ++;
						data->fa[cell]= (data->fdmod[cell] * 10000);
						data->SlopePtr[cell] = 1.3; // we will need slopes as well.
						if (counter <5) printf("value of fa at cell %d is %f\n", cell, data->fa[cell]);
					}
					else { (data->fa[cell]= 0);
						   (data->SlopePtr[cell] = 0); }
			}
	   }

	  data->dem[data->outletcellidx] = data->dem[data->outletcellidx] - 0.01 ;
	  return;
}

void bypassrouting(Data* data, Data* device, int iter)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int fullsize = ncell_x * ncell_y;

	data->Burnfile = "Burn/burnfd.asc";
	data->FAfile = "Burn/terraMFD30m.tif";
	data->fdmod = (int *) malloc(sizeof(int) * data->mapInfo.height * data->mapInfo.width); // use to read FA
    //readBurnFDfromFile(data); // read the FA
    transferFAtoFAgrid(data); //transfer readFA to data->fa

    checkCudaErrors(cudaMemcpy((void *) device->dem, data->dem, fullsize * sizeof(double), cudaMemcpyHostToDevice) );
    checkCudaErrors(cudaMemcpy((void *) device->SlopePtr, data->SlopePtr, fullsize * sizeof(double), cudaMemcpyHostToDevice) ); //new 30/01/16

	return;

}
