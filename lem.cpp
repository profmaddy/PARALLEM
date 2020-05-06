
//#include "lem.h" // this contains headers for all standard and add-in libraries
#include "lem.h"
#include "Data.h"
#include "io.h"
#include "memory.h"
#include "memory_dev.h"
#include "routeflow.h"
#include "newflow.h"
#include "mfd_simple.h"
#include "runoffweight.h"
#include "erosion.h"

int main(int argc, char* argv[]) {

	clock_t startrun;
	clock_t endrun;
	double Ptime;
	size_t freenow, total;
	errno_t err;

	Data data;
	Data device;

	int modelruntype;
	const char* pars = argv[1];
	if (!pars) pars = "default.par"; // use the default values
	createfilenamespace(&data);
	readinpar(&data, pars);
	data.demfile = "30mfamA.txt";
	//data.demfile = "pal25_bi20mod.asc";
	checkparread(&data);

	// output file for  data used in debugging
	//data.outlog = fopen(data.logfileOut, "w");
	if ((err = fopen_s(&data.outlog, data.logfileOut, "w")) != 0) fprintf(stderr, "cannot open file '%s': %s\n", data.logfileOut , strerror(err));

	fprintf(data.outlog,"Landscape file   =  %s\nClimate file     =  %s \nBedrock file   =  %s\n", data.demfile, data.clim_file, data.bedrockfile);
	fprintf(data.outlog, "logfile: %s\nMax_iterations:%d. \n", data.logfileOut, data.max_iterations);
	fprintf(data.outlog, "Summary file data in %s \n", data.outfilename);
	cudaMemGetInfo(&freenow, &total);
	fprintf(data.outlog,"Memory on CUDA card free at start of iteration: %zd total: %zd\n",freenow/1024,total/1024);
	
	startrun = clock();

	ReadClimateDataFromFile(&data); // read in appropriate number of data pairs, temp and rain and allocate sufficient memory

	readgrid(&data); // will read the ascii file it replaced ..
	//data.demfile = "testascii.txt";
	//write_double(&data, data.dem, data.demfile);

	createmask(&data);

	init_setoutletcell(&data);// set oulet cell for entire simulation i.e. nothing can be lower.

	//readgdalBedrockfromFile(&data, &device);
	createSoilTfromformula(&data);

	copyMask(&data, &device); // copy mask to device only needs to be done once

	SetClimateMatrices(&data, -10000) ; // setup the climate gradients for ppt and temp
	copylastclimate(&data, &device);// this copies the last climate matrices and frees host memory

	createProcessMatrices(&data); // allocate sufficient memory on host
	setProcessMatrices(&data);  // initialise values for attributes on host based upon parameter file inputs
	
	if (data.restart >0) {
		///retriveProcessMatrices(&data)	; uses GDAL
	} // call will also reset runiniter tp 0

	/// Perform the iterations
    printf("Performing the simulation\n");
    fprintf(data.outlog, "Performing the simulation\n");

    int start_pos_sim = data.start_iter - data.runiniter;
    int last_pos_sim = data.max_iterations;
    int first_iter = start_pos_sim;

	fflush(data.outlog);

    /************* START OF MAIN LOOP *********************/
	//for (int i = start_pos_sim; i < last_pos_sim; i++) {
		
	createDeviceSpace(&data, &device);
	modelruntype = 0;
		
    for (int i = -10; i < 0; i++) {
	
	printf("Iteration %d of %d :\n", i, last_pos_sim);

	fprintf(data.outlog,"\n***********************\n");
	fprintf(data.outlog,"Starting mode %s Iteration %d of %d\n", data.modelcode, i, last_pos_sim);
	fprintf(data.outlog,"***********************\n");

		if (i >= 0 ){
			SetClimateMatrices(&data,i);
		}
	
		setdevicespace_FD(&data, &device);
		cuFlowDirection(&data, &device, i); // calculate flow direction (MFD)
		cleardevicespace_FD(&data, &device);

	data.FDfile = "fd.txt";
	if (i == -1) write_int(&data, data.fd, data.FDfile);

		setdevicespace_FA(&data, &device);  // load matrices for runoffweight calculation
		computeRunOffWeights(&data, &device); // calculate runoff
		calcwater(&data);
		calcprops(&data); 
		correctmfdflow(&data, &device, i); // calculate flow accumulation (MFD)
		cleardevicespace_FA(&data, &device);

	data.FAfile = "fa.txt";
	if (i == -1) write_double(&data, data.fa, data.FAfile);

	    setdevicespace_Process(&data, &device);
	     erosionGPU(&data, &device, i);
	    cleardevicespace_Process(&data, &device);

		data.start_iter = -10;
		writeSummaryDataToFile(&data,  i); // note summary output is every iteration
		zerogrids(&data);

	//if (modelruntype ==0){copytolastclimate(&data, &device);} //copy climate matrices to lastclimte before erasing

		cudaMemGetInfo(&freenow, &total);
		fprintf(data.outlog,"Memory on CUDA card free at end of iteration: %zd total: %zd\n",freenow/1024,total/1024);

		fprintf(data.outlog,"Finished Iteration %d\n", i);
		fflush(data.outlog);

		//write_double(&data, data.dem, "demfile.asc");



	} // end of iterations ****************************************************

	clearDeviceSpace(&data, &device);
	free(data.dem);

	endrun = clock();
	Ptime = (double) (endrun-startrun)/ CLOCKS_PER_SEC ;
	fprintf(data.outlog, "Total simulation time of : %20.17lf \n\n", Ptime);
	printf("Total simulation time of : %20.17lf seconds \n\n", Ptime);

	cudaFree(device.mask);
	cudaFree(device.Flow_C);

	deleteProcessMatrices(&data);  //delete memory on host

	fflush(data.outlog);
	fclose(data.outlog);


}
