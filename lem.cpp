
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
	createDeviceSpace(&data, &device);

	/// Perform the iterations
    printf("Performing the simulation\n");
    fprintf(data.outlog, "Performing the simulation\n");
	fflush(data.outlog);

    int start_pos_sim = data.start_iter - data.runiniter;
    int last_pos_sim = data.max_iterations;

	data.start_iter = -1000;

	/************* START OF MAIN LOOP *********************/
	//for (int i = start_pos_sim; i < last_pos_sim; i++) {
    for (int i = data.start_iter; i < 0; i++) {
	
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

		setdevicespace_FA(&data, &device);  // load matrices for runoffweight calculation
		 computeRunOffWeights(&data, &device); // calculate runoff
		 calcwater(&data);
		 calcpropsSFD(&data, &device);
		 //calcprops(&data); 
		 correctmfdflow(&data, &device, i); // calculate flow accumulation (MFD)
		cleardevicespace_FA(&data, &device);

	    setdevicespace_Process(&data, &device);
	     erosionGPU(&data, &device, i);
		cleardevicespace_Process(&data, &device);
		
		writeSummaryDataToFile(&data,  i); // note summary output is every iteration
		zerogrids(&data);

	    //copytolastclimate(&data, &device); //copy climate matrices to lastclimte before erasing

		cudaMemGetInfo(&freenow, &total);
		fprintf(data.outlog,"Memory on CUDA card free at end of iteration: %zd total: %zd\n",freenow/1024,total/1024);
		printf("Memory on CUDA card free at end of iteration: %zd total: %zd\n", freenow / 1024, total / 1024);
		fprintf(data.outlog,"Finished Iteration %d\n", i);
		fflush(data.outlog);

		writegrids(&data, i);

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
