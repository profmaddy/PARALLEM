// Simple Process Kernels for calculating diffuse and concentrated erosion, returning them in erosion and incision matrices

#include "eroincidep.h"

__global__ void diff_erosion(int *mask, double *ePtr_d, double *FAPtr_d, double *SlopePtr_d, double *finesPtr_d, double *nutPtr_d,  double *soilTPtr_d, double *stonePtr_d, double *TotBPtr_d,
				                                     int soil_struct, int profile_permeability, int ncell_x, int ncell_y, double Qthreshold)
{
	/*!
	  Diffuse erosion by unconcentrated overland flow
	  erodibility is the k factor in the USLE.  Need to check whether a zero value is likely
	  The cell flow is the FA corrected for cell size (100. / cell_size) * (100. / cell_size). It is not clear whether this is really appropriate
	  This is based on a modified Musgrave approach. Calculation is based upon Mitchell and Bubenzer, 1980, Zhang et al., 2002
	  See Wainwright 2005 p.106
	  OF is runoff depth in mm and Slope is surface slope in m per m. according to Wainwright
	*/
	int icol = blockIdx.x * blockDim.x + threadIdx.x; //column
	int irow = blockIdx.y * blockDim.y + threadIdx.y; //row

	// simple out of bounds check
	if( icol >= ncell_x  || irow >= ncell_y )
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] == 0)
	{
		ePtr_d[self] = 0.0;
		return; // don't calculate if not in catchment(s) of interest
		}

	double erodibility = (2.77e-6 * pow (finesPtr_d[self], 1.14) * (12. - nutPtr_d[self])) + (0.043 * (soil_struct - 2.)) + (0.033 * (profile_permeability - 3.));
	// both soil_strict and profile_permeability are unused i.e. both are zero at all times

	if (erodibility < 0.0) erodibility = 0.0;
	double OF = ((FAPtr_d[self]) * 1000) / 25 ; // now in mm;

	if (OF > Qthreshold)            // limit to when concentrated erosion sets in
		OF = Qthreshold;

	double slope = (SlopePtr_d[self]);

	//slope constraint placed back on 3/10/18
	if (slope > 5.0) slope = 5.0;

	double answer = ((  erodibility * OF * OF * pow (slope, 1.67 ) *
						(pow ((110. - (stonePtr_d[self])), 0.877) / 61.7023))) *
								exp (-0.07 * (TotBPtr_d[self]) / 76.5) ;


	if (answer > (soilTPtr_d[self])) answer = soilTPtr_d[self];
	ePtr_d[self] = answer;

}

__global__ void conc_erosion( int *mask, double *inciPtr, double *FAPtr, double *SlopePtr, double *TotBPtr,
		                           double *soilTPtr, int ncell_x,int ncell_y, double cellsize, double flow_thresh,
		                                 double flow_A, double flow_B, double* flow_C, double maxincision )
{
	/*!
	 * Concentrated erosion by channelized flow [mm]. Simuluated using a modified stream-power approach
	 * see Eq [3] of Wainwright 1995 p. 106
	 * the concentrated erosion parameter k2 set at 0.00005 seems low - need to check
	 */
	int icol = blockIdx.x * blockDim.x + threadIdx.x; //column
	int irow = blockIdx.y * blockDim.y + threadIdx.y; //row

	// simple out of bounds check
	if( icol >= ncell_x  || irow >= ncell_y )
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] == 0)
		{
		inciPtr[self] = 0.0;
		return; // don't calculate if not in catchment(s) of interest
		}

	//double flow_C = 0.000001 ;
	double slope;
	slope = SlopePtr[self];
	if (slope > 0.5) slope = 0.5;
	double max_incision = maxincision;
			//slope * cellsize ; // changed for scaling 14/01/16

	//if (max_incision > maxincision) max_incision = 0.02;// from 0.0004

	double f_temp = (( FAPtr[self] ) * 1000 / 25)  - flow_thresh;
	//double f_temp =  (FAPtr[self] * 1000) - flow_thresh;
	
	if (f_temp < 0.0)
	{
		f_temp = 0.0;
		inciPtr[self] = 0.0;
	}
	else
	{
		//    but still account for vegetation effects
		if (slope == 0.0)
			inciPtr[self]= flow_A * pow (f_temp, flow_B) * 0.01 * exp (-0.07 * (TotBPtr[self]) / 76.5);
		else
			inciPtr[self] = flow_A * pow (f_temp, flow_B) * (slope) * exp (-0.07 * (TotBPtr[self]) / 76.5);

	}

	if ( inciPtr[self] > max_incision)	inciPtr[self] = max_incision;

	if ( inciPtr[self] > soilTPtr[self]) {

		double softrockinc = soilTPtr[self];
		double bedrockinc = ( 1.0- (softrockinc/inciPtr[self])) * 0.0001 * pow (f_temp, flow_B) * (slope) * exp (-0.07 * (TotBPtr[self]) / 76.5);
		inciPtr[self] = softrockinc + bedrockinc;
		if ( inciPtr[self] > max_incision)	inciPtr[self] = max_incision;
	}
}

__global__ void gelifluction( int *mask, double* geliPtr, double* temperature, double *eroPtr, double *soilTPtr, double *SlopePtr,
		                                                 int ncell_x, int ncell_y, double *FA, double FAMax, double Gslope, double kappa, double Gflowlim)
{
	/*! JW solifluxion component using version of Anderson (Culling) model Q=-k dz/dx where k=0.5*beta*f*zeta
	JW beta base value = 0.01-0.03 --> 0.02 (Anderson pers comm),
	JW f = 0.9236-0.0601*MAT-0.005MAT^2 (Matsuoka)
	JW zeta = 17.011-0.1745*MAT-0.0487MAT^2 (Grossi et al.)
	JW multiplying together gives fourth-order polynomial between -20�C<=MAT<=8.75�C, k set to 0.001 outside this range
    */

	int icol = blockIdx.x * blockDim.x + threadIdx.x; //column
	int irow = blockIdx.y * blockDim.y + threadIdx.y; //row

	// simple out of bounds check
	if( icol >= ncell_x  || irow >= ncell_y )
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] == 0)
		{
		geliPtr[self] = 0.0;
		return; // give nodatavalue if not in catchment(s) of interest
		}

	if (temperature[self] >= -20. && temperature[self] <= 8.87)
	{
		kappa = 0.15711 - 0.011835 * temperature[self] - 0.001195 * temperature[self] * temperature[self] + 0.000038 * temperature[self] * temperature[self] * temperature[self]  + 0.000002 * temperature[self] * temperature[self] * temperature[self] * temperature[self];
	}

	double slopethreshold;
	double scaler = 1 - (FA[self] / Gflowlim);
	if (FA[self] > Gflowlim) scaler = 0.0; // turn gelifluction off in channel cells with >60cm flow depth.

	slopethreshold = SlopePtr[self];
	if (slopethreshold > Gslope) slopethreshold = Gslope; // placed this constraint back in

	geliPtr[self] = kappa * slopethreshold * scaler; // * sfx_d;  //JW sfx is the scaling parameter passed via the parameter file or defaulting to 1. if not present

	if ((eroPtr[self] + geliPtr[self]) > soilTPtr[self])
	{
		geliPtr[self] = soilTPtr[self] - eroPtr[self];
	}

	if (geliPtr[self] < 0.0) geliPtr[self] = 0.0;
}

__global__ void surface_change (int *mask, double* dz_d, double* eroPtr_d, double* inciPtr_d, double* depoPtr_d, double* geliPtr_d, int rows, int cols)
{

  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

  int self = irow * cols + icol;
	if (mask[self] == 0)
		{
		dz_d[self] = 0.0;
		return;
		} // give nodatavalue if not in catchment(s) of interest

  dz_d[self] = (depoPtr_d[self] - eroPtr_d[self] - inciPtr_d[self] - geliPtr_d[self]) ;
}


__global__ void weathering (int *mask, double* temperature, double* rain, double *FAPtr_d, int type, double *weatherC_d, double *weatherP_d, double *soilTPtr_d,
	                                 double* finesPtr_d, double* stonePtr_d, double* soilMPtr_d, double* eroPtr_d, double* dz_d, double cell_size, int rows, int cols, double rescale)
{

  double surface_flow;
  double evap_trans;
  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

  int self = irow * cols + icol;

	if (mask[self] == 0)
		{
		weatherP_d[self] = 0.0;
		weatherC_d[self] = 0.0;
		return;
		} // don't calculate if not in catchment(s) of interest

	double PE_d = 1.44 * pow (temperature[self] + 6., 2);
	double P0_d = 7.7e-5;
	double P1_d = 2.3 ;

	// Based upon Dreybrodt as used by Kaufmann and Braun 2001 Terra Nova 13, 313-20.
	// Do these calculations once before entering kernel and transfer values into shared memory
    double t_abs = temperature[self] + 273.16;
    double t2 = t_abs * t_abs;
    double ltemp = log10 (t_abs);
    double DHA = -0.4883 + 8.074e-4 * temperature[self];
    double DHB = -0.3241 + 1.600e-4 * temperature[self];
    double sqI = sqrt (0.1);
    double denom1 = (-4. * DHA * sqI) / (1. + 5.e-8 * DHB * sqI);
    double gammCa_d = pow (10., denom1);
    double gammHCO3_d = pow (10., (-DHA * sqI) /  (1. + 5.4e-8 * DHB * sqI));
    double K1_d = pow (10., -356.3094 - 0.06091964 * t_abs + 21834.37 / t_abs + 126.8339 * ltemp -  1684915 / t2);
    double K2_d = pow (10., -107.8871 - 0.03252849 * t_abs +  5151.79 / t_abs + 38.92561 * ltemp - 563713.9 / t2);
    double KC_d = pow (10., -171.9065 - 0.077993 * t_abs +  2839.319 / t_abs + 71.595 * ltemp);
    double KH_d = pow (10., 108.3865 + 0.01985076 * t_abs - 6919.53 / t_abs - 40.45154 * ltemp +  669365.0 / t2);


  surface_flow = FAPtr_d[self] * 1000; // / 5;  // converted back to mm here
  if (rain[self] + surface_flow <= PE_d)
			evap_trans = rain[self] + surface_flow;
  else
			evap_trans = rain[self] + surface_flow - PE_d;


	switch (type)
	{
        case 0:  // White regression method (see Dreybrodt)
	   	    weatherC_d[self] = 6.3 + 0.049 * (rain[self] - evap_trans);
	   	    break;

        case 1:  // Dreybrodt as used by Kaufmann and Braun 2001 Terra Nova 13, 313-20.
			double PCO2;
			if (soilTPtr_d[self] > 0.1)
				PCO2 = pow (10., -2.5 + ((temperature[self] - 15.) / 10.));
			else
				PCO2 = pow (10., -3.5 + ((temperature[self] - 15.) / 10.));

			double Ca2eq = pow (PCO2 * ((K1_d * KC_d * KH_d) / (4. * K2_d * gammCa_d *  gammHCO3_d * gammHCO3_d)), (1. / 3.));

			//weatherC_d[self] = 0.0148148148 * Ca2eq * surface_flow / (cell_size * cell_size);
			weatherC_d[self] = 0.0148148148 * Ca2eq * (surface_flow / rescale);  //denominator now has same scalar as the original 100*100, changed from 10k on 30/11/15, see above 14/01/16

			if (weatherC_d[self] <= 0.0) weatherC_d[self] = 0.0;
			break;
	}
	if (weatherC_d[self]>0.001) weatherC_d[self]=0.001; // new constraint added 13/10/15

	weatherP_d[self] = P0_d * exp (-(P1_d) * soilTPtr_d[self]);

	// following makes the assumption of 2% solid residue from chemical
	//    weathering
	soilTPtr_d[self] = soilTPtr_d[self] + weatherC_d[self] * 0.02 + weatherP_d[self]; //soilthickness already adjusted in mfd_accum

	if (soilTPtr_d[self] <= 0.0)
	{
		soilTPtr_d[self] = 0.0;
        stonePtr_d[self] = 0.0;
        finesPtr_d[self] = 0.0;
        soilMPtr_d[self] = 0.02;
	}

	else
	{
		if (weatherP_d[self] > 0)
		//	stonePtr_d[self] += stonePtr_d[self] * (tanh (((eroPtr_d[self] / soilTPtr_d[self]) - 0.5) / 0.25) * -tanh (((weatherC_d[self] / weatherP_d[self]) - 1.) / 0.25));
		stonePtr_d[self] += stonePtr_d[self] * (tanh (((eroPtr_d[self] / soilTPtr_d[self])) ) * -tanh (((weatherC_d[self] / weatherP_d[self]))) );
		else
		//	stonePtr_d[self] += stonePtr_d[self] * tanh (((eroPtr_d[self] / soilTPtr_d[self]) - 0.5) / 0.25);
		stonePtr_d[self] += stonePtr_d[self] * tanh (((eroPtr_d[self] / soilTPtr_d[self]) ));

		if (stonePtr_d[self] < 0.01) stonePtr_d[self] = 0.01;
		else if (stonePtr_d[self] > 99.9) stonePtr_d[self] = 99.9;

		finesPtr_d[self] = 100. - stonePtr_d[self];

		soilMPtr_d[self] -= (1.e-3 * evap_trans / soilTPtr_d[self]); // 1.e-3 converts ET to m

		if (soilMPtr_d[self] < 0.02) soilMPtr_d[self] = 0.02;
        else if (soilMPtr_d[self] > 0.5) // assume constant porosity of 50% for the moment
		                   soilMPtr_d[self]  = 0.5;  // but don't feedback water sitting on surface at end of year yet
	}

}

void calc_diff_erosion(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;

	int block_ncell_y = 16;
	int block_ncell_x = 16;

	int full_size = ncell_x * ncell_y;

	dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
	dim3 dimBlock(block_ncell_x, block_ncell_y);

	diff_erosion<<< dimGrid, dimBlock >>>( device->mask, device->eroPtr, device->fa, device->SlopePtr, device->finesPtr,
	                                  device->nutPtr, device->soilTPtr, device->stonePtr, device->TotBPtr, data->soil_struct, data->profile_permeability,
	                                  ncell_x, ncell_y, data->Qthreshold);

	fprintf(data->outlog, "MOD: diff_erosion :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors( cudaMemcpy ( data->eroPtr,   device->eroPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	thrust::device_ptr<double> difftot_d = thrust::device_pointer_cast(device->eroPtr);
	cudaSetDevice(0);
	data->totE = thrust::reduce(difftot_d, difftot_d + full_size, (double) 0);
	fprintf(data->outlog, "total concentrated  from thrust is %10.8lf \n", data->totE);

	if (data->totE == NAN)
	{
		printf("data from totE is NaN \n");
		exit(0);
	}
}

void calc_conc_erosion(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	int block_ncell_y = 16;
	int block_ncell_x = 16;

	dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
	dim3 dimBlock(block_ncell_x, block_ncell_y);


	conc_erosion<<< dimGrid, dimBlock >>>( device->mask, device->inciPtr, device->fa, device->SlopePtr, device->TotBPtr,
			                                 device->soilTPtr, ncell_x, ncell_y, data->mapInfo.cellsize, data->Qthreshold,
			                                           data->Flow_A, data->Flow_B, data->Flow_C, data->max_incision);

	fprintf(data->outlog, "MOD: conc_erosion :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors( cudaMemcpy ( data->eroPtr,   device->eroPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	thrust::device_ptr<double> incitot_d = thrust::device_pointer_cast(device->inciPtr);
	cudaSetDevice(0);
	data->totI = thrust::reduce(incitot_d, incitot_d + full_size, (double) 0);
	fprintf(data->outlog, "total Incision from thrust is %10.8lf \n", data->totI);

	if (data->totI == NAN) {
		printf("data from totI is NaN \n");
		exit(0);
	}
}

void calc_gelifluction(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	int block_ncell_y = 16;
	int block_ncell_x = 16;

	dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
	dim3 dimBlock(block_ncell_x, block_ncell_y);

	const double sfxh = 1.0;
	const double* hsfx = &sfxh;

	if( !(cudaMemcpyToSymbol(sfx_d, hsfx,              sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess))
		printf ("failed \n");

	fprintf(data->outlog, "FA_max before gelifluction %f \n", data->FA_max);

	gelifluction<<< dimGrid, dimBlock >>>( device->mask, device->geliPtr, data->tempmat, device->eroPtr, device->soilTPtr, device->SlopePtr,
			                                                ncell_x, ncell_y, device->fa, data->FA_max, data->Gslope, data->kappa, data->Gflowlim);

	fprintf(data->outlog, "MOD: gelifluction :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors( cudaMemcpy ( data->geliPtr,   device->geliPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	thrust::device_ptr<double> gelitot_d = thrust::device_pointer_cast(device->geliPtr);
	cudaSetDevice(0);
	data->totG = thrust::reduce(gelitot_d, gelitot_d + full_size, (double) 0);
	fprintf(data->outlog, "total gelifluction from thrust is %10.8lf \n", data->totG);


	if (data->totG == NAN) {
		printf("data from totG is NaN \n");
		exit(0);
	}
}


void calc_dz(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	int block_ncell_y = 16;
	int block_ncell_x = 16;

	dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
	dim3 dimBlock(block_ncell_x, block_ncell_y);

	// Calculate dz due to erosion and depositional processes
	surface_change<<<dimGrid, dimBlock>>>(device->mask, device->dz, device->eroPtr, device->inciPtr, device->depoPtr, device->geliPtr, data->mapInfo.height, data->mapInfo.width);
	fprintf(data->outlog, "MOD: surface_change :%s\n", cudaGetErrorString(cudaGetLastError()));

	double max_dz, min_dz;
	thrust::device_ptr<double> max_dz_d = thrust::device_pointer_cast(device->dz);
	thrust::device_ptr<double> min_dz_d = thrust::device_pointer_cast(device->dz);
	max_dz = thrust::reduce(max_dz_d, max_dz_d + full_size, (double) 0, thrust::maximum<double>() );
	min_dz = thrust::reduce(min_dz_d, min_dz_d + full_size, (double) 0, thrust::minimum<double>() );
	printf("Maximum dz = %f, Minimum dz = %f \n", max_dz, min_dz);
}


void calc_weathering(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;

	int type = 1;
	int block_ncell_y = 16;
	int block_ncell_x = 16;

	int full_size = ncell_x * ncell_y;

	dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
	dim3 dimBlock(block_ncell_x, block_ncell_y);

	  double rescale = 10000 * (10000 / ( data->mapInfo.cellsize * data->mapInfo.cellsize ));
	  fprintf(data->outlog, "wC scaler = %f \n", rescale);

	weathering<<<dimGrid, dimBlock>>>(device->mask, data->tempmat, data->rainmat, device->fa, type, device->weatherC, device->weatherP, device->soilTPtr, device->finesPtr, device->stonePtr, device->soilMPtr,
		                              device->eroPtr, device->dz, data->mapInfo.cellsize, data->mapInfo.height, data->mapInfo.width, rescale);

	thrust::device_ptr<double> weatherC_d = thrust::device_pointer_cast(device->weatherC);
	data->totweatherC = thrust::reduce(weatherC_d, weatherC_d + full_size, (double) 0);
	fprintf(data->outlog, "total wC from thrust is %30.28lf \n", data->totweatherC);
	if (data->totweatherC == NAN) {
		printf("data from totweatherC is NaN \n");
		exit(0);
	}

	thrust::device_ptr<double> weatherP_d = thrust::device_pointer_cast(device->weatherP);
	data->totweatherP = thrust::reduce(weatherP_d, weatherP_d + full_size, (double) 0);
	fprintf(data->outlog, "total wP from thrust is %10.8lf \n", data->totweatherP);
	fprintf(data->outlog, "MOD: weathering :%s\n", cudaGetErrorString(cudaGetLastError()));
}

