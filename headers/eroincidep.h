#ifndef eroincidepH
#define eroincidepH	

//#include "lem.h"
#include "Data.h"
#include "Directions.h"
#include "device_constants.cuh"
#include "production.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_runtime_api.h"
#include "helper_cuda.h"

#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/fill.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>
#include "production.h"

void calc_diff_erosion(Data* data, Data* device);
void calc_conc_erosion(Data* data, Data* device);
void calc_gelifluction(Data* data, Data* device);
void calc_dz(Data* data, Data* device) ;
void calc_weathering(Data* data, Data* device);


#endif
