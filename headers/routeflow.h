#ifndef SFDH
#define SFDH

#include "lem.h"
#include "Data.h"
#include "Directions.h"

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

#include "newflow.h"


typedef int TYPE;

__constant__ const TYPE code[] = {64 , 128, 1, 2, 4, 8, 16, 32, 0};

void cuFlowDirection(Data* data, Data* device, int iter);

#endif
