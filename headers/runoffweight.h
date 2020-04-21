#include "lem.h"
#include "Data.h"
#include "config.h"
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



#ifndef RUNOFFWEIGHT
#define RUNOFFWEIGHT

int computeRunOffWeights(Data* data, Data* device);
void calcprops(Data* data);
void calcwater(Data* data);

#endif
