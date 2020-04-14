/*
 * memory
 *
 *  Created on: Jan 31, 2014
 *      Author: ndm12
 */

#ifndef MEMORY_
#define MEMORY_

#include "lem.h"

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

/*
struct is_not_zero
{
	__host__ __device__ bool operator() (double x)
	{
		return (x != 0.0);
	}
};

struct is_not_negative
{
	__host__ __device__ bool operator() (double x)
	{
		return (x >= 0.0);
	}
};
*/

int clearcontribAspace(Data* data);

int createfilenamespace(Data* data);

int createProcessMatrices(Data* data);

int deleteProcessMatrices(Data* data);

int createCatchmentSpace(Data* data, Catchment* Catchments);

void createSoilTfromformula(Data* data);

int createmask(Data* data);

#endif /* MEMORY_ */
