/*
 * mfd_accum.h
 *
 *  Created on: 6 Oct 2016
 *      Author: ndm12
 */

#ifndef MFD_ACCUM_H_
#define MFD_ACCUM_H_

//#include "lem.h"
#include "Data.h"
#include "Directions.h"
#include "config.h"
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

void sedmfdaccum(Data* data, Data* device);

#endif /* MFD_ACCUM_H_ */
