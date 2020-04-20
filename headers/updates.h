/*
 * updates.h
 *
 *  Created on: Feb 25, 2014
 *      Author: ndm12
 */
#include "lem.h"
#include "Data.h"

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

#ifndef UPDATES_H_
#define UPDATES_H_

void update_newSurface(Data* data, Data* device, int iter);
void update_nutrients(Data* data, Data* device);
void update_vegetation(Data* data, Data* device);

#endif /* UPDATES_H_ */
