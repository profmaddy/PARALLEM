/*
 * lem.h
 *
 *  Created on: Feb 18, 2014
 *      Author: ndm12
 */

#ifndef LEM_H_
#define LEM_H_

// Skeleton outline of CUDA enabled Landscape Evolution Model
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <float.h>
#include <iostream>

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>

#include "config.h"
#include "MapInfo.h"
#include "production.h"

#include "Data.h"
#include "io.h"
//#include "ioGDAL.h"

#include "memory.h"
#include "memory_dev.h"

//#include "bypass.h"
#include "routeflow.h" // gives access to flow direction routine
#include "runoffweight.h" // gives access to runoffweight calculations (losses)


//#include "headers/FA_SFD.h" //gives access to flow accumulation algorithm
#include "mfd_simple.h" // gives access to mfd flow accumulation
#include "erosion.h" // gives access to sediment erosion, transport, deposition model


//#include <cuda.h>
//#include "device_launch_parameters.h"
//#include "helper_cuda.h"
//#include <device_functions.h>

//#include "helper_timer.h"
//#include "helper_functions.h"
//#include "cuda_runtime_api.h"
//#include "cuda_runtime.h"
//#include <device_launch_parameters.h>

//#include <thrust/functional.h>
//#include <thrust/reduce.h>
//#include <thrust/copy.h>
//#include <thrust/count.h>
//#include <thrust/fill.h>
//#include <thrust/extrema.h>
//#include <thrust/execution_policy.h>
//#include "production.h"


#endif /* LEM_H_ */
