
#ifndef EROSION_H_
#define EROSION_H_

//#include "lem.h"
//#include "Data.h"
//#include "eroincidep.h"
//#include "config.h"
//#include "Directions.h"
//#include "production.h"
#include "eroincidep.h"
#include "mfd_accum.h"


void erosionGPU( Data* data, Data* device, int iter);


#endif /* EROSION_H_ */
