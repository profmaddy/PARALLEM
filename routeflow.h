#ifndef SFDH
#define SFDH

#include "lem.h"
#include "Data.h"


//#include "Directions.h"
#include "newflow.h"
//#include "config.h"
//#include "production.h"
//#include "newflow.h"

typedef int TYPE;

__constant__ const TYPE code[] = {64 , 128, 1, 2, 4, 8, 16, 32, 0};

void cuFlowDirection(Data* data, Data* device, int iter);

#endif
