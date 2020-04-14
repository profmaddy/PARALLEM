#ifndef eroincidepH
#define eroincidepH	

#include "lem.h"
#include "Data.h"
#include "Directions.h"
#include "device_constants.cuh"
#include "production.h"

void calc_diff_erosion(Data* data, Data* device);
void calc_conc_erosion(Data* data, Data* device);
void calc_gelifluction(Data* data, Data* device);
void calc_dz(Data* data, Data* device) ;
void calc_weathering(Data* data, Data* device);


#endif
