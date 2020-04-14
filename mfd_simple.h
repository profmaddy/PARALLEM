#ifndef MFDH
#define MFDH

#include "lem.h"
#include "Data.h"
#include "config.h"
#include "Directions.h"
#include "production.h"

void correctmfdflow(Data* data, Data* device, int iter);
int calcslopetotal(Data* data, Data* device);

#endif
