/*
 * ioGDAL.h
 *
 *  Created on: 30 Mar 2020
 *      Author: ndm12
 */

#ifndef IOGDAL_H_
#define IOGDAL_H_

#include "lem.h"
#include "MapInfo.h"
#include "production.h"

#include "gdal.h" // note this is a different header file to the c++ version
#include "cpl_port.h"
#include "gdal_version.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_core.h"
#include "ogr_srs_api.h"

int readgdalDEMfromFile(Data* data);

int readBurnFDfromFile(Data* data);

int readgdalBedrockfromFile(Data* data, Data* device);

int readgdalmaskfromFile(Data* data) ;

void savegrids(Data* data, Catchment* catchments,  int iteration);

int writegdalGRIDtoFile(Data* data, Catchment* catchments, char* fileName, int whichtype, int what) ;

int writeGRIDtoFile(Data* data, char* fileName, int whichtype, int what);

int retriveProcessMatrices(Data* data);

#endif /* IOGDAL_H_ */
