/*
 * io.h
 *
 *  Created on: 30 Mar 2020
 *      Author: ndm12
 */

#ifndef IO_H_
#define IO_H_

#include "lem.h"

struct header
{
    int ncols;
    int nrows;
    float xllcorner;
    float yllcorner;
    float cellsize;
    float nodata;
};

//io
void readinpar(Data* data, const char* file);
int checkparread(Data* data);

int readgrid(Data* data);
int write_float(Data *data, float *grid, char *filename);
int write_int(Data *data, int *grid, char *filename);
int write_double(Data *data, double *grid, char *filename);

int writeSummaryDataToFile(Data* data, int iteration);
int ReadClimateDataFromFile(Data* data);



#endif /* IO_H_ */
