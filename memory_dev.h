
#include "lem.h"
#include "Data.h"
#include "production.h"

#ifndef memory_devH
#define memory_devH

/*struct is_not_zero
{
	__host__ __device__ bool operator() (double x)
	{
	return (x != 0.0) ;
	}
};

struct is_not_negative
{
	__host__ __device__ bool operator() (double x)
	{
	return (x >= 0.0) ;
	}
};
*/


void setslopesspace(Data* data, Data* device);

void setdevicespace_FD(Data* data, Data* device);
void setdevicespace_FA(Data* data, Data* device);
void setdevicespace_Process(Data* data, Data* device);

void cleardevicespace_FD(Data* data, Data* device);
void cleardevicespace_FA(Data* data, Data* device);
void cleardevicespace_Process(Data* data, Data* device);

void createDeviceSpace(Data* data, Data* device);
int clearDeviceSpace(Data* data, Data* device);

int copyMask(Data* data, Data* device);
int copylastclimate(Data* data, Data* device);

int cleargrid(Data* data);

int zerogrids(Data* data);

int copytolastclimate(Data* data, Data* device);

#endif
