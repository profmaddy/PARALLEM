#ifndef newflowH
#define newflowH

#include "lem.h"
#include "Data.h"
#include "Directions.h"

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


typedef int TYPE;

class watershed_edge
{
    private:
        int node1;
        int node2;
        HASH_TYPE weight;

	public:
        watershed_edge();
        watershed_edge(int n1, int n2, HASH_TYPE w);
        ~watershed_edge();
        int get_node1();
        int get_node2();
        void set_node1(int n);
        void set_node2(int n);
		void set_weight(HASH_TYPE w);
        HASH_TYPE get_weight() const;
        bool operator<(const watershed_edge &other);
};

class union_set
{
    private:
        int *parent;
        int *rank;
        int size;
    public:
        union_set(int s);
        ~union_set();
        bool in_set(int x);
        void make_set(int x);
        int find_root(int x);
        void merge(int a, int b);
};

void DirFlow (Data *data);

__global__ void FlowDirs(int *mask, int *flatmask, double *zs, double *slopes, int *SFD, double *prop, int *mfd, int flowtype, int cell_size, int ncell_x, int ncell_y, int *dx, int *dy, int last,int *sinkcounter, int *flatcounter);

__global__ void flow_boundaryMFD( int *mask, double *zs, double *slopes, int *SFD, int *mfd, int ncell_x, int ncell_y);

__global__ void shortest_paths_plateaus_initMFD(int *mask, double *zs, int *fd, int* SFD, float* shortest_paths, int cell_size, int ncell_x, int ncell_y, double *lowHeight);

__global__ void route_plateausMFD(int *mask, double *zs, double*slopes, double*props, int *fd, int* SFD, float* shortest_paths,	int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag, double *lowHeight);

__global__ void init_watershedMFD(int *mask, int *SFD, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y);

__global__ void init_watershed_sinkMFD(int *mask, int *SFD, int *watershed_id, double *zs, int *dx, int *dy, int ncell_x, int ncell_y, int *changed);


__device__ double atomicMin(double* address, double val);

__global__ void final_watershedMFD(int *mask, int *aspects, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y);

__global__ void identify_watershedMFD(int *mask, int *aspects, int *watershed_id, int *dx, int *dy, int ncell_x, int ncell_y, int *changed);

__global__ void init_hash_tableMFD(HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size);

__global__ void identify_watershed_boundary_write_backMFD(int *mask, int *watershed_id, double *zs, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, int *dx, int *dy, int ncell_x, int ncell_y, int *is_succeeded, int* counter);

__global__ void raise_watershedMFD(int *mask, int *watershed_id, HASH_TYPE *raise, double *zs, int ncell_x, int ncell_y);

int init_and_sort_watershed_bound_labelsMFD(Data* data, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, watershed_edge *sorted_list, int sorted_list_size, int *max_watershed_p, int max_sink_id);



__global__ void slope_plateausMFD(int *mask, double *zs, float* shortest_paths, int ncell_x, int ncell_y, double *lowHeight, double *slopes, int cell_size);

__global__ void comp_shortest_paths_sinkMFD(int *mask, double *zs, int *aspects, float* shortest_paths, int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag);








#endif
