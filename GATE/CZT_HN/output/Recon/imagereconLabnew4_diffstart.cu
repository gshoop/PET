
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <list>
#include <cstdlib>
#include <sstream>
//#include <random>
#include <iomanip>
#include <limits>
#include <cuda_runtime.h>
#include "device_functions.h"
#include <sys/time.h>
using namespace std;

#define sharesize 110
#define blocksPerGrid 30
#define threadsperBlock 32
#define ThreshLineValue 0.001
#define SRTWO 1.414
#define reducsize 1024
#define UpholdVox 20
#define CylRadius 55.0
#define CylHeight 126.0
#define MU 0.0096

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

class TimingCPU {

    private:
        long cur_time_;

    public:

        TimingCPU(): cur_time_(0) {};

        ~TimingCPU() {};

        void StartCounter(){
			struct timeval time;
			if(gettimeofday( &time, 0 )) return;
			cur_time_ = 1000000 * time.tv_sec + time.tv_usec;
		}

        double GetCounter(){
			struct timeval time;
			if(gettimeofday( &time, 0 )) return -1;

			long cur_time = 1000000 * time.tv_sec + time.tv_usec;
			double sec = (cur_time - cur_time_) / 1000000.0;
			if(sec < 0) sec += 86400;
		    cur_time_ = cur_time;
		
		    return 1000.*sec; //wall clock time (ms)
		}

} timing; // CPU timer.

struct PrivateTimingGPU {
    cudaEvent_t     start;
    cudaEvent_t     stop;
};

class TimingGPU {
    private:
        PrivateTimingGPU *privateTimingGPU;

    public:

        TimingGPU() { privateTimingGPU = new PrivateTimingGPU; }

        ~TimingGPU() {}

        void StartCounter() {
            cudaEventCreate(&((*privateTimingGPU).start));
            cudaEventCreate(&((*privateTimingGPU).stop));
            cudaEventRecord((*privateTimingGPU).start,0);
        }

        float GetCounter(){
            float   time;
            cudaEventRecord((*privateTimingGPU).stop, 0);
            cudaEventSynchronize((*privateTimingGPU).stop);
            cudaEventElapsedTime(&time,(*privateTimingGPU).start,(*privateTimingGPU).stop);
            return time;
        }

} timinggpu; // TimingGPU class

class timevar{
	public:
		float txforward;
		float txbackward;
		float tyforward;
		float tybackward;
		float tzforward;
		float tzbackward;
		float tpostimageprocess;
		float memoryIO;
		float lorsorting;
		timevar(): txforward(0.0), txbackward(0.0), tyforward(0.0), tybackward(0.0), tzforward(0.0), tzbackward(0.0), tpostimageprocess(0.0), memoryIO(0.0), lorsorting(0.0) {}
		~timevar() {}
	
		void printvalue(){
			printf("%12s%12s%12s%12s%12s%12s%20s%12s%12s", "(ms) xf","yf", "zf","xb", "yb","zb","postprocess","memoryIO","lorsorting\n");
			printf("%12.1f%12.1f%12.1f%12.1f%12.1f%12.1f%20.1f%12.1f%12.1f\n", txforward, tyforward, tzforward, txbackward, tybackward, tzbackward, tpostimageprocess, memoryIO, lorsorting);
	}
	
} timeall;

//__constant__ int *imageindex;
//__constant__ int *lorindex;
//__constant__ float *info;
__device__ __constant__ float aves[1] , avep[1], aveunf[1], d_bsgm[3], d_rads[3], d_info[4], d_norm[2];
__device__ __constant__ int d_indr[3], d_imageindex[4], d_lorindex[3];

struct lor
{
	float x1;
	float y1;
	float z1;
	float x2;
	float y2;
	float z2;
	int mainaxis;	//0 for x, 1 for y, 2 for z
	float weight;
};

struct cudalor
{
	float *x1;
	float *y1;
	float *z1;
	float *x2;
	float *y2;
	float *z2;
	float *linevalue;
	float *weight;
} xlor, ylor, zlor, dev_xlor, dev_ylor, dev_zlor;

int numcal=0, numcal2=0;

vector<string> explode(string s, char c);
int gaussianblur(int nx, int ny, int nz);
int preplor(string fin, int senmap);
__global__ void calnewmatrix000(float *snmatrix, float *smatrix);
__global__ void calnewmatrix100(float *snmatrix, float *smatrix, float *normimage);
__global__ void calnewmatrix010(float *snmatrix, float *smatrix, float *poimage);
__global__ void calnewmatrix110(float *snmatrix, float *smatrix, float *normimage, float *poimage);
__global__ void calnewmatrix011(float *snmatrix, float *smatrix, float *poimage, float *bmatrix, float *allweight);
__global__ void calnewmatrix111(float *snmatrix, float *smatrix, float *normimage, float *poimage, float *bmatrix, float *allweight);
__global__ void calave(float *smatrix, float *gave);
__global__ void calavewithfilter(float *smatrix, float *gave, int *gnumave);
__global__ void gpublur(float *smatrix, float *bmatrix, float *allweight);
__global__ void calLogLike(float *xlinevalue, double *gloglike, const int lorindex);	
__global__ void calLogLikeS(float *smatrix, float *normimage, double *gloglike, const int msize, const int norma);
__global__ void calLogR(float *smatrix, float *poimage, double *gloglike, const int msize);


float a = 4.; //grid size for sensitivity map.
float bndry[3] = {200., 200., 216.};    //FOV x, y, z. Unit: mm.
int msize;
float torhw, torsgm;
float beta, sumpop, sumimp, sumb;
int rgl;	//indicator for regularization. 0: no regularization, 1: regularization
int blur;	//indicator for blurring. 0: no blurring, 1: blurring in regularization
int norma; //indicator for normalization. 0: no normalization, 1: normalization using data and generate sensitivity image, 2: normalization using sensitivity image
float ThreshNorm;	//threshold for normalization
float bsgm[3], bthresh=0.01;	//sigma x,y,z in blurring Gaussian function
float deltax = 0.0001;  //delta used in calculating derivative numerically
float rads[3];	//radius for Gaussian blurring
int indr[3];	//number of grid for the radius

float *smatrix;  //matrix for image
float *snmatrix; //matrix for image
float *poimage; //prior image
float *lastimage; //prior image
float *bmatrix;	//blurred image of smatrix
float *deri;    //derivative of regularization function
float *allweight;	//all the sumed weight in blurred image.
float *normimage, *dev_normimage;	//normalization image (sensitivity image)

int numline = 0;    //number of lines/lors in each input file.
int *nummainaxis;   //number of lines in each main axis (x,y,z)
int wgt = 0;	//whether include attenuation and normalization in system matrix.
	

//__global__ void testf(float *dev_test)
//{
//	//if(threadIdx.x == 3) dev_test[0] = threadIdx.x + 1;
//	dev_test[0] = powf(maxs[0],2);
//}

__global__ void xfpro( cudalor lor, float *smatrix ) 
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
	float a = d_info[0], torhw = d_info[1], torsgm2 = d_info[2];
	int lornum = d_lorindex[0];

	__shared__ float cache[sharesize][sharesize];
	int tid ;
	int cacheIndex1 , cacheIndex2;

	float linevalue = 0.;
	float ulen2, t, oy, oz;
	int mlyy,mhyy,mlzz,mhzz;
	float x1,x2,y1,y2,z1,z2,weight;

	int tilenum1 = (ny + sharesize - 1)/ sharesize, tilenum2 = (nz + sharesize - 1) / sharesize ;

	for(int i=0; i< nx; i++)
	{

        for(int tn1 = 0; tn1 < tilenum1; tn1++)	//for each tile of image. This is due to limited shared memory.
        {
		for(int tn2 = 0; tn2 < tilenum2; tn2++)
		{

        cacheIndex1 = threadIdx.x;

		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < ny )	//load a tile of image
        {
		cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < nz)
        {
            cache[cacheIndex1][cacheIndex2] = smatrix[i + ((sharesize * tn1) + cacheIndex1) * nx + ((sharesize * tn2) + cacheIndex2) * nx * ny];
            cacheIndex2 += blockDim.y;
        }
            cacheIndex1 += blockDim.x;
        }
        __syncthreads();

		tid = threadIdx.x + threadIdx.y * blockDim.x +  blockIdx.x * blockDim.x * blockDim.y;

		while(tid < lornum)
		{
			x1 = lor.x1[tid];
			y1 = lor.y1[tid];
			z1 = lor.z1[tid];
			x2 = lor.x2[tid];
			y2 = lor.y2[tid];
			z2 = lor.z2[tid];
			weight = lor.weight[tid];
			linevalue = 0.;
	
			ulen2 = powf(x1-x2,2) + powf(y1-y2,2) + powf(z1-z2,2);

			t = ( i * a - x1) / (x2 - x1);
			
			oy = y1 + t * (y2 - y1);
			oz = z1 + t * (z2 - z1);
			
			mlyy = max((int)truncf((oy - (SRTWO * torhw ))/a)+1, 0);
			mhyy = min((int)truncf((oy + (SRTWO * torhw ))/a), ny - 1);
			mlzz = max((int)truncf((oz - (SRTWO * torhw ))/a)+1, 0);
			mhzz = min((int)truncf((oz + (SRTWO * torhw ))/a), nz - 1);

			mlyy = max(mlyy, sharesize * tn1);
			mhyy = min(mhyy, sharesize * (tn1 + 1)-1);
			mlzz = max(mlzz, sharesize * tn2);
			mhzz = min(mhzz, sharesize * (tn2 + 1)-1);
			
			for(int ky = mlyy; ky <= mhyy; ky++)
			{
			    for(int kz = mlzz; kz <= mhzz; kz++)
			    {
					float dy = oy - ky*a, dz = oz - kz*a;
					float inner = dy * (y1-y2) + dz * (z1 - z2);
					float dst2 = dy * dy + dz * dz - inner * inner / ulen2;
					float maxdst2 = torhw * torhw;
			        //dst = sqrtf(pow(oy-ky*a, 2) + powf(oz-kz*a, 2) - powf((oy-ky*a) * (y1-y2) + (oz-kz*a) * (z1-z2), 2) / ulen2);
			        if(dst2 < maxdst2) linevalue += cache[ky -sharesize * tn1 ][ kz - sharesize * tn2] * expf(-dst2/(2.0f * torsgm2)) * weight;
			    }
			}
			lor.linevalue[tid] += linevalue;
			//finish one tile for one lor
			tid += blockDim.x * blockDim.y * gridDim.x;
		}
		__syncthreads();
		}
		}
	}


}

__global__ void xbpro( cudalor lor, float *snmatrix ) 
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
	float a = d_info[0], torhw = d_info[1], torsgm2 = d_info[2];
	int lornum = d_lorindex[0];

	__shared__ float cache[sharesize][sharesize];
	int tid ;
	int cacheIndex1 , cacheIndex2;

	float linevalue, rlinevalue;
	float ulen2, t, oy, oz;
	int mlyy,mhyy,mlzz,mhzz;
	float x1,x2,y1,y2,z1,z2,weight;

	int tilenum1 = (ny + sharesize - 1)/ sharesize, tilenum2 = (nz + sharesize - 1) / sharesize ;

	for(int i=0; i< nx; i++)
	{

        for(int tn1 = 0; tn1 < tilenum1; tn1++)	//for each tile of image. This is due to limited shared memory.
        {
		for(int tn2 = 0; tn2 < tilenum2; tn2++)
		{

        cacheIndex1 = threadIdx.x;

		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < ny )	//load a tile of image
        {
		cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < nz)
        {
            cache[cacheIndex1][cacheIndex2] = 0.0f;
            cacheIndex2 += blockDim.y;
        }
            cacheIndex1 += blockDim.x;
        }
        __syncthreads();

		tid = threadIdx.x + threadIdx.y * blockDim.x +  blockIdx.x * blockDim.x * blockDim.y;

		while(tid < lornum)
		{
			x1 = lor.x1[tid];
			y1 = lor.y1[tid];
			z1 = lor.z1[tid];
			x2 = lor.x2[tid];
			y2 = lor.y2[tid];
			z2 = lor.z2[tid];
			weight = lor.weight[tid];
	        linevalue = lor.linevalue[tid];
	
	        if(linevalue < ThreshLineValue) rlinevalue = 0.0f;
	        else rlinevalue = 1.0f / linevalue;
	
			ulen2 = powf(x1-x2,2) + powf(y1-y2,2) + powf(z1-z2,2);

			t = ( i * a - x1) / (x2 - x1);
			
			oy = y1 + t * (y2 - y1);
			oz = z1 + t * (z2 - z1);
			
			mlyy = max((int)truncf((oy - (SRTWO * torhw ))/a)+1, 0);
			mhyy = min((int)truncf((oy + (SRTWO * torhw ))/a), ny - 1);
			mlzz = max((int)truncf((oz - (SRTWO * torhw ))/a)+1, 0);
			mhzz = min((int)truncf((oz + (SRTWO * torhw ))/a), nz - 1);

			mlyy = max(mlyy, sharesize * tn1);
			mhyy = min(mhyy, sharesize * (tn1 + 1)-1);
			mlzz = max(mlzz, sharesize * tn2);
			mhzz = min(mhzz, sharesize * (tn2 + 1)-1);
			
			for(int ky = mlyy; ky <= mhyy; ky++)
			{
			    for(int kz = mlzz; kz <= mhzz; kz++)
			    {
					float dy = oy - ky*a, dz = oz - kz*a;
					float inner = dy * (y1-y2) + dz * (z1 - z2);
					float dst2 = dy * dy + dz * dz - inner * inner / ulen2;
					float maxdst2 = torhw * torhw;
			        //dst = sqrtf(pow(oy-ky*a, 2) + powf(oz-kz*a, 2) - powf((oy-ky*a) * (y1-y2) + (oz-kz*a) * (z1-z2), 2) / ulen2);
			        if(dst2 < maxdst2) atomicAdd(&cache[ky -sharesize * tn1 ][ kz - sharesize * tn2], expf(-dst2/(2.0f * torsgm2)) * rlinevalue * weight) ;
			    }
			}

			tid += blockDim.x * blockDim.y * gridDim.x;
		}
		__syncthreads();

		//write the tile of image to global memory
		cacheIndex1 = threadIdx.x;

		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < ny )	//load a tile of image
        {
        cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < nz)
        {
            atomicAdd(&snmatrix[i  + ((sharesize * tn1) + cacheIndex1) * nx + ((sharesize * tn2) + cacheIndex2) * nx * ny], cache[cacheIndex1][cacheIndex2]);
            cacheIndex2 += blockDim.y;
        }
            cacheIndex1 += blockDim.x;
        }
		__syncthreads();

		}
		}
	}


}


//add f to functions such as expf(), powf(), which stands for single precision. reduce some redundant variables. the trouble comes from cache[][], need to move cacheIndex2 increment out of while.some cores not launched due to tid limit.
__global__ void yfpro( cudalor lor, float *smatrix ) 
{

    int nx = d_imageindex[0], ny =d_imageindex[1], nz = d_imageindex[2];
    float a = d_info[0], torhw = d_info[1], torsgm2 = d_info[2];
	int lornum = d_lorindex[1];

	__shared__ float cache[sharesize][sharesize];
	int tid;
	int cacheIndex1, cacheIndex2;

    float linevalue = 0.;
    float ulen2, t, ox,oz;
    int mlxx,mhxx,mlzz,mhzz;
    float x1,x2,y1,y2,z1,z2,weight;

	int tilenum1 = (nx + sharesize - 1)/ sharesize, tilenum2 = (nz + sharesize - 1) / sharesize ;

//	testmatrix[0] = 20.0f;

    for(int i=0; i< ny; i++)
    {

        for(int tn1 = 0; tn1 < tilenum1; tn1++)	//for each tile of image. This is due to limited shared memory.
        {
		for(int tn2 = 0; tn2 < tilenum2; tn2++)
		{

        cacheIndex1 = threadIdx.x;

		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < nx )	//load a tile of image
        {
		cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < nz)
        {

            cache[cacheIndex1][cacheIndex2] = smatrix[((sharesize * tn1) + cacheIndex1) + i * nx + ((sharesize * tn2) + cacheIndex2)* nx * ny];
//			if(tn1 == 0 && tn2 == 0 && i == 0) testmatrix[cacheIndex1*sharesize + cacheIndex2] = threadIdx.y;
			cacheIndex2 += blockDim.y;

        }
            cacheIndex1 += blockDim.x;
        }
        __syncthreads();

		tid = threadIdx.x + threadIdx.y * blockDim.x +  blockIdx.x * blockDim.x * blockDim.y;

	    while(tid < lornum)
	    {
	        x1 = lor.x1[tid];
	        y1 = lor.y1[tid];
	        z1 = lor.z1[tid];
	        x2 = lor.x2[tid];
	        y2 = lor.y2[tid];
	        z2 = lor.z2[tid];
			weight = lor.weight[tid];
			linevalue = 0.;
	
	        ulen2 = powf(x1-x2,2) + powf(y1-y2,2) + powf(z1-z2,2);


            t = ( i * a - y1) / (y2 - y1);

            ox = x1 + t * (x2 - x1);
            oz = z1 + t * (z2 - z1);

			mlxx = max((int)truncf((ox - (SRTWO * torhw ))/a)+1, 0);
			mhxx = min((int)truncf((ox + (SRTWO * torhw ))/a), nx - 1);
			mlzz = max((int)truncf((oz - (SRTWO * torhw ))/a)+1, 0);
			mhzz = min((int)truncf((oz + (SRTWO * torhw ))/a), nz - 1);

			mlxx = max(mlxx, sharesize * tn1);
			mhxx = min(mhxx, sharesize * (tn1 + 1)-1);
			mlzz = max(mlzz, sharesize * tn2);
			mhzz = min(mhzz, sharesize * (tn2 + 1)-1);
			                                                                                                                                                                          
            for(int kx = mlxx; kx <= mhxx; kx++)
            {
                for(int kz = mlzz; kz <= mhzz; kz++)
                {
					float dx = ox - kx*a, dz = oz - kz*a;
					float inner = dx * (x1-x2) + dz * (z1 - z2);
					float dst2 = dx * dx + dz * dz - inner * inner / ulen2;
					float maxdst2 = torhw * torhw;
					//dst = sqrtf(powf(ox-kx*a, 2) + powf(oz-kz*a, 2) - powf((ox-kx*a) * (x1-x2) + (oz-kz*a) * (z1-z2), 2) / ulen2);
					if(dst2 < maxdst2) linevalue += cache[kx -sharesize * tn1 ][ kz - sharesize * tn2] * expf(-dst2/(2.0f * torsgm2)) * weight;
				}
            }

    		lor.linevalue[tid] += linevalue;
    		//finish one tile for one lor
	        tid += blockDim.x * blockDim.y * gridDim.x;
	    }
		//finish one tile for all lors
		__syncthreads();

        }
        }
		//finish all tiles in a slice
    }
	//finish all slices


}

__global__ void ybpro( cudalor lor, float *snmatrix )
{
    int nx = d_imageindex[0], ny =d_imageindex[1], nz = d_imageindex[2];
    float a = d_info[0], torhw = d_info[1], torsgm2 = d_info[2];
	int lornum = d_lorindex[1];

	__shared__ float cache[sharesize][sharesize];
	int tid ;
	int cacheIndex1 , cacheIndex2 ;

    float linevalue , rlinevalue;
    float ulen2, t, ox,oz;
    int mlxx,mhxx,mlzz,mhzz;
    float x1,x2,y1,y2,z1,z2,weight;

	int tilenum1 = (nx + sharesize - 1)/ sharesize, tilenum2 = (nz + sharesize - 1) / sharesize ;


    for(int i=0; i< ny; i++)
    {
        for(int tn1 = 0; tn1 < tilenum1; tn1++)	//for each tile of image. This is due to limited shared memory.
        {
		for(int tn2 = 0; tn2 < tilenum2; tn2++)
		{

		//initialize a tile of image in shared memory to zero.
        cacheIndex1 = threadIdx.x;
		
		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < nx )	
        {
		cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < nz)
        {
            cache[cacheIndex1][cacheIndex2] = 0.0f ;
            cacheIndex2 += blockDim.y;
        }
            cacheIndex1 += blockDim.x;
        }
        __syncthreads();

		tid = threadIdx.x + threadIdx.y * blockDim.x +  blockIdx.x * blockDim.x * blockDim.y;

		//calculate contribution of all lors to the tile of image
	    while(tid < lornum)
	    {
	        linevalue = lor.linevalue[tid];
	
	        if(linevalue < ThreshLineValue) rlinevalue = 0.0f;
	        else rlinevalue = 1.0f / linevalue;
	
	        x1 = lor.x1[tid];
	        y1 = lor.y1[tid];
	        z1 = lor.z1[tid];
	        x2 = lor.x2[tid];
	        y2 = lor.y2[tid];
	        z2 = lor.z2[tid];
			weight = lor.weight[tid];
	
	        ulen2 = powf(x1-x2,2) + powf(y1-y2,2) + powf(z1-z2,2);
	
            t = ( i * a - y1) / (y2 - y1);

            //if(t<0. || t> 1.) continue;

            ox = x1 + t * (x2 - x1);
            oz = z1 + t * (z2 - z1);

			mlxx = max((int)truncf((ox - (SRTWO * torhw ))/a)+1, 0);
			mhxx = min((int)truncf((ox + (SRTWO * torhw ))/a), nx - 1);
			mlzz = max((int)truncf((oz - (SRTWO * torhw ))/a)+1, 0);
			mhzz = min((int)truncf((oz + (SRTWO * torhw ))/a), nz - 1);

			mlxx = max(mlxx, sharesize * tn1);
			mhxx = min(mhxx, sharesize * (tn1 + 1)-1);
			mlzz = max(mlzz, sharesize * tn2);
			mhzz = min(mhzz, sharesize * (tn2 + 1)-1);

            for(int kx = mlxx; kx <= mhxx; kx++)
            {
                for(int kz = mlzz; kz <= mhzz; kz++)
                {
					float dx = ox - kx*a, dz = oz - kz*a;
					float inner = dx * (x1-x2) + dz * (z1 - z2);
					float dst2 = dx * dx + dz * dz - inner * inner / ulen2;
					float maxdst2 = torhw * torhw;
                    //dst = sqrtf(powf(ox-kx*a, 2) + powf(oz-kz*a, 2) - powf((ox-kx*a) * (x1-x2) + (oz-kz*a) * (z1-z2), 2) / ulen2);
                    if(dst2 < maxdst2) atomicAdd(&cache[kx -sharesize * tn1 ][ kz - sharesize * tn2], expf(-dst2/(2.0f * torsgm2)) * rlinevalue * weight) ;

                }
            }
	    
			tid += blockDim.x * blockDim.y * gridDim.x;
		}

        __syncthreads();

		//add the tile of image to global memory
		cacheIndex1 = threadIdx.x;
		
		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < nx )	//load a tile of image
		{
		cacheIndex2 = threadIdx.y;
		while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < nz)
		{
		    atomicAdd(&snmatrix[((sharesize * tn1)  + cacheIndex1)  + i * nx + ((sharesize * tn2) + cacheIndex2) * nx * ny ],cache[cacheIndex1][cacheIndex2]);
		    cacheIndex2 += blockDim.y;
		}
		    cacheIndex1 += blockDim.x;
		}
		__syncthreads();

        }
        }
    }


}


__global__ void zfpro( cudalor lor, float *smatrix)
{
    int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    float a = d_info[0], torhw = d_info[1], torsgm2 = d_info[2];
	int lornum = d_lorindex[2];

	__shared__ float cache[sharesize][sharesize];
	int tid ;
	int cacheIndex1 , cacheIndex2 ;

    float linevalue = 0.;
    float ulen2, t, ox, oy;
    int mlxx,mhxx,mlyy,mhyy;
    float x1,x2,y1,y2,z1,z2,weight;

	int tilenum1 = (nx + sharesize - 1)/ sharesize, tilenum2 = (ny + sharesize - 1) / sharesize ;

    for(int i=0; i< nz; i++)
    {

        for(int tn1 = 0; tn1 < tilenum1; tn1++)	//for each tile of image. This is due to limited shared memory.
        {
		for(int tn2 = 0; tn2 < tilenum2; tn2++)
		{

        cacheIndex1 = threadIdx.x;

		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < nx )	//load a tile of image
        {
		cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < ny)
        {
            cache[cacheIndex1][cacheIndex2] = smatrix[((sharesize * tn1) + cacheIndex1) + ((sharesize * tn2) + cacheIndex2) * nx + i * nx * ny];
            cacheIndex2 += blockDim.y;
        }
            cacheIndex1 += blockDim.x;
        }
      	__syncthreads();

		tid = threadIdx.x + threadIdx.y * blockDim.x +  blockIdx.x * blockDim.x * blockDim.y;

	    while(tid < lornum)
	    {
	        x1 = lor.x1[tid];
	        y1 = lor.y1[tid];
	        z1 = lor.z1[tid];
	        x2 = lor.x2[tid];
	        y2 = lor.y2[tid];
	        z2 = lor.z2[tid];
			weight = lor.weight[tid];
			linevalue = 0.;
	
	        ulen2 = powf(x1-x2,2) + powf(y1-y2,2) + powf(z1-z2,2);

            t = ( i * a - z1) / (z2 - z1);

            oy = y1 + t * (y2 - y1);
            ox = x1 + t * (x2 - x1);

            mlyy = max((int)truncf((oy - (SRTWO * torhw ))/a)+1, 0);
            mhyy = min((int)truncf((oy + (SRTWO * torhw ))/a), ny - 1);
            mlxx = max((int)truncf((ox - (SRTWO * torhw ))/a)+1, 0);
            mhxx = min((int)truncf((ox + (SRTWO * torhw ))/a), nx - 1);

			mlxx = max(mlxx, sharesize * tn1);
			mhxx = min(mhxx, sharesize * (tn1 + 1)-1);
			mlyy = max(mlyy, sharesize * tn2);
			mhyy = min(mhyy, sharesize * (tn2 + 1)-1);

            for(int kx = mlxx; kx <= mhxx; kx++)
            {
                for(int ky = mlyy; ky <= mhyy; ky++)
                {
					float dy = oy - ky*a, dx = ox - kx*a;
					float inner = dy * (y1-y2) + dx * (x1 - x2);
					float dst2 = dy * dy + dx * dx - inner * inner / ulen2;
					float maxdst2 = torhw * torhw;
                    //dst = sqrtf(powf(oy-ky*a, 2) + powf(ox-kx*a, 2) - powf((oy-ky*a) * (y1-y2) + (ox-kx*a) * (x1-x2), 2) / ulen2);
                    if(dst2 < maxdst2) linevalue += cache[kx -sharesize * tn1 ][ ky - sharesize * tn2] * expf(-dst2/(2.0f * torsgm2)) * weight;
                }
            }
	        lor.linevalue[tid] += linevalue;
    	    //finish one tile for one lor
	       	tid += blockDim.x * blockDim.y * gridDim.x;
		}

        __syncthreads();
        }
        }
    }

  
}


__global__ void zbpro( cudalor lor, float *snmatrix )
{
    int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    float a = d_info[0], torhw = d_info[1], torsgm2 = d_info[2];
	int lornum = d_lorindex[2];

	__shared__ float cache[sharesize][sharesize];
	int tid ;
	int cacheIndex1 , cacheIndex2 ;

    float linevalue, rlinevalue;
    float ulen2, t, ox, oy;
    int mlxx,mhxx,mlyy,mhyy;
    float x1,x2,y1,y2,z1,z2,weight;

	int tilenum1 = (nx + sharesize - 1)/ sharesize, tilenum2 = (ny + sharesize - 1) / sharesize ;

    for(int i=0; i< nz; i++)
    {

        for(int tn1 = 0; tn1 < tilenum1; tn1++)	//for each tile of image. This is due to limited shared memory.
        {
		for(int tn2 = 0; tn2 < tilenum2; tn2++)
		{

        cacheIndex1 = threadIdx.x;

		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < nx )	//load a tile of image
        {
		cacheIndex2 = threadIdx.y;
        while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < ny)
        {
            cache[cacheIndex1][cacheIndex2] = 0.0f;
            cacheIndex2 += blockDim.y;
        }
            cacheIndex1 += blockDim.x;
        }
      	__syncthreads();

		tid = threadIdx.x + threadIdx.y * blockDim.x +  blockIdx.x * blockDim.x * blockDim.y;

	    while(tid < lornum)
	    {
	        x1 = lor.x1[tid];
	        y1 = lor.y1[tid];
	        z1 = lor.z1[tid];
	        x2 = lor.x2[tid];
	        y2 = lor.y2[tid];
	        z2 = lor.z2[tid];
			weight = lor.weight[tid];
	        linevalue = lor.linevalue[tid];
	
	        if(linevalue < ThreshLineValue) rlinevalue = 0.0f;
	        else rlinevalue = 1.0f / linevalue;
	
	        ulen2 = powf(x1-x2,2) + powf(y1-y2,2) + powf(z1-z2,2);

            t = ( i * a - z1) / (z2 - z1);

            oy = y1 + t * (y2 - y1);
            ox = x1 + t * (x2 - x1);

            mlyy = max((int)truncf((oy - (SRTWO * torhw ))/a)+1, 0);
            mhyy = min((int)truncf((oy + (SRTWO * torhw ))/a), ny - 1);
            mlxx = max((int)truncf((ox - (SRTWO * torhw ))/a)+1, 0);
            mhxx = min((int)truncf((ox + (SRTWO * torhw ))/a), nx - 1);

			mlxx = max(mlxx, sharesize * tn1);
			mhxx = min(mhxx, sharesize * (tn1 + 1)-1);
			mlyy = max(mlyy, sharesize * tn2);
			mhyy = min(mhyy, sharesize * (tn2 + 1)-1);

            for(int kx = mlxx; kx <= mhxx; kx++)
            {
                for(int ky = mlyy; ky <= mhyy; ky++)
                {
					float dy = oy - ky*a, dx = ox - kx*a;
					float inner = dy * (y1-y2) + dx * (x1 - x2);
					float dst2 = dy * dy + dx * dx - inner * inner / ulen2;
					float maxdst2 = torhw * torhw;
                    //dst = sqrtf(powf(oy-ky*a, 2) + powf(ox-kx*a, 2) - powf((oy-ky*a) * (y1-y2) + (ox-kx*a) * (x1-x2), 2) / ulen2);
                    if(dst2 < maxdst2) atomicAdd(&cache[kx -sharesize * tn1 ][ ky - sharesize * tn2], expf(-dst2/(2.0f * torsgm2)) * rlinevalue * weight) ;
                }
            }

	       	tid += blockDim.x * blockDim.y * gridDim.x;
		}

        __syncthreads();

		//write the tile of image to global memory
		cacheIndex1 = threadIdx.x;
		
		while(cacheIndex1 < sharesize && ((sharesize * tn1) + cacheIndex1) < nx )	//load a tile of image
		{
		cacheIndex2 = threadIdx.y;
		while(cacheIndex2 < sharesize && ((sharesize * tn2) + cacheIndex2) < ny)
		{
		    atomicAdd(&snmatrix[((sharesize * tn1) + cacheIndex1) + ((sharesize * tn2) + cacheIndex2) * nx + i * nx * ny], cache[cacheIndex1][cacheIndex2]);
		    cacheIndex2 += blockDim.y;
		}
		    cacheIndex1 += blockDim.x;
		}
		
		__syncthreads();

        }
        }
    }

  
}


int main(int argc, char* argv[])
{
	vector<string> vpara;
	string line;
	stringstream ss;
	ifstream config ("configRecon.txt");
    int itenum; //number of iterations.
	int startDiff = 1;

	if (config.is_open())
	{
		while ( getline (config,line) )
		{
			vpara=explode(line,' ');
			if(vpara[0]=="FOV") {ss<<vpara[2];ss>>bndry[0];ss.clear();ss<<vpara[3];ss>>bndry[1];ss.clear();ss<<vpara[4];ss>>bndry[2];ss.clear();} else
			if(vpara[0]=="GridSize") {ss<<vpara[2];ss>>a;ss.clear();} else
			if(vpara[0]=="TorHalfWidth") {ss<<vpara[2];ss>>torhw;ss.clear();} else
			if(vpara[0]=="TorSigma") {ss<<vpara[2];ss>>torsgm;ss.clear();} else
			if(vpara[0]=="NumberOfIterations") {ss<<vpara[2];ss>>itenum;ss.clear();} else
			if(vpara[0]=="Regularization") {ss<<vpara[2];ss>>rgl;ss.clear();} else
			if(vpara[0]=="Normalization") {ss<<vpara[2];ss>>norma;ss.clear();} else
			if(vpara[0]=="ThreshNorm") {ss<<vpara[2];ss>>ThreshNorm;ss.clear();} else
			if(vpara[0]=="BetaR") {ss<<vpara[2];ss>>beta;ss.clear();} else
			if(vpara[0]=="BlurR") {ss<<vpara[2];ss>>blur;ss.clear();} else
			if(vpara[0]=="XsigmaRB") {ss<<vpara[2];ss>>bsgm[0];ss.clear();} else
			if(vpara[0]=="YsigmaRB") {ss<<vpara[2];ss>>bsgm[1];ss.clear();} else
			if(vpara[0]=="ZsigmaRB") {ss<<vpara[2];ss>>bsgm[2];ss.clear();} else
			if(vpara[0]=="Weight") {ss<<vpara[2];ss>>wgt;ss.clear();}

		}
		config.close();
	}
	else cout << "Unable to open config file"<<endl;

	cout<<"-------------------------------------------"<<endl;
	cout<<"Input parameters:"<<endl;
	cout<<"FOV: "<<bndry[0]<<" mm x "<<bndry[1]<<" mm x "<<bndry[2]<<" mm"<<endl;
	cout<<"Grid size: "<<a<<" mm"<<endl;
	cout<<"TOR half width: "<<torhw<<" mm"<<endl;
	cout<<"TOR sigma: "<<torsgm<<" mm"<<endl;
	cout<<"Number of iterations: "<<itenum<<endl;
	cout<<"Normalization?: "<<norma<<endl;
	if(norma != 0) cout<<"ThreshNorm: "<<ThreshNorm<<endl;
	cout<<"Regularization?: "<<rgl<<endl;
	cout<<"Include weight?: "<<wgt<<endl;
	if(rgl==1)
	{
		cout<<"Beta for regularization: "<<beta<<endl;
		cout<<"Blur?: "<<blur<<endl;
        if(blur==1)
        {
            cout<<"Xsigma for blur: "<<bsgm[0]<<" mm"<<endl;
            cout<<"Ysigma for blur: "<<bsgm[1]<<" mm"<<endl;
            cout<<"Zsigma for blur: "<<bsgm[2]<<" mm"<<endl;
			for(int i=0; i<3; i++) {rads[i] = bsgm[i] * sqrt(-2. * log (bthresh)); indr[i] = trunc(rads[i]/a);}
			cudaMemcpyToSymbol(d_bsgm, bsgm, 3 * sizeof(float), 0, cudaMemcpyHostToDevice);
			cudaMemcpyToSymbol(d_rads, rads, 3 * sizeof(float), 0, cudaMemcpyHostToDevice); 
			cudaMemcpyToSymbol(d_indr, indr, 3 * sizeof(int), 0, cudaMemcpyHostToDevice); 
        }
		
	}
	else beta = 0.;
	cout<<"-------------------------------------------"<<endl;




	//float stp;
	//vector<float> poimage;  //prior image
	float buff2;
	//float sumpop=0., sumimp=0.;
	ifstream finrgl;
	ifstream finlast;
    string rglname;
    string lastname;

	dim3 threads(threadsperBlock, threadsperBlock);
    nummainaxis = (int*) malloc(3 * sizeof(int));
	//cudalor dev_xlor, dev_ylor, dev_zlor;   //define variables for cudalor.

	msize = ceil(bndry[0] / a) * ceil(bndry[1] / a) *ceil( bndry[2] / a);
	int nx = ceil(bndry[0] / a);
	int ny = ceil(bndry[1] / a);
	int nz = ceil( bndry[2] / a);

//copy fundamental variables to cuda.
	int *temp_imageindex;
	temp_imageindex = (int*) malloc(4 * sizeof(int));
	temp_imageindex[0] = nx;
	temp_imageindex[1] = ny;
	temp_imageindex[2] = nz;
	temp_imageindex[3] = msize;

	float *temp_info;
	temp_info = (float*) malloc(4 * sizeof(float));
	temp_info[0] = a;
	temp_info[1] = torhw;
	temp_info[2] = pow(torsgm,2);	//for higher efficiency in gpu
	temp_info[3] = beta;
	
	//int *dev_imageindex, *dev_lorindex;
	//float *dev_info;
	//cudaMalloc((void**) &dev_imageindex, 4*sizeof(int) );
	//cudaMalloc((void**) &dev_lorindex, 3*sizeof(int) );
	//cudaMalloc((void**) &dev_info, 5*sizeof(float) );

//	cout<<"before copying to constant memory"<<endl;
	cudaMemcpyToSymbol(d_imageindex, temp_imageindex, 4 * sizeof(int), 0, cudaMemcpyHostToDevice);
	//cudaMemcpy(dev_imageindex, temp_imageindex, 4 * sizeof(int), cudaMemcpyHostToDevice);
//	cudaMemcpy(dev_lorindex, nummainaxis, 3 * sizeof(int), cudaMemcpyHostToDevice); //should be after taking the input data
	cudaMemcpyToSymbol(d_info, temp_info, 4 * sizeof(float), 0, cudaMemcpyHostToDevice);
	//cudaMemcpy(dev_info, temp_info,  5* sizeof(float), cudaMemcpyHostToDevice);	//should be after getting maxnorm and beta.

	//delete [] temp_imageindex;
	//delete [] temp_info;



	ifstream fnorm;	//use this if input is normalization image
	string filenorm;
	ofstream fnormout;
	float *lineval;	//used for initializing linevalue of cudalor. As cudaMemset deals with each byte, not each floating point number.
	float *normpara;
	normpara = (float*) malloc(2 * sizeof(float));
	float maxnorm;
	normpara[1] = ThreshNorm;

	if(norma == 1)	//generate sensitivity image from backprojecting all possible LORs or simulation data with weight being 1
	{
		filenorm = argv[4];
		cout<<"Sorting LORs for normalization and copying to device memory......"<<endl;
		preplor(filenorm, 1); //  read lors in the file, sort lors, copy to cud
        cout<<"Normalization: Number of LORs in each main axis (x,y,x): "<<nummainaxis[0]<<" "<<nummainaxis[1]<<" "<<nummainaxis[2]<<endl;
        //cudaMemcpy(dev_lorindex, nummainaxis, 3 * sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(d_lorindex, nummainaxis, 3 * sizeof(int), 0, cudaMemcpyHostToDevice);
		
		normimage = (float*) malloc(msize * sizeof(float));
		cudaMalloc((void**) &dev_normimage, msize*sizeof(float) );
		cudaMemset( dev_normimage, 0.0f, msize*sizeof(float));
	
		timinggpu.StartCounter();
		lineval = (float*) malloc(nummainaxis[0]*sizeof(float));
		for(int ii=0; ii<nummainaxis[0]; ii++) lineval[ii] = 1.0f;
		cudaMemcpy(dev_xlor.linevalue, lineval, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice);
		lineval = (float*) malloc(nummainaxis[1]*sizeof(float));
		for(int ii=0; ii<nummainaxis[1]; ii++) lineval[ii] = 1.0f;
		cudaMemcpy(dev_ylor.linevalue, lineval, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice);  
		lineval = (float*) malloc(nummainaxis[2]*sizeof(float));
		for(int ii=0; ii<nummainaxis[2]; ii++) lineval[ii] = 1.0f;
		cudaMemcpy(dev_zlor.linevalue, lineval, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice);  

		free(lineval);
		timeall.memoryIO += timinggpu.GetCounter();

		xbpro<<<blocksPerGrid, threads>>>(dev_xlor, dev_normimage);
		ybpro<<<blocksPerGrid, threads>>>(dev_ylor, dev_normimage);
		zbpro<<<blocksPerGrid, threads>>>(dev_zlor, dev_normimage);

		timinggpu.StartCounter();
		cudaMemcpy(normimage, dev_normimage, msize*sizeof(float), cudaMemcpyDeviceToHost);
		timeall.memoryIO += timinggpu.GetCounter();
		//cudaMemcpy(ylor.linevalue, dev_ylor.linevalue, nummainaxis[1]*sizeof(float), cudaMemcpyDeviceToHost);
		//cout<<ylor.linevalue[0]<<endl;
		
		fnormout.open("normImage", ios::out | ios::binary);
		maxnorm = 0.0;
		for(int iii=0; iii< msize; iii++)
		{
			fnormout.write( (char*)&normimage[iii], sizeof(float));
			if(maxnorm < normimage[iii]) maxnorm = normimage[iii];
		}
		fnormout.close();
		normpara[0] = maxnorm;
		cudaMemcpyToSymbol(d_norm, normpara,  2 * sizeof(float), 0, cudaMemcpyHostToDevice);
		cout<<"Finish creating normalization image."<<endl;

	}

	else if(norma == 2)	//read sensitivity image from input file.
	{
		filenorm = argv[4];
		cout<<"Reading normalization image......"<<endl;
		fnorm.open(filenorm.c_str(), ios::in | ios::binary);
		if (fnorm.is_open()){
			normimage = (float*) malloc(msize * sizeof(float));
			maxnorm = 0.0;
			for(int iii=0; iii< msize; iii++)
			{
				fnorm.read( (char*)&normimage[iii], sizeof(float));
				if(maxnorm < normimage[iii]) maxnorm = normimage[iii];
			}
		}
		else cout<<"Unable to open normImage file!!"<<endl;
		fnorm.close();
		normpara[0] = maxnorm;

		timinggpu.StartCounter();
		cudaMalloc((void**) &dev_normimage, msize*sizeof(float) );
		cudaMemcpy(dev_normimage, normimage, msize*sizeof(float), cudaMemcpyHostToDevice);
		timeall.memoryIO += timinggpu.GetCounter();

		cudaMemcpyToSymbol(d_norm, normpara, 2 *  sizeof(float), 0, cudaMemcpyHostToDevice);
		cout<<"Finish reading normalization image."<<endl;
	}


	else if(norma == 3)	//generate sensitivity image from backprojecting simulation data with different weight
	{
		filenorm = argv[4];
		cout<<"Sorting LORs for normalization and copying to device memory......"<<endl;
		preplor(filenorm, 1); //  read lors in the file, sort lors, copy to cud
        cout<<"Normalization: Number of LORs in each main axis (x,y,x): "<<nummainaxis[0]<<" "<<nummainaxis[1]<<" "<<nummainaxis[2]<<endl;
        //cudaMemcpy(dev_lorindex, nummainaxis, 3 * sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(d_lorindex, nummainaxis, 3 * sizeof(int), 0, cudaMemcpyHostToDevice);
		
		normimage = (float*) malloc(msize * sizeof(float));
		cudaMalloc((void**) &dev_normimage, msize*sizeof(float) );
		cudaMemset( dev_normimage, 0, msize*sizeof(float)); 

		cudaMemset( dev_xlor.linevalue, 0, nummainaxis[0]*sizeof(float));	// improved normalization method
		cudaMemset( dev_ylor.linevalue, 0, nummainaxis[1]*sizeof(float));
		cudaMemset( dev_zlor.linevalue, 0, nummainaxis[2]*sizeof(float));
		float *normphantom = (float*) malloc(msize * sizeof(float)), *dev_normphantom;
		cudaMalloc((void**) &dev_normphantom, msize * sizeof(float));
		for(int k=0; k< nz; k++){
			for(int j=0; j< ny; j++){
				for(int i=0; i< nx; i++){
					float cx = i * a - (bndry[0]/2 - 0.5 * a);
					float cy = j * a - (bndry[1]/2 - 0.5 * a);
					float cz = k * a - (bndry[2]/2 - 0.5 * a);
					float temp = 0.0f;
					if(cx * cx + cy * cy < 3025.0 && cz >= -63 && cz <= 63) temp = 0.01f;
					//if(cz >= -63 && cz <= 63) temp = 0.001f;
					normphantom[i + j * nx + k * nx * ny] = temp;

				}
			}
		}
		cudaMemcpy(dev_normphantom, normphantom, msize*sizeof(float), cudaMemcpyHostToDevice);
		
		xfpro<<<blocksPerGrid, threads>>>(dev_xlor, dev_normphantom);
		yfpro<<<blocksPerGrid, threads>>>(dev_ylor, dev_normphantom);
		zfpro<<<blocksPerGrid, threads>>>(dev_zlor, dev_normphantom);
		xbpro<<<blocksPerGrid, threads>>>(dev_xlor, dev_normimage);
		ybpro<<<blocksPerGrid, threads>>>(dev_ylor, dev_normimage);
		zbpro<<<blocksPerGrid, threads>>>(dev_zlor, dev_normimage);

		timinggpu.StartCounter();
		cudaMemcpy(normimage, dev_normimage, msize*sizeof(float), cudaMemcpyDeviceToHost);
		timeall.memoryIO += timinggpu.GetCounter();
		//cudaMemcpy(ylor.linevalue, dev_ylor.linevalue, nummainaxis[1]*sizeof(float), cudaMemcpyDeviceToHost);
		//cout<<ylor.linevalue[0]<<endl;
		
		fnormout.open("normImage", ios::out | ios::binary);
		maxnorm = 0.0;
		for(int iii=0; iii< msize; iii++)
		{
			fnormout.write( (char*)&normimage[iii], sizeof(float));
			if(maxnorm < normimage[iii]) maxnorm = normimage[iii];
		}
		fnormout.close();
		normpara[0] = maxnorm;
		cudaMemcpyToSymbol(d_norm, normpara,  2 * sizeof(float), 0, cudaMemcpyHostToDevice);
		cout<<"Finish creating normalization image."<<endl;

	}

	else if(norma == 0) {}

	else cout<<"Unkown indicator for normalization option!!"<<endl;




// open file that contains input lors for image reconstruction. This should be after normalization file read.
	//ifstream fin;
	string filein=argv[1];
	//fin.open(filein.c_str());
	cout<<"Sorting LORs and copying to device memory......"<<endl;
    preplor(filein, 0); //  read lors in the file, sort lors, copy to cuda
	cout<<"Finish sorting and copying."<<endl;
	//cudaMemcpy(dev_lorindex, nummainaxis, 3 * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_lorindex, nummainaxis, 3 * sizeof(int), 0, cudaMemcpyHostToDevice);
	//fin.close();


	timing.StartCounter();
	smatrix = (float*) malloc(msize * sizeof(float));
	snmatrix = (float*) malloc(msize * sizeof(float));
	poimage = (float*) malloc(msize * sizeof(float));
	lastimage = (float*) malloc(msize * sizeof(float));
	bmatrix = (float*) malloc(msize * sizeof(float));
    //deri = (float*) malloc(msize * sizeof(float));
	allweight = (float*) malloc(msize * sizeof(float));
	timeall.memoryIO += timing.GetCounter();

//	cout<<"after feeding data to cudalor"<<endl;
	if(rgl == 1){
		rglname=argv[3];
		finrgl.open(rglname.c_str(), ios::in | ios::binary);
	
		//read prior image into ram
        //sumpop = 0.0;
		if (finrgl.is_open()){
			for(int i=0; i< msize; i++)
			{
				finrgl.read((char *)&poimage[i],sizeof(buff2));
				//poimage[i] = buff2;
				//sumpop += buff2;
			}
		}
		else cout<<"Unable to open rgl image!!"<<endl;
		finrgl.close();

		lastname=argv[5];
		finlast.open(lastname.c_str(), ios::in | ios::binary);
		if (finlast.is_open()){
			for(int i=0; i< msize; i++)
			{
				finlast.read((char *)&lastimage[i],sizeof(buff2));
			}
		}
		else cout<<"Unable to open rgl image!!"<<endl;
		finlast.close();
	}


//re-define beta, beta_new = beta * A / 2. Also initialize smatrix.
	if(norma == 0) {
		beta = beta * float(numline)/msize / 2.0f;
		for(int i=0; i<msize; i++) smatrix[i] = float(numline)/msize;
	}
	else {
		float sumNormimage = 0.0f;
		for(int i=0; i< msize; i++) sumNormimage += normimage[i];
		sumNormimage = sumNormimage / maxnorm;
		beta = beta * float(numline)/sumNormimage / 2.0f;
		if(startDiff == true){	
			for(int i=0; i<msize; i++) {
				smatrix[i] = lastimage[i];
			}
		}
		else {
			for(int i=0; i<msize; i++) 
				smatrix[i] = float(numline)/sumNormimage;
		}
	}

	temp_info[3] = beta;
	cudaMemcpyToSymbol(d_info, temp_info, 4 * sizeof(float), 0, cudaMemcpyHostToDevice);
	
	//set aves for scaling. In this case, it remains unchanged across iterations.
	float allave = smatrix[0];
	cudaMemcpyToSymbol(aves, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);
	float hostAve = allave; //for storing the value of A.


//	cout<<"after cuda memory copy"<<endl;

	float *dev_smatrix, *dev_snmatrix, *dev_poimage, *dev_deri, *dev_bmatrix, *dev_allweight;

	timinggpu.StartCounter();
	cudaMalloc((void**) &dev_smatrix, msize*sizeof(float) );
	cudaMalloc((void**) &dev_snmatrix, msize*sizeof(float) );
	if(rgl == 1)
	{
		cudaMalloc((void**) &dev_poimage, msize*sizeof(float) );
		//cudaMalloc((void**) &dev_deri, msize*sizeof(float) );
		cudaMalloc((void**) &dev_bmatrix, msize*sizeof(float) );
		cudaMemcpy( dev_poimage, poimage, msize*sizeof(float), cudaMemcpyHostToDevice );
		if(blur == 1)
		{
			cudaMalloc((void**) &dev_bmatrix, msize*sizeof(float) );
			cudaMalloc((void**) &dev_allweight, msize*sizeof(float) );
		}
	}
	timeall.memoryIO += timinggpu.GetCounter();


	float rgimg, sumi, sumde2, sumde;
	//int idxy;
	//int ci,cj,ck,li,hi,lj,hj,lk,hk;

	//float senratio = 0.0;

//	cout<<"before iteration"<<endl;

//	float *test, *dev_test;
//	test = (float*) malloc( sizeof(float));
//	
//	cudaMalloc( (void**)&dev_test, sizeof(float) );

//	float *testmatrix, *dev_testmatrix;
//	cudaMalloc((void**) &dev_testmatrix, sharesize * sharesize *sizeof(float) );
//	cudaMemset( dev_testmatrix, 0.0f, sharesize * sharesize * sizeof(float));
//	testmatrix = (float *) malloc(sharesize * sharesize * sizeof(float));


	cout<<"Number of LORs in each main axis (x,y,x): "<<nummainaxis[0]<<" "<<nummainaxis[1]<<" "<<nummainaxis[2]<<endl;

	float *gave, *dev_gave;
    int *gnumave, *dev_gnumave;
	gave = (float*) malloc(blocksPerGrid * sizeof(float));
	cudaMalloc((void**) &dev_gave, blocksPerGrid*sizeof(float) );
	gnumave = (int*) malloc(blocksPerGrid * sizeof(int));
	cudaMalloc((void**) &dev_gnumave, blocksPerGrid*sizeof(int) );
    int allnumave = 0;

//if using regularization, calculate average value of voxels and substitute too large voxels value.
	ofstream fpriorout;
	float hostAvep;
    if(rgl == 1){
	allave = 0.0;
	calave<<<blocksPerGrid, reducsize>>>(dev_poimage, dev_gave);
	cudaMemcpy(gave, dev_gave, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);
	for(int jj=0; jj< blocksPerGrid; jj++)  allave += gave[jj];
	allave /= msize;
	cudaMemcpyToSymbol(avep, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);

	cout<<"Prior image average value, before filter: "<<allave<<endl;

	cudaMemcpyToSymbol(aveunf, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);
	allave = 0.0;
	allnumave = 0;
	calavewithfilter<<<blocksPerGrid, reducsize>>>(dev_poimage, dev_gave,dev_gnumave);
	gpuErrchk(cudaPeekAtLastError());gpuErrchk(cudaDeviceSynchronize());
	cudaMemcpy(gave, dev_gave, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(gnumave, dev_gnumave, blocksPerGrid*sizeof(int), cudaMemcpyDeviceToHost);
	for(int jj=0; jj< blocksPerGrid; jj++)  {allave += gave[jj]; allnumave += gnumave[jj];}
	allave /= allnumave;
	cudaMemcpyToSymbol(avep, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);
	hostAvep = allave;	//for storing A_P

	cout<<"Prior image average value, after filter: "<<allave<<" "<<allnumave<<endl;
	cudaMemcpy( poimage, dev_poimage, msize*sizeof(float), cudaMemcpyDeviceToHost );
		fpriorout.open("priorImage", ios::out | ios::binary);
		for(int iii=0; iii< msize; iii++)
		{
			fpriorout.write( (char*)&poimage[iii], sizeof(float));
		}
		fpriorout.close();
    }

	vector<double> logLikelihood, logR;	//value of objective functions in all iterations
//start iterations for image reconstruction
	for(int ij=0; ij<itenum; ij++){

		double templogLikelihood = 0.0;
		double templogR = 0.0;

		cout<<"Starting "<<ij<<" iteration."<<endl;

//        if(rgl == 1)
//        {
//            //sumimp = 0.;
//            //for(int jj=0; jj<msize; jj++) sumimp += smatrix[jj];
//			if(blur == 1) 
//			{
//				gaussianblur(nx, ny, nz);
//				sumb = 0.;
//				for(int jj=0; jj<msize; jj++) sumb += bmatrix[jj];
//			}
//		}

		timinggpu.StartCounter();
		cout<<"pre1 "<<smatrix[0]<<" "<<smatrix[100]<<endl;
		cout<<"pre2 "<<snmatrix[0]<<" "<<snmatrix[100]<<endl;
		cout<<"pre3 "<<lastimage[0]<<" "<<lastimage[100]<<endl;
		cudaMemcpy( dev_smatrix, smatrix, msize*sizeof(float), cudaMemcpyHostToDevice );
		cudaMemset( dev_snmatrix, 0, msize*sizeof(float));
		cudaMemset( dev_xlor.linevalue, 0, nummainaxis[0]*sizeof(float));
		cudaMemset( dev_ylor.linevalue, 0, nummainaxis[1]*sizeof(float)); 
		cudaMemset( dev_zlor.linevalue, 0, nummainaxis[2]*sizeof(float)); 
        //cudaMemset( dev_deri, 0.0f, msize*sizeof(float));
		timeall.memoryIO += timinggpu.GetCounter();


		timinggpu.StartCounter();
		xfpro<<<blocksPerGrid, threads>>>(dev_xlor, dev_smatrix);
		timeall.txforward += timinggpu.GetCounter();

		timinggpu.StartCounter();
		yfpro<<<blocksPerGrid, threads>>>(dev_ylor, dev_smatrix);
		timeall.tyforward += timinggpu.GetCounter();

		timinggpu.StartCounter();
		zfpro<<<blocksPerGrid, threads>>>(dev_zlor, dev_smatrix);
		timeall.tzforward += timinggpu.GetCounter();

		timinggpu.StartCounter();
		xbpro<<<blocksPerGrid, threads>>>(dev_xlor, dev_snmatrix);
		timeall.txbackward += timinggpu.GetCounter();

		timinggpu.StartCounter();
		ybpro<<<blocksPerGrid, threads>>>(dev_ylor, dev_snmatrix);
		timeall.tybackward += timinggpu.GetCounter();

		timinggpu.StartCounter();
		zbpro<<<blocksPerGrid, threads>>>(dev_zlor, dev_snmatrix);
        timeall.tzbackward += timinggpu.GetCounter();


	//	testf<<<blocksPerGrid, threads>>>(dev_test);
	//	testf<<<1,1>>>(dev_test);

	//	cudaMemcpy(snmatrix, dev_snmatrix, msize*sizeof(float), cudaMemcpyDeviceToHost);
	//	cudaMemcpy(ylor.linevalue, dev_ylor.x1, nummainaxis[1]*sizeof(float), cudaMemcpyDeviceToHost);
	//	cout<<ylor.linevalue[0]<<" "<<ylor.linevalue[1]<<" "<<ylor.linevalue[2]<<" "<<ylor.linevalue[3]<<" "<<ylor.linevalue[4]<<" "<<ylor.linevalue[5]<<endl;

	//	cudaMemcpy(xlor.linevalue, dev_xlor.linevalue, nummainaxis[0]*sizeof(float), cudaMemcpyDeviceToHost);
	//	string filein=argv[1];
	//	cout<<xlor.linevalue[0]<<" "<<xlor.linevalue[1]<<endl;
		
	//	cudaMemcpy(testmatrix, dev_testmatrix, sharesize * sharesize *sizeof(float), cudaMemcpyDeviceToHost);

	//	for(int ii=0; ii<sharesize; ii++)
	//	{
	//	for(int ij = 0; ij<sharesize; ij++)
	//	{cout<<testmatrix[ii*sharesize + ij]<<" ";}
	//	cout<<endl;
	//	}

	//	cudaMemcpy(test, dev_test, sizeof(float), cudaMemcpyDeviceToHost);
	//	cout<<"Test data is "<<test[0]<<endl;

		sumde2 = 0.0;

		//no regularization
		if(rgl == 0)
		{
			if(norma == 0){
				timinggpu.StartCounter();
				calnewmatrix000<<<blocksPerGrid, threads>>>(dev_snmatrix, dev_smatrix);//snmatrix[jj] = smatrix[jj] * snmatrix[jj];
				timeall.tpostimageprocess += timinggpu.GetCounter();
			}
			else{ 
				timinggpu.StartCounter();
				calnewmatrix100<<<blocksPerGrid, threads>>>(dev_snmatrix, dev_smatrix, dev_normimage);
				timeall.tpostimageprocess += timinggpu.GetCounter();
			}
		}
		
		//regularization, no blur
		else if(rgl == 1 && blur == 0)
		{
            //calculate average value without filter
			//allave = 0.0;
			//calave<<<blocksPerGrid, reducsize>>>(dev_smatrix, dev_gave);
			//cudaMemcpy(gave, dev_gave, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);
			//for(int jj=0; jj< blocksPerGrid; jj++)  allave += gave[jj];
			//allave /= msize;
			//cudaMemcpyToSymbol(aves, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);

            ////calculate average value with filter
			//cudaMemcpyToSymbol(aveunf, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);
			//allave = 0.0;
			//allnumave = 0;
			//calavewithfilter<<<blocksPerGrid, reducsize>>>(dev_smatrix, dev_gave,dev_gnumave);
			//cudaMemcpy(gave, dev_gave, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(gnumave, dev_gnumave, blocksPerGrid*sizeof(int), cudaMemcpyDeviceToHost);
			//for(int jj=0; jj< blocksPerGrid; jj++)  {allave += gave[jj]; allnumave += gnumave[jj];}
			//allave /= allnumave;
			//cudaMemcpyToSymbol(aves, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);



			if(norma == 0){
				timinggpu.StartCounter();
				calnewmatrix010<<<blocksPerGrid, threads>>>(dev_snmatrix, dev_smatrix, dev_poimage);//Error check: gpuErrchk(cudaPeekAtLastError());gpuErrchk(cudaDeviceSynchronize());}
				timeall.tpostimageprocess += timinggpu.GetCounter();
			}
			else{ 
				timinggpu.StartCounter();
				calnewmatrix110<<<blocksPerGrid, threads>>>(dev_snmatrix, dev_smatrix, dev_normimage, dev_poimage);
				timeall.tpostimageprocess += timinggpu.GetCounter();
			}
		}

		//regularizatin, blur
		else if(rgl == 1 && blur == 1)
		{
			cudaMemset( dev_bmatrix, 0, msize*sizeof(float));
			cudaMemset( dev_allweight, 0, msize*sizeof(float));
			gpublur<<<blocksPerGrid, threads>>>(dev_smatrix, dev_bmatrix, dev_allweight);	//blur image from last iteration

            //calculate average value for the blurred image without filter
			//allave = 0.0;
			//calave<<<blocksPerGrid, reducsize>>>(dev_bmatrix, dev_gave);
			//cudaMemcpy(gave, dev_gave, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);
			//for(int jj=0; jj< blocksPerGrid; jj++)  allave += gave[jj];
			//allave /= msize;
			//cudaMemcpyToSymbol(aves, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);

            //calculate average value with filter
			//cudaMemcpyToSymbol(aveunf, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);
			//allave = 0.0;
			//allnumave = 0;
			//calavewithfilter<<<blocksPerGrid, reducsize>>>(dev_bmatrix, dev_gave,dev_gnumave);
			//cudaMemcpy(gave, dev_gave, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(gnumave, dev_gnumave, blocksPerGrid*sizeof(int), cudaMemcpyDeviceToHost);
			//for(int jj=0; jj< blocksPerGrid; jj++)  {allave += gave[jj]; allnumave += gnumave[jj];}
			//allave /= allnumave;
			//cudaMemcpyToSymbol(aves, &allave, sizeof(float), 0, cudaMemcpyHostToDevice);

			//calculate new image for this iteration
			if(norma == 0){
				timinggpu.StartCounter();
				calnewmatrix011<<<blocksPerGrid, threads>>>(dev_snmatrix, dev_smatrix, dev_poimage, dev_bmatrix, dev_allweight);//Error check: gpuErrchk(cudaPeekAtLastError());gpuErrchk(cudaDeviceSynchronize());}
				timeall.tpostimageprocess += timinggpu.GetCounter();
			}
			else{
				timinggpu.StartCounter();
				calnewmatrix111<<<blocksPerGrid, threads>>>(dev_snmatrix, dev_smatrix, dev_normimage, dev_poimage, dev_bmatrix, dev_allweight);
				timeall.tpostimageprocess += timinggpu.GetCounter();
			}

		}

		else cout<<"Unknown indentifier for regularization or blur!!!"<<endl;

		timinggpu.StartCounter();
		cudaMemcpy(snmatrix, dev_snmatrix, msize*sizeof(float), cudaMemcpyDeviceToHost);
		timeall.memoryIO += timinggpu.GetCounter();

		cout<<"pre2 "<<snmatrix[0]<<" "<<snmatrix[100]<<endl;
		cout<<"Finish "<<ij<<" iteration."<<endl;

//write new image to file.
		ostringstream convert;
		convert<<(ij+1);
		ofstream fout;
		string fileout=argv[2];
		fileout.append(convert.str());
		fout.open(fileout.c_str(), ios::out | ios::binary);
		if (fout.is_open()){
			for(int iii=0; iii< msize; iii++)
			{
				fout.write( (char*)&snmatrix[iii], sizeof(float));
			}
		}
		else cout<<"Unable to write image to file!!"<<endl;

		fout.close();

		//calculate objective function values
		double *gloglike, *dev_gloglike;
		gloglike = (double*) malloc(blocksPerGrid * sizeof(double));
		cudaMalloc((void**) &dev_gloglike, blocksPerGrid*sizeof(double) );
		calLogLike<<<blocksPerGrid, reducsize>>>(dev_xlor.linevalue, dev_gloglike, nummainaxis[0]);
		cudaMemcpy(gloglike, dev_gloglike, blocksPerGrid*sizeof(double), cudaMemcpyDeviceToHost);
		for(int iobj = 0; iobj < blocksPerGrid; iobj++) templogLikelihood += gloglike[iobj];
	
		calLogLike<<<blocksPerGrid, reducsize>>>(dev_ylor.linevalue, dev_gloglike, nummainaxis[1]);
		cudaMemcpy(gloglike, dev_gloglike, blocksPerGrid*sizeof(double), cudaMemcpyDeviceToHost);
		for(int iobj = 0; iobj < blocksPerGrid; iobj++) templogLikelihood += gloglike[iobj];
	
		calLogLike<<<blocksPerGrid, reducsize>>>(dev_zlor.linevalue, dev_gloglike, nummainaxis[2]);
		cudaMemcpy(gloglike, dev_gloglike, blocksPerGrid*sizeof(double), cudaMemcpyDeviceToHost);
		for(int iobj = 0; iobj < blocksPerGrid; iobj++) templogLikelihood += gloglike[iobj];


		calLogLikeS<<<blocksPerGrid, reducsize>>>(dev_smatrix, dev_normimage, dev_gloglike, msize, norma);
		cudaMemcpy(gloglike, dev_gloglike, blocksPerGrid*sizeof(double), cudaMemcpyDeviceToHost);
		for(int iobj = 0; iobj < blocksPerGrid; iobj++) templogLikelihood += gloglike[iobj];


		if(rgl == 1 && blur == 0) calLogR<<<blocksPerGrid, reducsize>>>(dev_smatrix, dev_poimage, dev_gloglike, msize);
		else if(rgl == 1 && blur == 1) calLogR<<<blocksPerGrid, reducsize>>>(dev_bmatrix, dev_poimage, dev_gloglike, msize);
		if(rgl == 1) {
			cudaMemcpy(gloglike, dev_gloglike, blocksPerGrid*sizeof(double), cudaMemcpyDeviceToHost);
			for(int iobj = 0; iobj < blocksPerGrid; iobj++) templogR += gloglike[iobj];	
		}

		//cudaMemcpy(xlor.linevalue, dev_xlor.linevalue, nummainaxis[0]*sizeof(float), cudaMemcpyDeviceToHost);
		//cudaMemcpy(ylor.linevalue, dev_ylor.linevalue, nummainaxis[1]*sizeof(float), cudaMemcpyDeviceToHost);
		//cudaMemcpy(zlor.linevalue, dev_zlor.linevalue, nummainaxis[2]*sizeof(float), cudaMemcpyDeviceToHost);
		//for(int iobj = 0; iobj < nummainaxis[0]; iobj++) {
		//	if (xlor.linevalue[iobj] > ThreshLineValue ) templogLikelihood += log(xlor.linevalue[iobj]);
		//}
		//for(int iobj = 0; iobj < nummainaxis[1]; iobj++) {
		//	if (ylor.linevalue[iobj] > ThreshLineValue) templogLikelihood += log(ylor.linevalue[iobj]);
		//}
		//for(int iobj = 0; iobj < nummainaxis[2]; iobj++) {
		//	if (zlor.linevalue[iobj] > ThreshLineValue) templogLikelihood += log(zlor.linevalue[iobj]);
		//}

		//for(int iobj = 0; iobj < msize; iobj++){
		//	if(norm == 0) templogLikelihood -= smatrix[iobj];
		//	else {
		//		float tempSen = normimage[iobj] / maxnorm;
		//		if(tempSen < ThreshNorm) tempSen = ThreshNorm;
		//		templogLikelihood -= tempSen * smatrix[iobj];
		//	}
		//}
		//
		//if(rgl == 1 && blur == 0) for(int iobj = 0; iobj < msize; iobj++) templogR -= pow(smatrix[iobj]/hostAve - poimage[iobj]/hostAvep,2);
		//else if(rgl == 1 && blur == 1) {
		//	cudaMemcpy(bmatrix, dev_bmatrix, msize * sizeof(float), cudaMemcpyDeviceToHost);
		//	for(int iobj = 0; iobj < msize; iobj++) templogR -= pow(bmatrix[iobj]/hostAve - poimage[iobj]/hostAvep,2);
		//}

		templogR *= beta;
		
		logLikelihood.push_back(templogLikelihood);
		logR.push_back(templogR);

		//prepare for next iteration
		for(int iii=0; iii< msize; iii++)
		{
			smatrix[iii] = snmatrix[iii];
			snmatrix[iii] = 0.;
		}


	}

	ofstream fObjFunc ("ObjectiveFuncValue.txt");
	if(fObjFunc.is_open()){
		for (int i=0; i< itenum; i++) fObjFunc << i << " "<< logLikelihood[i] << " " << logR[i] << " " << logLikelihood[i] + logR[i] << endl;
	}
	else cout<< "Can not open ObjectiveFuncValue.txt!!" <<endl;
	fObjFunc.close();

	timeall.printvalue();	//print out timing information about cuda execution.

	//cout<< "Total number of voxel lor distance calculation is "<< numcal << endl;
	//cout<< "Total number of voxel lor kernel calculation is " << numcal2 << endl;
	cudaFree(dev_xlor.x1);
	cudaFree(dev_xlor.y1);
	cudaFree(dev_xlor.z1);
	cudaFree(dev_xlor.x2);
	cudaFree(dev_xlor.y2);
	cudaFree(dev_xlor.z2);
    cudaFree(dev_ylor.x1);
    cudaFree(dev_ylor.y1);
    cudaFree(dev_ylor.z1);
    cudaFree(dev_ylor.x2);
    cudaFree(dev_ylor.y2);
    cudaFree(dev_ylor.z2);
    cudaFree(dev_zlor.x1);
    cudaFree(dev_zlor.y1);
    cudaFree(dev_zlor.z1);
    cudaFree(dev_zlor.x2);
    cudaFree(dev_zlor.y2);
    cudaFree(dev_zlor.z2);
	cudaFree(dev_xlor.linevalue);
	cudaFree(dev_ylor.linevalue);
	cudaFree(dev_zlor.linevalue);
	cudaFree(dev_smatrix);
	cudaFree(dev_snmatrix);
	cudaFree(dev_poimage);
	cudaFree(dev_deri);
	cudaFree(dev_bmatrix);
    free(xlor.x1);
    free(xlor.y1);
    free(xlor.z1);
    free(xlor.x2);
    free(xlor.y2);
    free(xlor.z2);
    free(ylor.x1);
    free(ylor.y1);
    free(ylor.z1);
    free(ylor.x2);
    free(ylor.y2);
    free(ylor.z2);
    free(zlor.x1);
    free(zlor.y1);
    free(zlor.z1);
    free(zlor.x2);
    free(zlor.y2);
    free(zlor.z2);	
	free(xlor.linevalue);
	free(ylor.linevalue);
	free(zlor.linevalue);
	free(smatrix);
	free(snmatrix);
	free(bmatrix);
	free(poimage);
	free(deri);
	free(allweight);
	return 0;
}

vector<string> explode(string s, char c)
{
	string buff="";
	vector<string> v;
	char n;
	
	for(unsigned i=0; i<s.length(); ++i)
	{
		n=s.at(i);
		if(n != c) buff+=n; else
		if(n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if(buff != "") v.push_back(buff);
	
	return v;
}



int gaussianblur(int nx, int ny, int nz)
{
    float sumweight;

	int ci,cj,ck,li,hi,lj,hj,lk,hk;
	for(int i=0; i< msize;i++ )
	{
		bmatrix[i] = 0.;
		ci = (i % (nx*ny)) % nx;
		cj = (i % (nx*ny)) / nx;
		ck = i / (nx*ny);
		li = max(0, ci - indr[0]);
		lj = max(0, cj - indr[1]);
		lk = max(0, ck - indr[2]);
		hi = min(nx - 1, ci + indr[0]);
		hj = min(ny - 1, cj + indr[1]);
		hk = min(nz - 1, ck + indr[2]);
        sumweight = 0.0;

		for(int ii= li; ii<= hi; ii++)
		{
			for(int jj = lj; jj<= hj; jj++)
			{
				for(int kk = lk; kk<= hk; kk++)
				{
					if((pow((ii-ci)/rads[0],2) + pow((jj-cj)/rads[1],2) + pow((kk-ck)/rads[2],2)) * pow(a,2) < 1.0) 
					{
						bmatrix[i] += smatrix[ii + jj * nx + kk * nx * ny] * exp(-(pow(ii-ci,2)/(2. * pow(bsgm[0],2)) + pow(jj-cj,2)/(2. * pow(bsgm[1],2)) + pow(kk-ck,2)/(2. * pow(bsgm[2],2))) * pow(a,2));
						sumweight += exp(-(pow(ii-ci,2)/(2. * pow(bsgm[0],2)) + pow(jj-cj,2)/(2. * pow(bsgm[1],2)) + pow(kk-ck,2)/(2. * pow(bsgm[2],2))) * pow(a,2));
					}

				}
			}
		}
		allweight[i] = sumweight;
		bmatrix[i] = bmatrix[i] / sumweight;
		//for(int j=li; j<=hi; j++) bmatrix[i] += smatrix[j + cj * nx + ck * nx * ny] * exp(-pow(j-ci,2) * pow(a,2)/(2. * pow(bsgm,2)));
	}
	return 0;
}

//function that read lor from fin, sort lor, and copy lor to cuda
int preplor(string filein, int senmap)
{
	ifstream fin;
	fin.open(filein.c_str(), ios::in | ios::binary);

	vector<lor> alllor;		//matrix for all lor
	numline = 0;
	string line;
	
	nummainaxis[0] = 0;
	nummainaxis[1] = 0;
	nummainaxis[2] = 0;

	timing.StartCounter();
	if (fin.is_open()){

		while ( !fin.eof() )
		{
			numline += 1;
			float coordlor[6];
			fin.read((char*)coordlor, 6 * sizeof(float));
			lor bufflor;
			if(wgt == 1 && senmap == 0) {
				float aa,bb,cc, delta, ts1, ts2, ulen2;
				aa = pow(coordlor[0]-coordlor[3],2) + pow(coordlor[1]-coordlor[4],2);
				bb = 2. * (coordlor[0] * (coordlor[3] - coordlor[0]) + coordlor[1] * (coordlor[4] - coordlor[1]));
				cc = pow(coordlor[0], 2) + pow(coordlor[1], 2) - pow(CylRadius,2);
				ulen2 = pow(coordlor[0]-coordlor[3],2) + pow(coordlor[1]-coordlor[4],2) + pow(coordlor[2]-coordlor[5],2);
				delta = pow(bb,2) - 4.* aa * cc;
		
				if(delta <= 0.) bufflor.weight = 1.0;
				else
				{
					ts1 = (-bb + sqrt(delta)) / (2. * aa);
					ts2 = (-bb - sqrt(delta)) / (2. * aa);
					if(coordlor[5] == coordlor[2] && (coordlor[2] > CylHeight / 2.0 || coordlor[2] < -CylHeight / 2.0)) ts1 = ts2;
					else {
						float zz1 = coordlor[2] + ts1 * (coordlor[5] - coordlor[2]);
						if(zz1 > CylHeight / 2.0) ts1 = (CylHeight/2.0 - coordlor[2]) / (coordlor[5] - coordlor[2]);
						else if(zz1 < -CylHeight/2.0) ts1 = (-CylHeight/2.0 - coordlor[2]) / (coordlor[5] - coordlor[2]);
						float zz2 = coordlor[2] + ts2 * (coordlor[5] - coordlor[2]);
						if(zz2 > CylHeight / 2.0) ts2 = (CylHeight/2.0 - coordlor[2]) / (coordlor[5] - coordlor[2]);
						else if(zz2 < -CylHeight/2.0) ts2 = (-CylHeight/2.0 - coordlor[2]) / (coordlor[5] - coordlor[2]);
					}
			
					bufflor.weight = sqrt(ulen2) * abs(ts1 - ts2);
					bufflor.weight = exp(-(MU * bufflor.weight));

				}

			}
			else bufflor.weight = 1.0;
			bufflor.x1 = coordlor[0] - (-bndry[0]/2. + 0.5 * a);
			bufflor.y1 = coordlor[1] - (-bndry[1]/2. + 0.5 * a);
			bufflor.z1 = coordlor[2] - (-bndry[2]/2. + 0.5 * a);
			bufflor.x2 = coordlor[3] - (-bndry[0]/2. + 0.5 * a);
			bufflor.y2 = coordlor[4] - (-bndry[1]/2. + 0.5 * a);
			bufflor.z2 = coordlor[5] - (-bndry[2]/2. + 0.5 * a);
	
			if(abs(bufflor.x1-bufflor.x2) >= abs(bufflor.y1-bufflor.y2) && abs(bufflor.x1-bufflor.x2) >= abs(bufflor.z1-bufflor.z2)) {bufflor.mainaxis = 0; nummainaxis[0] += 1;}
			else if(abs(bufflor.y1-bufflor.y2) >= abs(bufflor.x1-bufflor.x2) && abs(bufflor.y1-bufflor.y2) >= abs(bufflor.z1-bufflor.z2)) {bufflor.mainaxis = 1; nummainaxis[1] += 1;}
			else if(abs(bufflor.z1-bufflor.z2) >= abs(bufflor.x1-bufflor.x2) && abs(bufflor.z1-bufflor.z2) >= abs(bufflor.y1-bufflor.y2)) {bufflor.mainaxis = 2; nummainaxis[2] += 1;}
			else cout<<"Cannot fing the main axis!!"<<endl;
	
			alllor.push_back(bufflor);
		}
	}
	else cout<<"Unable to open input lor file!!"<<endl;

	fin.close();
	timeall.lorsorting += timing.GetCounter();
//	cout<<"before lor malloc"<<endl;

	timing.StartCounter();
	xlor.x1 = (float*) malloc(nummainaxis[0] * sizeof(float)); 
	xlor.y1 = (float*) malloc(nummainaxis[0] * sizeof(float));
	xlor.z1 = (float*) malloc(nummainaxis[0] * sizeof(float));
	xlor.x2 = (float*) malloc(nummainaxis[0] * sizeof(float));
	xlor.y2 = (float*) malloc(nummainaxis[0] * sizeof(float));
	xlor.z2 = (float*) malloc(nummainaxis[0] * sizeof(float));
	xlor.linevalue = (float*) malloc(nummainaxis[0] * sizeof(float));
	xlor.weight = (float*) malloc(nummainaxis[0] * sizeof(float));

    ylor.x1 = (float*) malloc(nummainaxis[1] * sizeof(float)); 
    ylor.y1 = (float*) malloc(nummainaxis[1] * sizeof(float));
    ylor.z1 = (float*) malloc(nummainaxis[1] * sizeof(float));
    ylor.x2 = (float*) malloc(nummainaxis[1] * sizeof(float));
    ylor.y2 = (float*) malloc(nummainaxis[1] * sizeof(float));
    ylor.z2 = (float*) malloc(nummainaxis[1] * sizeof(float));
	ylor.linevalue = (float*) malloc(nummainaxis[1] * sizeof(float));
	ylor.weight = (float*) malloc(nummainaxis[1] * sizeof(float));

    zlor.x1 = (float*) malloc(nummainaxis[2] * sizeof(float)); 
    zlor.y1 = (float*) malloc(nummainaxis[2] * sizeof(float));
    zlor.z1 = (float*) malloc(nummainaxis[2] * sizeof(float));
    zlor.x2 = (float*) malloc(nummainaxis[2] * sizeof(float));
    zlor.y2 = (float*) malloc(nummainaxis[2] * sizeof(float));
    zlor.z2 = (float*) malloc(nummainaxis[2] * sizeof(float));
	zlor.linevalue = (float*) malloc(nummainaxis[2] * sizeof(float));
	zlor.weight = (float*) malloc(nummainaxis[2] * sizeof(float));
	timeall.memoryIO += timing.GetCounter();

//	cout<<"after lor malloc"<<endl;
	timing.StartCounter();
	int cma[3] = {0,0,0};	//structure of arrays. 
	for(int i=0; i< numline; i++)
	{
		lor bufflor = alllor[i];
		if(bufflor.mainaxis == 0) 
		{
			xlor.x1[cma[0]] = bufflor.x1;
			xlor.y1[cma[0]] = bufflor.y1;
			xlor.z1[cma[0]] = bufflor.z1;
			xlor.x2[cma[0]] = bufflor.x2;
			xlor.y2[cma[0]] = bufflor.y2;
			xlor.z2[cma[0]] = bufflor.z2;
			xlor.weight[cma[0]] = bufflor.weight;
			cma[0] += 1;
		}
		else if(bufflor.mainaxis == 1)
		{
            ylor.x1[cma[1]] = bufflor.x1;
            ylor.y1[cma[1]] = bufflor.y1;
            ylor.z1[cma[1]] = bufflor.z1;
            ylor.x2[cma[1]] = bufflor.x2;
            ylor.y2[cma[1]] = bufflor.y2;
            ylor.z2[cma[1]] = bufflor.z2;
			ylor.weight[cma[1]] = bufflor.weight;

            cma[1] += 1;
	
		}
		else if(bufflor.mainaxis == 2)
		{
            zlor.x1[cma[2]] = bufflor.x1;
            zlor.y1[cma[2]] = bufflor.y1;
            zlor.z1[cma[2]] = bufflor.z1;
            zlor.x2[cma[2]] = bufflor.x2;
            zlor.y2[cma[2]] = bufflor.y2;
            zlor.z2[cma[2]] = bufflor.z2;
			zlor.weight[cma[2]] = bufflor.weight;
            cma[2] += 1;
		}
	}
	if(cma[0] !=  nummainaxis[0] || cma[1] != nummainaxis[1] || cma[2] != nummainaxis[2]) cout<< "Something wrong with the number of lors for each main axis!!" <<endl;


	vector<lor>().swap(alllor);		//deallocate lor

	timeall.lorsorting += timing.GetCounter();

//	icout<<"before cudamalloc"<<endl;
	timinggpu.StartCounter();
    cudaMalloc((void**) &dev_xlor.x1, nummainaxis[0]*sizeof(float) );
	cudaMalloc((void**) &dev_xlor.y1, nummainaxis[0]*sizeof(float) ); 
	cudaMalloc((void**) &dev_xlor.z1, nummainaxis[0]*sizeof(float) ); 
	cudaMalloc((void**) &dev_xlor.x2, nummainaxis[0]*sizeof(float) ); 
	cudaMalloc((void**) &dev_xlor.y2, nummainaxis[0]*sizeof(float) ); 
	cudaMalloc((void**) &dev_xlor.z2, nummainaxis[0]*sizeof(float) ); 	
	cudaMalloc((void**) &dev_xlor.linevalue, nummainaxis[0]*sizeof(float) );
	cudaMalloc((void**) &dev_xlor.weight, nummainaxis[0]*sizeof(float) );

    cudaMalloc((void**) &dev_ylor.x1, nummainaxis[1]*sizeof(float) );
    cudaMalloc((void**) &dev_ylor.y1, nummainaxis[1]*sizeof(float) ); 
    cudaMalloc((void**) &dev_ylor.z1, nummainaxis[1]*sizeof(float) ); 
    cudaMalloc((void**) &dev_ylor.x2, nummainaxis[1]*sizeof(float) ); 
    cudaMalloc((void**) &dev_ylor.y2, nummainaxis[1]*sizeof(float) ); 
    cudaMalloc((void**) &dev_ylor.z2, nummainaxis[1]*sizeof(float) );
    cudaMalloc((void**) &dev_ylor.linevalue, nummainaxis[1]*sizeof(float) );
	cudaMalloc((void**) &dev_ylor.weight, nummainaxis[1]*sizeof(float) );	

    cudaMalloc((void**) &dev_zlor.x1, nummainaxis[2]*sizeof(float) );
    cudaMalloc((void**) &dev_zlor.y1, nummainaxis[2]*sizeof(float) ); 
    cudaMalloc((void**) &dev_zlor.z1, nummainaxis[2]*sizeof(float) ); 
    cudaMalloc((void**) &dev_zlor.x2, nummainaxis[2]*sizeof(float) ); 
    cudaMalloc((void**) &dev_zlor.y2, nummainaxis[2]*sizeof(float) ); 
    cudaMalloc((void**) &dev_zlor.z2, nummainaxis[2]*sizeof(float) ); 
	cudaMalloc((void**) &dev_zlor.linevalue, nummainaxis[2]*sizeof(float) );
	cudaMalloc((void**) &dev_zlor.weight, nummainaxis[2]*sizeof(float) );

//	cout<<"before cuda memory copy"<<endl;
	cudaMemcpy( dev_xlor.x1, xlor.x1, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_xlor.y1, xlor.y1, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_xlor.z1, xlor.z1, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_xlor.x2, xlor.x2, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_xlor.y2, xlor.y2, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_xlor.z2, xlor.z2, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_xlor.weight, xlor.weight, nummainaxis[0]*sizeof(float), cudaMemcpyHostToDevice );

    cudaMemcpy( dev_ylor.x1, ylor.x1, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_ylor.y1, ylor.y1, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_ylor.z1, ylor.z1, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_ylor.x2, ylor.x2, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_ylor.y2, ylor.y2, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_ylor.z2, ylor.z2, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_ylor.weight, ylor.weight, nummainaxis[1]*sizeof(float), cudaMemcpyHostToDevice );

    cudaMemcpy( dev_zlor.x1, zlor.x1, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_zlor.y1, zlor.y1, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_zlor.z1, zlor.z1, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_zlor.x2, zlor.x2, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_zlor.y2, zlor.y2, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_zlor.z2, zlor.z2, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_zlor.weight, zlor.weight, nummainaxis[2]*sizeof(float), cudaMemcpyHostToDevice );
	timeall.memoryIO += timinggpu.GetCounter();

    return 0;
}

//calculate snmatrix based on projection and previous value. No normalization, no regularization, no blur
__global__ void calnewmatrix000(float *snmatrix, float *smatrix)
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
	int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
	int jj;	//image index in 1D
	while(z < nz)
	{
		y = threadIdx.y;
		while(y < ny)
		{
			x = threadIdx.x;
			while(x < nx)
			{
			    jj = x + y * nx + z * nx * ny;
				snmatrix[jj] = snmatrix[jj] * smatrix[jj];
				x += blockDim.x;
			}
			y += blockDim.y;
		}
		z += gridDim.x;
	}
}

//Yes normalization, no regularization, no blur
__global__ void calnewmatrix100(float *snmatrix, float *smatrix, float *normimage)
{
    int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
	int jj;
	float senratio, maxnorm = d_norm[0], ThreshNorm = d_norm[1];
	while(z < nz)
	{
		y = threadIdx.y;
		while(y < ny)
		{
			x = threadIdx.x;
			while(x < nx)
			{
			    jj = x + y * nx + z * nx * ny;
				if(normimage[jj] / maxnorm < ThreshNorm) senratio = ThreshNorm;
				else senratio = normimage[jj] / maxnorm;
       			snmatrix[jj] = snmatrix[jj] * smatrix[jj] / senratio;
				x += blockDim.x;
			}
			y += blockDim.y;
		}
		z += gridDim.x;
	}
}       


//No normalization, yes regularization, no blur
__global__ void calnewmatrix010(float *snmatrix, float *smatrix, float *poimage)
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
    int jj; //image index in 1D
	float beta = d_info[3], aa, bb, cc, laves = aves[0], lavep = avep[0];
    while(z < nz)
    {
        y = threadIdx.y;
        while(y < ny)
        {
            x = threadIdx.x;
            while(x < nx)
            {
                jj = x + y * nx + z * nx * ny;
				aa = 2.0f * beta / powf(laves,2);
				bb = 1.0f - 2.0f * beta * poimage[jj] / (laves * lavep);
				cc = -snmatrix[jj] * smatrix[jj];
                snmatrix[jj] = (-bb + sqrtf(powf(bb,2) - 4.0f * aa * cc)) / (2.0f * aa);
                x += blockDim.x;
            }
            y += blockDim.y;
        }
        z += gridDim.x;
    }

}

//Yes normalization, yes regularization, no blur
__global__ void calnewmatrix110(float *snmatrix, float *smatrix, float *normimage, float *poimage)
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
    int jj; //image index in 1D
	float beta = d_info[3], aa, bb, cc, laves = aves[0], lavep = avep[0];
	float senratio, maxnorm = d_norm[0], ThreshNorm = d_norm[1];
    while(z < nz)
    {
        y = threadIdx.y;
        while(y < ny)
        {
            x = threadIdx.x;
            while(x < nx)
            {
                jj = x + y * nx + z * nx * ny;
				if(normimage[jj] / maxnorm < ThreshNorm) senratio = ThreshNorm;
				else senratio = normimage[jj] / maxnorm;
				aa = 2.0f * beta / powf(laves,2);
				bb = senratio - 2.0f * beta * poimage[jj] / (laves * lavep);
				cc = -snmatrix[jj] * smatrix[jj];
                snmatrix[jj] = (-bb + sqrtf(powf(bb,2) - 4.0f * aa * cc)) / (2.0f * aa);
                x += blockDim.x;
            }
            y += blockDim.y;
        }
        z += gridDim.x;
    }

}


//No normalization, ues regularization, yes blur.
__global__ void calnewmatrix011(float *snmatrix, float *smatrix, float *poimage, float *bmatrix, float *allweight)
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
    int jjj; //image index in 1D
	float beta = d_info[3], aa, bb, cc, laves = aves[0], lavep = avep[0], wi, a=d_info[0];
	int li,hi,lj,hj,lk,hk, idxy;
    while(z < nz)
    {
        y = threadIdx.y;
        while(y < ny)
        {
            x = threadIdx.x;
            while(x < nx)
            {
                jjj = x + y * nx + z * nx * ny;

				li = max(0, x - d_indr[0]);
				lj = max(0, y - d_indr[1]);
				lk = max(0, z - d_indr[2]);
				hi = min(nx - 1, x + d_indr[0]);
				hj = min(ny - 1, y + d_indr[1]);
				hk = min(nz - 1, z + d_indr[2]);

				aa = 0.0f;
				bb = 0.0f;
		
				for(int ii= li; ii<= hi; ii++)
				{
					for(int jj = lj; jj<= hj; jj++)
					{
						for(int kk = lk; kk<= hk; kk++)
						{
							if((powf((ii-x)/d_rads[0],2) + powf((jj-y)/d_rads[1],2) + powf((kk-z)/d_rads[2],2)) * powf(a,2) < 1.0) 
							{
								idxy = ii + jj * nx + kk * nx * ny;
                                wi = expf(-(powf(ii-x,2)/(2.0f * powf(d_bsgm[0],2)) + powf(jj-y,2)/(2.0f * powf(d_bsgm[1],2)) + powf(kk-z,2)/(2.0f * powf(d_bsgm[2],2))) * powf(a,2)) / allweight[idxy];
                                aa += wi;
                                bb += ((bmatrix[idxy] - smatrix[jjj] )/laves - poimage[idxy]/lavep) * wi;

							}
		
						}
					}
				}

				aa = aa * 2.0f * beta / powf(laves,2);
				bb = bb * 2.0f * beta / laves + 1.0f;
				cc = -snmatrix[jjj] * smatrix[jjj];
                snmatrix[jjj] = (-bb + sqrtf(powf(bb,2) - 4.0f * aa * cc)) / (2.0f * aa);
                x += blockDim.x;
            }
            y += blockDim.y;
        }
        z += gridDim.x;
    }

}

//Yes normalization, yes regularization, yes blur
__global__ void calnewmatrix111(float *snmatrix, float *smatrix, float *normimage, float *poimage, float *bmatrix, float *allweight)
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
    int jjj; //image index in 1D
	float beta = d_info[3], aa, bb, cc, laves = aves[0], lavep = avep[0], wi, a = d_info[0];
	int li,hi,lj,hj,lk,hk, idxy;
	float senratio, maxnorm = d_norm[0], ThreshNorm = d_norm[1];
    while(z < nz)
    {
        y = threadIdx.y;
        while(y < ny)
        {
            x = threadIdx.x;
            while(x < nx)
            {
                jjj = x + y * nx + z * nx * ny;
				if(normimage[jjj] / maxnorm < ThreshNorm) senratio = ThreshNorm;
				else senratio = normimage[jjj] / maxnorm;

				li = max(0, x - d_indr[0]);
				lj = max(0, y - d_indr[1]);
				lk = max(0, z - d_indr[2]);
				hi = min(nx - 1, x + d_indr[0]);
				hj = min(ny - 1, y + d_indr[1]);
				hk = min(nz - 1, z + d_indr[2]);

				aa = 0.0f;
				bb = 0.0f;
		
				for(int ii= li; ii<= hi; ii++)
				{
					for(int jj = lj; jj<= hj; jj++)
					{
						for(int kk = lk; kk<= hk; kk++)
						{
							if((powf((ii-x)/d_rads[0],2) + powf((jj-y)/d_rads[1],2) + powf((kk-z)/d_rads[2],2)) * powf(a,2) < 1.0) 
							{
								idxy = ii + jj * nx + kk * nx * ny;
                                wi = expf(-(powf(ii-x,2)/(2.0f * powf(d_bsgm[0],2)) + powf(jj-y,2)/(2.0f * powf(d_bsgm[1],2)) + powf(kk-z,2)/(2.0f * powf(d_bsgm[2],2))) * powf(a,2)) / allweight[idxy];
                                aa += wi;
                                bb += ((bmatrix[idxy] - smatrix[jjj])/laves - poimage[idxy]/lavep) * wi;

							}
		
						}
					}
				}

				aa = aa * 2.0f * beta / powf(laves,2);
				bb = bb * 2.0f * beta / laves + senratio;
				cc = -snmatrix[jjj] * smatrix[jjj];
                snmatrix[jjj] = (-bb + sqrtf(powf(bb,2) - 4.0f * aa * cc)) / (2.0f * aa);
                x += blockDim.x;
            }
            y += blockDim.y;
        }
        z += gridDim.x;
    }

}

//calculate average of voxel values
__global__ void calave(float *smatrix, float *gave)
{
	int msize = d_imageindex[3];
	int cacheindex = threadIdx.x, tid = threadIdx.x + blockIdx.x * blockDim.x;
   	__shared__ float buffave[reducsize];

	float buff = 0.0f;
	while(tid < msize)
	{
		buff += smatrix[tid];
		tid += blockDim.x * gridDim.x ;
	}
	buffave[cacheindex] = buff;
	__syncthreads();

	int i = blockDim.x / 2;
	while( i != 0)
	{
		if(cacheindex < i)  buffave[cacheindex] += buffave[cacheindex + i];
		__syncthreads();
		i /= 2;
	}

	if(cacheindex == 0) gave[blockIdx.x] = buffave[0];
}

//calculate average of voxel values, excluding voxel with too large values, based on average value calculated without filter. Also change the voxel value to average value to reduce artifacts, especially in image recon with normalization.
__global__ void calavewithfilter(float *smatrix, float *gave, int *gnumave)
{
	int msize = d_imageindex[3];
	int cacheindex = threadIdx.x, tid = threadIdx.x + blockIdx.x * blockDim.x;
   	__shared__ float buffave[reducsize];
	__shared__ int buffnumave[reducsize];
    float avepre = aveunf[0];	//average voxel value before filtering

	float buff = 0.0f;
	int buffnum = 0;
	while(tid < msize)
	{
		//if voxelvalue is less than UpholdVox * previous average value, then count it towards new average calculation. Otherwise, change the voxel value to average value to reduce artifacts, especially in image recon with normalization.
		if(smatrix[tid] < UpholdVox * avepre){
		buff += smatrix[tid];
		buffnum += 1;
		}
		else smatrix[tid]  = avepre;

		tid += blockDim.x * gridDim.x ;
	}
	buffave[cacheindex] = buff;
	buffnumave[cacheindex] = buffnum;
	__syncthreads();

	int i = blockDim.x / 2;
	while( i != 0)
	{
		if(cacheindex < i)  {buffave[cacheindex] += buffave[cacheindex + i]; buffnumave[cacheindex] += buffnumave[cacheindex + i];}
		__syncthreads();
		i /= 2;
	}

	if(cacheindex == 0) {gave[blockIdx.x] = buffave[0]; gnumave[blockIdx.x] = buffnumave[0];}
}


//Gaussian blur to image
__global__ void gpublur(float *smatrix, float *bmatrix, float *allweight)
{
	int nx = d_imageindex[0], ny = d_imageindex[1], nz = d_imageindex[2];
    int x = threadIdx.x, y = threadIdx.y, z = blockIdx.x;
    int i; //image index in 1D
    float sumweight, sumval, a = d_info[0];
	int li,hi,lj,hj,lk,hk;
    float	wi;

    while(z < nz)
    {
        y = threadIdx.y;
        while(y < ny)
        {
            x = threadIdx.x;
            while(x < nx)
            {
                i = x + y * nx + z * nx * ny;
				li = max(0, x - d_indr[0]);
				lj = max(0, y - d_indr[1]);
				lk = max(0, z - d_indr[2]);
				hi = min(nx - 1, x + d_indr[0]);
				hj = min(ny - 1, y + d_indr[1]);
				hk = min(nz - 1, z + d_indr[2]);
		        sumweight = 0.0f;
		        sumval = 0.0f;
		
				for(int ii= li; ii<= hi; ii++)
				{
					for(int jj = lj; jj<= hj; jj++)
					{
						for(int kk = lk; kk<= hk; kk++)
						{
							if((powf((ii-x)/d_rads[0],2) + powf((jj-y)/d_rads[1],2) + powf((kk-z)/d_rads[2],2)) * powf(a,2) < 1.0) 
							{
								wi = expf(-(powf(ii-x,2)/(2.0f * powf(d_bsgm[0],2)) + powf(jj-y,2)/(2.0f * powf(d_bsgm[1],2)) + powf(kk-z,2)/(2.0f * powf(d_bsgm[2],2))) * powf(a,2));
								sumval += smatrix[ii + jj * nx + kk * nx * ny] * wi;
								sumweight += wi;
							}
		
						}
					}
				}
				allweight[i] = sumweight;
				bmatrix[i] = sumval / sumweight;
                x += blockDim.x;
            }
            y += blockDim.y;
        }
        z += gridDim.x;
    }


}

// For calculating loglikelihood function value
__global__ void calLogLike(float *xlinevalue, double *gloglike, const int lorindex)
{
	//int msize = d_imageindex[3];
	int cacheindex = threadIdx.x, tid = threadIdx.x + blockIdx.x * blockDim.x;
   	__shared__ double buffave[reducsize];

	double buff = 0.0;
	while(tid < lorindex)
	{
		if (xlinevalue[tid] > ThreshLineValue ) buff += logf(xlinevalue[tid]);
		tid += blockDim.x * gridDim.x ;
	}
	buffave[cacheindex] = buff;
	__syncthreads();

	int i = blockDim.x / 2;
	while( i != 0)
	{
		if(cacheindex < i)  buffave[cacheindex] += buffave[cacheindex + i];
		__syncthreads();
		i /= 2;
	}

	if(cacheindex == 0) gloglike[blockIdx.x] = buffave[0];
}

// For calculating loglikelihood function value 
__global__ void calLogLikeS(float *smatrix, float *normimage, double *gloglike, const int msize, const int norma)
{
	//int msize = d_imageindex[3];
	int cacheindex = threadIdx.x, tid = threadIdx.x + blockIdx.x * blockDim.x;
	float maxnorm = d_norm[0], ThreshNorm = d_norm[1];
   	__shared__ double buffave[reducsize];

	double buff = 0.0;
	while(tid < msize)
	{
		if(norma == 0) buff -= smatrix[tid];
		else {
			float tempSen = normimage[tid] / maxnorm;
			if(tempSen < ThreshNorm) tempSen = ThreshNorm;
			buff -= tempSen * smatrix[tid];
		}
		tid += blockDim.x * gridDim.x ;
	}
	buffave[cacheindex] = buff;
	__syncthreads();

	int i = blockDim.x / 2;
	while( i != 0)
	{
		if(cacheindex < i)  buffave[cacheindex] += buffave[cacheindex + i];
		__syncthreads();
		i /= 2;
	}

	if(cacheindex == 0) gloglike[blockIdx.x] = buffave[0];
}

// For calculating loglikelihood function value 
__global__ void calLogR(float *smatrix, float *poimage, double *gloglike, const int msize)
{
	//int msize = d_imageindex[3];
	int cacheindex = threadIdx.x, tid = threadIdx.x + blockIdx.x * blockDim.x;
   	__shared__ double buffave[reducsize];

	double buff = 0.0;
	while(tid < msize)
	{
		buff -= powf(smatrix[tid]/aves[0] - poimage[tid]/avep[0],2);
		tid += blockDim.x * gridDim.x ;
	}
	buffave[cacheindex] = buff;
	__syncthreads();

	int i = blockDim.x / 2;
	while( i != 0)
	{
		if(cacheindex < i)  buffave[cacheindex] += buffave[cacheindex + i];
		__syncthreads();
		i /= 2;
	}

	if(cacheindex == 0) gloglike[blockIdx.x] = buffave[0];
}
