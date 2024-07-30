typedef float real;
// typedef double real;
typedef {data_type} data_type;	// np.int32 -> int, np.float32 -> float
#define c 3e8	// the speed of light

__global__ void bp(data_type* d_data, real* d_laserpoints, real* d_detectpoints, data_type* d_vol,
	float timeRes, int N_voxel_s, int N_voxel_t, int N_bin, float x_range0, float x_range1,
	float y_range0, float y_range1, float z_range0, float z_range1){
	// gridsize: (N_ld*N_ld, 1, 1); blocksize: (min(N_bin, 1024), 1, 1)
	// d_data: (N_ld, N_ld, N_bin)
	// d_laserpoints, d_detectpoints: (N_ld*N_ld, 2)
	// d_vol: (N_voxel_s, N_voxel_s, N_voxel_t)
	int bx = blockIdx.x;	// 0 ~ N_ld * N_ld - 1
	int tx = threadIdx.x;	// 0 ~ 1023
	// int tid = tx + bx * blockDim.x;	// not correct

	// one block deal with one pair of lp and dp, shares the following:
	__shared__ real s_lpx;
	__shared__ real s_lpy;
	__shared__ real s_lpz;
	__shared__ real s_dpx;
	__shared__ real s_dpy;
	__shared__ real s_dpz;
	extern __shared__ data_type s_data[];	// dynamic shared memory, length: N_bin

	// all blocks intialize d_vol to zeros, no need after python np.zeros() initialization
	// for (int k = 0; tid + blockDim.x * gridDim.x < N_voxel_s * N_voxel_s * N_voxel_t; k++){
	// 	d_vol[tid + blockDim.x * gridDim.x] = 0;
	// }

	// transport from d_data to s_data
	for (int k = 0; tx + k * blockDim.x < N_bin; k++){
		s_data[tx + k * blockDim.x] = d_data[tx + k * blockDim.x + bx * N_bin];
	}

	if (tx == 0){
		s_lpx = d_laserpoints[bx * 3 + 0];
		s_lpy = d_laserpoints[bx * 3 + 1];
		s_lpz = d_laserpoints[bx * 3 + 2];
		s_dpx = d_detectpoints[bx * 3 + 0];
		s_dpy = d_detectpoints[bx * 3 + 1];
		s_dpz = d_detectpoints[bx * 3 + 2];
	}
	__syncthreads();

	// backprojection, one thread deal with one voxel (vx, vy, vz)
	for (int k = 0; tx + k * blockDim.x < N_voxel_s * N_voxel_s * N_voxel_t; k++){
		int vx = (tx + k * blockDim.x) / N_voxel_s / N_voxel_t;
		int vy = (tx + k * blockDim.x - vx * N_voxel_s * N_voxel_t) / N_voxel_t;
		int vz = tx + k * blockDim.x - vx * N_voxel_s * N_voxel_t - vy * N_voxel_t;
		
		real x = x_range0 + (x_range1 - x_range0) / (N_voxel_s - 1) * real(vx);
		real y = y_range0 + (y_range1 - y_range0) / (N_voxel_s - 1) * real(vy);
		real z = z_range0 + (z_range1 - z_range0) / (N_voxel_t - 1) * real(vz);
		
		real t1 = sqrt((x - s_lpx) * (x - s_lpx) + (y - s_lpy) * (y - s_lpy) + (z - s_lpz) * (z - s_lpz));
		real t2 = sqrt((x - s_dpx) * (x - s_dpx) + (y - s_dpy) * (y - s_dpy) + (z - s_dpz) * (z - s_dpz));
		real d = t1 + t2;

		int bin = int(d / timeRes / c + 0.5);
		if (bin >= 10 && bin < N_bin){
			data_type num = s_data[bin];	// good
			// data_type num = data_type(float(s_data[bin]) * t1 * t1 * t2 * t2);	// bad
			// data_type num = data_type(float(s_data[bin]) * t1 * t2);	// good
			atomicAdd(&d_vol[tx + k * blockDim.x], num);
		}
	}
}
