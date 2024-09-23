import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
from scipy import signal
import pycuda.driver as cuda
import pycuda.tools
import pycuda.autoinit
from pycuda.compiler import SourceModule
import sys
from utils.ParseArgs import parse_args
from utils.write_ply import write_ply
import time
from utils.show_result import show_vol

##### define and load parameters #####
write_3D_point_cloud = 1

# choose one config file
args = parse_args('./config/config_0527_tuan.ini')							### Fig. 3 (d)
# args = parse_args('./config/config_0527_wall.ini')							### Fig. 3 (g)
# args = parse_args('./config/CC-SOCR/config_bunny_full.ini')
# args = parse_args('./config/CC-SOCR/config_bunny_sticks.ini')					### Fig. 4
# args = parse_args('./config/CC-SOCR/config_bunny_window.ini')
# args = parse_args('./config/CC-SOCR/config_non_planar_letters_NT_full.ini')	### Fig. 4
# args = parse_args('./config/CC-SOCR/config_statue_NLOS.ini')					### Fig. 4
# args = parse_args('./config/CC-SOCR/config_statue_sticks.ini')
# args = parse_args('./config/CC-SOCR/config_figure_4_random.ini')				### Fig. 4
# args = parse_args('./config/3D-RSD/config_data_0.ini')						### Fig. 5
# args = parse_args('./config/3D-RSD/config_data_2.ini')						### Fig. 5
# args = parse_args('./config/3D-RSD/config_data_6.ini')						### Fig. 5
# args = parse_args('./config/3D-RSD/config_data_10.ini')						### Fig. 5

N_ld = args['N_ld']
N_bin = args['N_bin']
N_voxel_s = args['N_voxel_s']
N_voxel_t = args['N_voxel_t']
timeRes = args['timeRes']
# data_type = np.int32
data_type = np.float32

# read data
data = si.loadmat(args['path_data'])['data']

# read setup
setup = si.loadmat(args['path_setup'])

if data.shape != (N_ld, N_ld, N_bin):
	sys.exit('Data size error!')


##### filtering #####
wall_size = args['width'] * 2
sampling_coeff = 1.2	##### 1.2

s_lamda_limit = wall_size / (N_ld - 1)	# sample space on the wall

virtual_wavelength = sampling_coeff * (s_lamda_limit * 2)	# virtual wavelength (cm)
cycles = 5	# number of wave cycles in the wavelet, typically 4~6
s_z = timeRes * 3e8
samples = int(np.round(cycles * virtual_wavelength / s_z))
sigma = args['sigma']	##### 0.3

start0 = time.perf_counter()
sin_wave = np.sin(2*np.pi*(cycles * np.linspace(1,samples,samples))/samples)
cos_wave = np.cos(2*np.pi*(cycles * np.linspace(1,samples,samples))/samples)
window = signal.windows.gaussian(samples, std = (samples-1)*sigma/2)
virtual_wave_sin = (sin_wave * window).reshape(1, 1, -1)
virtual_wave_cos = (cos_wave * window).reshape(1, 1, -1)
filterd_data_sin = signal.convolve(data, virtual_wave_sin, 'same')
filterd_data_cos = signal.convolve(data, virtual_wave_cos, 'same')
filterd_data_sin = np.asarray(filterd_data_sin, dtype = data_type)
filterd_data_cos = np.asarray(filterd_data_cos, dtype = data_type)
end0 = time.perf_counter()

# plt.plot(sin_wave)
# plt.plot(cos_wave)
# plt.show()

# plt.plot(window)
# plt.show()

# plt.plot(virtual_wave_sin.reshape(-1))
# plt.plot(virtual_wave_cos.reshape(-1))
# plt.show()
# sys.exit()


##### back-projection #####
laserpoints = np.array(setup['laserpoints'], dtype = data_type)
detectpoints = np.array(setup['detectpoints'], dtype = data_type)
if laserpoints.shape != detectpoints.shape or laserpoints.shape != (N_ld, N_ld, 3):
	sys.exit('Setup size error!')
if write_3D_point_cloud == 1:
	write_ply(laserpoints.reshape(-1, 3), './laserpoints.ply')
	write_ply(detectpoints.reshape(-1, 3), './detectpoints.ply')
# sys.exit()

with open('./utils/kernel_BP.cu', 'r') as f:
	source = f.read()
if data_type == np.int32:
	source = source.replace('{data_type}', 'int')
elif data_type == np.float32:
	source = source.replace('{data_type}', 'float')
mod = SourceModule(source)
bp_cuda = mod.get_function('bp')

def BP(data, laserpoints, detectpoints, bp_cuda):
	vol = np.zeros((N_voxel_s, N_voxel_s, N_voxel_t), dtype = data_type)

	d_data = cuda.mem_alloc(data.nbytes)
	d_laserpoints = cuda.mem_alloc(laserpoints.nbytes)
	d_detectpoints = cuda.mem_alloc(detectpoints.nbytes)
	d_vol = cuda.mem_alloc(vol.nbytes)

	cuda.memcpy_htod(d_data, data.reshape(N_ld*N_ld, N_bin))
	cuda.memcpy_htod(d_laserpoints, laserpoints.reshape(N_ld*N_ld, 3))
	cuda.memcpy_htod(d_detectpoints, detectpoints.reshape(N_ld*N_ld, 3))
	grid = (N_ld*N_ld, 1, 1)
	block = (1024, 1, 1)
	shared_memory_size = N_bin * 4	# np.float32 -> 4; np.float64 -> 8

	bp_cuda(d_data, d_laserpoints, d_detectpoints, d_vol, np.float32(timeRes), 
		np.int32(N_voxel_s), np.int32(N_voxel_t), np.int32(N_bin), 
		np.float32(args['x_range0']), np.float32(args['x_range1']), np.float32(args['y_range0']),
		np.float32(args['y_range1']), np.float32(args['z_range0']), np.float32(args['z_range1']), 
		block = block, grid = grid, shared = shared_memory_size)
	cuda.memcpy_dtoh(vol, d_vol)

	d_data.free()
	d_laserpoints.free()
	d_detectpoints.free()
	d_vol.free()

	return vol

# vol_tmp = BP(data, laserpoints, detectpoints, bp_cuda)
# si.savemat('./data/vol_tmp.mat', {'vol_tmp': vol_tmp})

start = time.perf_counter()
vol_sin = BP(filterd_data_sin, laserpoints, detectpoints, bp_cuda)
vol_cos	= BP(filterd_data_cos, laserpoints, detectpoints, bp_cuda)
vol = np.sqrt(vol_sin**2 + vol_cos**2)
end = time.perf_counter()
print('filtering time: {:.4f}s, BP time: {:.4f}s'.format(end0 - start0, end - start))
print('total time: {:.4f}s'.format(end - start + end0 - start0))

si.savemat('./data_tmp/vol.mat', {'vol': vol})
show_vol(vol, args['th'])
