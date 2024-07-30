import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
import sys

dtype = np.float32
bz = 9.6e-3
bx = 2/63

C_d = si.loadmat('C_d.mat')['C_d']	# (1229, 3)
C_i = si.loadmat('C_i.mat')['C_i']
N = int(np.ceil(np.sqrt(C_d.shape[0])))
print(N)
# sys.exit()

laserpoints = np.zeros((N * N, 3), dtype = dtype)
detectpoints = np.zeros((N * N, 3), dtype = dtype)

laserpoints[:C_i.shape[0], 0] = C_i[:, 1] * bx
laserpoints[:C_i.shape[0], 1] = C_i[:, 2] * bx
laserpoints[:C_i.shape[0], 2] = C_i[:, 0] * bz
detectpoints[:C_i.shape[0], 0] = C_d[:, 1] * bx
detectpoints[:C_i.shape[0], 1] = C_d[:, 2] * bx
detectpoints[:C_i.shape[0], 2] = C_d[:, 0] * bz

laserpoints[C_i.shape[0]:, 0] = C_i[0, 1] * bx
laserpoints[C_i.shape[0]:, 1] = C_i[0, 2] * bx
laserpoints[C_i.shape[0]:, 2] = C_i[0, 0] * bx
detectpoints[C_d.shape[0]:, 0] = C_d[0, 1] * bx
detectpoints[C_d.shape[0]:, 1] = C_d[0, 2] * bx
detectpoints[C_d.shape[0]:, 2] = C_d[0, 0] * bx

laserpoints = laserpoints.reshape(N, N, -1)
detectpoints = detectpoints.reshape(N, N, -1)

si.savemat('setup.mat', {'laserpoints': laserpoints, 'detectpoints': detectpoints})


sub_Sig = si.loadmat('sub_Sig.mat')['sub_Sig']	# (1229, 512)
print(sub_Sig.shape)
# sys.exit()
data = np.zeros((N * N, sub_Sig.shape[1]), dtype = dtype)
data[:sub_Sig.shape[0], :] = sub_Sig
data = data.reshape(N, N, -1)

si.savemat('data.mat', {'data': data})
