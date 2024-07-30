import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
import sys

N = 64

path = 'data_tuan.mat'
data = si.loadmat(path)['counts'][:4096, :]
data = data.reshape(N, N, -1)
for i in range(1, N, 2):
	data[i, :] = np.flip(data[i, :], 0)
print(data.shape)

peak1 = si.loadmat('peak1.mat')['peak1']
peak2 = si.loadmat('peak2.mat')['peak2']
peak = (peak1 + peak2) / 2

for i in range(N):
	for j in range(N):
		data[i, j] = np.roll(data[i, j], - int(np.round(peak[i, j])))

data = data[:, :, :512]
data[:, :, :40] = 0

plt.figure('hist_tot')
hist_tot = np.sum(np.sum(data, 0), 0)
plt.plot(hist_tot)
plt.show()

si.savemat('data.mat', {'data': data})
