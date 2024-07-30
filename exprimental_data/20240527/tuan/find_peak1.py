import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
import sys
import scipy.signal as signal

N = 64
show = 0

path = 'data_laser.mat'
data = si.loadmat(path)['counts'][:4096,:]
# print(data.shape)
# sys.exit()
data = data.reshape(N, N, -1)
for i in range(1, N, 2):
	data[i, :] = np.flip(data[i, :], 0)
print(data.shape)

# hist = data[63, 62]
# plt.plot(hist)
# plt.show()
# sys.exit()

peak1 = np.argmax(data[:, :, :300], 2)
# peak1[40, 46] = 90

plt.imshow(peak1, cmap = 'jet')
plt.show()

si.savemat('peak1.mat', {'peak1': peak1})
