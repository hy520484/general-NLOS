import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
import sys
import write_ply

dtype = np.float32
N = 190
tt = '10'

path = 'scene_{}_raw.mat'.format(tt)
a = si.loadmat(path)
# print(a.keys())

timeRes = a['ts'][0][0]
print('timeRes: ', timeRes)
data = a['total_rect_data']	# (36100, 1024)
laserpoints = a['camera_pos']	# (36100, 3)
detectpoints = a['laser_pos_base']	# (1, 1, 3)
# print(detectpoints)

# hist_tot = np.sum(data, 0)
# plt.plot(hist_tot)
# plt.show()
sys.exit()

write_ply.write_ply(laserpoints, 'test.ply')

detectpoints = np.squeeze(detectpoints, 0)
detectpoints = detectpoints.repeat(36100, 0)	# (36100, 0, 0)
print(detectpoints.shape)

laserpoints = np.array(laserpoints.reshape(N, N, -1), dtype = dtype)
detectpoints = np.array(detectpoints.reshape(N, N, -1), dtype = dtype)
data = np.array(data.reshape(N, N, -1), dtype = dtype)

si.savemat('setup_{}.mat'.format(tt), {'laserpoints': laserpoints, 'detectpoints': detectpoints})
si.savemat('data_{}.mat'.format(tt), {'data': data})
