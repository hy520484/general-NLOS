import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt

data = si.loadmat('C_d.mat')['C_d']	
print(data.shape)	# (1229, 3)

x = data[:, 0]#.reshape(64, 64)
y = data[:, 1]#.reshape(64, 64)
z = data[:, 2]#.reshape(64, 64)

# plt.figure('xyz')
# plt.subplot(131)
# plt.title('x')
# plt.imshow(x, cmap = 'jet')
# plt.subplot(132)
# plt.title('y')
# plt.imshow(y, cmap = 'jet')
# plt.subplot(133)
# plt.title('z')
# plt.imshow(z, cmap = 'jet')
# plt.show()

# fig = plt.figure('3D scatter')
# ax = fig.add_subplot(111, projection = '3d')
# ax.scatter3D(x, y, z)
# ax.set_xlabel('x')
# ax.set_ylabel('z')
# ax.set_zlabel('y')
# plt.show()

