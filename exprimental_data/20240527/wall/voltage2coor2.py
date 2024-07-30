import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
import sys
import write_ply

N = 128
timeRes = 32e-12
c = 3e8
show = 0

# relation between galvanometer's rotate (mechanical) angle and voltage
scale_X = - 0.5	# V/degree
scale_Y = 0.5	# V/degree
# delta_p = 443	# absolute distance = (p + delta_p) * timeRes * c / 2
delta_p = 436

# a = si.loadmat('../../scanning points/voltages.mat')
# print(a.keys())
XY_real = si.loadmat('v_detect.mat')['XY_real2']	# unit: V, shape: (2, 4096), 0 for X, 1 for Y
XY_real = XY_real.reshape(2, N, N)
# XY_real = np.flip(XY_real, 1)
for i in range(1, N, 2):
	XY_real[0, i, :] = np.flip(XY_real[0, i, :], 0)
	XY_real[1, i, :] = np.flip(XY_real[1, i, :], 0)

if show == 1:
	plt.figure('voltages')
	plt.subplot(121)
	plt.title('voltages_X')
	plt.imshow(XY_real[0], cmap = 'jet')
	plt.subplot(122)
	plt.title('voltages_Y')
	plt.imshow(XY_real[1], cmap = 'jet')
	# X_real = XY_real.reshape(2, 64, 64)[0]
	# Y_real = XY_real.reshape(2, 64, 64)[1]
	# plt.scatter(X_real, Y_real)
	# plt.gca().set_aspect(1)
	# plt.show()

# p = si.loadmat('startposition_{}.mat'.format('4ms'))['startposition'].reshape(-1)
p = si.loadmat('peak2.mat')['peak2']
# p = si.loadmat('peak_shangxia_fuza.mat')['peak']
# p = p.reshape(N, N)	# (64, 64)
# for i in range(1, N, 2):
# 	p[i, :] = np.flip(p[i, :])
# # print(p.shape)
if show == 1:
	plt.figure('p')
	plt.imshow(p, cmap = 'jet')
# plt.show()
# sys.exit()

# ori (0, 0, 0): galvanometer
angle_X = XY_real[0] / scale_X * 2	# degree
angle_Y = XY_real[1] / scale_Y * 2	# degree
distance = (p + delta_p) * timeRes * c / 2	# (64, 64)
# plt.title('distance')
# plt.imshow(distance, cmap = 'jet')
# plt.show()
# sys.exit()
points = np.zeros((N, N, 3))

# calculate the 3D coordinates, origin located on the galvanometer
# points[:, :, 0] = distance * np.cos(angle_Y/180*np.pi) * np.sin(angle_X/180*np.pi)	# x
# points[:, :, 1] = distance * np.sin(angle_Y/180*np.pi)	# y
# points[:, :, 2] = - distance * np.cos(angle_Y/180*np.pi)	# z

gamma = np.sqrt(angle_X**2 + angle_Y**2)	# (64, 64)
points[:, :, 0] = distance * np.sin(gamma/180*np.pi) / gamma * angle_X	# x
points[:, :, 1] = distance * np.sin(gamma/180*np.pi) / gamma * angle_Y	# y
points[:, :, 2] = - distance * np.cos(gamma/180*np.pi)	# z

p_left = points[31, 0, :]
p_right = points[31, 63,:]
# print(p_left, p_right)
# sys.exit()

if show == 1:
	plt.figure('points xyz')
	plt.subplot(131)
	plt.title('x')
	plt.imshow(points[:, :, 0], cmap = 'jet')
	plt.subplot(132)
	plt.title('y')
	plt.imshow(points[:, :, 1], cmap = 'jet')
	plt.subplot(133)
	plt.title('z')
	plt.imshow(points[:, :, 2], cmap = 'jet')
# plt.show()
# sys.exit()

points = points.reshape(-1, 3)
if show == 1:
	fig = plt.figure('3D scatter')
	ax = fig.add_subplot(111, projection = '3d')
	ax.scatter3D(points[:, 0], - points[:, 2], points[:, 1])
	ax.set_xlabel('x')
	ax.set_ylabel('z')
	ax.set_zlabel('y')
	plt.show()
	sys.exit()

# write_ply.write_ply(points, 'test.ply')
# print('write to tmp.ply')

# coordinate transformation
transformation = si.loadmat('transformation.mat')
# new_ori = (p_left + p_right) / 2
new_ori = transformation['new_ori'].reshape(-1)
angle_along_Yaxis = np.arctan((p_right[2] - p_left[2]) / (p_right[0] - p_left[0]))	# rad
beta = angle_along_Yaxis
print('rotate angle:\t', angle_along_Yaxis / np.pi * 180)

# translation
points[:, 0] -= new_ori[0]
points[:, 1] -= new_ori[1]
points[:, 2] -= new_ori[2]

# rotate
# R = np.array([[np.cos(beta), 0, np.sin(beta)],\
# 			[0, 1, 0],\
# 			[- np.sin(beta), 0, np.cos(beta)]])
R = transformation['R']
points = np.matmul(R, points.transpose()).transpose()	# (4096, 3)
print('z_max - z_min:\t', np.max(points[:, 2]) - np.min(points[:, 2]))
# sys.exit()

write_ply.write_ply(points, 'test.ply')

points = points.reshape(N, N, -1)
if show == 1:
	plt.figure('points xyz')
	plt.subplot(131)
	plt.title('x')
	plt.imshow(points[:, :, 0], cmap = 'jet')
	plt.subplot(132)
	plt.title('y')
	plt.imshow(points[:, :, 1], cmap = 'jet')
	plt.subplot(133)
	plt.title('z')
	plt.imshow(points[:, :, 2], cmap = 'jet')
	plt.show()

laserpoints = np.array(points, dtype = np.float32)
detectpoints = np.array(points, dtype = np.float32)
# detectpoints[:, :, 1] += 0.1

# plt.figure('laserpoints')
# plt.subplot(121)
# plt.title('x')
# plt.imshow(laserpoints[:, :, 0], cmap = 'jet')
# plt.subplot(122)
# plt.title('y')
# plt.imshow(laserpoints[:, :, 1], cmap = 'jet')
# plt.show()

si.savemat('detectpoints.mat', {'detectpoints': detectpoints})
