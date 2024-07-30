import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt
import sys

def show_vol(vol, th):
	axis = 0
	# cmap = 'gray'
	cmap = 'hot'

	# vol_phasor = si.loadmat('./data_tmp/vol_phasor.mat')
	# vol_sin = vol_phasor['vol_sin']
	# vol_cos = vol_phasor['vol_cos']
	print(vol.shape)
	# print(vol_sin.dtype)

	def simple_denoise0(im, th):
		im[im < np.max(im) * th] = 0
		return im

	def simple_denoise(im, th):
		mm = np.min(im)
		MM = np.max(im)
		im[im < mm + (MM - mm) * th] = mm
		return im

	# vol = vol_sin
	# vol = vol_cos
	# print(np.min(vol), np.max(vol))
	# vol = np.sqrt(vol_sin**2 + vol_cos**2)
	vol = np.flip(vol, 0)
	vol = np.flip(vol, 1)
	vol = vol.transpose(1, 0, 2)
	# print(vol.dtype)

	simple_denoise(vol, th)

	im0 = np.max(vol, 2)
	im1 = np.max(vol, 1)
	im2 = np.max(vol, 0)

	plt.figure('result', figsize = (15, 4.5))
	plt.subplot(131)
	plt.title('front')
	if axis == 0:
		plt.axis('off')
	plt.imshow(im0, cmap = cmap)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.subplot(132)
	plt.title('side')
	if axis == 0:
		plt.axis('off')
	plt.imshow(im1, cmap = cmap)
	plt.xlabel('z')
	plt.ylabel('y')
	plt.subplot(133)
	plt.title('top')
	if axis == 0:
		plt.axis('off')
	plt.imshow(im2, cmap = cmap)
	plt.xlabel('z')
	plt.ylabel('x')
	# if axis == 0:
	# 	plt.savefig('result.png')
	plt.show()


def show(vol_sin, vol_cos, th):
	axis = 0
	# cmap = 'gray'
	cmap = 'hot'

	# vol_phasor = si.loadmat('./data_tmp/vol_phasor.mat')
	# vol_sin = vol_phasor['vol_sin']
	# vol_cos = vol_phasor['vol_cos']
	print(vol_sin.shape)
	# print(vol_sin.dtype)

	def simple_denoise0(im, th):
		im[im < np.max(im) * th] = 0
		return im

	def simple_denoise(im, th):
		mm = np.min(im)
		MM = np.max(im)
		im[im < mm + (MM - mm) * th] = mm
		return im

	# vol = vol_sin
	# vol = vol_cos
	# print(np.min(vol), np.max(vol))
	vol = np.sqrt(vol_sin**2 + vol_cos**2)
	vol = np.flip(vol, 0)
	vol = np.flip(vol, 1)
	vol = vol.transpose(1, 0, 2)
	# print(vol.dtype)

	simple_denoise(vol, th)

	im0 = np.max(vol, 2)
	im1 = np.max(vol, 1)
	im2 = np.max(vol, 0)

	plt.figure('result', figsize = (15, 4.5))
	plt.subplot(131)
	plt.title('front', fontsize=300)
	if axis == 0:
		plt.axis('off')
	plt.imshow(im0, cmap = cmap)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.subplot(132)
	plt.title('side')
	if axis == 0:
		plt.axis('off')
	plt.imshow(im1, cmap = cmap)
	plt.xlabel('z')
	plt.ylabel('y')
	plt.subplot(133)
	plt.title('top')
	if axis == 0:
		plt.axis('off')
	plt.imshow(im2, cmap = cmap)
	plt.xlabel('z')
	plt.ylabel('x')
	if axis == 0:
		plt.savefig('phasor_result.png')
	plt.show()
