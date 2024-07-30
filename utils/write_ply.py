import numpy as np
import sys

def write_ply(points, path):
	f = open(path, 'w')
	header='''ply
	format ascii 1.0
	comment by HY
	element vertex '''
	f.write(header)
	f.write(str(points.shape[0]))
	header='''
	property float x
	property float y
	property float z
	end_header
	'''
	f.write(header)

	for i in range(int(points.shape[0])):
		f.write(str(points[i, 0]) + '\t' + str(points[i, 1]) + '\t' + str(points[i, 2]) + '\n')
	f.close()
