import os
import sys
import configparser
from configparser import ConfigParser, ExtendedInterpolation

def parse_args(config_path):
	if not os.path.exists(config_path):
		print('No such file: {}'.format(config_path))
		sys.exit('No File Error')

	args = {}
	config = ConfigParser(interpolation = ExtendedInterpolation())
	config.read(config_path)

	args['path_data'] = config.get('params', 'path_data')
	args['path_setup'] = config.get('params', 'path_setup')

	args['width'] = float(eval(config.get('params', 'width')))
	
	args['N_ld'] = int(eval(config.get('params', 'N_ld')))
	args['N_bin'] = config.getint('params', 'N_bin')
	args['timeRes'] = config.getfloat('params', 'timeRes')

	args['N_voxel_s'] = config.getint('params', 'N_voxel_s')
	args['N_voxel_t'] = config.getint('params', 'N_voxel_t')

	args['x_range0'] = config.getfloat('params', 'x_range0')
	args['x_range1'] = config.getfloat('params', 'x_range1')
	args['y_range0'] = config.getfloat('params', 'y_range0')
	args['y_range1'] = config.getfloat('params', 'y_range1')
	args['z_range0'] = config.getfloat('params', 'z_range0')
	args['z_range1'] = config.getfloat('params', 'z_range1')
	try:
		args['th'] = config.getfloat('params', 'th')
	except:
		args['th'] = 0
	try:
		args['sigma'] = config.getfloat('params', 'sigma')
	except:
		args['sigma'] = 0.23

	return args
