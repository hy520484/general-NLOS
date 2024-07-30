import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt

laserpoints = si.loadmat('laserpoints.mat')['laserpoints']
detectpoints = si.loadmat('detectpoints.mat')['detectpoints']
si.savemat('setup.mat', {'laserpoints': laserpoints, 'detectpoints': detectpoints})
