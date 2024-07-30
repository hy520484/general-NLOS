import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt

data = si.loadmat('sub_Sig.mat')['sub_Sig']
print(data.shape)	# (100, 512)

hist_tot = np.sum(data, 0)
plt.plot(hist_tot)
plt.show()
