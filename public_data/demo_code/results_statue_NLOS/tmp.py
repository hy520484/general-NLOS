import numpy as np
import scipy.io as si
import matplotlib.pyplot as plt

# a = si.loadmat('alb.mat')
# vol = a['alb']
# print(vol.shape)

# # vol = u[:, :, :, 2]
# im = np.sum(vol, 0)
# im = np.flip(im, 0)
# im = np.rot90(im)


# plt.imshow(im, cmap = 'hot')
# plt.axis('off')
# plt.show()

ttt = si.loadmat('ttt.mat')['ttt']
plt.imshow(ttt, cmap = 'hot')
plt.axis('off')
plt.show()
