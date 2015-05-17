""" ------------------------------
bior2_6
cwt
  by BEAR, 05/04/14
------------------------------ """

import scipy as sp
import numpy as np
from scipy.signal import convolve
#import pywt

_scale_max = 1024
_scale_max = int(2**(sp.ceil(sp.log2(_scale_max))))
tmp = np.loadtxt('bior2.6_1024.txt')
_x_bior2_6   = tmp[:,0]
_psi_bior2_6 = tmp[:,1]
#_, _psi_bior2_6, _, _, _x_bior2_6 = pywt.Wavelet('bior2.6').wavefun(sp.log2(_scale_max))

def bior2_6(length, width):
    length = int(length)
    width = int(width)
    i = sp.arange(0, 13*width)
    u = _psi_bior2_6[_scale_max*i/width]/sp.sqrt(width)
    n = int(abs((length-width*13)/2))
    if length > width*13:
        u = sp.concatenate((u,sp.zeros(length-width*13)), axis=0)
        u = sp.roll(u, n)
    elif length < width*13:
        u = u[n:n+length]
 
    return u
    

def cwt(x, scales, wname, bplot=False):
    coefs = sp.zeros((len(scales), len(x)))
    for i in range(0, len(scales)):
        if wname == 'bior2.6':
            length = min(13*scales[i], len(x))
            wavelet = bior2_6
        coefs[i-1, :] = convolve(x, wavelet(length, i), mode='same')
    
    if bplot:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        plt.ion()
        fig = plt.figure(num=None, figsize=(14,5), dpi=100, facecolor='w', edgecolor='k')
        plt.clf()
        gs = gridspec.GridSpec(3, 1)
        ax1 = fig.add_subplot(gs[0,0])
        ax1.plot(x,'b-')
    
        ax2 = fig.add_subplot(gs[1:,0])
        im = ax2.imshow(coefs[::-1,:], extent=[0, len(x), scales[0], scales[-1]], aspect='auto', cmap='jet')
        ax2.invert_yaxis()
        ax2.set_xlabel('t')
        ax2.set_ylabel('scale')
        l, b, w, h = ax2.get_position().bounds
        cax = fig.add_axes([l+w+0.01, b, 0.02, h])
        plt.colorbar(im, cax=cax)
        plt.suptitle('cwt by python')
        plt.draw()
        plt.show(block=True)
            
    return coefs
