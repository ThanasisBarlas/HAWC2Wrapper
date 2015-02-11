# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 14:58:39 2014

@author: tlbl
"""
import numpy as np
from scipy.interpolate import pchip


def sharp_curves_interpolation(x, y, x_i, x_max):

    for ixx, xx in enumerate(x):
        if (xx-x_max)>0:
            break
    for ixx_i, xx_i in enumerate(x_i):
        if (xx_i-x_max)>0:
            break
    ppf_l = np.polyfit(x[:ixx], y[:ixx], 2)
    y_i = np.zeros(len(x_i))
    y_i[:ixx_i+1] = np.polyval(ppf_l, x_i[:ixx_i+1])

    x_l = np.array([x_i[ixx_i]]+x[ixx:].tolist()) 
    y_l = np.array([y_i[ixx_i]]+y[ixx:].tolist())

    y_i[ixx_i:] = np.interp(x_i[ixx_i:], x_l, y_l)

    return y_i


if __name__ == "__main__":
    import matplotlib.pylab as plt

  
    a =np.array([[  4.      ,  -5.61355 ],
              [  6.      ,  -4.30368 ],
           [  8.      ,  -2.89396 ],
           [ 10.      ,  -0.636955],
           [ 11.      ,   0.625755],
           [ 12.      ,  -0.743788],
           [ 16.      ,  -4.07054 ],
           [ 20.      ,  -5.657   ],
           [ 25.      ,  -6.97827 ]])
    """
    a =np.array([[  4.      ,  -5.61355 ],
              [  6.      ,  -4.30368 ],
           [  8.      ,  -2.89396 ],
           [ 10.      ,  -0.636955],
           [ 11.      ,  -0.425755],
           [ 12.      ,  -0.743788],
           [ 16.      ,  -4.07054 ],
           [ 20.      ,  -5.657   ],
           [ 25.      ,  -6.97827 ]])
    
    a = np.array([[  6.00000000e+00,   4.95695264e+05],
       [  8.00000000e+00,   7.93284671e+05],
       [  1.00000000e+01,   1.23950730e+06],
       [  1.10000000e+01,   1.49980383e+06],
       [  1.20000000e+01,   1.28264278e+06],
       [  1.40000000e+01,   9.75256053e+05]])

           
    
    
    a =np.array([[  4.      ,  -5.61355 ],
              [  6.      ,  -4.30368 ],
           [  8.      ,  -2.89396 ],
           [ 10.      ,  -0.636955],
           [ 12.      ,  -0.743788],
           [ 16.      ,  -4.07054 ],
           [ 20.      ,  -5.657   ],
           [ 25.      ,  -6.97827 ]])
           
    a =np.array([[  4.      ,  -5.61355 ],
              [  6.      ,  -4.30368 ],
           [  8.      ,  -2.89396 ],
           [ 10.      ,  -0.836955],
           [ 12.      ,  -0.743788],
           [ 16.      ,  -4.07054 ],
           [ 20.      ,  -5.657   ],
           [ 25.      ,  -6.97827 ]])
    
    """
    x_i = np.linspace(a[0, 0], a[-1, 0], 200)
    y_i = sharp_curves_interpolation(a[:, 0], a[:, 1], x_i, 11.5)
    plt.figure(10)
    plt.plot(a[:, 0], a[:, 1], 'o')
    plt.plot(x_i, y_i, '-')
    plt.show()
