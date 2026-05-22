# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots LBL cross-sections from corr_k netcdf files
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':
    if (len(sys.argv) > 2):
        filename1 = sys.argv[1]
        filename2 = sys.argv[2]
    elif (len(sys.argv) > 1):
        filename1 = sys.argv[1]
    else:
        raise RuntimeError('please enter a corr_k netcdf format LBL file')

    xsc1 = Dataset(filename1)
    n_nu = len(xsc1.dimensions['nu'])
    n_pt_pair = len(xsc1.dimensions['pt_pair'])
    p_calc = xsc1.variables['p_calc'][:]
    t_calc = xsc1.variables['t_calc'][:]
    nu = xsc1.variables['nu'][:]
    wl=1.0e9/nu
    kabs = xsc1.variables['kabs'][:,:]

    fig=plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(111)
    ax1.plot(wl, kabs[0,:], color='green')
    ax1.set_xlabel('Wavelength (nm)')
#    ax1.plot(nu, kabs[475,:], color='green')
#    ax1.set_xlabel('Wavenumber (m$^{-1}$)')
    ax1.set_ylabel('Absorption (m$^{2}$ kg$^{-1}$)')
    ax1.set_yscale('symlog')

#    leg=plt.legend()
#    leg.set_draggable(True)
    plt.tight_layout()
    plt.show()
