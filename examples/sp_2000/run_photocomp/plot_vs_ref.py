# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots a mean profile of the variable in the supplied file
against height (calculated from pressure assuming isothermal atmos.)
If channel dimension exists, spectrum at boundaries is also plotted.
'''

CONST = 287. * 250. / (9.80665 * 1000.0) #

from netCDF4 import Dataset
import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if (len(sys.argv) > 4):
        filename4 = sys.argv[4]
        name4 = filename4[filename4.find('.') + 1:]
        dgs4 = Dataset(filename4)
        var4 = dgs4.variables[name4][:]
    if (len(sys.argv) > 3):
        filename3 = sys.argv[3]
        name3 = filename3[filename3.find('.') + 1:]
        dgs3 = Dataset(filename3)
        var3 = dgs3.variables[name3][:]
    if (len(sys.argv) > 2):
        filename2 = sys.argv[2]
        name2 = filename2[filename2.find('.') + 1:]
        dgs2 = Dataset(filename2)
        var2 = dgs2.variables[name2][:]
    if (len(sys.argv) > 1):
        filename = sys.argv[1]
        name = filename[filename.find('.') + 1:]
        dgs = Dataset(filename)
        var = dgs.variables[name][:]
    else:
        raise RuntimeError('please enter a file name')

lon = dgs.variables['lon'][:]
lat = dgs.variables['lat'][:]
p = dgs.variables['plev'][:]

n_lon = len(lon)
n_lat = len(lat)
layers= len(p)
e_factor=2.718

if ('ph_rate' in name):
    e_factor=1.0e2
    xlim = 0.412

#if (len(sys.argv) > 2):
#    var = var - var2

try:
    width     = dgs.variables['bandwidth'][:]
    wl_short  = dgs.variables['wl_short'][:]
    wl_long   = dgs.variables['wl_long'][:]
    n_channel = len(width)
except:
    vmean = np.zeros(layers)
    for i in np.arange(layers):
        vmean[i] = np.sum(var[i, :, :]) / (n_lon * n_lat)
    ax1 = plt.figure().add_subplot(111)
    ax1.plot(vmean, p, label=name)
    ax1.invert_yaxis()
    ax1.set_yscale('log')
    ax1.set_title(dgs.variables[name].title)
    if ('ph_rate' in name):
        ax1.set_xlabel('J (s$^{-1}$)')
        ax1.set_xscale('log')
    elif (name == 'aflx'):
        ax1.set_xlabel('Actinic Flux (W m$^{-2}$)')
    elif (name == 'hrts'):
        ax1.set_xlabel('Heating rate (K day$^{-1}$)')
    else:
        ax1.set_xlabel('Flux (W m$^{-2}$)')
    ax1.set_ylabel('Pressure (Pa)')
    if (len(sys.argv) > 2):
        for i in np.arange(layers):
            vmean[i] = np.sum(var2[i, :, :]) / (n_lon * n_lat)
        ax1.plot(vmean, p, linestyle='dashed', label='Ref')
    plt.legend()
else:
    width2     = dgs2.variables['bandwidth'][:]
    wl_short2  = dgs2.variables['wl_short'][:]
    wl_long2   = dgs2.variables['wl_long'][:]
    n_channel2 = len(width2)
    fig=plt.figure()
    vmean = np.zeros(layers)
    vtop = np.sum(var[:, 0, :, :]) / (n_lon * n_lat)
    for i in np.arange(layers):
        vmean[i] = np.sum(var[:, i, :, :]) / (n_lon * n_lat)
        if (vmean[i] > vtop/e_factor):
          e_layer = i
    e_layer = 40
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(vmean, p)
    ax1.invert_yaxis()
    ax1.set_yscale('log')
    ax1.set_title(dgs.variables[name].title)
    ax1.set_xlabel(name)
    ax1.set_ylabel('Pressure (Pa)')
    e_height = p[e_layer]
    ax1.plot([min(vmean),max(vmean)], [e_height,e_height], color='green')
    ax1.set_xscale('log')
    if (len(sys.argv) > 2):
        for i in np.arange(layers):
            vmean[i] = np.sum(var2[:, i, :, :]) / (n_lon * n_lat)
        ax1.plot(vmean, p, linestyle='dashed', label='Ref')
    wn = 0.5e-2/wl_short + 0.5e-2/wl_long
    wl = 0.5e6*(wl_short + wl_long)
    wl2 = 0.5e6*(wl_short2 + wl_long2)
    toa_spec = np.zeros(n_channel)
    surf_spec = np.zeros(n_channel)
    mid_spec = np.zeros(n_channel)
    mid_spec2 = np.zeros(n_channel2)
    for ch in range(0, n_channel):
        toa_spec[ch]  = np.sum(var[ch,0, :,:])/(width[ch]*n_lon*n_lat)
        surf_spec[ch] = np.sum(var[ch,14,:,:])/(width[ch]*n_lon*n_lat)
        mid_spec[ch]  = np.sum(var[ch,e_layer,:,:])/(width[ch]*n_lon*n_lat)
    for ch in range(0, n_channel2):
        mid_spec2[ch]  = np.sum(var2[ch,e_layer, :,:])/(width2[ch]*n_lon*n_lat)
#    ax2.plot(wl, toa_spec, color='blue', label='TOA')
#    ax2.set_xlabel('Wavelength (micron)')
    ax2.plot(wl, mid_spec, color='blue', label='Mid atmos')
    ax2.plot(wl2, mid_spec2, color='green', label='Ref')
    if ('ph_rate' in name):
        ax2.set_ylabel('J (s$^{-1}$ m$^{-1}$)')
        ax2.set_yscale('symlog')
        ax2.set_title('Mid atmosphere spectra')
    elif (name == 'aflx'):
        ax2.set_ylabel('Actinic Flux (W m$^{-2} m$^{-1})')
        ax2.set_yscale('symlog')
        ax2.set_title('Mid atmosphere spectra')
    elif (name == 'hrts'):
        ax2.set_ylabel('Heating rate (K day$^{-1} m$^{-1})')
        ax2.set_title('Mid atmosphere spectra')
    else:
        ax2.set_title('Mid atmosphere spectrum')
        ax2.set_ylabel('Flux (Wm-2m-1)')
        ax2.set_xscale('symlog')
    ax2.set_xlim(left=0.1,right=5.0)
    ax2.set_xlabel('Wavelength (micron)')
    plt.legend()

plt.tight_layout()
plt.show()
