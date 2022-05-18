#Mass Budget Python Scrip



#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset

# Retrieving Data from the nc files
fname = '/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc'
data = Dataset(fname, 'r')

#Defining Variables
time = data['time']
rho_gas = data['rho'][:,:,:,0]
pressure = data['press'][:,:,:,0]
mass_flux = data['v1rho']
print(mass_flux.shape)
ntime, nx1, nx2 = rho_gas.shape

rho_0 = rho_gas[0,:,0]
density_variations = rho_gas - rho_0.reshape((1,nx1,1))
#Plotting Density Variations

pavg = mean(pressure, axis = (0,2))
Pressure, Time = meshgrid(pavg,time)


print(Pressure.shape, Time.shape, mass_flux.shape)

fig,ax=subplots(1,1)

contour_plot_density_var = ax.contourf(Time, Pressure, mean(density_variations, axis =2))
title(r'$\frac{\delta \rho}{\delta t}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), min(pavg))
ax.set_yscale('log')
fig.colorbar(contour_plot_density_var)
savefig('contour_plot_density.png')

fig, ax = subplots(1,1)
contour_plot_mass_flux = ax.contourf(Time, Pressure, mass_flux)
title(r'$\frac{\delta \rho}{\delta t}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), min(pavg))
ax.set_yscale('log')
fig.colorbar(contour_plot_mass_flux)
savefig('contour_plot_mass_flux.png')

show()
