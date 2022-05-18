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
mass_flux = data['v1rho'][:,:]
ntime, nx1, nx2 = rho_gas.shape

rho_0 = rho_gas[0,:,0]
density_variations = rho_gas - rho_0.reshape((1,nx1,1))


horizontally_averaged_density_var = mean(density_variations, axis =2)


heat_rate= horizontally_averaged_density_var + mass_flux

#Plotting Density Variations

pavg = mean(pressure, axis = (0,2))
Pressure, Time = meshgrid(pavg,time)

fig,ax=subplots(1,1)

contour_plot_density_var = ax.contourf(Time, Pressure,horizontally_averaged_density_var)
title(r'$\frac{\delta \rho}{\delta t}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), min(pavg))
ax.set_yscale('log')
fig.colorbar(contour_plot_density_var)
savefig('contour_plot_density.png')

#Plotting the Vertical Mass Flux

fig,ax= subplots(1,1)
contour_plot_mass_flux = ax.contourf(Time, Pressure, mass_flux)
title(r'$\frac{\delta(w \rho)}{\delta z}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), min(pavg))
ax.set_yscale('log')
fig.colorbar(contour_plot_mass_flux)
savefig('contour_plot_mass_flux.png')


#Plotting the Heating Rate due to condensation? 
# unsure what the term is called without the Latent heat

fig,ax= subplots(1,1)
contour_plot_heat_rate = ax.contourf(Time, Pressure, heat_rate)
title(r'$-<\rho \dot{q_c}>$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), min(pavg))
ax.set_yscale('log')
fig.colorbar(contour_plot_heat_rate)
savefig('contour_plot_heat_rate.png')



show()
