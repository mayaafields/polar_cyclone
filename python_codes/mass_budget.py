#Mass Budget Python Script



#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset

# Retrieving Data from the nc files
fname = '/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc'
data = Dataset(fname, 'r')

#Defining Variables
time = data['time']
rho = data['rho'][:,:,:,0]
pres = data['press'][:,:,:,0]
mass_flux = data['v1rho'][:,:]
ntime, nx1, nx2 = rho.shape

rho_0 = rho[0,:,0]
density_variations = rho - rho_0.reshape((1,nx1,1))


rhoa_avg = mean(density_variations, axis =2)


#heat_rate= rhoa_avg + mass_flux
# this formula is incorrect 

#Plotting Density Variations

p_avg = mean(pres, axis = (0,2))
pres, time = meshgrid(p_avg,time)

fig,ax=subplots(1,1)

h = ax.contourf(time, pres,rhoa_avg)
title(r'$\frac{\delta \rho}{\delta t}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_density.png', bbox_inches = 'tight')

#Plotting the Vertical Mass Flux

fig,ax= subplots(1,1)
h = ax.contourf(time, pres, mass_flux)
title(r'$\frac{\delta(w \rho)}{\delta z}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_mass_flux.png', bbox_inches = 'tight')


#Plotting the Heating Rate due to condensation? 
# unsure what the term is called without the Latent heat
# formula above is incorrect. Will replot this when I update my understanding. 


#fig,ax= subplots(1,1)
#h = ax.contourf(Time, Pressure, heat_rate)
#title(r'$-<\rho \dot{q_c}>$')
#ax.set_xlabel(r'Time(s)',fontsize = 12)
#ax.set_ylabel(r'Pressure',fontsize = 12)
#ax.set_ylim(max(p_avg), min(p_avg))
#ax.set_yscale('log')
#fig.colorbar(h)
#savefig('contour_plot_heat_rate.png', bbox_inches = 'tight')



show()
