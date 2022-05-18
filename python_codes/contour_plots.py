#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset

# Retrieving Data from the nc files
fname = '/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc'
data = Dataset(fname, 'r')

# Defining Variables
time = data['time']
pressure =data['press'][:,:,:,0]
rho_gas = data['rho'][:,:,:,0]
temp = data['temp'][:,:,:,0]

horiz_v1 = data['vel1'][:,:,:,0]
horiz_v2 = data['vel2'][:,:,:,0]
kinetic_energy = (1/2)*(horiz_v1**2+horiz_v2**2)

H2Ov = data['vapor1'][:,:,:,0]
H2Oc = data['H2Op1'][:,:,:,0]
H2Op = data['H2Op2'][:,:,:,0]

NH3v = data['vapor2'][:,:,:,0]
NH3c = data['NH3p1'][:,:,:,0]
NH3p = data['NH3p2'][:,:,:,0]

dry_air = rho_gas*H2Ov + rho_gas*NH3v
H2O_total = rho_gas*H2Ov + H2Oc + H2Op
NH3_total = rho_gas*NH3v + NH3c + NH3p

# Finding total mass
dry_air_t = sum(dry_air, axis = (1,2))
H2Ot = sum(H2O_total, axis = (1,2))
NH3t = sum(NH3_total, axis = (1,2))


# Deviations of Temperature and Concentrations
ntime, nx1, nx2 = temp.shape

t_0 = temp[0,:,0]
temp_variations = temp - t_0.reshape((1,nx1,1))

# NH3 Variations

NH3avg = mean(NH3_total, axis = (0,2))
NH3_variations = NH3_total - NH3avg.reshape((1,nx1,1))

NH3_0 = NH3_total[0,:,0]
NH3_variations_from_initial = NH3_total - NH3_0.reshape((1,nx1,1))


#H2O Variations

H2O_0 = H2O_total[0,:,0]
H2O_variations_from_initial = H2O_total - H2O_0.reshape((1,nx1,1))

H2Oavg = mean(H2O_total, axis = (0,2))
H2O_variations = H2O_total - H2Oavg.reshape((1,nx1,1))

# Plotting 2D contour plots for ammonia and water vapor

pavg = mean(pressure, axis = (0,2))
Pressure, Time = meshgrid(pavg,time)

fig,ax=subplots(1,1)

contour_plot_temp_var = ax.contourf(Time, Pressure, mean(temp_variations, axis =2))
title(r'Temperature Variations')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), 10**6)
ax.set_yscale('log')
fig.colorbar(contour_plot_temp_var)
savefig('contour_plot_temp_deep.png')

fig,ax=subplots(1,1)

contour_plot_H2O = ax.contourf(Time, Pressure, sum(H2O_variations_from_initial, axis = 2))
title(r'$H_2O$ Concentration')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), 10**6)
ax.set_yscale('log')
fig.colorbar(contour_plot_H2O)
savefig('contour_plot_H2O_deep.png')

fig,ax=subplots(1,1)

contour_plot_NH3 = ax.contourf(Time, Pressure, sum(NH3_variations_from_initial, axis = 2))
title(r'$NH_3$ Concentration')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), 10**6)
ax.set_yscale('log')
fig.colorbar(contour_plot_NH3)
savefig('contour_plot_NH3_deep.png')

fig,ax=subplots(1,1)
contour_plot_KE = ax.contourf(Time, Pressure, kinetic_energy[:,:,0])
ax.set_title(r'Mean Horizontal Zonal Velocity')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(max(pavg), min(pavg))
ax.set_yscale('log')
fig.colorbar(contour_plot_KE)
savefig('contour_plot_KE.png')


#Gradient Plots with respect to time

grad_NH3 = gradient(NH3_variations_from_initial)
print(grad_NH3)

show()
