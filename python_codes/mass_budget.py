#Mass Budget Python Script

#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset
from scipy.ndimage.filters import uniform_filter1d

# Retrieving Data from the nc files
fname = '/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc'
data = Dataset(fname, 'r')

#Defining Variables
time = data['time']
rho = data['rho'][:,:,:,0]
pres = data['press'][:,:,:,0]
mass_flux = data['v1rho'][:,:]
ntime, nx1, nx2 = rho.shape

rho_avg = mean(rho, axis = (2))
p_avg = mean(pres, axis = (0,2))

#mental check to find the time with no convection
deltime1 = time[280]
deltime2 = time[500]
print(deltime1, deltime2)

#Defining time with no convection
time_nc = time[280:500] # time inbetween convection events

#Creating density @ one pressure level to test method
rho_lev1 = rho_avg[:,0] #170 bar
rho_lev1_nc = rho_lev1[280:500]

# Plotting 1D density versus time 
fig,ax=subplots(1,1)

#Full timeseries 
plot(time, rho_lev1, label = 'level 1 (170 bar)')
title(r'Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('density_timeseries.png',bbox_inches= 'tight')

#Cut timeseries inbetween convection events
fig,ax=subplots(1,1)

plot(time_nc, rho_lev1_nc, label = 'level 1 (170 bar)')
title(r'Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('noconvec_density_timeseries.png',bbox_inches= 'tight')

#Linear Fits
beta = np.zeros(48)
for i in range(nx1):
    rho_level = rho_avg[:,i]
    rho_level_nc = rho_level[280:500]
    linfit = polyfit(time_nc, rho_level_nc,1)
    beta[i] = linfit[0] #beta is an array of d\rho/dt values for each pressure level

#d\rho/dt plot 
fig,ax=subplots(1,1)
plot(beta,p_avg)
title(r'$\frac{\delta \rho}{\delta t}$')
ax.set_xlabel(r'$\frac{\delta \rho}{\delta t}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg),min(p_avg))
savefig('noconvec_density_variation.png',bbox_inches= 'tight')

#Coefficents from linear fit
lev1_linfit = polyfit(time_nc,rho_lev1_nc,1)

#Lines
y1 = lev1_linfit[0]*time_nc + lev1_linfit[1]

fig,ax = subplots(1,1)

#Linear fit of just one pressure level
plot(time_nc, y1, label = 'level 1 (170 bar)')
title(r'Linear Fits of Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('linfits_of_density_timeseries.png',bbox_inches= 'tight')

#Plotting the Vertical Mass Flux
mass_flux_nc = mass_flux[280:500,:]
mass_flux_avg = mean(mass_flux_nc, axis = 0)
print(mass_flux_avg)

fig,ax= subplots(1,1)
plot(mass_flux_avg, p_avg)
title(r'$(w \rho)$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg),min(p_avg))
savefig('mass_flux_nc.png', bbox_inches = 'tight')

#Comparison to d\rho/dt
fig,ax= subplots(1,1)
plot(mass_flux_avg,p_avg, label = r'$(w \rho)$')
plot(beta,p_avg,label = r'$ \frac{\delta \rho}{\delta t}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg),min(p_avg))
legend()
savefig('mf_rhovar_nc.png', bbox_inches = 'tight')



#Moving Average of Mass Flux
N = 5
moving_average_mf = uniform_filter1d(mass_flux_nc, size=N)
print(moving_average_mf.shape)
moving_average_mf_gradz = mean(moving_average_mf, axis =0)


fig,ax= subplots(1,1)
plot(moving_average_mf_gradz,p_avg, label = r'$(w \rho)$')
plot(beta,p_avg,label = r'$ \frac{\delta \rho}{\delta t}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg),min(p_avg))
legend()
savefig('mfmovingavg_rhovar_nc.png', bbox_inches = 'tight')


show()
