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

rho_avg = mean(rho, axis = (2))

#mental check to find the time with no convection
deltime1 = time[280]
deltime2 = time[500]
print(deltime1, deltime2)

time_nc = time[280:500] # time inbetween convection events

rho_lev1 = rho_avg[:,0] #170 bar
rho_lev1_nc = rho_lev1[280:500]

# Plotting 1D density versus time 
fig,ax=subplots(1,1)

plot(time, rho_lev1, label = 'level 1 (170 bar)')
title(r'Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('density_timeseries.png',bbox_inches= 'tight')

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
    beta[i] = linfit[0]

z = arange(0,48,1)

fig,ax=subplots(1,1)
plot(beta,z)
title(r'$\frac{\delta \rho}{\delta t}$')
ax.set_xlabel(r'$\frac{\delta \rho}{\delta t}$',fontsize = 12)
ax.set_ylabel('z',fontsize = 12)
savefig('noconvec_density_variation.png',bbox_inches= 'tight')

#Coefficents
lev1_linfit = polyfit(time_nc,rho_lev1_nc,1)

#Lines
y1 = lev1_linfit[0]*time_nc + lev1_linfit[1]

fig,ax = subplots(1,1)

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
plot(mass_flux_avg,z)
title(r'$(w \rho)$')
ax.set_ylabel('z',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
savefig('mass_flux_nc.png', bbox_inches = 'tight')


fig,ax= subplots(1,1)
plot(mass_flux_avg,z, label = r'$(w \rho)$')
plot(beta,z ,label = r'$ \frac{\delta \rho}{\delta t}$')
ax.set_ylabel('z',fontsize = 12)
#ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
legend()
savefig('mf_rhovar_nc.png', bbox_inches = 'tight')


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
