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

#Creating Levels
#Seem to be getting results where the 1 bar and 170 bar density levels seems similar. This doesn't make sense to me since the 2nd level is different. The 

rho_lev1 = rho_avg[:,0] #170 bar
rho_lev1_nc = rho_lev1[280:500]
rho_lev2 = rho_avg[:,24] #85 bar
rho_lev2_nc = rho_lev2[280:500]
rho_lev3 = rho_avg[:,47] #1 bar
rho_lev3_nc = rho_lev3[280:500]

# Plotting 1D density versus time 
fig,ax=subplots(1,1)

plot(time, rho_lev1, label = 'level 1 (170 bar)')
#plot(time, rho_lev2, label = 'level 2 (85 bar)')
#plot(time, rho_lev3, label = 'level 3 (1 bar)')
title(r'Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('density_timeseries.png',bbox_inches= 'tight')


fig,ax=subplots(1,1)

plot(time_nc, rho_lev1_nc, label = 'level 1 (170 bar)')
#plot(time_nc, rho_lev2_nc, label ='level 2 (85 bar)')
#plot(time_nc, rho_lev3_nc, label ='level 3 (1 bar)')
title(r'Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('noconvec_density_timeseries.png',bbox_inches= 'tight')

#Linear Fits
print(rho_lev3_nc.shape)

#Coefficents
lev1_linfit = polyfit(time_nc,rho_lev1_nc,1)
lev2_linfit = polyfit(time_nc,rho_lev2_nc,1)
lev3_linfit = polyfit(time_nc,rho_lev3_nc,1)

#Lines
y1 = lev1_linfit[0]*time_nc + lev1_linfit[1]
y2 = lev1_linfit[0]*time_nc + lev2_linfit[1]
y3 = lev3_linfit[0]*time_nc + lev3_linfit[1]


fig,ax = subplots(1,1)

plot(time_nc, y1, label = 'level 1 (170 bar)')
#plot(time_nc, y2, label ='level 2 (85 bar)')
#plot(time_nc, y3, label ='level 3 (1 bar)')
title(r'Linear Fits of Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
savefig('linfits_of_density_timeseries.png',bbox_inches= 'tight')

##Old Method

#Plotting Density Variations

#p_avg = mean(pres, axis = (0,2))
#pres, time = meshgrid(p_avg,time)

#fig,ax=subplots(1,1)

#h = ax.contourf(time, pres,rhoa_avg)
#title(r'$\frac{\delta \rho}{\delta t}$')
#ax.set_xlabel(r'Time(s)',fontsize = 12)
#ax.set_ylabel(r'Pressure',fontsize = 12)
#ax.set_ylim(max(p_avg), 10**6)
#ax.set_yscale('log')
#fig.colorbar(h)
#savefig('contour_plot_density_deep.png', bbox_inches = 'tight')

#Plotting the Vertical Mass Flux
mass_flux_lev1 =mass_flux[:,0] #170 bar
mass_flux_lev1_nc = mass_flux_lev1[280:500]
mass_flux_lev2 =mass_flux[:,24] #85 bar
mass_flux_lev2_nc = mass_flux_lev1[280:500]
mass_flux_lev3 =mass_flux[:,47] #1 bar
mass_flux_lev3_nc = mass_flux_lev1[280:500]


fig,ax= subplots(1,1)
plot(time_nc, mass_flux_lev1_nc, label = 'level 1 (170 bar)')
plot(time_nc, mass_flux_lev2_nc, label ='level 2 (85 bar)')
plot(time_nc, mass_flux_lev3_nc, label ='level 3 (1 bar)')
title(r'$\frac{\delta(w \rho)}{\delta z}$')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
savefig('mass_flux_nc.png', bbox_inches = 'tight')


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
