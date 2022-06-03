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
x1 = data['x1']
x2 = data['x2']
rho_avg = mean(rho, axis = (2))
pres_avg = mean(pres, axis = (0,2))

#Time Periods determined by visually identifying convection anomalies
#Period before 1st convection event
t1s = 0
t1e = 270

#Period inbetween convection events
t2s = 280
t2e = 500

#Period after 2nd convection event
t3s = 510
t3e = 800

#Defining time with no convection
t1 = time[t1s:t1e]
t2 = time[t2s:t2e]
t3 = time[t3s:t3e]

#Linear Fits
beta1 = np.zeros(48)
beta2 = np.zeros(48)
beta3 = np.zeros(48)
for i in range(nx1):
    rho = rho_avg[:,i]
    rho1 = rho[t1s:t1e]
    rho2 = rho[t2s:t2e]
    rho3 = rho[t3s:t3e]
    linfit1 = polyfit(t1,rho1,1)
    linfit2 = polyfit(t2,rho2,1)
    linfit3 = polyfit(t3,rho3,1)
    #beta is an array of d\rho/dt values for each pressure level
    beta1[i] = linfit1[0]
    beta2[i] = linfit2[0]
    beta3[i] = linfit3[0]

#d\rho/dt plot 
fig,ax=subplots(1,1)
plot(beta1,pres_avg, label = r'$\frac{d(\rho)_1}{dt}$')
plot(beta2,pres_avg, label = r'$\frac{d(\rho)_2}{dt}$')
plot(beta3,pres_avg, label = r'$\frac{d(\rho)_3}{dt}$')
title(r'$\frac{d\rho}{dt}$')
ax.set_xlabel(r'$\frac{d\rho}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('drhodt.png',bbox_inches= 'tight')

#Defining the Vertical Mass Flux per time period
mass_flux1 = mass_flux[t1s:t1e,:]
mass_flux1_avg = mean(mass_flux1, axis = 0)

mass_flux2 = mass_flux[t2s:t2e,:]
mass_flux2_avg = mean(mass_flux2, axis = 0)

mass_flux3 = mass_flux[t3s:t3e,:]
mass_flux3_avg = mean(mass_flux3, axis = 0)

#plotting the vertical mass flux
fig,ax= subplots(1,1)
plot(mass_flux1_avg, pres_avg, label = r'$(w\rho)_1$')
plot(mass_flux2_avg, pres_avg, label = r'$(w\rho)_2$')
plot(mass_flux3_avg, pres_avg, label = r'$(w\rho)_3$')
title(r'$(w \rho)$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
savefig('mass_flux.png', bbox_inches = 'tight')

#Moving Average of Mass Flux
N = 5 #5 point moving average

mass_flux1 = uniform_filter1d(mass_flux1, size=N)
mass_flux1_avg = mean(mass_flux1, axis =0)
dwrhodz1 = (mass_flux1_avg[1:] - mass_flux1_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux2 = uniform_filter1d(mass_flux2, size=N)
mass_flux2_avg = mean(mass_flux2, axis =0)
dwrhodz2 = (mass_flux2_avg[1:] - mass_flux2_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux3 = uniform_filter1d(mass_flux3, size=N)
mass_flux3_avg = mean(mass_flux3, axis =0)
dwrhodz3 = (mass_flux3_avg[1:] - mass_flux3_avg[:-1])/(x1[1:] - x1[:-1])

pres_avg = sqrt(pres_avg[:1]*pres_avg[:-1])

fig,ax= subplots(1,1)
plot(-1*dwrhodz1,pres_avg,label = r'$\frac{-d(w \rho)_1}{dz}$')
plot(-1*dwrhodz2,pres_avg,label = r'$\frac{-d(w \rho)_2}{dz}$')
plot(-1*dwrhodz3,pres_avg,label = r'$\frac{-d(w \rho)_3}{dz}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
ax.set_xlim(-0.5e-9,0.5e-9)
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg.png', bbox_inches = 'tight')


show()
