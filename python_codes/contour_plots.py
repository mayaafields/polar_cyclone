#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset

# Retrieving Data from the nc files
fname = '/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc'
data = Dataset(fname, 'r')

# Defining Variables
time = data['time']
pres =data['press'][:,:,:,0]
rho = data['rho'][:,:,:,0]
temp = data['temp'][:,:,:,0]
print(temp.shape)
v1 = data['vel1'][:,:,:,0]
v2 = data['vel2'][:,:,:,0]
kinetic_energy = (1/2)*(v1**2+v2**2)

qH2Ov = data['vapor1'][:,:,:,0]
qH2Oc = data['H2Op1'][:,:,:,0]
qH2Op = data['H2Op2'][:,:,:,0]

qNH3v = data['vapor2'][:,:,:,0]
qNH3c = data['NH3p1'][:,:,:,0]
qNH3p = data['NH3p2'][:,:,:,0]

rhoNH3v = rho*qNH3v
rhoH2Ov = rho*qH2Ov

dry_air = rho - (rhoNH3v + rhoH2Ov)
qH2O = rhoH2Ov + qH2Oc + qH2Op
qNH3 = rhoNH3v + qNH3c + qNH3p
print(qH2O.shape)

# Finding total mass
dry_air_t = sum(dry_air, axis = (1,2))
qH2Ot = sum(qH2O, axis = (1,2))
qNH3t = sum(qNH3, axis = (1,2))

# Deviations of Temperature and Concentrations
ntime, nx1, nx2 = temp.shape

t_0 = temp[0,:,0]
tempa = temp - t_0.reshape((1,nx1,1))
tempa_avg = mean(tempa, axis =2) 

# NH3 Variations

#anomaly from horizontal average
qNH3_avg = mean(qNH3, axis = (0,2))
qNH3a1 = qNH3 - qNH3_avg.reshape((1,nx1,1))

#anomaly from initial qNH3
qNH3_0 = qNH3[0,:,0]
qNH3a2 = qNH3 - qNH3_0.reshape((1,nx1,1))
print(qNH3a2.shape)

#H2O Variations

#anomaly from initial qH20
qH2O_0 = qH2O[0,:,0]
qH2Oa2 = qH2O - qH2O_0.reshape((1,nx1,1))

#anomaly from horizontal average 
qH2O_avg = mean(qH2O, axis = (0,2))
qH2Oa1 = qH2O - qH2O_avg.reshape((1,nx1,1))
print(qH2Oa1.shape)
# Plotting 2D contour plots for ammonia and water vapor

p_avg = mean(pres, axis = (0,2))
Pressure, Time = meshgrid(p_avg,time)
print(Pressure.shape)
fig,ax=subplots(1,1)

h = ax.contourf(Time, Pressure, tempa_avg)
ax.set_title('Temperature Variations')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), 10**6)
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_tempa_deep.png')

fig,ax=subplots(1,1)

h = ax.contourf(Time, Pressure, sum(qH2Oa2, axis = 2))
ax.set_title(r'$H_2O$ Concentration')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), 10**6)
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_qH2Oa2_deep.png')

fig,ax=subplots(1,1)

h = ax.contourf(Time, Pressure, sum(qNH3a2, axis = 2))
ax.set_title(r'$NH_3$ Concentration')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), 10**6)
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_qNH3a2_deep.png')

fig,ax=subplots(1,1)
h = ax.contourf(Time, Pressure, kinetic_energy[:,:,0])
ax.set_title('Mean Horizontal Zonal Velocity')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_KE.png')

show()
