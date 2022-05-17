#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset
import numpy as np

'''Retrieving Data from the nc files'''
data = Dataset('/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc')

'''Defining Variables'''
time = data['time']
pressure =data['press']
rho_gas = data['rho'][:,:,:,0]
temp = data['temp']
print(temp.shape)
#mass = rho_gas*volume

horiz_v1 = data['vel1'][:,:,:,0]
horiz_v2 = data['vel2'][:,:,:,0]

print(horiz_v1.shape)

kinetic_energy = (1/2)*(horiz_v1**2+horiz_v2**2)
print(kinetic_energy)
print(kinetic_energy.shape)

#print(time.shape)
#print(pressure.shape)

H2Ov = data['vapor1'][:,:,:,0]
H2Oc = data['H2Op1'][:,:,:,0]
H2Op = data['H2Op2'][:,:,:,0]

#print(H2Ov.shape)
#print(H2Oc.shape)
#print(H2Op.shape)

NH3v = data['vapor2'][:,:,:,0]
NH3c = data['NH3p1'][:,:,:,0]
NH3p = data['NH3p2'][:,:,:,0]

dry_air = rho_gas*H2Ov + rho_gas*NH3v
H2O_total = rho_gas*H2Ov + H2Oc + H2Op
NH3_total = rho_gas*NH3v + NH3c + NH3p

'''Finding total mass'''
dry_air_t = sum(dry_air, axis = (1,2))
H2Ot = sum(H2O_total, axis = (1,2))
NH3t = sum(NH3_total, axis = (1,2))


'''Deviations of Temperature and Concentrations'''

'''temperature'''

t_0 = np.zeros((len(temp[0,0,:,0])))

for i in range(len(temp[0,0,:,0])): 
    t_0 = temp[0,0,i,0]
print(t_0)

temp_variations = np.zeros((len(temp[:,0,0,0]),len(temp[0,:,0,0]),len(temp[0,0,:,0])))
for i in range(len(temp[:,0,0,0])):
    for j in range(len(temp[0,:,0,0])):
        for k in range(len(temp[0,0,:,0])):
            temp_variations[i,j,k] = temp[i,j,k,0] - t_0

print(temp_variations)
print(np.mean(temp_variations))
temp_small_changes = temp_variations - np.mean(temp_variations)
print(temp_small_changes)

'''NH3 Variations'''
mean_NH3 = np.mean(NH3_total[:,:,0])
NH3_variations = np.zeros((len(NH3_total[:,0,0]),len(NH3_total[0,:,0])))
for i in range(len(NH3_total[:,0,0])):
    for j in range(len(NH3_total[0,:,0])):
            NH3_variations[i,j] = NH3_total[i,j,0] - mean_NH3

'''H2O Variations'''
mean_H2O = np.mean(H2O_total[:,:,0])
H2O_variations = np.zeros((len(H2O_total[:,0,0]),len(H2O_total[0,:,0])))
for i in range(len(H2O_total[:,0,0])):
    for j in range(len(H2O_total[0,:,0])):
            H2O_variations[i,j] = H2O_total[i,j,0] - mean_H2O



'''Plotting 2D contour plots for ammonia and water vapor'''

Pressure, Time = np.meshgrid(pressure[0,:,0,0],time)

fig,ax=plt.subplots(1,1)

contour_plot_temp_var = ax.contourf(Time, Pressure, temp_variations[:,:])
plt.title(r'Temperature Variations')
plt.xlabel(r'Time(s)',fontsize = 12)
plt.ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(np.max(Pressure), np.min(Pressure))
plt.yscale('log')
fig.colorbar(contour_plot_temp_var)
plt.savefig('contour_plot_temp.png')

fig,ax=plt.subplots(1,1)

contour_plot_H2O = ax.contourf(Time, Pressure, H2O_total[:,:,0])
plt.title(r'$H_2O$ Concentration')
plt.xlabel(r'Time(s)',fontsize = 12)
plt.ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(np.max(Pressure), np.min(Pressure))
plt.yscale('log')
fig.colorbar(contour_plot_H2O)
plt.savefig('contour_plot_H2O.png')

fig,ax=plt.subplots(1,1)

contour_plot_NH3 = ax.contourf(Time, Pressure, NH3_total[:,:,0])
plt.title(r'$NH_3$ Concentration')
plt.xlabel(r'Time(s)',fontsize = 12)
plt.ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(np.max(Pressure), np.min(Pressure))
plt.yscale('log')
fig.colorbar(contour_plot_NH3)
plt.savefig('contour_plot_NH3.png')


fig,ax=plt.subplots(1,1)
contour_plot_KE = ax.contourf(Time, Pressure, kinetic_energy[:,:,0])
plt.title(r'Mean Horizontal Zonal Velocity')
plt.xlabel(r'Time(s)',fontsize = 12)
plt.ylabel(r'Pressure',fontsize = 12)
ax.set_ylim(np.max(Pressure), np.min(Pressure))
plt.yscale('log')
fig.colorbar(contour_plot_KE)
plt.savefig('contour_plot_KE.png')
show()
