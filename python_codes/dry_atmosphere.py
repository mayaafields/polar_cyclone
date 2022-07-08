#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset
from scipy.ndimage.filters import uniform_filter1d
from scipy import signal
# Retrieving Data from the nc files


fname = '/data1/mafields/jupiter_modelruns_data/240622_1000_64256_c8_1e5_dry/jupiter-main.nc'
data = Dataset(fname, 'r')

# Defining Variables
time = data['time']
pres =data['press'][:,:,:,0]
rho = data['rho'][:,:,:,0]
temp = data['temp'][:,:,:,0]
theta = data['theta'][:,:,:,0]
v1 = data['vel1'][:,:,:,0]
v2 = data['vel2'][:,:,:,0]
kinetic_energy = (1/2)*(v1**2+v2**2)
x1 = data['x1']
print(x1)
x2 = data['x2']
x1f = data['x1f']
x2f = data['x2f']
grav_acc1= 24.79

# Finding total mass
dry_air_t = sum(rho, axis = (1,2))

# Deviations of Temperature
ntime, nx1, nx2 = temp.shape

t_0 = temp[0,:,0]
tempa = temp - t_0.reshape((1,nx1,1))
tempa_avg = mean(tempa, axis = (2))
theta_avg2d = mean(theta, axis = (2))
theta1 = theta_avg2d[0]
theta2 = theta_avg2d[1000]
theta3 = theta_avg2d[2000]

tempa_avg_1 = tempa_avg[0]
tempa_avg_2 = tempa_avg[1000]
tempa_avg_3 = tempa_avg[2000]

temp = mean(temp, axis = (2))
temp1 = temp[0]
temp2 = temp[1000]
temp3 = temp[2000]

rho_0 = rho[0,:,0]
rhoa = rho - rho_0.reshape((1,nx1,1))
rhoa_avg = mean(rhoa, axis = (2))

#Moving Average of potential temperature
n  = 5 #5 point moving average

theta_avg = mean(theta, axis = (0,2))

theta_avg = uniform_filter1d(theta_avg, size=n)
dthetadz = (theta_avg[1:] - theta_avg[:-1])/(x1[1:] - x1[:-1])

theta_avg1 = np.array(theta_avg)
theta_avg1 = signal.resample(theta_avg1, 63)

print(theta_avg.shape)

## Brunt Vaisala Frequency is N^2 = (gravity/theta)*(d theta/dz)'
N_squared = (grav_acc1/theta_avg1)*(dthetadz)
N = sqrt(N_squared)
H = 27000

nh = N*H

# Plotting 2D contour plots for ammonia and water vapor

p_avg = mean(pres, axis = (0,2))
p_avg2 = sqrt(p_avg[1:]*p_avg[:-1])
pres, time = meshgrid(p_avg,time)

fig,ax= subplots(1,1)
plot(dthetadz,p_avg2,label = r'$\frac{d(\theta)_1}{dz}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg2),min(p_avg2))
#axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dthetadz_mvavg.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
plot(N_squared,p_avg2,label = r'$N^2 = \frac{g}{\theta}\frac{d(\theta)}{dz}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg2),min(p_avg2))
#axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('brunt_vaisala.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
plot(nh,p_avg2,label = r'$NH =H\cdot \sqrt{\frac{g}{\theta}\frac{d(\theta)}{dz}}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(p_avg2),min(p_avg2))
#axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('NH.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
h = ax.contourf(time, pres, tempa_avg)
ax.set_title('Temperature Variations')
ax.axhline(y = 2.7e4, linestyle = '--', color = 'k', label = '2.7e4 Pascals')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
legend(loc = 'lower right')
savefig('contour_plot_temp_dryatm.png', bbox_inches = 'tight')

fig,ax = subplots(1,1)

h = ax.contourf(time, pres, theta_avg2d)
ax.set_title('Potential Temperature')
#ax.axhline(y = 2.7e4, linestyle = '--', color = 'k', label = '2.7e4 Pascals')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
legend(loc = 'lower right')
savefig('contour_plot_theta_dryatm.png', bbox_inches = 'tight')


fig,ax=subplots(1,1)

h = ax.contourf(time, pres, rhoa_avg)
ax.set_title('Density Variations')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_density_dryatm.png', bbox_inches = 'tight')

fig,ax=subplots(1,1)

plot(tempa_avg_1, p_avg, label = 'time = 0s')
plot(tempa_avg_2, p_avg, label = 'time = .5e8')
plot(tempa_avg_3, p_avg, label = 'time = 1e8')
ax.axhline(y = 2.7e4, linestyle = '--', color = 'k', label = '2.7e4 Pascals')
ax.set_title('Temperature Variation Profiles')
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
legend()
savefig('pres_tempv_dryatm.png', bbox_inches = 'tight')

fig,ax=subplots(1,1)

plot(temp1, p_avg, label = 'time = 0s')
plot(temp2, p_avg, label = 'time = .5e8')
plot(temp3, p_avg, label = 'time = 1e8')
ax.axhline(y = 2.7e4, linestyle = '--', color = 'k', label = '2.7e4 Pascals')
ax.set_title('Temperature Profiles')
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
legend()
savefig('pres_temp_dryatm.png', bbox_inches = 'tight')

fig,ax=subplots(1,1)

plot(theta1, p_avg, label = 'time = 0s')
plot(theta2, p_avg, label = 'time = .5e8')
plot(theta3, p_avg, label = 'time = 1e8')
ax.axhline(y = 2.7e4, linestyle = '--', color = 'k', label = '2.7e4 Pascals')
ax.set_title('Potential Temperature Profiles')
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
legend()
savefig('pres_theta_dryatm.png', bbox_inches = 'tight')



fig,ax=subplots(1,1)
h = ax.contourf(time, pres, kinetic_energy[:,:,0])
ax.set_title('Mean Horizontal Zonal Velocity')
ax.set_xlabel('Time(s)',fontsize = 12)
ax.set_ylabel('Pressure',fontsize = 12)
ax.set_ylim(max(p_avg), min(p_avg))
ax.set_yscale('log')
fig.colorbar(h)
savefig('contour_plot_KE_01.png',bbox_inches = 'tight')

show()

