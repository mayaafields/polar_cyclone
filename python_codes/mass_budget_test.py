#Mass Budget Python Script

#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset
from scipy.ndimage.filters import uniform_filter1d
from scipy import signal

# Retrieving Data from the nc files
#fname = '/home/mafields/CylconeDiagnostics/jupiter_runs/jupiter-2d01-main.nc'
fname = '/home/mafields/CylconeDiagnostics/data/jupiter2d-long-main.nc'
data = Dataset(fname, 'r')

#Defining Variables
time = data['time']
rho = data['rho'][:,:,:,0]
pres = data['press'][:,:,:,0]
mass_flux = data['v1rho'][:,:]
w = data['vel1'][:,:,:,0]
ntime, nx1, nx2 = rho.shape
x1 = data['x1']
x2 = data['x2']
rho_avg = mean(rho, axis = (2))
pres_avg = mean(pres, axis = (0,2))

#Solving For the qNH3 and qH20 for the advective equation
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

mass_flux_qH2O = w*qH2O
mass_flux_qNH3 = w*qNH3

mass_flux_qH2O_avg = mean(mass_flux_qH2O, axis = (2))
mass_flux_qNH3_avg = mean(mass_flux_qNH3, axis = (2))
qH2O_avg = mean(qH2O, axis = (2))
qNH3_avg = mean(qNH3, axis = (2))

print(mass_flux_qH2O_avg.shape)

qH2O_lev1 = qH2O_avg[:,0]
qNH3_lev1 = qNH3_avg[:,0]

fig,ax = subplots(1,1)
plot(time, qH2O_lev1, label = 'level 1 (170 bar) H2O')
legend()
close()

fig,ax = subplots(1,1)
plot(time, qNH3_lev1, label = 'level 1 (170 bar) NH3') 
legend()
close()

#Looking at density spikes in the atmosphere (possible convection events)
rho_lev1 = rho_avg[:,0]
rho_lev2 = rho_avg[:,24]
rho_lev3 = rho_avg[:,47]
fig,ax=subplots(1,1)
plot(time, rho_lev1, label = 'level 1 (170 bar)')
#plot(time, rho_lev2, label ='level 2 (85 bar)')
#plot(time, rho_lev3, label ='level 3 (1 bar)')
title(r'Density versus Time')
ax.set_xlabel(r'Time(s)',fontsize = 12)
ax.set_ylabel(r'$\rho$',fontsize = 12)
legend()
close()

#Time Periods determined by visually identifying convection anomalies
#Period before 1st convection event
t1s = 0
t1e = 201

#Period inbetween convection events
t2s = 280
t2e = 500

#Period after 2nd convection event
t3s = 510
t3e = 800

#1st Anomaly
t4s = t1e
t4e = t2s

#2nd Anomaly
t5s = t2e
t5e = t3s

#Defining time with no anomalies
t1 = time[t1s:t1e]
t2 = time[t2s:t2e]
t3 = time[t3s:t3e]

#Defining time with anomalies
t4 = time[t4s:t4e]
t5 = time[t5s:t5e]

#Linear Fits

#Density trends
beta1 = np.zeros(len(x1))
beta2 = np.zeros(len(x1))
beta3 = np.zeros(len(x1))
beta4 = np.zeros(len(x1))
beta5 = np.zeros(len(x1))

#H2O trends
alpha1 = np.zeros(len(x1))
alpha2 = np.zeros(len(x1))
alpha3 = np.zeros(len(x1))
alpha4 = np.zeros(len(x1))
alpha5 = np.zeros(len(x1))

#NH3 trends
eta1 = np.zeros(len(x1))
eta2 = np.zeros(len(x1))
eta3 = np.zeros(len(x1))
eta4 = np.zeros(len(x1))
eta5 = np.zeros(len(x1))
for i in range(nx1):
    rho = rho_avg[:,i]
    rho1 = rho[t1s:t1e]
    rho2 = rho[t2s:t2e]
    rho3 = rho[t3s:t3e]
    rho4 = rho[t4s:t4e]
    rho5 = rho[t5s:t5e]

    h = qH2O_avg[:,i]
    H2O1 = h[t1s:t1e]
    H2O2 = h[t2s:t2e]
    H2O3 = h[t3s:t3e]
    H2O4 = h[t4s:t4e]
    H2O5 = h[t5s:t5e]

    n = qNH3_avg[:,i]
    NH3_1 = n[t1s:t1e]
    NH3_2 = n[t2s:t2e]
    NH3_3 = n[t3s:t3e]
    NH3_4 = n[t4s:t4e]
    NH3_5 = n[t5s:t5e]

    linfit_rho1 = polyfit(t1,rho1,1)
    linfit_rho2 = polyfit(t2,rho2,1)
    linfit_rho3 = polyfit(t3,rho3,1)
    linfit_rho4 = polyfit(t4,rho4,1)
    linfit_rho5 = polyfit(t5,rho5,1)

    linfit_H2O1 = polyfit(t1,H2O1,1)
    linfit_H2O2 = polyfit(t2,H2O2,1)
    linfit_H2O3 = polyfit(t3,H2O3,1)
    linfit_H2O4 = polyfit(t4,H2O4,1)
    linfit_H2O5 = polyfit(t5,H2O5,1)

    linfit_NH3_1 = polyfit(t1,NH3_1,1)
    linfit_NH3_2 = polyfit(t2,NH3_2,1)
    linfit_NH3_3 = polyfit(t3,NH3_3,1)
    linfit_NH3_4 = polyfit(t4,NH3_4,1)
    linfit_NH3_5 = polyfit(t5,NH3_5,1)
    #beta is an array of d\rho/dt values for each pressure level
    beta1[i] = linfit_rho1[0]
    beta2[i] = linfit_rho2[0]
    beta3[i] = linfit_rho3[0]
    beta4[i] = linfit_rho4[0]
    beta5[i] = linfit_rho5[0]

    alpha1[i] = linfit_H2O1[0]
    alpha2[i] = linfit_H2O2[0]
    alpha3[i] = linfit_H2O3[0]
    alpha4[i] = linfit_H2O4[0]
    alpha5[i] = linfit_H2O5[0]
    
    eta1[i] = linfit_NH3_1[0]
    eta2[i] = linfit_NH3_2[0]
    eta3[i] = linfit_NH3_3[0]
    eta4[i] = linfit_NH3_4[0]
    eta5[i] = linfit_NH3_5[0]

#d\rho/dt plot 
fig,ax=subplots(1,1)
plot(beta1,pres_avg, label = r'$\frac{d(\rho)_1}{dt}$')
plot(beta2,pres_avg, label = r'$\frac{d(\rho)_2}{dt}$')
plot(beta3,pres_avg, label = r'$\frac{d(\rho)_3}{dt}$')
ax.set_xlabel(r'$\frac{d\rho}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('drhodt.png',bbox_inches= 'tight')
close()

#d qH2O/dt plot
fig,ax=subplots(1,1)
plot(alpha1,pres_avg, label = r'$\frac{d(qH_2O)_1}{dt}$')
plot(alpha2,pres_avg, label = r'$\frac{d(qH_2O)_2}{dt}$')
plot(alpha3,pres_avg, label = r'$\frac{d(qH_2O)_3}{dt}$')
ax.set_xlabel(r'$\frac{d(qH_2O)}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('dh2Odt.png',bbox_inches= 'tight')

#d qNH3/dt plot
fig,ax=subplots(1,1)
plot(eta1,pres_avg, label = r'$\frac{d(qNH_3)_1}{dt}$')
plot(eta2,pres_avg, label = r'$\frac{d(qNH_3)_2}{dt}$')
plot(eta3,pres_avg, label = r'$\frac{d(qNH_3)_3}{dt}$')
ax.set_xlabel(r'$\frac{d(qNH_3)}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('dNH3dt.png',bbox_inches= 'tight')

#drho/dt anomaly plot
fig,ax=subplots(1,1)
plot(beta4,pres_avg, label = r'$\frac{d(\rho)_4}{dt}$ 1st Anomaly')
plot(beta5,pres_avg, label = r'$\frac{d(\rho)_5}{dt}$ 2nd Anomaly')
ax.set_xlabel(r'$\frac{d(\rho)}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
#ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('drhodt_anomaly.png',bbox_inches= 'tight')
close()

#dH2O/dt anomaly plot
fig,ax=subplots(1,1)
plot(alpha4,pres_avg, label = r'$\frac{d(qH_2O)_4}{dt}$ 1st Anomaly')
plot(alpha5,pres_avg, label = r'$\frac{d(qH_2O)_5}{dt}$ 2nd Anomaly')
ax.set_xlabel(r'$\frac{d(qH_2O)}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
#ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('dH2Odt_anomaly.png',bbox_inches= 'tight')

#dNH3/dt anomaly plot
fig,ax=subplots(1,1)
plot(eta4,pres_avg, label = r'$\frac{d(qNH_3)_4}{dt}$ 1st Anomaly')
plot(eta5,pres_avg, label = r'$\frac{d(qNH_3)_5}{dt}$ 2nd Anomaly')
ax.set_xlabel(r'$\frac{d(qNH_3)}{dt}$',fontsize = 12)
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
#ax.set_xlim(-0.5e-9,0.5e-9)
legend()
savefig('dNH3dt_anomaly.png',bbox_inches= 'tight')

#Defining the Vertical Mass Flux per time period
mass_flux1 = mass_flux[t1s:t1e,:]
mass_flux1_avg = mean(mass_flux1, axis = 0)

mass_flux2 = mass_flux[t2s:t2e,:]
mass_flux2_avg = mean(mass_flux2, axis = 0)

mass_flux3 = mass_flux[t3s:t3e,:]
mass_flux3_avg = mean(mass_flux3, axis = 0)

mass_flux4 = mass_flux[t4s:t4e,:]
mass_flux4_avg = mean(mass_flux4, axis = 0)

mass_flux5 = mass_flux[t5s:t5e,:]
mass_flux5_avg = mean(mass_flux5, axis = 0)


#Defining the Vertical H2O Mass Flux per time period
mass_flux1_qH2O = mass_flux_qH2O_avg[t1s:t1e,:]
mass_flux1_qH2O_avg = mean(mass_flux1_qH2O, axis = 0)

mass_flux2_qH2O = mass_flux_qH2O_avg[t2s:t2e,:]
mass_flux2_qH2O_avg = mean(mass_flux2_qH2O, axis = 0)

mass_flux3_qH2O = mass_flux_qH2O_avg[t3s:t3e,:]
mass_flux3_qH2O_avg = mean(mass_flux3_qH2O, axis = 0)

mass_flux4_qH2O = mass_flux_qH2O_avg[t4s:t4e,:]
mass_flux4_qH2O_avg = mean(mass_flux4_qH2O, axis = 0)

mass_flux5_qH2O = mass_flux_qH2O_avg[t5s:t5e,:]
mass_flux5_qH2O_avg = mean(mass_flux5_qH2O, axis = 0)

#Defining the Vertical NH3 Mass Flux per time period
mass_flux1_qNH3 = mass_flux_qNH3_avg[t1s:t1e,:]
mass_flux1_qNH3_avg = mean(mass_flux1_qNH3, axis = 0)

mass_flux2_qNH3 = mass_flux_qNH3_avg[t2s:t2e,:]
mass_flux2_qNH3_avg = mean(mass_flux2_qNH3, axis = 0)

mass_flux3_qNH3 = mass_flux_qNH3_avg[t3s:t3e,:]
mass_flux3_qNH3_avg = mean(mass_flux3_qNH3, axis = 0)

mass_flux4_qNH3 = mass_flux_qNH3_avg[t4s:t4e,:]
mass_flux4_qNH3_avg = mean(mass_flux4_qNH3, axis = 0)

mass_flux5_qNH3 = mass_flux_qNH3_avg[t5s:t5e,:]
mass_flux5_qNH3_avg = mean(mass_flux5_qNH3, axis = 0)

##Mass Flux (w\rho)
#plotting the vertical mass flux
fig,ax= subplots(1,1)
plot(mass_flux1_avg, pres_avg, label = r'$(w\rho)_1$')
plot(mass_flux2_avg, pres_avg, label = r'$(w\rho)_2$')
plot(mass_flux3_avg, pres_avg, label = r'$(w\rho)_3$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux.png', bbox_inches = 'tight')
close()

#plotting the vertical mass flux anomalies
fig,ax= subplots(1,1)
plot(mass_flux4_avg, pres_avg, label = r'$(w\rho)_4$ 1st Anomaly')
plot(mass_flux5_avg, pres_avg, label = r'$(w\rho)_5$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_anomaly.png', bbox_inches = 'tight')
close()

##Mass Flux (w*qH2O)
#plotting the vertical mass flux
fig,ax= subplots(1,1)
plot(mass_flux1_qH2O_avg, pres_avg, label = r'$(w\rho)_1$')
plot(mass_flux2_qH2O_avg, pres_avg, label = r'$(w\rho)_2$')
plot(mass_flux3_qH2O_avg, pres_avg, label = r'$(w\rho)_3$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qH_2O$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_H2O.png', bbox_inches = 'tight')
close()

#plotting the vertical mass flux anomalies
fig,ax= subplots(1,1)
plot(mass_flux4_qH2O_avg, pres_avg, label = r'$(w\rho)_4$ 1st Anomaly')
plot(mass_flux5_qH2O_avg, pres_avg, label = r'$(w\rho)_5$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qH_2O$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_qH2O_anomaly.png', bbox_inches = 'tight')
close()

##Mass Flux (w*qH2O)
#plotting the vertical mass flux
fig,ax= subplots(1,1)
plot(mass_flux1_qNH3_avg, pres_avg, label = r'$(w\rho)_1$')
plot(mass_flux2_qNH3_avg, pres_avg, label = r'$(w\rho)_2$')
plot(mass_flux3_qNH3_avg, pres_avg, label = r'$(w\rho)_3$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qNH_3$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_NH3.png', bbox_inches = 'tight')
close()

#plotting the vertical mass flux anomalies
fig,ax= subplots(1,1)
plot(mass_flux4_qNH3_avg, pres_avg, label = r'$(w\rho)_4$ 1st Anomaly')
plot(mass_flux5_qNH3_avg, pres_avg, label = r'$(w\rho)_5$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qNH_3$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_qNH3_anomaly.png', bbox_inches = 'tight')
close()

#Moving Average of Mass Flux
N = 5 #5 point moving average

mass_flux1_avg = uniform_filter1d(mass_flux1_avg, size=N)
dwrhodz1 = (mass_flux1_avg[1:] - mass_flux1_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux2_avg = uniform_filter1d(mass_flux2_avg, size=N)
dwrhodz2 = (mass_flux2_avg[1:] - mass_flux2_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux3_avg = uniform_filter1d(mass_flux3_avg, size=N)
dwrhodz3 = (mass_flux3_avg[1:] - mass_flux3_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux4_avg = uniform_filter1d(mass_flux4_avg, size=N)
dwrhodz4 = (mass_flux4_avg[1:] - mass_flux4_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux5_avg = uniform_filter1d(mass_flux5_avg, size=N)
dwrhodz5 = (mass_flux5_avg[1:] - mass_flux5_avg[:-1])/(x1[1:] - x1[:-1])

#Moving Average of Mass Flux

mass_flux1_qH2O_avg = uniform_filter1d(mass_flux1_qH2O_avg, size=N)
dwqH2Odz1 = (mass_flux1_qH2O_avg[1:] - mass_flux1_qH2O_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux2_qH2O_avg = uniform_filter1d(mass_flux2_qH2O_avg, size=N)
dwqH2Odz2 = (mass_flux2_qH2O_avg[1:] - mass_flux2_qH2O_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux3_qH2O_avg = uniform_filter1d(mass_flux3_qH2O_avg, size=N)
dwqH2Odz3 = (mass_flux3_qH2O_avg[1:] - mass_flux3_qH2O_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux4_qH2O_avg = uniform_filter1d(mass_flux4_qH2O_avg, size=N)
dwqH2Odz4 = (mass_flux4_qH2O_avg[1:] - mass_flux4_qH2O_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux5_qH2O_avg = uniform_filter1d(mass_flux5_qH2O_avg, size=N)
dwqH2Odz5 = (mass_flux5_qH2O_avg[1:] - mass_flux5_qH2O_avg[:-1])/(x1[1:] - x1[:-1])

#Moving Average of Mass Flux

mass_flux1_qNH3_avg = uniform_filter1d(mass_flux1_qNH3_avg, size=N)
dwqNH3dz1 = (mass_flux1_qNH3_avg[1:] - mass_flux1_qNH3_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux2_qNH3_avg = uniform_filter1d(mass_flux2_qNH3_avg, size=N)
dwqNH3dz2 = (mass_flux2_qNH3_avg[1:] - mass_flux2_qNH3_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux3_qNH3_avg = uniform_filter1d(mass_flux3_qNH3_avg, size=N)
dwqNH3dz3 = (mass_flux3_qNH3_avg[1:] - mass_flux3_qNH3_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux4_qNH3_avg = uniform_filter1d(mass_flux4_qNH3_avg, size=N)
dwqNH3dz4 = (mass_flux4_qNH3_avg[1:] - mass_flux4_qNH3_avg[:-1])/(x1[1:] - x1[:-1])

mass_flux5_qNH3_avg = uniform_filter1d(mass_flux5_qNH3_avg, size=N)
dwqNH3dz5 = (mass_flux5_qNH3_avg[1:] - mass_flux5_qNH3_avg[:-1])/(x1[1:] - x1[:-1])

#plotting the vertical mass flux
fig,ax= subplots(1,1)
plot(mass_flux1_avg, pres_avg, label = r'$(w\rho)_1$')
plot(mass_flux2_avg, pres_avg, label = r'$(w\rho)_2$')
plot(mass_flux3_avg, pres_avg, label = r'$(w\rho)_3$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_mvavg.png', bbox_inches = 'tight')
close()

#plotting the vertical mass flux Anomalies
fig,ax= subplots(1,1)
plot(mass_flux4_avg, pres_avg, label = r'$(w\rho)_4$ 1st Anomaly')
plot(mass_flux5_avg, pres_avg, label = r'$(w\rho)_5$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_mvavg_anomaly.png', bbox_inches = 'tight')
close()

#plotting the vertical mass flux H2O
fig,ax= subplots(1,1)
plot(mass_flux1_qH2O_avg, pres_avg, label = r'$(w\cdot qH_2O)_1$')
plot(mass_flux2_qH2O_avg, pres_avg, label = r'$(w\cdot qH_2O)_2$')
plot(mass_flux3_qH2O_avg, pres_avg, label = r'$(w\cdot qH_2O)_3$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qH_2O$($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_qH2O_mvavg.png', bbox_inches = 'tight')


#plotting the vertical H2O mass flux Anomalies
fig,ax= subplots(1,1)
plot(mass_flux4_qH2O_avg, pres_avg, label = r'$(w\cdot qH_2O)_4$ 1st Anomaly')
plot(mass_flux5_qH2O_avg, pres_avg, label = r'$(w\cdot qH_2O)_5$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qH_2O$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
savefig('mass_flux_qH2Omvavg.png', bbox_inches = 'tight')
legend()

#plotting the vertical mass flux NH3
fig,ax= subplots(1,1)
plot(mass_flux1_qNH3_avg, pres_avg, label = r'$(w\cdot qNH_3)_1$')
plot(mass_flux2_qNH3_avg, pres_avg, label = r'$(w\cdot qNH_3)_2$')
plot(mass_flux3_qNH3_avg, pres_avg, label = r'$(w\cdot qNH_3)_3$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qNH_3$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
legend()
savefig('mass_flux_qNH3_mvavg.png', bbox_inches = 'tight')


#plotting the vertical NH3 mass flux Anomalies
fig,ax= subplots(1,1)
plot(mass_flux4_qNH3_avg, pres_avg, label = r'$(w \cdot qNH_3)_4$ 1st Anomaly')
plot(mass_flux5_qNH3_avg, pres_avg, label = r'$(w\cdot qNH_3)_5$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_xlabel(r'Mass Flux $qNH_3$ ($\frac{kg}{s \cdot m^2}$)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
savefig('mass_flux_qNH3_mvavg.png', bbox_inches = 'tight')
legend()

pres_avg2 = sqrt(pres_avg[1:]*pres_avg[:-1])

#lines where clouds occur 
y1 = 3e6 
y2 = 3e5
y3 = 6e4
y4 = 6e3

fig,ax= subplots(1,1)
plot(-1*dwrhodz1,pres_avg2,label = r'$\frac{-d(w \rho)_1}{dz}$')
plot(-1*dwrhodz2,pres_avg2,label = r'$\frac{-d(w \rho)_2}{dz}$')
plot(-1*dwrhodz3,pres_avg2,label = r'$\frac{-d(w \rho)_3}{dz}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg2),min(pres_avg2))

axhline(y1, linestyle = '-', color = 'black')
axhline(y2, linestyle = '-', color = 'black')
ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

axhline(y3, linestyle = '-', color = 'black')
axhline(y4, linestyle = '-', color = 'black')
ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg.png', bbox_inches = 'tight')
close()

fig,ax= subplots(1,1)
plot(-1*dwqH2Odz1,pres_avg2,label = r'$\frac{-d(w\cdot qH_2O)_1}{dz}$')
plot(-1*dwqH2Odz2,pres_avg2,label = r'$\frac{-d(w\cdot qH_2O)_2}{dz}$')
plot(-1*dwqH2Odz3,pres_avg2,label = r'$\frac{-d(w\cdot qH_2O)_3}{dz}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg2),min(pres_avg2))

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwqH2Odz_mvavg.png', bbox_inches = 'tight')


fig,ax= subplots(1,1)
plot(-1*dwqNH3dz1,pres_avg2,label = r'$\frac{-d(w\cdot qNH_3)_1}{dz}$')
plot(-1*dwqNH3dz2,pres_avg2,label = r'$\frac{-d(w\cdot qNH_3)_2}{dz}$')
plot(-1*dwqNH3dz3,pres_avg2,label = r'$\frac{-d(w\cdot qNH_3)_3}{dz}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg2),min(pres_avg2))

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwqNH3dz_mvavg.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
plot(-1*dwrhodz4,pres_avg2,label = r'$\frac{-d(w \rho)_4}{dz}$ 1st Anomaly')
plot(-1*dwrhodz5,pres_avg2,label = r'$\frac{-d(w \rho)_5}{dz}$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg2),min(pres_avg2))

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg_anomaly.png', bbox_inches = 'tight')
close()

fig,ax= subplots(1,1)
plot(-1*dwqH2Odz4,pres_avg2,label = r'$\frac{-d(w \cdot qH_2O)_4}{dz}$ 1st Anomaly')
plot(-1*dwqH2Odz5,pres_avg2,label = r'$\frac{-d(w \cdot qH_2O)_5}{dz}$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg2),min(pres_avg2))

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg_qH2O_anomaly.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
plot(-1*dwqNH3dz4,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_4}{dz}$ 1st Anomaly')
plot(-1*dwqNH3dz5,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_5}{dz}$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg2),min(pres_avg2))

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg_qNH3_anomaly.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
plot(-1*dwrhodz1,pres_avg2,label = r'$\frac{-d(w \rho)_1}{dz}$')
plot(-1*dwrhodz2,pres_avg2,label = r'$\frac{-d(w \rho)_2}{dz}$')
plot(-1*dwrhodz3,pres_avg2,label = r'$\frac{-d(w \rho)_3}{dz}$')

axhline(y1, linestyle = '-', color = 'black')
axhline(y2, linestyle = '-', color = 'black')
ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

axhline(y3, linestyle = '-', color = 'black')
axhline(y4, linestyle = '-', color = 'black')
ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

plot(beta1,pres_avg, label = r'$\frac{d(\rho)_1}{dt}$')
plot(beta2,pres_avg, label = r'$\frac{d(\rho)_2}{dt}$')
plot(beta3,pres_avg, label = r'$\frac{d(\rho)_3}{dt}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg_compare.png', bbox_inches = 'tight')
close()

fig,ax= subplots(1,1)
plot(-1*dwqH2Odz1,pres_avg2,label = r'$\frac{-d(w \cdot H_2O)_1}{dz}$')
plot(-1*dwqH2Odz2,pres_avg2,label = r'$\frac{-d(w \cdot H_2O)_2}{dz}$')
plot(-1*dwqH2Odz3,pres_avg2,label = r'$\frac{-d(w \cdot H_2O)_3}{dz}$')

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

plot(alpha1,pres_avg, label = r'$\frac{d(qH_2O)_1}{dt}$')
plot(alpha2,pres_avg, label = r'$\frac{d(qH_2O)_2}{dt}$')
plot(alpha3,pres_avg, label = r'$\frac{d(qH_2O)_3}{dt}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwqH2Odt_mvavg_compare.png', bbox_inches = 'tight')
close()

fig,ax= subplots(1,1)
plot(-1*dwqNH3dz1,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_1}{dz}$')
plot(-1*dwqNH3dz2,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_2}{dz}$')
plot(-1*dwqNH3dz3,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_3}{dz}$')

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

plot(eta1,pres_avg, label = r'$\frac{d(qNH_3)_1}{dt}$')
plot(eta2,pres_avg, label = r'$\frac{d(qNH_3)_2}{dt}$')
plot(eta3,pres_avg, label = r'$\frac{d(qNH_3)_3}{dt}$')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwqNH3dt_mvavg_compare.png', bbox_inches = 'tight')
close()

fig,ax= subplots(1,1)
plot(-1*dwrhodz4,pres_avg2,label = r'$\frac{-d(w \rho)_4}{dz}$ 1st Anomaly')
plot(-1*dwrhodz5,pres_avg2,label = r'$\frac{-d(w \rho)_5}{dz}$ 2nd Anomaly')

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

plot(beta4,pres_avg, label = r'$\frac{d(\rho)_4}{dt}$ 1st Anomaly')
plot(beta5,pres_avg, label = r'$\frac{d(\rho)_5}{dt}$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwrhodt_mvavg_compare_anomaly.png', bbox_inches = 'tight')
close()

fig,ax= subplots(1,1)
plot(-1*dwqH2Odz4,pres_avg2,label = r'$\frac{-d(w \dot H_2O)_4}{dz}$ 1st Anomaly')
plot(-1*dwqH2Odz5,pres_avg2,label = r'$\frac{-d(w \dot H_2O)_5}{dz}$ 2nd Anomaly')

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

plot(alpha4,pres_avg, label = r'$\frac{d(qH_2O)_4}{dt}$ 1st Anomaly')
plot(alpha5,pres_avg, label = r'$\frac{d(qH_2O)_5}{dt}$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwqH2Odt_mvavg_compare_anomaly.png', bbox_inches = 'tight')

fig,ax= subplots(1,1)
plot(-1*dwqNH3dz4,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_4}{dz}$ 1st Anomaly')
plot(-1*dwqNH3dz5,pres_avg2,label = r'$\frac{-d(w \cdot qNH_3)_5}{dz}$ 2nd Anomaly')

#axhline(y1, linestyle = '-', color = 'black')
#axhline(y2, linestyle = '-', color = 'black')
#ax.axhspan(y1, y2, facecolor='yellow', alpha=0.4)

#axhline(y3, linestyle = '-', color = 'black')
#axhline(y4, linestyle = '-', color = 'black')
#ax.axhspan(y3, y4, facecolor='yellow', alpha=0.4)

plot(eta4,pres_avg, label = r'$\frac{d(qNH_3)_4}{dt}$ 1st Anomaly')
plot(eta5,pres_avg, label = r'$\frac{d(qNH_3)_5}{dt}$ 2nd Anomaly')
ax.set_ylabel('Pressure (Pascals)',fontsize = 12)
ax.set_yscale('log')
ax.set_ylim(max(pres_avg),min(pres_avg))
axvline(x = 0, linestyle = '--', color = 'black')
legend()
savefig('dwqNH3dt_mvavg_compare_anomaly.png', bbox_inches = 'tight')


show();
