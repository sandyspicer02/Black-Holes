# This program reads in master.txt which includes all gas particle ID's, scale factors, and radii for all snapshots and plots the particle trajectories for each

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import read_snap
import cosmolopy.distance as cd
import cosmolopy.constants as cc

cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.674}
cosmo = cd.set_omega_k_0(cosmo)
agefunc = cd.quick_age_function(**cosmo)

# Specify galaxy
gal = 858

# Specify path
path = '/Volumes/Happy/Choi16_Fiducial/BH_trace/m0'+str(gal).zfill(3)+'/'

# Read in files
gasparticleID = np.loadtxt('/Users/sandyspicer02/Research/m0'+str(gal).zfill(3)+'/r_vs_t/master.txt',dtype = 'uint',usecols=(0,))
scalef = np.loadtxt('/Users/sandyspicer02/Research/m0'+str(gal).zfill(3)+'/r_vs_t/master.txt',usecols=(1,))
radius = np.loadtxt('/Users/sandyspicer02/Research/m0'+str(gal).zfill(3)+'/r_vs_t/master.txt',usecols=(2,))
cosmictime = np.loadtxt('/Users/sandyspicer02/Research/m0'+str(gal).zfill(3)+'/r_vs_t/master.txt',usecols=(3,))
swallowarray = np.loadtxt(path+'blackhole_details/Swallow_shortlist.txt',dtype='uint', skiprows=1,usecols=(1,2))
accretion = np.loadtxt(path+'blackhole_details/Swallow_shortlist.txt',skiprows=1,usecols=(0,))
rvir = np.loadtxt(path+'Rvir.txt',skiprows=5)
#time = np.loadtxt('/Users/brennan/Documents/winds_analysis/snap_z_time.txt', skiprows=5,usecols=(3,))


# Transpose matrices
trans_rvir = np.transpose(rvir)
redshift = trans_rvir[0]
virial_radius = trans_rvir[1]

trans_swallowarray = np.transpose(swallowarray)
swallow_ID = trans_swallowarray[1]
trans_accretion = np.transpose(accretion)

# Swallow time
age_ = agefunc(redshift)
time = age_/cc.Gyr_s

for g in range(0,len(swallow_ID)):
    match_radius = read_snap.check_in(radius,gasparticleID, swallow_ID[0+g])
    match_scalef = read_snap.check_in(scalef,gasparticleID, swallow_ID[0+g])
    match_ID = read_snap.check_in(gasparticleID, gasparticleID, swallow_ID[0+g])
    match_time = read_snap.check_in(cosmictime, gasparticleID, swallow_ID[0+g])
   # print match_time
    
    # Simulation time
    age = agefunc((1/match_scalef)-1)
    t = age/cc.Gyr_s

    # Accretion time
    A_time = agefunc((1/trans_accretion[g])-1)
    a_time = A_time/cc.Gyr_s
    #print a_time

    test = agefunc(0.69)
    print test/cc.Gyr_s
    
    for f in range(0,len(match_radius)):
        if len(match_radius)==1:
            plt.plot(match_time,match_radius,'ro',label='Trajectory')
            plt.plot(time,virial_radius,'b--',label='Rvir')
            plt.plot(time,0.1*virial_radius,'k--',label='0.1*Rvir')
            plt.axvline(x=a_time,color = 'g', linestyle = 'dashed',label='Accretion')
            plt.xlim(1,8)
            plt.xlabel('Time (Gyr)')
            plt.ylabel('Radius (kpc)')
            plt.title('Gas Particle '+str(swallow_ID[g]))
            plt.legend(['Trajectory','Rvir','0.1*Rvir','Accretion'],loc='upper right',prop={'size':10},numpoints=1)

            plt.savefig('./m0'+str(gal).zfill(3)+'/r_vs_t/'+str(swallow_ID[g])+'_trajectory.png')
            plt.close()
                
        else:
            plt.plot(match_time,match_radius,'r-',label='Trajectory')
            plt.plot(time,virial_radius,'b--',label='Rvir')
            plt.plot(time,0.1*virial_radius,'k--',label='0.1*Rvir')
            plt.axvline(x=a_time,color = 'g', linestyle = 'dashed',label='Accretion')
            plt.xlim(1,8)
            plt.xlabel('Time (Gyr)')
            plt.ylabel('Radius (kpc)')
            plt.title('Gas Particle '+str(swallow_ID[g]))
            plt.legend(['Trajectory','Rvir','0.1*Rvir','Accretion'],loc='upper right',prop={'size':10})
            
            plt.savefig('./m0'+str(gal).zfill(3)+'/r_vs_t/'+str(swallow_ID[g])+'_trajectory.png')
            plt.close()
