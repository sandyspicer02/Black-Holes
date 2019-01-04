# This program calculates the trajectory of all swallowed gas particles for a given galaxy and outputs the gas ID, scale factor, and radii as a text file.

import pygad
import pygad.plotting
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import read_snap 

# Set which galaxy you want to look at (these have names like m0053, m0163, etc. For gal, you can just use the number, 53, 163, etc.) at which snapshot                              
gal = 858
snapmin = 19
snap = range(snapmin,95) # This creates a list of snaps from the snapmin to the final snapshot of 94                                                                                      
# This path directs to a folder on mhaas which contains information about the galaxy                                                                                                 
path = '/Volumes/Happy/Choi16_Fiducial/BH_trace/m0'+str(gal).zfill(3)+'/'

# Read in the files with these values as arrays                                                                                                                                      
centerarray = np.loadtxt(path+'Center_trace.txt',usecols=(0,1,2))
bharray = np.loadtxt(path+'BHSFactivity/bhlist10.txt', skiprows=1,dtype = 'uint',usecols=(2,))
swallowarray = np.loadtxt(path+'blackhole_details/Swallow_shortlist.txt',dtype='uint', skiprows=1,usecols=(1,2))

# This loops over all of the snapshots                                                                                                                                              

for i in snap:

# This reads in the snapshot, converts it into physical units, and assigns a center to the galaxy                                                                                    
    s = pygad.Snap('/Volumes/Happy/Choi16_Fiducial/BH_trace/m0'+str(gal).zfill(3)+'/snap_m0'+str(gal).zfill(3)+'_sf_x_2x_'+str(i).zfill(3), load_double_prec=True)
    s.to_physical_units()
    center = pygad.analysis.shrinking_sphere(s.stars, center = [s.boxsize/2]*3, R=s.boxsize)
    pygad.Translation(-center).apply(s)

    dx = []
    dy = []
    dz = []

    gas_id = s.gas['ID']
    gasID = np.array(gas_id)                                                                                                                                                         
    gasid = [id.astype(np.uint64) for id in gasID]
    gasid = np.array(gasid)

    bh_id = s.bh['ID']
    bhID = np.array(bh_id)

    bhid = [id.astype(np.uint64) for id in bhID]
    bhid = np.array(bhid)

    trans_swallowarray = np.transpose(swallowarray)
    swallow_ID = trans_swallowarray[1]

    trans_bharray = np.transpose(bharray)
    BHID = trans_bharray
    #print BHID

    # Need scale factor
    z = s.redshift
    a = 1/(z+1)

# If you want to check if a specific ID is in the simulation file

#    print gasid[gasid == 1729382282680115925]

    gasparticleID = read_snap.check_in(gasid, gasid, swallow_ID)
    gasparticle = read_snap.check_in(s.gas['pos'], gasid, swallow_ID)
    blackhole = read_snap.check_in(s.bh['pos'], bhid, BHID)
    
# Need to transpose the matrix to get x, y and z values
            
    gas_particle = np.array(gasparticle)
    black_hole = np.array(blackhole)
    trans_gas_particle = np.transpose(gas_particle)
    trans_bh_particle = np.transpose(black_hole)

    for g in range(0,len(swallow_ID)):
        for h in range(0,len(gasparticleID)):
            if gasparticleID[h]==swallow_ID[0+g]:
                
              #  if g==0:
              #      label = 'Gas ID', 'Scale Factor', 'Radius', 'Cosmic Time (Gyr)'
               
                x1 = trans_gas_particle[0][h]
                y1 = trans_gas_particle[1][h]
                z1 = trans_gas_particle[2][h]
    
                x2 = trans_bh_particle[0]
                y2 = trans_bh_particle[1]
                z2 = trans_bh_particle[2]
                
                dx = x1-x2
                dy = y1-y2
                dz = z1-z2


                r = np.sqrt(dx**2 + dy**2 + dz**2)
               # print r[0]

# Used this to create a master.txt file with all gas particle ID's, scale factors, and radii ('a' appends the values each time so this only needs to be run once and then commented out
                file = open('./m0'+str(gal).zfill(3)+'/r_vs_t/master2.txt','a')
                print >> file,gasparticleID[h],a,r[0],s.cosmic_time()
          
