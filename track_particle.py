# This program plots the gas map and black hole particles on top of the gas for each galaxy. It also cross-checks the particle ID's in the simulation file with the particle ID's in the bhlist10.txt file. This cross-check prints out what black hole particles are within 10% of the virial radius. Gas particles eaten by different black holes are denoted in different colors

# For a more thorough introduction to pygad, check out the documentation at https://bitbucket.org/broett/pygad/downloads/, especially the quickstart jupyter notebook

# First import pygad, numpy (which has nice functions for manipulating arrays of data), and matplotlib, which allows us to make figures

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
bharray = np.loadtxt(path+'BHSFactivity/bhlist10.txt', skiprows=1)
swallowarray = np.loadtxt(path+'blackhole_details/Swallow_shortlist.txt',dtype = 'int', skiprows=1, usecols=(1,2))

# This loops over all of the snapshots

for i in snap:

# This reads in the snapshot, converts it into physical units, and assigns a center to the galaxy
    s = pygad.Snap('/Volumes/Happy/Choi16_Fiducial/BH_trace/m0'+str(gal).zfill(3)+'/snap_m0'+str(gal).zfill(3)+'_sf_x_2x_'+str(i).zfill(3), load_double_prec=True)
    s.to_physical_units()
    center = pygad.analysis.shrinking_sphere(s.stars, center = [s.boxsize/2]*3, R=s.boxsize)
    pygad.Translation(-center).apply(s)
                                                                   
# This plots the gas map with an extent and color bar that can be modified. I also assign a variable bh_pos to equal the positions of the black holes. I transpose the matrix so my x, y and z values can be accessed from the same row. I also change the data type of the gas and BH ID's so that the full value can be read in
                           
    fig, ax, im, cbar =  pygad.plotting.image(s.gas, extent = '200 kpc',clim=[10**4.0,10**7.0])
    gas_id = s.gas['ID']
    gasID = np.array(gas_id)
    gasid = [id.astype(np.int64) for id in gasID]
    gasid = np.array(gasid)

    gas_pos = s.gas['pos']
    bh_pos = s.bh['pos']
    bh_ID = s.bh['ID']

    bhid = [id.astype(np.int64) for id in bh_ID]
    bhid = np.array(bhid)

    bh_mass = s.bh['mass']
    max_bh_mass = max(bh_mass)
    print 'The most massive black hole has a mass of', max_bh_mass
    ax.text(0.75, 0.06, 'z = %7.2f' % s.redshift, color='white', fontsize=10, transform=ax.transAxes)
    ax.text(0.75, 0.03, 't = %7.2f Gyr' % s.cosmic_time(), color='white', fontsize=10, transform=ax.transAxes)
   # print s.cosmic_time()
# This is my cross-check which uses a function from read_snap.py. The first entry is what is being printed out while the second and third entries are what you are comparing. I then print the particle ID's from the simulations that match bhlist10.txt for each snapshot
    
    trans_swallowarray = np.transpose(swallowarray)
    swallow_gasID = trans_swallowarray[1]
    match = read_snap.check_in(s.bh['ID'], s.bh['ID'], bharray[i][2:6])    
    print 'These are the matches for snapshot ',i, ':', match

# This is one way of matching ID's without using check_in function

    for bh in range(0,len(bh_pos)):
        if bhid[bh] == swallowarray[0][0]:
            for j in range(0,len(swallow_gasID)):
                match1 = swallow_gasID[i]

                for k in range(0,len(gasid)):
                    match2 = gasid[j]

                    if match1 == match2:
                        print swallow_gasID[j]
                        print gasid[k]
                        print('Swallow ID: %i is Simulation ID: %i'%(j,k))

# Matching ID's using check_in function from read_snap.py

    for a in range(0,len(swallowarray)):
        if swallowarray[a][0]==swallowarray[0][0]:
            particle1 = swallowarray[a][1]
            match_gasID1 = read_snap.check_in(gasid, gasid, particle1)
            if len(match_gasID1)>0:
                match_gaspos1 = read_snap.check_in(gas_pos, gasid, particle1)
                gaspos1 = np.array(match_gaspos1)
                
                for p in range(0, len(match_gaspos1)):
                    circle1 = plt.Circle((gaspos1[p][0],gaspos1[p][1]), 5, color = 'k', fill = False)                                                                           
                    ax.add_artist(circle1)
        
        else:
            particle2 = swallowarray[a][1]
            match_gasID2 = read_snap.check_in(gasid, gasid, particle2)
            if len(match_gasID2)>0:
                print match_gasID2
                match_gaspos2 = read_snap.check_in(gas_pos, gasid, particle2)
                gaspos2 = np.array(match_gaspos2)

                for q in range(0, len(match_gaspos2)):
                    circle2 = plt.Circle((gaspos2[q][0],gaspos2[q][1]), 5, color = 'w', fill = False)
                    ax.add_artist(circle2)
                

# This loops over the bh_pos array and plots a circle where each black hole exists on top of the gas map. Need to add artist so the circles show up on top. This is a multi-dimensional array so [b] tells us what row to look in, [0] tells us the x values, and [1] tells us the y values for each black hole position.

    for b in range(0, len(bh_pos)):
        if bhid[b]==swallowarray[0][0]:
            circle3 = plt.Circle((bh_pos[b][0],bh_pos[b][1]), 8, color = 'b', fill = False)
            ax.add_artist(circle3)

        if bhid[b]==trans_swallowarray[0][4]:
            circle4 = plt.Circle((bh_pos[b][0],bh_pos[b][1]), 8, color = 'g', fill = False)
            ax.add_artist(circle4)
           
        if bh_mass[b]==max_bh_mass:
            circle5 = plt.Circle((bh_pos[b][0],bh_pos[b][1]), 8, color = 'r', fill = False)
            ax.add_artist(circle5)


    plt.savefig('./m0'+str(gal).zfill(3)+'/gas_trace/gas_trace_m0'+str(gal).zfill(3)+'_snap_'+str(i).zfill(3)+'.png')
    plt.close()                                                             










# ax.set_title(r'm0'+str(gal).zfill(3)+'              $M_{vir}=$'+str(virial_mass)+' x $10^{12}  M_{\odot}$')

#If you want the radii of all the particles that are NOT in special_ids_z0, then use

#read_snap.check_notin(s.gas['r'], s.gas['ID'], special_ids_z0)

