# This program plots the star map and black hole particles on top of the gas for each galaxy. It also cross-checks the particle ID's in the simulation file with the particle ID's in the bhlist10.txt file. This cross-check prints out what black hole particles are within 10% of the virial radius for each corresponding snapshot.

# For a more thorough introduction to pygad, check out the documentation at https://bitbucket.org/broett/pygad/downloads/, especially the quickstart jupyter notebook

# First import pygad, numpy (which has nice functions for manipulating arrays of data), and matplotlib, which allows us to make figures

import pygad
import pygad.plotting
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import read_snap #ignore for now
import difflib

# Set which galaxy you want to look at (these have names like m0053, m0163, etc. For gal, you can just use the number, 53, 163, etc.) at which snapshot

gal = 858
snapmin = 19
snap = range(45,47) # This creates a list of snaps from the snapmin to the final snapshot of 94

# This path directs to a folder on mhaas which contains information about the galaxy

path = '/Volumes/Happy/Choi16_Fiducial/BH_trace/m0'+str(gal).zfill(3)+'/'

# Read in the files with these values as arrays 

centerarray = np.loadtxt(path+'Center_trace.txt',usecols=(0,1,2))
bharray = np.loadtxt(path+'BHSFactivity/bhlist10.txt', skiprows=1)
swallowarray = np.loadtxt(path+'blackhole_details/Swallow_shortlist.txt', skiprows=1)

# This loops over all of the snapshots

for i in snap:

# This reads in the snapshot, converts it into physical units, and assigns a center to the galaxy
    s = pygad.Snap('/Volumes/Happy/Choi16_Fiducial/BH_trace/m0'+str(gal).zfill(3)+'/snap_m0'+str(gal).zfill(3)+'_sf_x_2x_'+str(i).zfill(3), load_double_prec=True)
    s.to_physical_units()
    center = pygad.analysis.shrinking_sphere(s.stars, center = [s.boxsize/2]*3, R=s.boxsize)
    pygad.Translation(-center).apply(s)
                                                                   
# This plots the gas map with an extent and color bar that can be modified. I also assign a variable bh_pos to equal the positions of the black holes

    fig, ax, im, cbar =  pygad.plotting.image(s.stars, extent = '200 kpc',clim=[0,12], Npx=150)
    bh_pos = s.bh['pos']
    bh_mass = s.bh['mass']
    max_bh_mass = max(bh_mass)
    x = np.where(max_bh_mass)
    print 'The most massive black hole has a mass of', max_bh_mass
    ax.text(0.80, 0.05, 'z = %7.2f' % s.redshift, color='white', fontsize=10, transform=ax.transAxes) 

# This is my cross-check which uses a function from read_snap.py. The first entry is what is being printed out while the second and third entries are what you are comparing. I then print the particle ID's from the simulations that match bhlist10.txt for each snapshot

    match = read_snap.check_in(s.bh['ID'], s.bh['ID'], bharray[i][2:6])    
    print 'These are the matches for snapshot ',i, ':', match

# Need to transpose matrix so all data for the same thing is in the same row for python to read in
# Defining rows
    trans_swallowarray = np.transpose(swallowarray)

    scale_factor = trans_swallowarray[0]    
    bh_ID = trans_swallowarray[1]
    gas_ID = trans_swallowarray[2]
    BH_mass = trans_swallowarray[3]

    trans_bharray = np.transpose(bharray)

    scalef = trans_bharray[0]
    BHIDS = trans_bharray[2]

    track_gas_ID = []

    for ID in range(0,len(BHIDS)):
        for id in range(0,len(bh_ID)):
            print np.where(BHIDS[ID]==bh_ID[id])
                
            

# This loops over the bh_pos array and plots a circle where each black hole exists on top of the gas map. Need to add artist so the circles show up on top. This is a multi-dimensional array so [b] tells us what row to look in, [1] tells us the x values, and [0] tells us the y values for each black hole position (simulation files created this way)

    for b in range(0, len(bh_pos)):
        for m in range(0,len(bh_mass)):
            if bh_mass[b]==max_bh_mass:
                circle = plt.Circle((bh_pos[b][1],bh_pos[b][0]), 5, color = 'w', fill = False)
                ax.add_artist(circle)
        
    plt.savefig('./m0'+str(gal).zfill(3)+'/BH_trace/smd/BH_smd_trace_m0'+str(gal).zfill(3)+'_snap_'+str(i).zfill(3)+'.png')
    plt.close()                                                             



# ax.set_title(r'm0'+str(gal).zfill(3)+'              $M_{vir}=$'+str(virial_mass)+' x $10^{12}  M_{\odot}$')

#If you want the radii of all the particles that are NOT in special_ids_z0, then use

#read_snap.check_notin(s.gas['r'], s.gas['ID'], special_ids_z0)

