import numpy as np

path = '/Volumes/Happy/Choi16_Fiducial/BH_trace/m0858/'

bharray = np.loadtxt(path+'BHSFactivity/bhlist10.txt', skiprows=1)
swallowarray = np.loadtxt(path+'blackhole_details/Swallow_shortlist.txt', skiprows=1)

# Need to transpose matrix so all data for the same thing is in the same row for python to read in                                                                                                                                                                                                                                                         
trans_swallowarray = np.transpose(swallowarray)

scale_factor = trans_swallowarray[0]
bh_ID = trans_swallowarray[1]
gas_ID = trans_swallowarray[2]
BH_mass = trans_swallowarray[3]

trans_bharray = np.transpose(bharray)

scalef = trans_bharray[0]
BHIDS = trans_bharray[2]

track_gas_ID = []

#print np.where(BHIDS==bh_ID)
for ID in range(0,len(BHIDS)):
    for id in range(0,len(bh_ID)):
        if BHIDS[ID] == bh_ID[id]:
            track_gas_ID.append(gas_ID[id])
            print track_gas_ID
