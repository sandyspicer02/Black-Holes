import pygad
import numpy as np

def read_snap(gal, snap, snapmin, bh, centerarray, rvirarray, centertype):

    if centertype == 'rockstar':
        padding = 95 - len(centerarray)
        center = centerarray[snap - padding]
        rvir = rvirarray[snap - padding]

    if centertype == 'gtrace':
        center = centerarray[snap - snapmin + 1]
        rvir=rvirarray[snap - snapmin + 1]

    s = pygad.Snap('/Volumes/Happy/Choi16_Fiducial/Fiducial_'+str(bh)+'/m0'+str(gal).zfill(3)+'/snap_m0'+str(gal).zfill(3)+'_sf_x_2x_'+str(snap).zfill(3), load_double_prec=True)

    s.to_physical_units()

    pygad.Translation(-center).apply(s)

    vel_center = pygad.analysis.mass_weighted_mean(s.stars[pygad.BallMask('1 kpc', center=pygad.UnitArr([0.0, 0.0, 0.0], 'kpc'))], 'vel')

    s['vel'] -= vel_center

    return s

def check_in(array, checkarray, checklist):
    return array[np.nonzero(np.in1d(checkarray, checklist))]

def check_notin(array, checkarray, checklist):
    return array[np.nonzero(np.in1d(checkarray, checklist)==0)]


        
        
    
    
