import meep as mp
from meep import mpb
import numpy as np
import datetime
import math
import sys
import matplotlib.pyplot as plt
from IPython.display import Video
import warnings


from cavitysimulations.geometry.waveguide import *
from cavitysimulations.geometry.cavity import *
from cavitysimulations.visualization import *
from cavitysimulations.utilities.utilities import *


def main():
    geom = a_tapered_cavity()
    
    boundary_layers = get_boundary_layer(sim2d=False)
    
    fcen = 1/1.54
    df = 0.1
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), 
                         component=mp.Hz, 
                         center=mp.Vector3())]
    
    symmetries = [mp.Mirror(mp.X,+1), mp.Mirror(mp.Y,-1), mp.Mirror(mp.Z,+1)]
    
    sim = mp.Simulation(resolution=30, 
                        cell_size=mp.Vector3(20, 8, 8), 
                        geometry=geom, 
                        boundary_layers=boundary_layers, 
                        sources=sources,
                        symmetries=symmetries,
                        progress_interval=100,)
    
    h = mp.Harminv(mp.Hz, mp.Vector3(0, 0, 0), fcen, df)
    time_after_source = 500
    # Don't output eps anymore to save disk space
    sim.run(mp.after_sources(h),
            until_after_sources=time_after_source)
    
    visualize_sim_cell(sim)
    
    print("Modal Volume: {}".format(sim.modal_volume_in_box()))
    
    
if __name__ == '__main__':
    main()
