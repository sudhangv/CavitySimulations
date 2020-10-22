
import meep as mp
import argparse
import math
import numpy as np 
import matplotlib.pyplot as plt
import h5py
from meep import mpb
from math import sqrt, pi
import os

'''
CONTENTS:
--------
(1)convert_freq_to_Thz
(2)get_freq_Thz
(3)get_gamma_from_Thz
(4)get_freqs
(5)add_substrate
(6)visualise_geometry
(7)do_the_sweep
(8)get_index
(9)get_excitation_mode_from_string
(10)get_boundary_layer
(11)index_to_material
(12)get_value_from_index(index, param = 'a')
(13)get_gamma(freq, a, f_target = 1 / (1.54) * 1.01)
'''
def convert_freq_to_Thz(freq, a = 0):
    freq = np.array(freq)
    
#----------- Converts a MEEP frequency ( units of 2pi*c/a) to frequency (not angular velocity) in units of Thz--------#
    if a != 0:
        return ( freq * 3 / a * (10**2))
    else:
        return ( freq * 3 * 10**2)

def get_freq_Thz(hx , hy , a , w , h = 0.22, substrate = False, output_epsilon = False, mode = "zEyO"):
    
    freq = get_freqs(hx = hx , hy = hy , a = a , w = w , h = h , 
                     substrate = substrate, output_epsilon = output_epsilon, mode = mode) 
    
    return convert_freq_to_Thz(freq, a)

def get_gamma_from_Thz(band_edge_f, check_freq):
    if (check_freq < band_edge_f[1])  and (check_freq > band_edge_f[0]):
        
        f_mid = (band_edge_f[0] + band_edge_f[1])/2            
        diff = band_edge_f[0] - band_edge_f[1]
        delta = 1 - (check_freq/ f_mid)

        gamma = math.sqrt(abs(( 0.5 * diff/ f_mid ) ** 2 - delta**2 ))
    
    else:
        gamma = 0
        
    return gamma               


def get_freqs(hx , hy , a , w , h = 0.22, substrate = False, output_epsilon = False ,mode = "zEyO", num_bands = 2):
    
    # h = 0.23    # for manually setting waveguide height 
    res = 20
    #mode = "zEyO"
    resolution = res  # pixels/a, taken from simpetus example
    
    print(" h = " + str(h) + ", SUBSTRATE = " + str(substrate) + ", mode = " + str(mode))
    
    a = round(a,3)        # units of um
    h = round(h, 3)         # units of um
    w = round(w, 3)         # units of um
    hx = round(hx, 3)
    hy = round(hy, 3)
    
    h = h/a          # units of "a"       
    w = w/a          # units of "a"
    hx = hx/a        # units of "a"
    hy = hy/a        # units of "a"
    
    cell_x = 1
    cell_y = 4
    cell_z = 4


    nSi = 3.45
    Si = mp.Medium(index=nSi)

    geometry_lattice = mp.Lattice(size=mp.Vector3(cell_x , cell_y , cell_z)) # dimensions of lattice taken from simpetus example

    geometry = [ mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf, w ,h ), material=Si),
             mp.Ellipsoid(material=mp.air,
             center=mp.Vector3(),
             size=mp.Vector3(hx,hy,mp.inf)) ]
    
    if substrate:
        geometry = add_substrate(geom = geometry, waveguide_height = h, substrate_height = cell_z/2 - h/2 )  # substrate height normalized with a 


    k_points = [mp.Vector3(0.5, 0, 0)]
    
    num_bands = num_bands 

    ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
                        geometry=geometry,
                        k_points=k_points,
                        resolution=resolution,
                        num_bands=num_bands)

    if mode == "te":
        
        ms.run_te() # running for all modes and extracting parities
        
    if mode == "zEyO":
        
        ms.run_yodd_zeven()
    
    if mode == "yO":
        
        ms.run_yodd()
        
    if mode == "yE":
        
        ms.run_yeven()
        
#     if output_epsilon:
#         visualise_geometry(ms = ms, x = None , y = None, z = (4 * res)/ 2, value = None)
    if output_epsilon:
        return ms.get_epsilon()
        
    
        
#         with h5py.File('epsilon.hdf5', 'w') as f:
#             arr = ms.get_epsilon
#             dset = f.create_dataset("epsilon", data = arr)
    
    return ms.freqs

def add_substrate(geom=None, waveguide_height = 0.22, substrate_height = 5, material=mp.Medium(index=1.44)):
    """
    Creates a (by default SiO2) substrate.
    If the unlike case occurs that we need to change the center, we could again make everything relative to the center
    of the substrate.
    """
    if geom is None:
        geom = []

    _center = mp.Vector3(0, 0, -waveguide_height/2 - substrate_height/2)

    geom.append(mp.Block(material= material, 
                         center=_center, 
                         size= mp.Vector3(mp.inf, 
                                          mp.inf, 
                                          substrate_height)))
                
    return geom

def visualise_geometry( ms, x = None, y = None, z = None, value = None, filename = "epsilon.h5"):
    ms.output_epsilon()
    if x is not None:
        plane = 'x'
    
    elif y is not None:
        plane = 'y'
        
    else:
        plane = 'z'
    
    os.system('h5topng -z -Zc dkbluered -a yarg {} -d epsilon.xx'.format(filename))
    
def do_the_sweep():
    
    del_a = 0.001
    del_hy = 0.025
    del_hx = 0.025 
    del_w = 0.05

    a_min = 0.25
    a_max = 0.45        # upper limit of the sweep of a 

    w_min = 0.65         #  lower limit of w 
    w_max = 0.75        #  upper limit of w 

    hx_min = 0.05        # lower limit of the sweep of a
    hy_min = 0.1        #  lower limit of hy 

    hx_max = a_max - 0.07
    hy_max = w_max - 0.1  

    wvg_height = 0.22
    # hx_max = 0.4
    # hy_max = 0.4

    f_target = 1/1.54
    f_target_Thz = convert_freq_to_Thz(f_target) * 1.01

    parameters = []

    with h5py.File('sweep_data230.hdf5', 'w') as f:


        gamma_max = 0        # arbitrary small value
        mirror_strength = []   

        j = len(np.arange(w_min, w_max , del_w))
        k = len(np.arange(a_min , a_max, del_a))
        l = len(np.arange(hy_min, hy_max , del_hy))
        m = len(np.arange(hx_min, hx_max, del_hx ))

        dset = f.create_dataset("data", (j,k,l,m))
        dset[:,:,:,:] =  np.zeros((j,k,l,m))

        index = []
        index_count = 0
          # stores parameters for optimal value of gamma


        #---------------------------#
        #       WIDTH LOOP          
        #---------------------------#
        for w in np.arange(w_min, w_max , del_w):

            #w = round(w,3)

            freq1_Thz = convert_freq_to_Thz(get_freqs(hx_min, hy_min, a_max, w, wvg_height), a_max)  # getting lowest possible frequencies for width w 
            print(" -------------------- w loop ----------------------")
            if ((freq1_Thz[0]  > f_target_Thz)):

                continue
            else:

            #---------------------------#
            #   LATTICE CONSTANT LOOP          
            #---------------------------#

                for a in np.arange(a_min , a_max, del_a):
                    print("------------------------a loop-------------------------")
                    #a = round(a,3)

                    freq2_Thz = convert_freq_to_Thz(get_freqs(hx_min, hy_min, a, w, wvg_height), a)  # getting lowest possible frequences for (w,a)

                    if ( freq2_Thz[0] > f_target_Thz):
                        continue

                    #---------------------------#
                    #       HY LOOP          
                    #---------------------------#

                    for hy in np.arange(hy_min, w - 0.1 , del_hy):
                        print("-----------------------hy loop----------------------")
                        #hy = round(hy,3)

                        freq3_Thz = convert_freq_to_Thz(get_freqs(hx_min, hy, a, w, wvg_height), a)  # getting lowest possible frequences for (w,a, hy)

                        if ( (freq3_Thz[0]  > f_target_Thz) ):
                            continue

                        #---------------------------#
                        #       HX LOOP          
                        #---------------------------#

                        for hx in np.arange(hx_min, a - 0.07, del_hx ):
                            count  = 0
                            print(" ---------------- hx loop ---------------------")
        #                         a = round(a,3)
        #                         hy = round(hy,3)
        #                         hx = round(hx,3)
        #                         w = round(w,3)

                            freq4_Thz = convert_freq_to_Thz(get_freqs(hx, hy, a, w, wvg_height), a )

                            if ( freq4_Thz[0] > f_target_Thz):  # if w_target is outside the bandgap for 2 consecutive runs, break outside the loop
                                count = count + 1

                                if count == 2 :
                                    break

                            else:

                                if (f_target_Thz < freq4_Thz[1])  and (f_target_Thz > freq4_Thz[0]):  # final check to see that the target frequency is in the bandgap

                                    print(" ------------------- new gamma ------------------- at hx = {}, hy = {}, a = {}, w = {}".format(hx,hy,a,w)) 

                                    f_mid = (freq4_Thz[0] + freq4_Thz[1])/2            
                                    diff = freq4_Thz[0] - freq4_Thz[1]
                                    delta = 1 - (f_target_Thz/ f_mid)

                                    gamma =  math.sqrt(abs(( 0.5 * diff/ f_mid ) ** 2 - delta**2 ))

                                    mirror_strength.append(gamma)
                                    index_count = index_count + 1
                                    index.append(index_count)

                                    dset[int((w - w_min) / del_w + 0.1), 
                                         int((a - a_min) / del_a + 0.1), 
                                         int((hy - hy_min) / del_hy + 0.1), 
                                         int((hx - hx_min) / del_hx + 0.1) ] = round(gamma,4)

                                    if gamma > gamma_max:
                                        print(" ------------------- new gamma max ------------------- at hx = {}, hy = {}, a = {}, w = {}, gamma = {}".format(hx,hy,a,w, gamma))
                                        gamma_max = gamma
                                        parameters.append( (round(hx, 4), 
                                                           round(hy, 4), 
                                                           round(a,  4), 
                                                           round(w,  4),
                                                           freq4_Thz,
                                                           round(gamma_max,5)))


                                else:
                                    continue

        with open("parameters230.txt", "w") as file1: 
            for parameter in parameters:
        # Writing data to a file 
                file1.write("hx = {}, hy = {}, a = {}, w = {}, gamma = {}".format(*parameter)) 
                file1.write("\n")

    
def get_index(w, a, hy, hx):
    
        #------------- DEFAULTS -------------------------------------------#
        #  del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025
        #  w_max ,a_max = 0.7, 0.45
        #  w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05
        #------------------------------------------------------------------#
        
        del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025
        w_max ,a_max = 0.7, 0.45
        w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05
        
        index_a = int((a - a_min)/del_a + 0.1)
        index_w = int((w - w_min) / del_w + 0.1)
        index_hy = int((hy - hy_min) / del_hy + 0.1)
        index_hx = int((hx - hx_min) / del_hx + 0.1)
        
        return index_w, index_a, index_hy, index_hx

def get_excitation_mode_from_string(mode_string):
    if mode_string == 'Ex':
        excitation_mode = mp.Ex
    elif mode_string == 'Ey':
        excitation_mode = mp.Ey
    elif mode_string == 'Ez':
        excitation_mode = mp.Ez
    elif mode_string == 'Hx':
        excitation_mode = mp.Hx
    elif mode_string == 'Hy':
        excitation_mode = mp.Hy
    elif mode_string == 'Hz':
        excitation_mode = mp.Hz
    else:
        warnings.warn("Mode string not understood. mp.Hz set as excitation mode", UserWarning)
        excitation_mode = mp.Hz
    return excitation_mode


# Boundary layers
def get_boundary_layer(dpml=1, sim2d=False):
    """
    Get boundary layer of thickness dpml.
    """
    if sim2d:
        boundary_layers = [mp.PML(dpml, direction=mp.X),
                           mp.PML(dpml, direction=mp.Y)]
    else:
        boundary_layers = [mp.PML(dpml, direction=mp.X),
                           mp.PML(dpml, direction=mp.Y),
                           mp.PML(dpml, direction=mp.Z)]
    return boundary_layers

def index_to_material(element):
    if isinstance(element, mp.Medium):
        return element
    else:
        return mp.Medium(index=element)

def get_value_from_index(index, param = 'a'):
    del_a = 0.001
    del_hy = 0.025
    del_hx = 0.025 
    del_w = 0.05
    
    a_min = 0.25
    a_max = 0.45        # upper limit of the sweep of a 

    w_min = 0.65         #  lower limit of w 
    w_max = 0.75        #  upper limit of w 

    hx_min = 0.05        # lower limit of the sweep of a
    hy_min = 0.1        #  lower limit of hy 
    
    if param == 'a':
        return (a_min + del_a * index)
    if param == 'hy':
        return (hy_min + del_hy * index)
    if param == 'hx':
        return (hx_min + del_hx * index)
    if param == 'w':
        return (w_min + del_w * index)
 
def get_gamma(freq, a, f_target = 1 / (1.54) * 1.01):
    '''
    f_target is in terms of 1/lambda * 1.01
    freq is in terms of 2pi * c/ a
    
    RETURNS: 
    
    The mirror strength for the input tuple (freq) of the dielectric and air band edge frequencies with
    f_target
    '''
    import math 
    
    freq[0] = convert_freq_to_Thz(freq[0], a)
    freq[1] = convert_freq_to_Thz(freq[1], a)
    f_target = convert_freq_to_Thz(1/1.54)* 1.01
     


    if f_target < freq[0] or f_target > freq[1] :  # if f_target outside bandgap
        return 0
    
    w_mid = (freq[0] + freq[1])/2            
    diff = freq[0] - freq[1]

    delta = 1 - (w_target/ (w_mid))

    gamma =  math.sqrt((( 0.5 * diff/ w_mid ) ** 2 - delta**2 ))
    
    return gamma
