import meep as mp
import warnings
import numpy as np
import math


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
    
def convert_freq_to_Thz(freq, a = 0):
    '''
    Converts a MEEP frequency ( units of 2pi*c/a or 1/lambda_air) to frequency (not angular velocity) in units of Thz
    freq ------> freq_in_THz
    '''
    
    freq = np.array(freq)
    

    if a != 0:
        return ( freq * 3 / a * (10**2))
    else:
        return ( freq * 3 * 10**2)
    
    
def get_gamma_from_Thz(band_edge_f, check_freq):
    if (check_freq < band_edge_f[1])  and (check_freq > band_edge_f[0]):
        
        f_mid = (band_edge_f[0] + band_edge_f[1])/2            
        diff = band_edge_f[0] - band_edge_f[1]
        delta = 1 - (check_freq/ f_mid)

        gamma = math.sqrt(abs(( 0.5 * diff/ f_mid ) ** 2 - delta**2 ))
    
    else:
        gamma = 0
        
    return gamma