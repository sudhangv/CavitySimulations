import meep as mp
import numpy as np
import warnings
import os

from .lattice import OneDLattice
from .waveguide import *
from ..utilities.utilities import *

del_w, del_a, del_hy, del_hx = 0.05, 0.001, 0.025, 0.025
w_max ,a_max = 0.7, 0.45
w_min, a_min, hy_min, hx_min = 0.65, 0.25, 0.1, 0.05

def get_perturb_param(to_perturb, freq_data, gamma_data, w,  a,   hy,  hx , 
                      target_f_Thz, f_perturb_lower_Thz, f_perturb_upper_Thz, tol_Thz = 1):
    
    '''
    Here upper_param and lower_param correspond to the perturbed parameters that result in resonant 
    frequencies of f_perturb_upper_Thz and f_perturb_lower_Thz respectively '''
    
    del_a, del_hy, del_hx = 0.001,0.025 ,0.025
    lower_param = -1
    upper_param = -1
    index_w, index_a, index_hy, index_hx = get_index(w = w, a = a, hy = hy, hx = hx)
    
    if to_perturb == 'hx':
        
        run = True
        i = 0
        
        while run:
            i = i + 1
            #----------------- Checking by increasing the parameter to be inspected -----------------------#
            freq_to_check_Thz = data_freq[index_w, index_a, index_hy, index_hx + i][0] 
                                             
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = hx + (i * del_hx)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = hx + (i * del_hx)
            
            #----------------- Checking by decreasing the parameter to be inspected -----------------------#
            freq_to_check_Thz = data_freq[index_w, index_a, index_hy, index_hx - i][0] 
                                                    
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = hx + (-i * del_hx)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = hx + (-i * del_hx)
                
            if upper_param != -1 and lower_param != -1:
                run = False
    
    if to_perturb == 'hy':
        
        run = True
        i = 0
        
        while run:
            i = i + 1
            #----------------- Checking by increasing the parameter to be inspected -----------------------#
            freq_to_check_Thz = data_freq[index_w, index_a, index_hy  + i, index_hx][0]
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = hy + (i * del_hy)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = hy + (i * del_hy)
            
            #----------------- Checking by decreasing the parameter to be inspected -----------------------#
            freq_to_check_Thz = data_freq[index_w, index_a, index_hy - i ,index_hx][0]
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = hy + (-i * del_hy)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = hy + (-i * del_hy)
            
            if upper_param != -1 and lower_param != -1:
                run = False
                
        if to_perturb == 'a':
        
        run = True
        i = 0
        
        while run:
            i = i + 1
            #----------------- Checking by increasing the parameter to be inspected -----------------------#
            freq_to_check_Thz = freq = data_freq[index_w, index_a + i, index_hy , index_hx][0]
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = hy + (i * del_hy)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = hy + (i * del_hy)
            
            #----------------- Checking by decreasing the parameter to be inspected -----------------------#
            freq_to_check_Thz = data_freq[index_w, index_a + i, index_hy  ,index_hx][0]
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = a + (-i * del_a)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = a + (-i * del_a)
            
            if upper_param != -1 and lower_param != -1:
                run = False
            
    return lower_param, upper_param

def get_mirror_param(to_perturb, freq_data, gamma_data, w,  a_mirror,   hy,  hx , 
                      target_f_Thz, f_perturb_lower_Thz, f_perturb_upper_Thz, steps = 20):
    '''
    This function takes the new central segment parameters (hx, hy, a_mirror ,w ) and 
    the new resonant frequencies (f_perturb_lower_Thz,f_perturb_upper_Thz) and returns
    the parameters for the mirror segment that maximize the mirror strength
    
    steps: (a_mirror +- steps * del_a) specify the range within which the new maxima for gamma will be searched
    '''
    
    index_w, index_a_mirror, index_hy, index_hx = get_index(w = w, a = a_mirror, hy = hy, hx = hx)
    
    if to_perturb == 'hx':
        
        freq_mirror = freq_data[index_w, index_a_mirror , index_hy, index_hx]
        gamma_mirror = get_gamma(band_edge_f = freq_plus_Thz, check_freq = f_perturb_lower_Thz)
        
        gamma_max_upper = gamma_mirror # upper denotes that the quantity is relevant for the increased resonant frequency           
        gamma_max_lower = gamma_mirror # lower denotes that the quantity is relevant for the decreased resonant frequency
        
        for i in range( steps + 1):
            
            # The plus(minus) here denote that the quantity refers to increasing(decreasing) a_mirror 
            # upper and lower retain their meaning from the previous comments
            
            freq_plus_Thz = freq_data[index_w, index_a_mirror + i, index_hy, index_hx]  
            freq_minus_Thz = freq_data[index_w, index_a_mirror - i, index_hy, index_hx]
            
            gamma_plus_lower = get_gamma(band_edge_f = freq_plus_Thz, check_freq = f_perturb_lower_Thz)
            gamma_minus_lower = get_gamma(band_edge_f = freq_minus_Thz, check_freq = f_perturb_lower_Thz)
            
            gamma_plus_upper = get_gamma(band_edge_f = freq_plus_Thz, check_freq = f_perturb_upper_Thz)
            gamma_minus_upper = get_gamma(band_edge_f = freq_minus_Thz, check_freq = f_perturb_upper_Thz)
            
            if gamma_plus_lower > gamma_max_lower:
                gamma_max_lower = gamma_plus_lower
                index_lower = i
            
            if gamma_minus_lower > gamma_max_lower:
                gamma_max_lower = gamma_minus_lower
                index_lower = i
            
            if gamma_plus_upper > gamma_max_upper:
                gamma_max_upper = gamma_plus_upper
                index_upper = i
                
            if gamma_minus_upper > gamma_max_upper:
                gamma_max_upper = gamma_minus_upper
                index_upper = i
    
        a_mirror_new_lower = a_mirror + index_lower * del_a
        a_mirror_new_upper = a_mirror + index_upper * del_a
        
        return a_mirror_new_lower,a_mirror_new_upper
        
def _a_poly_tapering(geom=None, n_segments=20, material_holes=mp.vacuum):
    
    filename = "/bandstructure_data/<filename>"
    
    hf = h5py.File(filename, 'r')
    gamma_data = np.array( hf.get("gamma"))
    freq_data = np.array( hf.get("freq"))
    hf.close()
    
    if geom is None:
         geom = []
            
    material_holes = index_to_material(material_holes)
    
    #--------------- These are the parameters for DESIGNED geometry, which we want to perturb ------------
               
    hx = 0.143
    hy = 0.315
    w = 0.65
    a_cen = 0.313
    a_mirror = 0.361
    Lx = 20
    _cavity = OneDLattice(Lx = Lx)
    _n_taper = 10
    h  = 0.25
    
    # ------------------------------ PERTURBATION HERE -------------------------------------------- #
    
    to_perturb = "hx"            # one of hx, hy and a
    perturb_range = 0.04         # edges of the wavelength window will be (target_lambda +- perturb_range) 
    tol_Thz = 1                  # tolerance in Thz to select the perturbed segment parameters
    
    target_wvl = 1.54            # vaccum wavelength ( in um ) of the unperturbed cavity design
    target_f = 1/target_wvl
    target_f_Thz =  convert_freq_to_Thz(target_f)
    
    f_perturb_lower = 1 / (target_wvl + perturb_range )           # target_f - perturbation
    f_perturb_upper = 1 / (target_wvl - perturb_range )           # target_f + perturbation
    
    f_perturb_lower_Thz =  convert_freq_to_Thz(f_perturb_lower)
    f_perturb_upper_Thz =  convert_freq_to_Thz(f_perturb_upper)
    
    lower_param, upper_param = get_perturb_param(to_perturb = to_perturb , 
                                                 freq_data=  freq_data , gamma_data= gamma_data,
                                                 w = w,  a = a,   hy = hy,  hx = hx , 
                                                 target_f_Thz= target_f_Thz, 
                                                 f_perturb_lower_Thz= f_perturb_lower_Thz, 
                                                 f_perturb_upper_Thz= f_perturb_upper_Thz, tol_Thz = tol_Thz)
    
    a_mirror_new_lower,a_mirror_new_upper = get_mirror_param(to_perturb, freq_data, gamma_data, 
                                                             w,  a_mirror,   hy,  hx , 
                                                             target_f_Thz, f_perturb_lower_Thz, 
                                                             f_perturb_upper_Thz, steps = 20)
    
    # Here upper_param and lower_param correspond to the perturbed parameters that result in resonant 
    
    # frequencies of f_perturb_upper_Thz and f_perturb_lower_Thz respectively 
    
    # a_mirror_new_lower,a_mirror_new_upper refer to the new values of a_mirror that maximise gamma_mirror for
    # the the two ends of the perturbation window
    
    #----------------------------------------------------------------------------------------------------#
    
    _cavity.polynomial_elliptical_hole_taper(_n_taper, hx, hy, w, a_cen, a_mirror )
    _cavity.apply_poly_spacing()
    
    print("--------------------------------------------------------------------------------------------------------")
    print(" Poly Tapering : hx = {}, hy = {}, w = {}, h= {}, a_cen  = {},  a_mirror = {}, n_taper = {}, Lx = {}".format(hx, hy,w, h, a_cen,a_mirror,_n_taper,Lx))
    print("--------------------------------------------------------------------------------------------------------")

    # cavity holes
    for x, y, z, hx, hy in _cavity.coordinates:
        # holes are completely filled with tuning material:
        geom.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))

        geom.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(-x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))

    length = 2 * max(_cavity.coordinates[:, 0])

    return geom, length

def _a_pow_tapering(geom=None, n_segments=20, material_holes=mp.vacuum):
  
    if geom is None:
        geom = []
        
    material_holes = index_to_material(material_holes)
    hx = 0.143
    hy = 0.315                                                                                                              
    w = 0.65
    a_cen = 0.303
    a_mirror = 0.346
    Lx = 20
    h = 0.25
    _cavity = OneDLattice(Lx = Lx)
    _n_taper = 10
    _cavity.pow_degree_a_taper(_n_taper, 
                               hx = hx, 
                               hy = hy, 
                               w = w, 
                               a_center = a_cen,
                               a_mirror = a_mirror,
                               pow = 2)
    
    _cavity.apply_pow_spacing()
    
    print("--------------------------------------------------------------------------------------------------------")
    print(" Pow Tapering : hx = {}, hy = {}, w = {}, h = {}, a_cen  = {},  a_mirror = {}, n_taper = {}, Lx = {}".format(hx, hy, w, h, a_cen, a_mirror, _n_taper, Lx))    
    print("--------------------------------------------------------------------------------------------------------") 

    # cavity holes
    for x, y, z, hx, hy in _cavity.coordinates:
        # holes are completely filled with tuning material:
        geom.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))

        geom.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(-x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))

    length = 2 * max(_cavity.coordinates[:, 0])

    return geom, length

def _a_normal_tapering(geom=None, n_segments=20, material_holes=mp.vacuum):
     
    if geom is None:
        geom = []
    material_holes = index_to_material(material_holes)

    _cavity = OneDLattice(Lx = n_segments)
    _cavity.normal_spacing(a = 0.303, hx = 0.143, hy = 0.315)
    
    for x, y, z, hx, hy in _cavity.coordinates:
        # holes are completely filled with tuning material:
        geom.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))
        geom.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(-x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))
    length = 10                             
    return geom, length

def a_pow_tapered_cavity(geom = None, n_segments=20, waveguide_parameters= None, substrate_parameters=None):
    """
    Returns the geometry objects for a the air holes of 1D phc cavity with tapered lattice constants.
    """
    if geom is None:
        geom = []

    if waveguide_parameters is None:
        waveguide_parameters = {}

    if substrate_parameters is None:
        substrate_parameters = {}

    geom = add_waveguide_1d(geom=geom)

    geom, _ = _a_pow_tapering(geom=geom, n_segments=n_segments)

    # geom = add_substrate(geom=geom, **substrate_parameters)

    return geom

def a_normal_cavity(geom = None, n_segments=20, waveguide_parameters= None, substrate_parameters=None):    
    
    if geom is None:
        geom = []

    if waveguide_parameters is None:
        waveguide_parameters = {}

    if substrate_parameters is None:
        substrate_parameters = {}

    geom = add_waveguide_1d(geom=geom)

    geom, _ = _a_normal_tapering(geom=geom, n_segments=n_segments)

    # geom = add_substrate(geom=geom, **substrate_parameters)                                                           

    return geom
                                                                                        
                                                                                                                                                               
    
               

def a_poly_tapered_cavity(geom=None, n_segments=20, waveguide_parameters=None, substrate_parameters=None):
    if geom is None:
        geom = []

    if waveguide_parameters is None:
        waveguide_parameters = {}

    if substrate_parameters is None:
        substrate_parameters = {}

    geom = add_waveguide_1d(geom=geom, **waveguide_parameters)

    geom, _ = _a_poly_tapering(geom=geom, n_segments=n_segments)

    geom = add_substrate(geom=geom, **substrate_parameters)

    return geom
                                                                                                      
