import meep as mp
import numpy as np
import warnings
import os
import h5py

from lattice import OneDLattice
#from waveguide import *
from utilities import *
from sweep_util import *



def get_perturb_param(to_perturb, w,  a,   hy,  hx , h,
                      target_f_Thz, f_perturb_lower_Thz, f_perturb_upper_Thz, tol_Thz = 1,
                      substrate = False, mode = "yO"):
    
    '''
    Here upper_param and lower_param correspond to the perturbed parameters that result in resonant 
    frequencies of f_perturb_upper_Thz and f_perturb_lower_Thz respectively 
    '''
    del_a, del_hy, del_hx = 0.004, 0.004, 0.004  # increment values for each of the parameters
    
    step_w, step_a, step_hy, step_hx = 0.002, 0.002, 0.002, 0.002
    lower_param = -1             # flags to indicate whether a suitable parameter has been found
    upper_param = -1             # flags to indicate whether a suitable parameter has been found
    
    #  index_w, index_a, index_hy, index_hx = get_index(w = w, a = a, hy = hy, hx = hx)
    
    if to_perturb == 'hx':
        
        run = True
        i = 0
        
        while run:
            i = i + 1
            #----------------- Checking by increasing the parameter to be inspected -----------------------#
            
            freq_to_check_Thz = get_freq_Thz(w = w, a = a, hy = hy, hx = hx + i*del_hx, h = h ,
                                             substrate = substrate, mode = mode)[0]
                                             
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                
                upper_param = hx + (i * del_hx)
                # if the resonant frequency for the new parameter lies within a tolerance from the peturbed resonant freqeuncy
                # we have found the desired perturbed parameter
                                                            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                
                lower_param = hx + (i * del_hx)
            
            #----------------- Checking by decreasing the parameter to be inspected -----------------------#
 
            freq_to_check_Thz = get_freq_Thz(w = w, a = a, hy = hy, hx = hx - i*del_hx, h =h,
                                             substrate = substrate, mode = mode)[0]
                                                    
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

            freq_to_check_Thz = get_freq_Thz(w =w , a = a, hy = hy + i*del_hy, hx = hx, h = h,
                                            substrate = substrate, mode = mode)
            
            if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                upper_param = hy + (i * del_hy)
            
            if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                lower_param = hy + (i * del_hy)
            
            #----------------- Checking by decreasing the parameter to be inspected -----------------------#
            
            freq_to_check_Thz = get_freq_Thz(w =w , a = a, hy = hy - i*del_hy, hx = hx, h = h,
                                            substrate = substrate, mode = mode)
            
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
                
             
                freq_to_check_Thz = get_freq_Thz(w = w, a = a + i* del_a, hy = hy , hx = hx, h= h,
                                                substrate = substrate, mode =mode)
                
                if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                    upper_param = hy + (i * del_hy)

                if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                    lower_param = hy + (i * del_hy)

                #----------------- Checking by decreasing the parameter to be inspected -----------------------#
                
              
                freq_to_check_Thz = get_freq_Thz(w = w, a = a - i* del_a, hy = hy , hx = hx, h =h,
                                                substrate = substrate, mode = mode)
                
                if freq_to_check_Thz < (f_perturb_upper_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_upper_Thz - tol_Thz):
                    upper_param = a + (-i * del_a)

                if freq_to_check_Thz < (f_perturb_lower_Thz + tol_Thz) and freq_to_check_Thz > (f_perturb_lower_Thz - tol_Thz):
                    lower_param = a + (-i * del_a)

                if upper_param != -1 and lower_param != -1:
                    run = False
            
    return lower_param, upper_param

def get_mirror_param(to_perturb, w,  a_mirror,   hy,  hx , h, substrate, 
                     target_f_Thz, f_perturb_lower_Thz, f_perturb_upper_Thz, 
                     mode, step_size, n_step, lower_param, upper_param ):
    '''
    This function takes the new central segment parameters (hx, hy, a_mirror ,w ) and 
    the new resonant frequencies (f_perturb_lower_Thz,f_perturb_upper_Thz) and returns
    the parameters for the mirror segment that maximize the mirror strength
    
    steps: (a_mirror +- steps * del_a) specify the range within which the new maxima for gamma will be searched
    '''
    
    lower_param = lower_param
    upper_param = upper_param
        
    index_w, index_a_mirror, index_hy, index_hx = get_index(w = w, a = a_mirror, hy = hy, hx = hx)
    
    step_size = step_size
    n_step = n_step
    
    if to_perturb == 'hx':
        hx_lower = lower_param
        hy_lower = hy
        hx_upper = upper_param
        hy_upper = hy
    
    if to_perturb == 'hy':
        hy_lower = lower_param
        hx_lower = hx
        hy_upper = upper_param
        hx_upper = hx   
    
    if to_perturb == 'a':
        hy_lower = hy
        hy_upper = hy
        hx_lower = hx
        hx_upper = hx
        
    freq_mirror_lower = get_freq_Thz(w = w, a = a_mirror, hy =  hy_lower, hx = hx_lower, h = h, substrate = substrate,  mode = mode)  
    freq_mirror_upper = get_freq_Thz(w = w, a = a_mirror, hy =  hy_upper, hx = hx_upper, h = h, substrate = substrate,  mode = mode)  
         
    gamma_mirror_upper = get_gamma_from_Thz(band_edge_f = freq_mirror_upper, check_freq = f_perturb_upper_Thz)
    gamma_mirror_lower = get_gamma_from_Thz(band_edge_f = freq_mirror_lower, check_freq = f_perturb_lower_Thz)

    # freq_mirror : Previous mirror strength of the mirror segment
    # gamma_mirror_upper(lower) denote gamma for the mirror segment at the upper(lower) ends of the resonant 
    # frequency perturbation

    gamma_max_upper = gamma_mirror_upper # upper denotes that the quantity is relevant for the increased resonant frequency           
    gamma_max_lower = gamma_mirror_lower # lower denotes that the quantity is relevant for the decreased resonant frequency
    print("gamma_max_lower :" + str(gamma_max_lower) )
    print("gamma_max_upper :" + str(gamma_max_upper) )
    index_lower = 0
    index_upper = 0
    
    for i in range( n_step + 1):

        # The plus(minus) here denote that the quantity refers to increasing(decreasing) a_mirror 
        # upper and lower retain their meaning from the previous comments

        freq_plus_upper_Thz = get_freq_Thz( w = w, a = a_mirror + i* step_size, hy = hy_upper, hx = hx_upper, 
                                           h =h, substrate = substrate,  mode = mode)

        freq_minus_upper_Thz = get_freq_Thz( w = w, a = a_mirror - i* step_size, hy = hy_upper, hx = hx_upper, 
                                           h = h, substrate = substrate,  mode = mode)

        freq_plus_lower_Thz = get_freq_Thz( w = w, a = a_mirror + i* step_size, hy = hy_lower, hx = hx_lower, 
                                           h =h, substrate = substrate,  mode = mode)

        freq_minus_lower_Thz = get_freq_Thz( w = w, a = a_mirror - i* step_size, hy = hy_lower, hx = hx_lower, 
                                           h = h, substrate = substrate,  mode = mode)

        gamma_plus_upper = get_gamma_from_Thz(band_edge_f = freq_plus_upper_Thz, check_freq = f_perturb_upper_Thz)
        #print("gamma_plus_upper :" + str(gamma_plus_upper) )
        gamma_minus_upper = get_gamma_from_Thz(band_edge_f = freq_minus_upper_Thz, check_freq = f_perturb_upper_Thz)
        #print("gamma_minus_upper :" + str(gamma_minus_upper) )
        gamma_plus_lower = get_gamma_from_Thz(band_edge_f = freq_plus_lower_Thz, check_freq = f_perturb_lower_Thz)
        #print("gamma_plus_lower :" + str(gamma_plus_lower) )
        gamma_minus_lower = get_gamma_from_Thz(band_edge_f = freq_minus_lower_Thz, check_freq = f_perturb_lower_Thz)
        #print("gamma_minus_lower :" + str(gamma_minus_lower) )


        if gamma_plus_lower > gamma_max_lower:
            gamma_max_lower = gamma_plus_lower
            index_lower = i
            # index relative to the orignal a_mirror
        if gamma_minus_lower > gamma_max_lower:
            gamma_max_lower = gamma_minus_lower
            index_lower = -i 

        if gamma_plus_upper > gamma_max_upper:
            gamma_max_upper = gamma_plus_upper
            index_upper = i

        if gamma_minus_upper > gamma_max_upper:
            
            gamma_max_upper = gamma_minus_upper
            index_upper = -i

    a_mirror_new_lower = a_mirror + index_lower * step_size
    a_mirror_new_upper = a_mirror + index_upper * step_size
    print("FOUNDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD")
    return (a_mirror_new_lower, a_mirror_new_upper)
        
def _a_poly_tapering(geom=None, n_segments=20, material_holes=mp.vacuum):
    
    filename = "bandstructure_data/sweep_data.hdf5"
    
    hf = h5py.File(filename, 'r')
    gamma_data = np.array( hf.get("gamma"))
    freq_data = np.array( hf.get("freq"))
    hf.close()
    
    if geom is None:
         geom = []
            
    material_holes = index_to_material(material_holes)
    
    #--------------- These are the parameters for DESIGNED geometry, which we want to perturb ------------ #
               
    
    hx = 0.225
    hy = 0.4
    a_mirror = 0.414
    a_cen = 0.385
    w = 0.65
    h  = 0.19
    mode = "zEyO"
    
    substrate = True
    
    Lx = 20
    _n_taper = 10
    
    _cavity = OneDLattice(Lx = Lx)
    
    # ------------------------------ PERTURBATION HERE -------------------------------------------- #
    # Essentially how this works is that at the end of all the repetitive code, you get 4 values:
    # (1) lower_param : the value of the parameter that gives a resonant frequency slightly lower than the earlier one
    # (2) upper_param : the value of the parameter that gives a resonant frequency slightly higher than the earlier one
    # (3) a_mirror_new_lower : the new a_mirror that maximizes gamma for positive delta_freq and new hy/hx 
    # (4) a_mirror_new_upper : the new a_mirror that maximizes gamma for negative delta_freq and new hy/x
    
    to_perturb = "hx"            # one of hx, hy and a
    perturb_range = 0.04         # edges of the wavelength (in um) window will be (target_lambda +- perturb_range) 
    tol_Thz = 1                 # tolerance in Thz to select the perturbed segment parameters
    
    target_wvl = 1.54            # vaccum wavelength ( in um ) of the unperturbed cavity design
    target_f = 1/target_wvl
    target_f_Thz =  convert_freq_to_Thz(target_f)
    
    f_perturb_lower = 1 / (target_wvl + perturb_range )           # target_f - perturbation
    f_perturb_upper = 1 / (target_wvl - perturb_range )           # target_f + perturbation
    
    f_perturb_lower_Thz =  convert_freq_to_Thz(f_perturb_lower)
    f_perturb_upper_Thz =  convert_freq_to_Thz(f_perturb_upper)
    
    lower_param, upper_param = get_perturb_param(to_perturb = to_perturb , 
                                                 w = w,  a = a_cen,   hy = hy,  hx = hx , h =h, 
                                                 target_f_Thz= target_f_Thz, substrate = substrate,
                                                 f_perturb_lower_Thz= f_perturb_lower_Thz, 
                                                 f_perturb_upper_Thz= f_perturb_upper_Thz, tol_Thz = tol_Thz,
                                                 mode = mode)


    step_size = 0.004            # size of each of the n_steps
    n_step = 7                   # number of steps to search for the optimal gamma for the mirror segment
    
    a_mirror_new_lower, a_mirror_new_upper = get_mirror_param(to_perturb = to_perturb,
                                                             w = w,  a_mirror = a_mirror,   
                                                             hy = hy,  hx = hx , h = h,
                                                             target_f_Thz = target_f_Thz, 
                                                             f_perturb_lower_Thz = f_perturb_lower_Thz, 
                                                             f_perturb_upper_Thz = f_perturb_upper_Thz, 
                                                             mode = mode, substrate = substrate,
                                                             step_size = step_size, n_step = n_step,
                                                             lower_param = lower_param, upper_param = upper_param)

    
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
    a_cen = 0.36
    a_mirror = 0.414
    Lx = 20
    h = 0.19
    
    substrate = False
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
                                                                                                      
