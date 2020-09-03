import meep as mp
import numpy as np
import warnings
import os

from .lattice import OneDLattice
from .waveguide import *
from ..utilities.utilities import *


def _a_poly_tapering(geom=None, n_segments=20, material_holes=mp.vacuum):
    if geom is None:
         geom = []
    material_holes = index_to_material(material_holes)
    hx = 0.143
    hy = 0.315
    w = 0.65
    a_cen = 0.313
    a_mirror = 0.361
    Lx = 20
    _cavity = OneDLattice(Lx = Lx)
    _n_taper = 10
    h  = 0.25
    _cavity.polynomial_elliptical_hole_taper(_n_taper, hx, hy, w, a_cen, a_mirror )
    _cavity.apply_poly_spacing()
    print("--------------------------------------------------------------------------------------------------------")
    print(" Poly Tapering : hx = {}, hy = {}, w = {}, h= {}, a_cen  = {},  a_mirror = {}, n_taper = {}, Lx = {}".format(hx, hy,w, h, a_cen,a_mirror,_n_taper,Lx))
    print("--------------------------------------------------------------------------------------------------------")
#print(_cavity.coordinates)
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
#print(_cavity.coordinates)
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
                                                                                                      
