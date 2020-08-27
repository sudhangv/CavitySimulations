import meep as mp
import numpy as np
import warnings
import os

from .lattice import Lattice
from .waveguide import *
from ..utilities.utilities import *


def quadratic_radius_1d_symmetric(n_segments=20, a=.33, hx=.21, hy=.21, h_rel=.7, z_center=0,
                                  d_tuning=0, material_holes=mp.vacuum):
    """
    Returns the geometry objects for a the air holes of 1D phc cavity with quadraticly tapered radii.

    TODO: change between Ellipsoids and Cylinders if hx==hy.
    """

    material_holes = index_to_material(material_holes)

    cavity = Lattice(n_segments)
    cavity.set_z(z_center)
    cavity.quadratic_hole_taper(h_rel)

    cavity_holes = []

    # cavity holes
    for x, y, z, hx, hy in cavity.coordinates:
        # holes are completely filled with tuning material:
        cavity_holes.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))

        cavity_holes.append(mp.Ellipsoid(material=material_holes,
                                         center=mp.Vector3(-x, y, z),
                                         size=mp.Vector3(hx, hy, mp.inf)))

    length = 2 * max(cavity.coordinates[:, 0])

    return cavity_holes, length


def _a_tapering(geom=None, n_segments=20, material_holes=mp.vacuum):
    """
    Returns the geometry objects for a the air holes of 1D phc cavity with tapered lattice constants.
    TODO; allow change of resonant frequency to be able to make sweeps.
    """
    if geom is None:
        geom = []
    material_holes = index_to_material(material_holes)

    print(os.getcwd())

    _cavity = Lattice(Lx=n_segments)
    _cavity.polynomial_elliptical_hole_taper(20, 0.3, 0.5, 0.7)
    _cavity.apply_poly_spacing()

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


def a_tapered_cavity(geom=None, n_segments=20, waveguide_parameters=None, substrate_parameters=None):
    """
    Returns the geometry objects for a the air holes of 1D phc cavity with tapered lattice constants.
    """
    if geom is None:
        geom = []

    if waveguide_parameters is None:
        waveguide_parameters = {}

    if substrate_parameters is None:
        substrate_parameters = {}

    geom = add_waveguide_1d(geom=geom, **waveguide_parameters)

    geom, _ = _a_tapering(geom=geom, n_segments=n_segments)

    # geom = add_substrate(geom=geom, **substrate_parameters)

    return geom
