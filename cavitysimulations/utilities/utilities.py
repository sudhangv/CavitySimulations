import meep as mp
import warnings


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