def add_substrate(geom=None, waveguide_height=0.22, substrate_height=5, material=mp.Medium(index=1.44)):
    """
    Creates a (by default SiO2) substrate.
    If the unlike case occurs that we need to change the center, we could again make everything relative to the center
    of the substrate.
    """
    if geom is None:
        geom = []

    _center = mp.Vector3(0, 0, -waveguide_height/2-substrate_height/2)

    geom.append(mp.Block(material=index_to_material(material),
                         center=_center,
                         size=mp.Vector3(mp.inf, mp.inf, substrate_height)))

    return geom
