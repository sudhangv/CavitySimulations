import meep as mp
import numpy as np
import warnings

# IMPORT COSTUM CLASSES
import Simulations.lib.Lattice as Lattice
# import Simulations.lib.SemData as SemData
import Simulations.lib.costum_os_utilities as c_os


def main(args=None):
    fcen = 1 / 1.54
    df = 0.1
    wvg_width = .7
    wvg_height = .22
    dpml = 1
    dpad = 2
    sim2d = False
    time_after_source = 500
    eps_only = False
    subz = True
    resolution = 30
    _nSi = 3.45
    mode = get_excitation_mode("Ey")

    # we want to investigate tuning. Change this to add a different material instead of air:
    cover_material = mp.Medium(index=1.11)

    if args is not None:
        hx = args.diameter
        hy = args.diameter
        h_rel = args.h_rel
        n = args.n
        _nSi = args.n_si

    nSi = mp.Medium(index=_nSi)

    wvg = waveguide_1d(wvg_width=wvg_width, wvg_height=wvg_height, material=nSi, )
    cav, cav_length = quadratic_radius_1d_symmetric(n_segments=n, hx=hx, hy=hy, h_rel=h_rel,
                                                    material_tuning=cover_material)

    sx = dpml + dpad + cav_length + dpad + dpml
    sy = dpml + dpad + wvg_width + dpad + dpml
    sz = dpml + dpad + wvg_height + dpad + dpml

    subs_height = (sz - wvg_height) / 2
    subs = substrate(substrate_height=subs_height, center=mp.Vector3(0, 0, -sz / 2 + subs_height / 2))

    geometry = wvg + cav

    symmetries = [mp.Mirror(mp.X, +1), mp.Mirror(mp.Y, -1), mp.Mirror(mp.Z, +1)]

    if subz:
        geometry = wvg + cav + subs
        symmetries = [mp.Mirror(mp.X, +1), mp.Mirror(mp.Y, -1)]

    boundary = get_boundary_layer(dpml, sim2d=sim2d)

    if sim2d:
        sz = 0

    cell_size = mp.Vector3(sx, sy, sz)

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=mode,
                         center=mp.Vector3())]

    sim = mp.Simulation(cell_size=cell_size,
                        geometry=geometry,
                        sources=sources,
                        boundary_layers=boundary,
                        symmetries=symmetries,
                        resolution=resolution,
                        progress_interval=100,
                        default_material=cover_material)

    f = plt.figure(figsize=(10, 10))
    sim.plot2D(ax=f.gca(), output_plane=mp.Volume(center=mp.Vector3(0, 0, 0),
                                                  size=mp.Vector3(sx, sy, 0)))
    f.savefig("xy_plane.pdf", format="PDF")

    f = plt.figure(figsize=(10, 10))
    sim.plot2D(ax=f.gca(), output_plane=mp.Volume(center=mp.Vector3(0, 0, 0),
                                                  size=mp.Vector3(sx, 0, sz)))

    f.savefig("xz_plane.pdf", format="PDF")

    f = plt.figure(figsize=(10, 10))
    sim.plot2D(ax=f.gca(), output_plane=mp.Volume(center=mp.Vector3(0, 0, 0),
                                                  size=mp.Vector3(0, sy, sz)))

    f.savefig("yz_plane.pdf", format="PDF")

    h = mp.Harminv(mode, mp.Vector3(0, 0, 0), fcen, df)

    # h_displaced = mp.Harminv(mode, mp.Vector3(0.05, 0.1, 0.02), fcen, df)

    if eps_only:
        sim.run(mp.at_beginning(mp.output_epsilon),
                until=0)

    else:
        # Don't output eps anymore to save disk space
        sim.run(mp.after_sources(h),
                # mp.after_sources(h_displaced),
                until_after_sources=time_after_source)