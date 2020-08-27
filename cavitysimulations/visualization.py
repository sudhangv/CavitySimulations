import matplotlib.pyplot as plt
import meep as mp
plt.switch_backend('agg')

def plotter_simulation_cell(ax, sim, _output_plane=None, **kwargs):

    cell_size = sim.cell_size

    if _output_plane == "x":
        _slice = mp.Vector3(0, cell_size[1], cell_size[2])
        output_plane = mp.Volume(center=mp.Vector3(0, 0, 0), size=_slice)

    elif _output_plane == "y":
        _slice = mp.Vector3(cell_size[0], 0, cell_size[2])
        output_plane = mp.Volume(center=mp.Vector3(0, 0, 0), size=_slice)

    elif _output_plane == "z":
        _slice = mp.Vector3(cell_size[0], cell_size[1], 0)
        output_plane = mp.Volume(center=mp.Vector3(0, 0, 0), size=_slice)

    elif _output_plane == "2d":
        _slice = mp.Vector3(cell_size[0], cell_size[1])
        output_plane = mp.Volume(center=mp.Vector3(0, 0), size=_slice)

    else:
        output_plane = _output_plane

    sim.plot2D(ax=ax, output_plane=output_plane, **kwargs)


def visualize_sim_cell(sim):
    try:
        fig, ax = plt.subplots(figsize=(10, 10))
        plotter_simulation_cell(ax=ax, sim=sim, _output_plane="z")
        fig.savefig("xy_plane.pdf", format="PDF")

        fig, ax = plt.subplots(figsize=(10, 10))
        plotter_simulation_cell(ax=ax, sim=sim, _output_plane="y")
        fig.savefig("xz_plane.pdf", format="PDF")

        fig, ax = plt.subplots(figsize=(10, 10))
        plotter_simulation_cell(ax=ax, sim=sim, _output_plane="x")
        fig.savefig("yz_plane.pdf", format="PDF")
    except ValueError:
        print("Plotting failed. Try plotting a 2D structure.")
        fig, ax = plt.subplots(figsize=(10, 10))
        plotter_simulation_cell(ax=ax, sim=sim, _output_plane="2d")
        fig.savefig("yz_plane.pdf", format="PDF")
