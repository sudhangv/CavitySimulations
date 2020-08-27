import meep as mp
import argparse

def main(args):

    resolution = 30                       # pixels/um

    a_start = args.a_start                # starting periodicity
    a_end = args.a_end                    # ending periodicity
    s_cav = args.s_cav                    # cavity length
    r = args.r                            # hole radius  (units of a)
    h = args.hh                           # waveguide height
    w = args.w                            # waveguide width

    dair = 1.00                           # air padding
    dpml = 1.00                           # PML thickness

    Ndef = args.Ndef                      # number of defect periods
    a_taper = mp.interpolate(Ndef, [a_start,a_end])
    dgap = a_end-2*r*a_end

    Nwvg = args.Nwvg                      # number of waveguide periods
    sx = 2*(Nwvg*a_start+sum(a_taper))-dgap+s_cav
    sy = dpml+dair+w+dair+dpml
    sz = dpml+dair+h+dair+dpml

    cell_size = mp.Vector3(sx,sy,sz)
    boundary_layers = [mp.PML(dpml)]

    nSi = 3.45
    Si = mp.Medium(index=nSi)

    geometry = [mp.Block(material=Si, center=mp.Vector3(), size=mp.Vector3(mp.inf,w,h))]

    for mm in range(Nwvg):
        geometry.append(mp.Cylinder(material=mp.air, radius=r*a_start, height=mp.inf,
                                    center=mp.Vector3(-0.5*sx+0.5*a_start+mm*a_start,0,0)))
        geometry.append(mp.Cylinder(material=mp.air, radius=r*a_start, height=mp.inf,
                                    center=mp.Vector3(+0.5*sx-0.5*a_start-mm*a_start,0,0)))

    for mm in range(Ndef+2):
        geometry.append(mp.Cylinder(material=mp.air, radius=r*a_taper[mm], height=mp.inf,
                                    center=mp.Vector3(-0.5*sx+Nwvg*a_start+(sum(a_taper[:mm]) if mm>0 else 0)+0.5*a_taper[mm],0,0)))
        geometry.append(mp.Cylinder(material=mp.air, radius=r*a_taper[mm], height=mp.inf,
                                    center=mp.Vector3(+0.5*sx-Nwvg*a_start-(sum(a_taper[:mm]) if mm>0 else 0)-0.5*a_taper[mm],0,0)))

    lambda_min = 1.46        # minimum source wavelength
    lambda_max = 1.66        # maximum source wavelength
    fmin = 1/lambda_max
    fmax = 1/lambda_min
    fcen = 0.5*(fmin+fmax)
    df = fmax-fmin

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey, center=mp.Vector3())]

    symmetries = [mp.Mirror(mp.X,+1), mp.Mirror(mp.Y,-1), mp.Mirror(mp.Z,+1)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=geometry,
                        sources=sources,
                        dimensions=3,
                        symmetries=symmetries)

    #sim.run(mp.in_volume(mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx,sy,0)), mp.at_end(mp.output_epsilon, mp.output_efield_y)),
    #        mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(), fcen, df)),
    #        until_after_sources=500)
    sim.run(mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(), fcen, df)),
            until_after_sources=500)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a_start', type=float, default=0.43, help='starting periodicity (default: 0.43 um)')
    parser.add_argument('-a_end', type=float, default=0.33, help='ending periodicity (default: 0.33 um)')
    parser.add_argument('-s_cav', type=float, default=0.146, help='cavity length (default: 0.145 um)')
    parser.add_argument('-r', type=float, default=0.28, help='hole radius (default: 0.28 um)')
    parser.add_argument('-hh', type=float, default=0.22, help='waveguide height (default: 0.22 um)')
    parser.add_argument('-w', type=float, default=0.50, help='waveguide width (default: 0.50 um)')
    parser.add_argument('-Ndef', type=int, default=3, help='number of defect periods (default: 3)')
    parser.add_argument('-Nwvg', type=int, default=8, help='number of waveguide periods (default: 8)')
    args = parser.parse_args()
    main(args)
