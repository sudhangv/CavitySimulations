def get_gamma(freq, a, f_target = 1 / (1.54) * 1.01):
    '''
    w_target is in terms of 1/lambda * 1.01
    freq is in terms of 2pi * c/ a
    
    RETURNS: 
    
    The mirror strength for the input tuple (freq) of the dielectric and air band edge frequencies with
    f_target
    '''
    import math 
    
    freq[0] = convert_freq_to_Thz(freq[0], a)
    freq[1] = convert_freq_to_Thz(freq[1], a)
    f_target = convert_freq_to_Thz(1/1.54)* 1.01
     


    if f_target < freq[0] or f_target > freq[1] :  # if f_target outside bandgap
        return 0
    
    w_mid = (freq[0] + freq[1])/2            
    diff = freq[0] - freq[1]

    delta = 1 - (w_target/ (w_mid))

    gamma =  math.sqrt((( 0.5 * diff/ w_mid ) ** 2 - delta**2 ))
    
    return gamma

def get_freqs(hx = 0.24 , hy = 0.24, a = 0.33, w = 0.7, h = 0.22  ):
    '''
    Returns tuple of dielectric and air band edge frequencies for input parameters
    '''
    
    import meep as mp
    from meep import mpb

    res = 20
    mode = "zEyO"
    resolution = res  # pixels/a, taken from simpetus example
    
    a = round(a,3)        # units of um
    h = round(h, 3)         # units of um
    w = round(w, 3)         # units of um
    hx = round(hx, 3)
    hy = round(hy, 3)
    
    h = h/a          # units of "a"       
    w = w/a          # units of "a"
    hx = hx/a        # units of "a"
    hy = hy/a        # units of "a"


    nSi = 3.45
    Si = mp.Medium(index=nSi)

    geometry_lattice = mp.Lattice(size=mp.Vector3(1,4,4)) # dimensions of lattice taken from simpetus example

    geometry = [ mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf,w ,h ), material=Si),
             mp.Ellipsoid(material=mp.air,
             center=mp.Vector3(),
             size=mp.Vector3(hx,hy,mp.inf)) ]

    k_points = [mp.Vector3(0.5, 0, 0)]
    num_bands = 2 # from simpetus example

    ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
                        geometry=geometry,
                        k_points=k_points,
                        resolution=resolution,
                        num_bands=num_bands)

    if mode == "te":
        
        ms.run_te() # running for all modes and extracting parities
        
    if mode == "zEyO":
        
        ms.run_yodd_zeven()
    
    return ms.freqs

def convert_freq_to_Thz(freq, a = 0):
    
#----------- Converts a MEEP frequency ( units of 2pi*c/a) to frequency (not angular velocity) in units of Thz--------#
    if a != 0:
        return ( freq * 3 / a * (10**2))
    else:
        return ( freq * 3 * 10**2 )
    
def get_freq_in_Thz(hx, hy, a, w, h = 0.22):
    return convert_freq_to_Thz(np.array(get_freqs(hx, hy, a, w , h)), a )

def get_freqs_interpolate(hx = 0.24 , hy = 0.24, a = 0.33, wy = 0.7, h = 0.22  ):
    '''
    Useless 
    '''
    
    import meep as mp
    from meep import mpb

    mode = "zEyO"
    resolution = 20  # pixels/a, taken from simpetus example
    
    a = round(a, 3)        # units of um
    h = round(h, 3)         # units of um
    w = round(wy, 3)         # units of um
    hx = round(hx, 3)
    hy = round(hy, 3)
    
    h = h/a          # units of "a"       
    w = w/a          # units of "a"
    hx = hx/a        # units of "a"
    hy = hy/a        # units of "a"


    nSi = 3.45
    Si = mp.Medium(index=nSi)

    geometry_lattice = mp.Lattice(size=mp.Vector3(1,4,4)) # dimensions of lattice taken from simpetus example

    geometry = [ mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf,w ,h ), material=Si),
             mp.Ellipsoid(material=mp.air,
             center=mp.Vector3(),
             size=mp.Vector3(hx,hy,mp.inf)) ]
    
    num_k = 20  # from simpetus example, no. of k_points to evaluate the eigen frequency at
    k_points = mp.interpolate(num_k, [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0)])
    
    num_bands = 2 # from simpetus example

    ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
                        geometry=geometry,
                        k_points=k_points,
                        resolution=resolution,
                        num_bands=num_bands)

    if mode == "te":
        
        ms.run_te() # running for all modes and extracting parities
        
    if mode == "zEyO":
        
        ms.run_yodd_zeven()
    
    return ms.freqs

def quick_plot(x, y, title = "Unspecified", xlabel = "x", ylabel = "y"):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.plot(x, y, 'x')