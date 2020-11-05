import math
import numpy as np 
import matplotlib.pyplot as plt
import h5py
from meep import mpb
from math import sqrt, pi

from sweep_util import *

'''
---------------- FORMAT -------------------------------------------------------
The h5py file created from this script contains two datasets : 
(1) dset_gamma --> stores the mirror strength for parameters [w, a, hy, hx]
(2) dset_freq_lower  --> stores the bandedge frequencies for parameters [w, a, hy, hx]
(3) dset_freq_upper  --> stores the bandedge frequencies for parameters [w, a, hy, hx]
-------------------------------------------------------------------------------
'''

data_file = "perturb_sub_yO_220.hdf5"             # Name of the file where the data will be stored
param_file = "perturb_sub_yO_220_param.txt" 


# Name of the file where the parameters of interest will be stored     
SUBSTRATE = True
wvg_height = 0.22
mode = "yO"

#-----------------DEFAULTS-----------------------#
#     del_a = 0.001      
#     del_hy = 0.025
#     del_hx = 0.025 
#     del_w = 0.05
#------------------------------------------------#

#--------------- Increments --------------------#
del_a = 0.001
del_hy = 0.025
del_hx = 0.025
del_w = 0.05

#------------------ Ranges ---------------------#
a_min = 0.25
a_max = 0.45        # upper limit of the sweep of a 

w_min = 0.65         #  lower limit of w 
w_max = 0.70        #  upper limit of w 

hx_min = 0.125        # lower limit of the sweep of a
hy_min = 0.200        #  lower limit of hy 

hx_max = a_max - 0.07
hy_max = w_max - 0.1  
#------------------ Geometry Characteristics ---------------------#
perturb_range = 0.05
target_wvl = 1.54        # vaccum wavelength ( in um ) of the unperturbed cavity design

f_perturb_lower = 1 / (target_wvl + perturb_range )           # target_f - perturbation
f_perturb_upper = 1 / (target_wvl - perturb_range )           # target_f + perturbation

f_target = 1/target_wvl
f_target_Thz =  convert_freq_to_Thz(f_target) * 1.01

f_perturb_lower_Thz =  convert_freq_to_Thz(f_perturb_lower) * 1.01
f_perturb_upper_Thz =  convert_freq_to_Thz(f_perturb_upper) * 1.01

parameters = []

e = 0.0001

run_count = 0

with h5py.File(data_file, 'w') as f:

    #dt = h5py.special_dtype(vlen=np.float32)
    gamma_max = 0        # arbitrary small value
    mirror_strength = []   
    
    j = len(np.arange(w_min, w_max + e , del_w))
    k = len(np.arange(a_min , a_max + e, del_a))
    l = len(np.arange(hy_min, hy_max + e , del_hy))
    m = len(np.arange(hx_min, hx_max + e, del_hx )) 
    
    #breakpoint()
    
    dset_gamma = f.create_dataset("gamma", (j,k,l,m))
    dset_gamma.attrs['wvg_height'] = wvg_height
    dset_gamma.attrs['substrate'] = SUBSTRATE
    dset_gamma.attrs["mode"] = "yO"
    
    dset_gamma[:,:,:,:] =  np.zeros((j,k,l,m))

    #dset_freq = f.create_dataset("freq", (j,k,l,m), dtype=dt)
    
    dset_freq_lower =  f.create_dataset("freq_lower", (j,k,l,m))
    dset_freq_upper =  f.create_dataset("freq_upper", (j,k,l,m))
    
    dset_freq_lower[:,:,:,:] = np.zeros((j,k,l,m))
    dset_freq_upper[:,:,:,:] = np.zeros((j,k,l,m))


    index = []
    index_count = 0
      # stores parameters for optimal value of gamma


    #---------------------------#
    #       WIDTH LOOP          
    #---------------------------#
    for w in np.arange(w_min, w_max  + e , del_w):

        #w = round(w,3)

        freq1_Thz = convert_freq_to_Thz(get_freqs(hx = hx_min, hy = hy_min, a = a_max, w = w, 
                                                  h = wvg_height,  # getting lowest possible frequencies for width w 
                                                  substrate = SUBSTRATE, mode = mode), 
                                        a_max)
        run_count += 1
       # breakpoint()
#         dset_freq[int((w - w_min) / del_w + 0.1), 
#                   int((a_max - a_min) / del_a + 0.1), 
#                   int((hy_min - hy_min) / del_hy + 0.1), 
#                   int((hx_min - hx_min) / del_hx + 0.1)] = np.array((freq1_Thz[0], freq1_Thz[1]))
        
        dset_freq_lower[int((w - w_min) / del_w + 0.1), 
                        int((a_max - a_min) / del_a + 0.1), 
                        int((hy_min - hy_min) / del_hy + 0.1), 
                        int((hx_min - hx_min) / del_hx + 0.1)] = freq1_Thz[0]
        
        
        dset_freq_upper[int((w - w_min) / del_w + 0.1), 
                        int((a_max - a_min) / del_a + 0.1), 
                        int((hy_min - hy_min) / del_hy + 0.1), 
                        int((hx_min - hx_min) / del_hx + 0.1)] = freq1_Thz[1]

        print(" -------------------- w: {} loop ----------------------".format(w))
        #if ((freq1_Thz[0]  > f_target_Thz)):
        if ( freq1_Thz[0] > f_perturb_upper_Thz):

            continue
        else:

        #---------------------------#
        #   LATTICE CONSTANT LOOP          
        #---------------------------#

            for a in np.arange(a_min , a_max + e, del_a):
                print("------------------------a: {} loop-------------------------".format(a))
                #a = round(a,3)

                freq2_Thz = convert_freq_to_Thz(get_freqs(hx = hx_min, hy = hy_min, a = a, w = w,  h = wvg_height,
                                                          substrate = SUBSTRATE, mode = mode), 
                                                a)  # getting lowest possible frequences for (w,a)
                run_count += 1

                dset_freq_lower[int((w - w_min) / del_w + 0.1), 
                               int((a - a_min) / del_a + 0.1), 
                               int((hy_min - hy_min) / del_hy + 0.1), 
                               int((hx_min - hx_min) / del_hx + 0.1)] = freq2_Thz[0]
        
                dset_freq_upper[int((w - w_min) / del_w + 0.1), 
                               int((a - a_min) / del_a + 0.1), 
                               int((hy_min - hy_min) / del_hy + 0.1), 
                               int((hx_min - hx_min) / del_hx + 0.1)] = freq2_Thz[1]

                
                #if ( freq2_Thz[0] > f_target_Thz):
                if (freq2_Thz[0] > f_perturb_upper_Thz):
                    continue

                #---------------------------#
                #       HY LOOP          
                #---------------------------#

                for hy in np.arange(hy_min, w - 0.1 + e , del_hy):
                    print("-----------------------hy {} loop----------------------".format(hy))
                    #hy = round(hy,3)

                    freq3_Thz = convert_freq_to_Thz(get_freqs(hx = hx_min, hy = hy, a = a, w = w, h = wvg_height, 
                                                              substrate = SUBSTRATE, mode = mode)
                                                    , a)  # getting lowest possible frequences for (w,a, hy)
                    run_count += 1
#                     dset_freq[int((w - w_min) / del_w + 0.1), 
#                               int((a - a_min) / del_a + 0.1), 
#                               int((hy - hy_min) / del_hy + 0.1), 
#                               int((hx_min - hx_min) / del_hx + 0.1)] = np.array((freq3_Thz[0], freq3_Thz[1]))
                    
                    dset_freq_lower[int((w - w_min) / del_w + 0.1), 
                                   int((a - a_min) / del_a + 0.1), 
                                   int((hy - hy_min) / del_hy + 0.1), 
                                   int((hx_min - hx_min) / del_hx + 0.1)] = freq3_Thz[0]
        
                    dset_freq_upper[int((w - w_min) / del_w + 0.1), 
                                   int((a - a_min) / del_a + 0.1), 
                                   int((hy - hy_min) / del_hy + 0.1), 
                                   int((hx_min - hx_min) / del_hx + 0.1)] = freq3_Thz[1]
            
                    #if ( (freq3_Thz[0]  > f_target_Thz) ):
                    if ( freq3_Thz[0] > f_perturb_upper_Thz):
                        continue

                    #---------------------------#
                    #       HX LOOP          
                    #---------------------------#

                    for hx in np.arange(hx_min, a - 0.07 + e, del_hx ):
                        count  = 0
                        print(" ---------------- hx : {} loop ---------------------".format(hx))
                        a = round(a,4)
                        hy = round(hy,4)
                        hx = round(hx,4)
                        w = round(w,4)

                        freq4_Thz = convert_freq_to_Thz(get_freqs(hx = hx, hy = hy, a = a, w = w, h = wvg_height, 
                                                                  substrate = SUBSTRATE, mode = mode)
                                                        , a )
                        run_count += 1
#                         dset_freq[int((w - w_min) / del_w + 0.1), 
#                                   int((a - a_min) / del_a + 0.1), 
#                                   int((hy - hy_min) / del_hy + 0.1), 
#                                   int((hx - hx_min) / del_hx + 0.1)] = np.array((freq4_Thz[0], freq4_Thz[1]))

                        dset_freq_lower[int((w - w_min) / del_w + 0.1), 
                                       int((a - a_min) / del_a + 0.1), 
                                       int((hy - hy_min) / del_hy + 0.1), 
                                       int((hx - hx_min) / del_hx + 0.1)] = freq4_Thz[0]
        
                        dset_freq_upper[int((w - w_min) / del_w + 0.1), 
                                       int((a - a_min) / del_a + 0.1), 
                                       int((hy - hy_min) / del_hy + 0.1), 
                                       int((hx - hx_min) / del_hx + 0.1)] = freq4_Thz[1]
            
                        #if ( freq4_Thz[0] > f_target_Thz):                   
                        if ( freq4_Thz[0] > f_perturb_upper_Thz):
                            
                            count = count + 1  # if w_target is outside the bandgap for 2 consecutive runs, break outside the loop
                            if count == 3 :    
                                break

                        #else:

    #                             if (f_target_Thz < freq4_Thz[1])  and (f_target_Thz > freq4_Thz[0]):  # final check to see that the target frequency is in the bandgap
                        #if (f_perturb_lower_Thz < freq4_Thz[1])  and (f_perturb_upper_Thz > freq4_Thz[0]):
                        
                        
                        if (f_target_Thz < freq4_Thz[1])  and (f_target_Thz > freq4_Thz[0]):
                            
                                print(" ------------------- new gamma ------------------- at hx = {}, hy = {}, a = {}, w = {}".format(hx,hy,a,w))
                                
                                f_mid = (freq4_Thz[0] + freq4_Thz[1])/2            
                                diff = freq4_Thz[0] - freq4_Thz[1]
                                delta = 1 - (f_target_Thz/ f_mid)

                                gamma =  math.sqrt(abs(( 0.5 * diff/ f_mid ) ** 2 - delta**2 ))

                                mirror_strength.append(gamma)
                                index_count = index_count + 1
                                index.append(index_count)

                                dset_gamma[int((w - w_min) / del_w + 0.1), 
                                             int((a - a_min) / del_a + 0.1), 
                                             int((hy - hy_min) / del_hy + 0.1), 
                                             int((hx - hx_min) / del_hx + 0.1) ] = round(gamma,4)

#                         dset_freq[int((w - w_min) / del_w + 0.1), 
#                              int((a - a_min) / del_a + 0.1), 
#                              int((hy - hy_min) / del_hy + 0.1), 
#                              int((hx - hx_min) / del_hx + 0.1)] = np.array((freq4_Thz[0], freq4_Thz[1])) 

                                if gamma > gamma_max:
        
                                    print(" ------------------- new gamma max ------------------- at hx = {}, hy = {}, a = {}, w = {}, gamma = {}".format(hx,hy,a,w, gamma))
                                    gamma_max = gamma
                
        #                             parameters.append((round(hx, 4), 
        #                                                round(hy, 4), 
        #                                                round(a,  4), 
        #                                                round(w,  4), 
        #                                                freq4_Thz , 
        #                                                round(gamma_max,5))) 
        
                                    parameters.append((round(hx, 4), 
                                                        round(hy, 4), 
                                                        round(a,  4), 
                                                        round(w,  4), 
                                                        round(gamma_max,5),
                                                        round(freq4_Thz[0], 4),
                                                        round(freq4_Thz[1], 4))
                                                      )



#                             else:
#                                 continue

    with open(param_file, "w") as file1: 
        file1.write("wvg_height = {} , mode = {}, SUBSRATE = {}".format(wvg_height, mode, SUBSTRATE))
        file1.write("\n")
        for parameter in parameters:
    # Writing data to a file 
            file1.write("hx = {}, hy = {}, a = {}, w = {}, gamma = {}, freqs = {} , {} ".format(*parameter, )) 
            file1.write("\n")
        file1.write(f'run_count : {run_count}')

