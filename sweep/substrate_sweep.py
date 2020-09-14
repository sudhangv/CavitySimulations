#def do_the_sweep():
import math
import numpy as np 
import matplotlib.pyplot as plt
import h5py
from meep import mpb
from math import sqrt, pi

from sweep_util import *


data_file = "check.hdf5"      # Name of the file where the data will be stored
param_file = "param.txt"      # Name of the file where the parameters of interest will be stored     

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
w_max = 0.75        #  upper limit of w 

hx_min = 0.05        # lower limit of the sweep of a
hy_min = 0.1        #  lower limit of hy 

hx_max = a_max - 0.07
hy_max = w_max - 0.1  
#------------------ Geometry Characteristics ---------------------#

wvg_height = 0.22

SUBSTRATE = True


f_target = 1/1.54
f_target_Thz = convert_freq_to_Thz(f_target) * 1.01

parameters = []

with h5py.File(data_file, 'w') as f:

    dt = h5py.special_dtype(vlen=np.float32)
    gamma_max = 0        # arbitrary small value
    mirror_strength = []   

    j = len(np.arange(w_min, w_max , del_w))
    k = len(np.arange(a_min , a_max, del_a))
    l = len(np.arange(hy_min, hy_max , del_hy))
    m = len(np.arange(hx_min, hx_max, del_hx ))

    dset_gamma = f.create_dataset("gamma", (j,k,l,m))
    dset_gamma[:,:,:,:] =  np.zeros((j,k,l,m))

    dset_freq = f.create_dataset("freq", (j,k,l,m), dtype=dt)
    #dset_freq[:,:,:,:] =  np.zeros((j,k,l,m))


    index = []
    index_count = 0
      # stores parameters for optimal value of gamma


    #---------------------------#
    #       WIDTH LOOP          
    #---------------------------#
    for w in np.arange(w_min, w_max , del_w):

        #w = round(w,3)

        freq1_Thz = convert_freq_to_Thz(get_freqs(hx_min, hy_min, a_max, w, wvg_height, substrate = SUBSTRATE), a_max)  # getting lowest possible frequencies for width w 
        print(" -------------------- w: {} loop ----------------------".format(w))
        if ((freq1_Thz[0]  > f_target_Thz)):

            continue
        else:

        #---------------------------#
        #   LATTICE CONSTANT LOOP          
        #---------------------------#

            for a in np.arange(a_min , a_max, del_a):
                print("------------------------a: {} loop-------------------------".format(a))
                #a = round(a,3)

                freq2_Thz = convert_freq_to_Thz(get_freqs(hx_min, hy_min, a, w, wvg_height,substrate = SUBSTRATE), a)  # getting lowest possible frequences for (w,a)

                if ( freq2_Thz[0] > f_target_Thz):
                    continue

                #---------------------------#
                #       HY LOOP          
                #---------------------------#

                for hy in np.arange(hy_min, w - 0.1 , del_hy):
                    print("-----------------------hy {} loop----------------------".format(hy))
                    #hy = round(hy,3)

                    freq3_Thz = convert_freq_to_Thz(get_freqs(hx_min, hy, a, w, wvg_height, substrate = SUBSTRATE), a)  # getting lowest possible frequences for (w,a, hy)

                    if ( (freq3_Thz[0]  > f_target_Thz) ):
                        continue

                    #---------------------------#
                    #       HX LOOP          
                    #---------------------------#

                    for hx in np.arange(hx_min, a - 0.07, del_hx ):
                        count  = 0
                        print(" ---------------- hx : {} loop ---------------------".format(hx))
    #                         a = round(a,3)
    #                         hy = round(hy,3)
    #                         hx = round(hx,3)
    #                         w = round(w,3)

                        freq4_Thz = convert_freq_to_Thz(get_freqs(hx, hy, a, w, wvg_height, substrate = SUBSTRATE), a )

                        if ( freq4_Thz[0] > f_target_Thz):  # if w_target is outside the bandgap for 2 consecutive runs, break outside the loop
                            count = count + 1

                            if count == 2 :
                                break

                        else:

                            if (f_target_Thz < freq4_Thz[1])  and (f_target_Thz > freq4_Thz[0]):  # final check to see that the target frequency is in the bandgap

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

                                dset_freq[int((w - w_min) / del_w + 0.1), 
                                     int((a - a_min) / del_a + 0.1), 
                                     int((hy - hy_min) / del_hy + 0.1), 
                                     int((hx - hx_min) / del_hx + 0.1)] = np.array((freq4_Thz[0], freq4_Thz[1])) 
                                if gamma > gamma_max:
                                    print(" ------------------- new gamma max ------------------- at hx = {}, hy = {}, a = {}, w = {}, gamma = {}".format(hx,hy,a,w, gamma))
                                    gamma_max = gamma
                                    parameters.append((round(hx, 4), 
                                                       round(hy, 4), 
                                                       round(a,  4), 
                                                       round(w,  4), 
                                                       freq4_Thz , 
                                                       round(gamma_max,5))) 


                            else:
                                continue

    with open(param_file, "w") as file1: 
        for parameter in parameters:
    # Writing data to a file 
            file1.write("hx = {}, hy = {}, a = {}, w = {}, freqs = {}, gamma = {}".format(*parameter)) 
            file1.write("\n")
